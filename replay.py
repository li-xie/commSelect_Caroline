import numpy as np
import pickle
import commSelect
import copy

# Replay Cycle 1. Load the random number generator and selected Adults from Cycle 0
# Reproduce the selected Adults from Cycle 0, then replay Cycle 1.
# The Newborns and Adults from the replay of Cycle 1 are newbDataAll_replay and
# adultDataAll_replay, respectively. They are identical to those from the original simulation.
    
num_wells = 10; 
parallel = False;
fname = "Results"

# time params
dt = 0.05;
cycle_duration = 17;
# number of integration steps for mature cycle, positive int
nsteps = int(np.ceil(cycle_duration/dt));

#mutable params initial conditions, non-negative float
fp_init = 0.1;
K_MR_init = 1.;
K_MB_init = 1./3 * 500;
b_Mmax_init = 0.7 / 1.2;
b_Hmax_init = 0.3 / 1.2;
K_HR_init = 1.;

#mutable params upper bounds, non-negative float
fp_bound = 1.;
K_MR_bound = 1./3;
K_MB_bound = 100./3;
K_HR_bound = 0.2;
b_Mmax_bound = 0.7;
b_Hmax_bound = 0.3;

#note that the Monod constants are saved as their reciprocals.
#This is because the mutation function biases towards deleterious mutations
#and each trait needs to be more advantageous at higher values to
#preserve this effect. A higher Monod constant represents slower growth at
#low resource.
traits_manu_bound = np.array([fp_bound, b_Mmax_bound, 1./K_MR_bound, 1./K_MB_bound]);
traits_help_bound = np.array([b_Hmax_bound, 1./K_HR_bound]);
traits_bound = [traits_manu_bound, traits_help_bound];

#mutable params lower bounds, non-negative float
traits_manu_lowerbound = np.array([0., 0., 0., 0.]);
traits_help_lowerbound = np.array([0., 0.]);
traits_lowerbound = [traits_manu_lowerbound, traits_help_lowerbound];

#constant growth + metabolite params, any type
c_BM = 1./3;
c_RM = 10**-4;
c_RH = 10**-4;
r_P = 1.;

#death rates, non-negative float
d_M = 3.5 * 10**-3;
d_H = 1.5 * 10**-3;

#metabolite initial conc, non-negative float
R_init = 1.;

#cell initial pop size, non-negative int
M_init = 60;
H_init = 40;

#mutation params, non-negative float
p_mut = 2 * 10**-3;  # mutation rate
frac_null = 0.5; # fraction of mutations that are null (set trait to 0.)
sp0 = 0.05; # 'positive' mutation factor
sn0 = 0.067; # 'negative' mutation factor
mut_params = [p_mut, frac_null, sp0, sn0];

#reproduction params
nD = 100; # dilution factor for fixed fold pipetting, positive int
BM_target = 100.; # target biomass for regular pipetting, positive float
top_tier = 1.; #top percent of adults chosen for reproduction, float >= 1.
newborns_per_adult = int(np.floor(num_wells / top_tier)); # max newborns per adult, positive int

#system of diffeqs that describe community dynamics
#arg t: timestep (required for solve_ivp)
#arg y: state vector (required for solve_ivp)
#arg manuCell: manufacturer CellType
#arg helpCell: helper CellType
#returns: dydt
def RBPFG_prime(t, y, manuCell, helpCell):
    
    #biomass vectors
    #length depends on how many mutants of each type there are
    Bio_M = manuCell.biomass(); #vector (multiple variants of M)
    Bio_H = helpCell.biomass(); #vector (multiple variants of H)
    
    #trait vectors with same length as Bio_M
    fp = manuCell.traits[:,0];
    b_Mmax = manuCell.traits[:,1];
    K_MR = manuCell.traits[:,2];
    K_MB = manuCell.traits[:,3];
    
    #trait vectors with same length as Bio_H
    b_Hmax = helpCell.traits[:,0];
    K_HR = helpCell.traits[:,1];
    
    #preventing division error
    K_HR[K_HR < 1e-9] = 1e-9;
    K_MR[K_MR < 1e-9] = 1e-9;
    K_MB[K_MB < 1e-9] = 1e-9;
    
    #metabolites (P is not needed for diffeqs)
    R = y[0];
    B = y[1];
        
    # birth rate coefficients 
    # calculated separately for each mutant based on their traits
    #vector with same length as Bio_H
    b_Hcoef = R * K_HR / (R * K_HR + 1); # one value for each variant
    
    #vectors with same length as Bio_M
    R_M = R * K_MR;
    B_M = B * K_MB;
    b_Mcoef = (R_M * B_M) / (R_M + B_M) * (1. / (1 + R_M) + 1. / (1 + B_M));
    
    #rate of change of metabolites
    #find the production/consumption rate of each metabolite for each mutant
    #by multiplying vectors elementwise
    #then use sum() to add up all fluxes and find overall rate of change
    R_prime = -sum(b_Hmax * b_Hcoef * Bio_H) * c_RH \
              -sum(b_Mmax * b_Mcoef * Bio_M) * c_RM;
    B_prime = sum(b_Hmax * b_Hcoef * Bio_H) \
             -sum(b_Mmax * b_Mcoef * Bio_M) * c_BM;
    P_prime = sum(b_Mmax * b_Mcoef * fp * Bio_M) * r_P;
     
    #cell birth rates per biomass for each mutant
    #vector with same length as Bio_M
    M_prime_coef = (1 - fp) * b_Mmax * b_Mcoef;
    #vector with same length as Bio_H
    H_prime_coef = b_Hmax * b_Hcoef;
    
    metabol_prime = [R_prime, B_prime, P_prime];  # order is same as before
    bio_prime = np.append(M_prime_coef, H_prime_coef).tolist();
    
    #dydt (must be list, not numpy array)
    dydt = metabol_prime + bio_prime; #metabolites first, then cells

    return dydt;

#community function = Pfinal
#arg data: 3D array of time series data of all adults
#axis 0: community index
#axis 1: data type (M total biomass, H total biomass, R, B, P)
#axis 2: time
def commFunc(data):
    return data[:,-1,-1]; #P conc at the end of maturation
    #for each community
        
if __name__ == "__main__":
    # initialize community structure
    ic_traits_manu = np.array([[fp_init, b_Mmax_init, 1./K_MR_init, 1./K_MB_init]]);
    ic_L_manu = np.array([1.]); #cell length
    ic_N_manu = np.array([M_init]); #number of cells
    ic_traits_help = np.array([[b_Hmax_init, 1./K_HR_init]]);
    ic_L_help = np.array([1.]);
    ic_N_help = np.array([H_init]);
    N0 = ic_L_manu * ic_N_manu + ic_L_help * ic_N_help; #total number of cells
    
    #args for construction of CellType: traits, L, N, death rate
    manu = commSelect.CellType.CellType(ic_traits_manu, ic_L_manu, ic_N_manu, d_M * dt);
    manu.traits_bound = traits_manu_bound;
    manu.traits_lowerbound = traits_manu_lowerbound;
    
    helper = commSelect.CellType.CellType(ic_traits_help, ic_L_help, ic_N_help, d_H * dt);
    helper.traits_bound = traits_help_bound;
    helper.traits_lowerbound = traits_help_lowerbound;
    
    #cell types and metabolites
    newbData = [manu, helper]; # cell types order: M, H
    ic_metabol = [R_init, 0., 0.]; # metabolites order: R, B, P
    
    # Up to this point, the codes are the same as in HMSixPhenoTypes.py.
    
    # Replay the simulation from reproducing selected Adults in Cycle 0. 
    # The replayed results from Cycle 1 should be the same as the saved results.
    # newbDataAll_replay = newbDataAll_simu
    # adultDataAll_replay = adultDataAll_simu
    c = 1 # Cycle to be replayed
    # load rng_main from Cycle 0, used for generating random number seeds and reproduction
    with open(fname + f'/rng_main{c-1}.pkl', 'rb') as f:
        rng_main=pickle.load(f)
    rng_seeds = rng_main.integers(0, 2**32, num_wells)
    with open(fname + f'/adultDataSel{c-1}.pkl', 'rb') as f:
        adultSel=pickle.load(f)
    with open(fname + f'/newborns{c}.pkl', 'rb') as f:
        newbDataAll_simu=pickle.load(f)
    with open(fname + f'/adults{c}.pkl', 'rb') as f:
        adultDataAll_simu=pickle.load(f)
        
    #re-initialize newborns
    newbDataAll = [];
    for i in range(num_wells):
        newbDataAll += [list(map(commSelect.CellType.CellType.copy, newbData))];
    P_sorted = np.arange(9,-1,-1,dtype=int)
    newbDataAll, win_inds = commSelect.reproduce.reproduce(adultSel, newbDataAll, P_sorted, \
                            num_wells, BM_target, newborns_per_adult, rng_main, False);
    newbDataAll_replay = copy.deepcopy(newbDataAll)
    rng_seeds = rng_main.integers(0, 2**32, num_wells)
    adultDataAll_replay, data = commSelect.mature.mature(newbDataAll, ic_metabol, \
                            RBPFG_prime, nsteps, dt, mut_params, rng_seeds, parallel);