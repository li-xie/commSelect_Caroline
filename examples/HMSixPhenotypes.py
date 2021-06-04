#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:43:14 2020

@author: carolc24
"""

import numpy as np
import commSelect
import os
import time as tm

num_wells = 10; 
num_cycles = 1;
top_tier = 1.; #top percent of adults chosen for reproduction
parallel = False;
run_id = "default"

# time params
dt = 0.05;
cycle_duration = 17;
nsteps = int(np.ceil(cycle_duration/dt));
time = np.arange(0, cycle_duration, dt);

#mutable params
fp_init = 0.1;
K_MR_init = 1.;
K_MB_init = 1./3 * 500;
b_Mmax_init = 0.7 / 1.2;
b_Hmax_init = 0.3 / 1.2;
K_HR_init = 1.;

fp_bound = 1.;
K_MR_bound = 1./3;
K_MB_bound = 100./3;
K_HR_bound = 0.2;
b_Mmax_bound = 0.7;
b_Hmax_bound = 0.8;

#note that the Monod constants are saved as their reciprocals.
#This is because the mutation function biases towards deleterious mutations
#and each trait needs to be more advantageous at higher values to
#preserve this effect. A higher Monod constant represents slower growth at
#low resource.
traits_manu_bound = np.array([fp_bound, b_Mmax_bound, 1./K_MR_bound, 1./K_MB_bound]);
traits_help_bound = np.array([b_Hmax_bound, 1./K_HR_bound]);
traits_bound = [traits_manu_bound, traits_help_bound];

traits_manu_lowerbound = np.array([0., 0., 0., 0.]);
traits_help_lowerbound = np.array([0., 0.]);
traits_lowerbound = [traits_manu_lowerbound, traits_help_lowerbound];

#constant params
c_BM = 1./3;
c_RM = 10**-4;
c_RH = 10**-4;
r_P = 1.;
d_M = 3.5 * 10**-3;
d_H = 1.5 * 10**-3;
R_init = 10.;
M_init = 60;
H_init = 40;

#mutation params (required)
p_mut = 2 * 10**-3;  # mutation rate
frac_null = 0.5; # fraction of mutations that are null
sp0 = 0.05; # 'positive' mutation factor
sn0 = 0.067; # 'negative' mutation factor
mut_params = [p_mut, frac_null, sp0, sn0];

nD = 100; # dilution factor for fixed fold pipetting
BM_target = 100.; # target biomass for regular pipetting
newborns_per_adult = int(np.floor(num_wells / top_tier));

#system of diffeqs that describe community dynamics
#arg t: timestep (required for solve_ivp)
#arg y: state vector (required for solve_ivp)
#arg manuCell: manufacturer CellType
#arg helpCell: helper CellType
#returns: dydt
def RBPFG_prime(t, y, manuCell, helpCell):
    
    Bio_M = manuCell.biomass(); #vector (multiple variants of M)
    Bio_H = helpCell.biomass(); #vector (multiple variants of H)
    
    #making sure they are 1D vectors
    if np.ndim(Bio_M) > 1:
        Bio_M = Bio_M[:,0];
    if np.ndim(Bio_H) > 1:
        Bio_H = Bio_H[:,0];
    
    #trait vectors of length len(Bio_M)
    fp = manuCell.traits[:,0];
    b_Mmax = manuCell.traits[:,1];
    K_MR = manuCell.traits[:,2];
    K_MB = manuCell.traits[:,3];
    
    #trait vectors of length len(Bio_H)
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
    b_Hcoef = R * K_HR / (R * K_HR + 1); # one value for each variant
    R_M = R * K_MR;
    B_M = B * K_MB;
    b_Mcoef = (R_M * B_M) / (R_M + B_M) * (1. / (1 + R_M) + 1. / (1 + B_M));
    
    #rate of change of metabolites
    #sum up all the metabolite production and consumption by all variants
    R_prime = -sum(b_Hmax * b_Hcoef * c_RH * np.transpose(Bio_H)) \
              -sum(b_Mmax * b_Mcoef * c_RM * np.transpose(Bio_M));
    B_prime = sum(b_Hmax * b_Hcoef * np.transpose(Bio_H)) \
             -sum(b_Mmax * b_Mcoef * c_BM * np.transpose(Bio_M));
    P_prime = sum(b_Mmax * b_Mcoef * r_P * fp * np.transpose(Bio_M));
     
    #cell birth rates per biomass
    M_prime_coef = (1 - fp) * b_Mmax * b_Mcoef;
    H_prime_coef = b_Hmax * b_Hcoef;
    
    metabol_prime = [R_prime, B_prime, P_prime];  # order is same as before
    bio_prime = np.append(M_prime_coef, H_prime_coef).tolist();
    
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
    
    #initialize output dir
    init_seed = np.random.randint(1e7);
    np.random.seed(init_seed);
    run_id += str(init_seed);
    os.mkdir(run_id);
    
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
    
    # make (num_wells) copies of CellTypes
    newbDataAll = [];
    for i in range(num_wells):
        newbDataAll += [list(map(commSelect.CellType.CellType.copy, newbData))];
 
    print("first newborns initialized");
    
    #cycle through
    for c in range(num_cycles):
        print("cycle %d..." % c);
            
        #save newborn data
        if (c % 10 == 0):
            commSelect.save.saveNewborns(newbDataAll, run_id + "/newb_%d.txt" % c);
        
        #mature
        #args for mature function: 
        #newbDataAll (list of newborn communities)
        #ic_metabol (initial metabolite concs)
        #RBPFG_prime (diffeq method)
        #nsteps (num timesteps per maturation phase)
        #dt (timestep length)
        #mut_params (mutation parameters)
        #parallel (boolean, True if using parallel processing)
            
        #returns:
        #adultDataAll (list of adult communities)
        #data (time series data of biomass and metabolites for each community)
        #data structure: [community index, 'data type', time]
        #'data type' is total biomass of each cell type and conc of each metabolite
        #biomass comes first (M, H) then metabolites (R, B, P)
        #so for example data[:,1,-1] is the final H biomass of each community
        start = tm.time();
        adultDataAll, data = commSelect.mature.mature(newbDataAll, ic_metabol, \
                                RBPFG_prime, nsteps, dt, mut_params, parallel);
        end = tm.time();
        print("Cycle completed in %.2f secs" % (end - start));
            
        #save adult data
        if (c % 10 == 0):
            commSelect.save.saveAdults(adultDataAll, data, run_id + "/adult_%d.txt" % c);       

        #evaluate community function for each adult
        P = commFunc(data);
        P_sorted = np.argsort(P);
        
        #re-initialize newborns
        newbDataAll = [];
        for i in range(num_wells):
            newbDataAll += [list(map(commSelect.CellType.CellType.copy, newbData))];
            
        #reproduce based on community function
        #arg adultDataAll: adults to reproduce
        #arg newbDataAll: newborns to fill with cells
        #P_sorted: indices of adults sorted from lowest to highest P
        #num_wells: number of wells
        #BM_target: target biomass per newborn
        #newborns_per_adult: max number of newborns to be made with 1 adult
        #cs: True if cell sorting, False if pipetting
        newbDataAll = commSelect.reproduce.reproduce(adultDataAll, newbDataAll, P_sorted, \
                                num_wells, BM_target, newborns_per_adult, False);

