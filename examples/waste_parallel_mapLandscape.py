#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 21:10:54 2021

@author: carolc24
"""


import numpy as np
import commSelect
import matplotlib.pyplot as plt

#experiment parameters
num_cycles = 1;
num_wells = 100;
top_tier = 1; # top dog
parallel = False; # no parallel processing

#integration params
dt = 0.05;
cycle_duration = 15;
nsteps = int(np.ceil(cycle_duration/dt));
time = np.arange(0, cycle_duration, dt);

#mutation params (required)
p_mut = 0 * 10**-3;  # mutation rate
frac_null = 0.5; # fraction of mutations that are null
sp0 = 0.05; # 'positive' mutation factor
sn0 = 0.05; # 'negative' mutation factor
mut_params = [p_mut, frac_null, sp0, sn0];

#reproduction params
nD = 100; # dilution factor for fixed fold pipetting
BM_target = 100.; # target biomass for regular pipetting
newborns_per_adult = int(np.floor(num_wells / top_tier)); # max newborns from each adult

#system params
W0 = 200; #initial waste
v = 1;
R0 = 1*v; #initial rsrc
N0 = 100; #initial num cells

g10 = 1; # max growth rate
g20 = 1;
c1 = 1e-4; # rsrc consumption rate
c2 = 1e-4;
f1 = 0.6; # fp
f2 = 0.6;
r2 = 5; # waste deg rate
K_R1 = 2; # rsrc affinity
K_R2 = 2;
K_W1 = 0.6*W0; # waste affinity
K_W2 = 0.1*W0;
a = 7; # power for species 1 growth eq
W2 = 0.1*W0; # waste toxicity for species 2
d1 = 1e-3; # death rate
d2 = 1e-3;

#function to integrate
#y: [R(t), W(t), ln(S1(t)), ln(S2(t))]
#low: S1 tracker
#high: S2 tracker
def dydt(t, y, low, high):
    # state variables
    R = y[0];
    W = y[1];
    S1 = low.biomass(); # this is a vector of variable length
    S2 = high.biomass(); # this is a vector of variable length
    
    g1 = g10 * R/(R+K_R1); 
    g2 = g20 * R/(R+K_R2);
    
    f1 = low.traits[:,0]; # f1 of each genotype stored in tracker
    f2 = high.traits[:,0];
    
    dS1 = (1-f1)*g1; # S1 birth rate per capita
    dS2 = (1-f2)*g2; # S2 birth rate per capita
    #use sum function because S1 and S2 are vectors
    dR = -np.sum(c1*g1*S1) - np.sum(c2*g2*S2); #dR/dt
    dW = -np.sum(f1*g1*W**a/(W**a+K_W1**a)*S1) \
         -np.sum(r2*f2*g2*W/(W+K_W2)*np.exp(-W/W2)*S2);  #dW/dt
    
    #dydt must be a list so dS1 and dS2 must be converted back to list form
    dy = [dR, dW] + np.append(dS1, dS2).tolist();
    return dy;

# initialize community structure
ic_L1 = np.array([1.]); #initial cell length
ic_L2 = np.array([1.]); #initial cell length

ic_f1_range = np.linspace(0.1,0.9,num=10);
phi1_range = np.linspace(0.1,0.9,num=10);

# make 100 communities
# scan over 10 initial phi1 values x 10 initial f1 values
newbDataAll = [];
for i in range(num_wells):
    #make S1 CellType
    low = commSelect.CellType.CellType(np.atleast_2d(ic_f1_range[int(i/10)]), ic_L1, \
                    int(phi1_range[i % 10] * N0), d1 * dt);
    low.traits_bound = np.array([1.]); # mutation bounds
    low.traits_lowerbound = np.array([0.]);
    
    #make S2 CellType
    #args for construction of CellType: traits, L, N, death rate
    high = commSelect.CellType.CellType(np.atleast_2d(ic_f1_range[int(i/10)]), ic_L2, \
                    N0 - int(phi1_range[i % 10] * N0), d2 * dt);
    high.traits_bound = np.array([1.]); # mutation bounds
    high.traits_lowerbound = np.array([0.]);
    
    newbData = [low, high];
    newbDataAll += [newbData];

ic_metabol = [R0, W0]; # metabolite types order: R, W

for c in range(num_cycles):
    #mature all communities      
    #args for mature function: 
    #newbDataAll (list of newborn communities)
    #ic_metabol (initial metabolite concs)
    #RPF_prime (diffeq method)
    #nsteps (num timesteps per maturation phase)
    #dt (timestep length)
    #mut_params (mutation parameters)
    #parallel (boolean, True if using parallel processing)
        
    #returns:
    #adultDataAll (list of adult communities)
    #data (time series data of biomass and metabolites for each community)
    #data structure: [community index, 'data type', time]
    #'data type' is total biomass of each cell type and conc of each metabolite
    #biomass comes first (S1, S2) then metabolites (R, W)
    #so for example data[:,1,-1] is the final S2 biomass of each community

    adultDataAll, data = commSelect.mature.mature(newbDataAll, ic_metabol, \
                            dydt, nsteps, dt, mut_params, parallel);
        
    
    # remaining waste of each community
    Wf = data[:,3,-1];
    P = W0 - Wf; # community function: how much waste was degraded
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
    #cs=False: use pipetting method
    newbDataAll = commSelect.reproduce.reproduce(adultDataAll, newbDataAll, P_sorted, \
                            num_wells, BM_target, newborns_per_adult, False);
        
        
P_array = np.reshape(P,(10,10));

plt.contour(ic_f1_range, phi1_range, P_array.T);
plt.xlabel("fp")
plt.ylabel("frac S1")