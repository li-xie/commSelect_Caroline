#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 12:42:45 2020

@author: carolc24
"""

import numpy as np
import commSelect.CellType
import scipy.integrate

#simulate growth+mutation in one community
#input:
#newbData: newborn community, list of CellType objects
#initMetabolData: list, initial conditions for metabolites
#ivpFunc: function, net growth+production rate function to integrate
#nsteps: int, number of integration steps
#dt: float, size of integration step
#mut_params: list of size 4, mutation parameters
#output:
#adultData: mature community, list of CellType objects
#data: time series data of species biomass + metabolite concentration
def simulateOneWell(newbData, initMetabolData, ivpFunc, nsteps, dt, mut_params, rng_seed):  
    rng = np.random.default_rng(rng_seed)    
    #initialize metabolite time series data structure
    metabolData = np.append(np.transpose([initMetabolData]), \
                            np.zeros((len(initMetabolData),nsteps)), axis=1);
        
    n_metabols = len(initMetabolData);
    
    #initial biomass data for each cell type
    bioInitIter = map(commSelect.CellType.CellType.biomassSum, newbData);
    bioInit = [_ for _ in bioInitIter];
    
    #initialize species biomass time series data structure
    bioData = np.append(np.transpose([bioInit]), np.zeros((len(bioInit),nsteps)), axis=1);
    
    #for tracking number of genotypes per species
    n_genos_curr = [];
    for n in newbData:
        n_genos_curr += [n.n_genos];
     
    #mut parameters
    p_mut = mut_params[0];
    frac_null = mut_params[1];
    sp0 = mut_params[2];
    sn0 = mut_params[3];
    
    #simulate
    for tstep in range(1, nsteps + 1):
        
        #for ode solver
        now = np.append(metabolData[:,tstep - 1], np.zeros(sum(n_genos_curr)));
        tspan = [0., dt];
        ode_args = (newbData);
        soln = scipy.integrate.solve_ivp(ivpFunc, tspan, now, args=ode_args);
        
        #update rsrc and byproduct
        metabolData[:,tstep] = soln.y[:n_metabols,-1]
        
        #grow, divide, die, mutate, and mutation housekeeping
        count = n_metabols;
        for n in range(len(newbData)):
            #grow all genotypes in species n
            newbData[n].grow(np.exp(soln.y[count:(count+newbData[n].n_genos),-1]));
            #update tracker for next species
            count += newbData[n].n_genos;
            #divide cells of length >= 2, track daughter cells for mutation
            pot_mut_index = newbData[n].divide();
            #random cell death
            newbData[n].death(rng);
            #mutate random daughter cells
            newbData[n].mutateTraits(pot_mut_index, p_mut, frac_null, sp0, sn0, rng);
            #bound mutant cell genotypes
            newbData[n].boundTraits();
            #update number of genotypes
            n_genos_curr[n] = newbData[n].n_genos;
                                                         
        #update biomass
        bioIter = map(commSelect.CellType.CellType.biomassSum, newbData);
        bio = [_ for _ in bioIter];
        bioData[:,tstep] = bio;
    
    #compile data from end of run        
    data = np.append(bioData, metabolData, axis=0);
    return [newbData, data];
