#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 15:21:04 2020

@author: carolc24
"""

import commSelect
import numpy as np

#lets list of communities grow for given number of steps
#input:
#newbDataAll: list of communities (community is a list of CellTypes)
#ic_metabol: list of initial concentrations for each metabolite
#ivpFunc: function to integrate
#nsteps: int, number of integration steps
#dt: float, size of each integration step
#mut_params: list with size 4, mutation parameters
#parallel: bool, True if running in parallel
#verbose: print id of well when it's being simulated
#output:
#adultDataAll: list of adult communities (community is a list of CellTypes)
#data: time series data on biomass of each CellType + conc of each metabolite

def mature(newbDataAll, ic_metabol, ivpFunc, nsteps, dt, mut_params, parallel, verbose=False):

    num_wells = len(newbDataAll);
    if (parallel):
        from mpi4py.futures import MPIPoolExecutor;
        from itertools import repeat;
        
        with MPIPoolExecutor() as executor:
            dataIter = executor.map(commSelect.simulateOneWell.simulateOneWell, newbDataAll, \
                       repeat(ic_metabol), repeat(ivpFunc), repeat(nsteps), \
                       repeat(dt), repeat(mut_params));
        adultDataAll = [_ for _ in dataIter];
    else:
        adultDataAll = [0.] * num_wells;
        for i in range(num_wells):
            if (verbose): print(i);
            adultDataAll[i] = commSelect.simulateOneWell.simulateOneWell( \
                                newbDataAll[i], ic_metabol, \
                                ivpFunc, nsteps, dt, mut_params);
    
    #for saving adult data to files
    data = [[0.]] * num_wells;
    
    for i in range(num_wells):
        data[i] = adultDataAll[i][-1]; # biomass+metabolite time series data
        adultDataAll[i] = adultDataAll[i][:-1][0]; #CellType objects
        
    return [adultDataAll, np.array(data)];