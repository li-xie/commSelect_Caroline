#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 15:21:04 2020

@author: carolc24
"""

import commSelect

def mature(newbDataAll, ic_metabol, ivpFunc, nsteps, dt, mut_params, parallel):

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
            print(i);
            adultDataAll[i] = commSelect.simulateOneWell.simulateOneWell( \
                                newbDataAll[i], ic_metabol, \
                                ivpFunc, nsteps, dt, mut_params);
    
    #for saving adult data to files
    data = [[0.]] * num_wells;
    
    for i in range(num_wells):
        data[i] = adultDataAll[i][-1]; # biomass+metabolite time series data
        adultDataAll[i] = adultDataAll[i][:-1][0]; #CellType objects
        
    return [adultDataAll, data];