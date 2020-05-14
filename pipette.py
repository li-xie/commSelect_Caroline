#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 14:54:59 2020

@author: carolc24
"""


import commSelect.CellType
import numpy as np

def pipette(newbDataAll, winnerData, target_inds, nD):
    
    if (len(target_inds) > nD):
        print("Error: number of targets must be less than dilution factor!")
        return;
    
    empty_wells = [];
    for c in range(len(winnerData)):
    #distribute manu cells
        newb_N_mat = [];
        for i in winnerData[c].N:
            newb_N_mat.append(np.random.multinomial(i, np.ones(nD) * 1./nD));
        
        if len(newb_N_mat) > 0:    
            newb_N_mat = np.array(newb_N_mat);
        else:
            newb_N_mat = np.zeros((0,nD));
        
        extras = range(len(target_inds), nD);     
        #seed newborns
        for i in range(len(target_inds)):
            j = i;
            while sum(newb_N_mat[:,j]) <= 0 and len(extras) > 0:
                j = extras[0];
                extras = extras[1:];
                
            if sum(newb_N_mat[:,j]) > 0:
                nb_traits_temp = winnerData[c].traits[np.nonzero(newb_N_mat[:,j]),:];
                nb_L_temp = winnerData[c].L[np.nonzero(newb_N_mat[:,j])];
                nb_N_temp = newb_N_mat[np.nonzero(newb_N_mat[:,j]),j][0];
            else:
                empty_wells += [target_inds[i]];            
                        
            newbDataAll[target_inds[i]][c].traits = np.array(nb_traits_temp)[0]; #traits
            newbDataAll[target_inds[i]][c].L = np.array(nb_L_temp); #L
            newbDataAll[target_inds[i]][c].N = np.array(nb_N_temp); #N
            newbDataAll[target_inds[i]][c].n_genos = len(nb_L_temp);

    return newbDataAll, empty_wells;   