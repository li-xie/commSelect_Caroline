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
    
    #iterate over each cell type
    for c in range(len(winnerData)):

        #separate adult into nD partitions
        newb_N_mat = [];
        for i in winnerData[c].N:
            newb_N_mat.append(np.random.multinomial(i, np.ones(nD) * 1./nD));
        
        if len(newb_N_mat) > 0:    
            newb_N_mat = np.array(newb_N_mat);
        else:
            newb_N_mat = np.zeros((1,nD));
   
        #seed newborns
        for j in range(len(target_inds)):

                
            if sum(newb_N_mat[:,j]) > 0:
                nb_traits_temp = winnerData[c].traits[np.nonzero(newb_N_mat[:,j]),:][0];
                nb_L_temp = winnerData[c].L[np.nonzero(newb_N_mat[:,j])];
                nb_N_temp = newb_N_mat[np.nonzero(newb_N_mat[:,j]),j][0];
            else:
         
                nb_traits_temp = winnerData[c].traits;
                nb_L_temp = winnerData[c].L;
                nb_N_temp = winnerData[c].N * 0;
            
            newbDataAll[target_inds[j]][c].traits = np.array(nb_traits_temp); #traits
            newbDataAll[target_inds[j]][c].L = np.array(nb_L_temp); #L
            newbDataAll[target_inds[j]][c].N = np.array(nb_N_temp); #N
            newbDataAll[target_inds[j]][c].n_genos = len(nb_L_temp);

    return newbDataAll;   