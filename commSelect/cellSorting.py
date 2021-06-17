#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 14:53:58 2020

@author: carolc24
"""


import commSelect.CellType
import numpy as np

#sort cells from mature community into newborns (fix phi and biomass)
#input:
#newbDataAll: newborn community data structure
#winnerData: mature community to be reproduced
#target_inds: list of indices of newborns to be populated
#nD: dilution factor for mature community
#N0: target total number of cells
#output:
#newbDataAll: newborn community data structure
def cellSorting(newbDataAll, winnerData, target_inds, nD, N0):
    
    #calculate target fraction for each cell type
    bm = list(map(commSelect.CellType.CellType.biomassSum, winnerData));
    phi = bm / sum(bm);
    bm_targets = phi * N0;
    
    #iterate over each cell type
    for c in range(len(winnerData)):
        
        #biomass target
        bmt = bm_targets[c];
        
        newb_N_mat = [];
        for i in winnerData[c].N:
            newb_N_mat.append(np.random.multinomial(i, np.ones(nD) * 1./nD));
            
        if len(newb_N_mat) > 0:    
            newb_N_mat = np.array(newb_N_mat);
        else:
            newb_N_mat = np.zeros((1,nD));    
          
        #to hold extra cells
        donor = np.zeros((len(winnerData[c].N)));
        
        #tax overfilled wells and give their cells to donor pool
        for i in range(nD):
            
            while (sum(newb_N_mat[:,i] * winnerData[c].L) > bmt):
                ne_inds = newb_N_mat[:,i] > 0;
                temp1 = newb_N_mat[ne_inds,i];
                temp2 = donor[ne_inds];
                drawOneCell(temp1, temp2);
                newb_N_mat[ne_inds,i] = temp1;
                donor[ne_inds] = temp2;
                
        #fill underfilled wells with cells from donor pool
        for i in range(nD):
            while sum(newb_N_mat[:,i] * winnerData[c].L) < bmt:
                ne_inds = donor > 0;
                temp1 = newb_N_mat[ne_inds,i];
                temp2 = donor[ne_inds];
                indx = drawOneCell(temp2, temp1);
                newb_N_mat[ne_inds,i] = temp1;
                donor[ne_inds] = temp2;
                
            #remove the last cell added so biomass is just under target
            if sum(winnerData[c].L * newb_N_mat[:,i]) > bmt:
                temp1[indx] -= 1;
                temp2[indx] += 1;
                newb_N_mat[ne_inds,i] = temp1;
                donor[ne_inds] = temp2;
                
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
            
    
#take one cell randomly from i vector and give it to o vector    
#i: donor vector of cell counts
#o: recipient vector of cell counts
#indx: index of cell that moved
def drawOneCell(i, o):
    csN = np.cumsum(i);
    indx = np.nonzero(csN >= np.random.randint(csN[-1]))[0][0];
    i[indx] -= 1;
    o[indx] += 1;
    return indx;