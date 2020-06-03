#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 12:19:37 2020

@author: carolc24
"""

import numpy as np
import commSelect.CellType

def reproduce(adultDataAll, newbDataAll, P_sorted, num_wells, BM_target, newborns_per_adult):
    
    finished_pipetting = False;
    num_wells_filled = 0;
    win_inds = [];
    
    for i in range(num_wells):
        windex = int(P_sorted[-(i+1)]); #index of adult being reproduced
        win_inds += [windex];
        winnerData = adultDataAll[windex];
        
        #dilution factor
        nD = int(np.floor(sum(list(map(commSelect.CellType.CellType.biomassSum, winnerData)))/BM_target));

        #if biomass of winner is small then reproduce all of it in 1 newborn
        if nD < 1:
            print("small adult kept in slot %d" % num_wells_filled)
            newbDataAll[num_wells_filled] = winnerData;
            num_wells_filled += 1;
        else:
            #separate adult into nD newborns
            target_inds = np.asarray(range(num_wells_filled, num_wells_filled + nD));
            
            #or the max number of newborns per adult if nD is larger
            if (len(target_inds) > newborns_per_adult):
                target_inds = np.asarray(range(num_wells_filled, num_wells_filled + newborns_per_adult));
            
            #or however many newborns still need to be filled
            if (target_inds[-1] >= num_wells-1):
                target_inds = np.asarray(range(num_wells_filled, num_wells));
                finished_pipetting = True;
                
            num_wells_filled += len(target_inds);
    
            print("filling: " + str(target_inds));
            
            newbDataAll = commSelect.pipette.pipette(newbDataAll, winnerData, target_inds, nD);
            
        if (finished_pipetting):
            break;
    return newbDataAll;
