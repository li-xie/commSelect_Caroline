#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 12:19:37 2020

@author: carolc24
"""

import numpy as np
import commSelect.CellType

#create newborn communites from mature communities
#input:
#adultDataAll: list of mature communities (community is list of CellTypes)
#newbDataAll: initialized newborn community list
#P_sorted: list, indices of communities sorted by function from lowest to highest
#num_wells: int, number of communities
#BM_target: float, target biomass for newborns
#newborns_per_adult: int, max number of newborns to be made from 1 adult
#cs: bool, True if cell sorting, False if pipetting
#output:
#newbDataAll: populated newborn community list
def reproduce(adultDataAll, newbDataAll, P_sorted, num_wells, BM_target, newborns_per_adult, cs):
    
    finished_pipetting = False;
    num_wells_filled = 0;
    win_inds = [];
    
    #reproduce each adult from highest to lowest function until newborns full
    for i in range(num_wells):
        windex = int(P_sorted[-(i+1)]); #index of adult being reproduced
        win_inds += [windex];
        winnerData = adultDataAll[windex];
        
        #dilution factor (max number of newborns that can be populated with this adult)
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
            
            #or however many newborns still need to be filled if # available newborns < nD
            if (target_inds[-1] >= num_wells-1):
                target_inds = np.asarray(range(num_wells_filled, num_wells));
                finished_pipetting = True;
                
            num_wells_filled += len(target_inds);
    
            print("filling: " + str(target_inds));
            
            if cs: # cell sorting method
                newbDataAll = commSelect.cellSorting.cellSorting(newbDataAll, winnerData, target_inds, nD, BM_target);
            else: # pipetting method
                newbDataAll = commSelect.pipette.pipette(newbDataAll, winnerData, target_inds, nD);
            
        if (finished_pipetting):
            break;
    return newbDataAll;
