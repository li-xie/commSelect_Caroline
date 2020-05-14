#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 12:19:37 2020

@author: carolc24
"""

import numpy as np
import commSelect.CellType

def reproduce(adultDataAll, newbDataAll, P, num_wells, BM_target, newborns_per_adult):
    
    finished_pipetting = False;
    num_wells_filled = 0;
    win_inds = [];
    empty_wells = [];
    P_sorted = np.array([np.sort(P), np.argsort(P)]);
    
    for i in range(num_wells):
        windex = int(P_sorted[1,-(i+1)]);
        win_inds += [windex];
        winnerData = adultDataAll[windex];
        print("winner P: %.2f" % P_sorted[0,-(i+1)]);
        
        nD = int(np.floor(sum(list(map(commSelect.CellType.CellType.biomassSum, winnerData)))/BM_target));
                
        target_inds = np.asarray(range(num_wells_filled, num_wells_filled + nD));
        if (len(target_inds) > newborns_per_adult):
            target_inds = np.asarray(range(num_wells_filled, num_wells_filled + newborns_per_adult));
        if (target_inds[-1] >= num_wells-1):
            target_inds = np.asarray(range(num_wells_filled, num_wells));
            finished_pipetting = True;
        num_wells_filled += max(len(target_inds) - len(empty_wells), 0);
        if (len(empty_wells) > 0 and len(empty_wells) <= len(target_inds)):
            target_inds[-len(empty_wells):] = empty_wells;
        elif (len(empty_wells) > len(target_inds)):
            target_inds = empty_wells;

        print("filling: " + str(target_inds));
            
        newbDataAll, empty_wells = commSelect.pipette.pipette(newbDataAll, winnerData, target_inds, nD);
                
        if (finished_pipetting):
            return newbDataAll;
