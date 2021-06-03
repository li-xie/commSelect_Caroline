#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 13:10:18 2020

@author: carolc24
"""

import commSelect.CellType
import numpy as np

#save data on newborns to file
def saveNewborns(newbDataAll, fname):
    data_final = [0.] * len(newbDataAll);
    for i in range(len(newbDataAll)):
        traits_data = np.concatenate(list(map( \
            commSelect.CellType.CellType.meanTraits, newbDataAll[i])));
        bio_data = np.array(list(map( \
            commSelect.CellType.CellType.biomassSum, newbDataAll[i])));
        data_final[i] = np.append(traits_data, bio_data);
        
    np.savetxt(fname, np.array(data_final));
    
#save data on adults to file
def saveAdults(adultDataAll, add_data, fname):
    data_final = [0.] * len(adultDataAll);
    for i in range(len(adultDataAll)):
        #final average traits
        traits_data = np.concatenate(list(map( \
            commSelect.CellType.CellType.meanTraits, adultDataAll[i])))
        #final average traits, biomass, and metabolite conc
        data_final[i] = np.append(traits_data, add_data[i][:,-1]);
        
    np.savetxt(fname, np.array(data_final));