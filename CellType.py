#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:35:07 2020

@author: carolc24
"""

import commSelect.randFcns

import numpy as np
# contains information on all variants of this type
class CellType:
    
    #create new cell type
    def __init__(self, traits, L, N, d):
        if len(traits) == len(N) and len(N) == len(L):
            self.traits = traits.copy();
            self.N = N.copy();
            self.L = L.copy();
            self.n_genos = len(traits);
            self.death_rate = d;
        elif len(N) == 1 and len(L) == 1:
            self.traits = traits.copy();
            self.N = np.repeat(N.copy(), len(traits));
            self.L = np.repeat(L.copy(), len(traits));
            self.n_genos = len(traits);
            self.death_rate = d;
        else:
            print("Error: invalid input (check sizes of vectors)");
            return;
        self.traits_bound = self.traits.copy();
        self.inviable_index = [];
        
    def copy(self):
        other = CellType(self.traits, self.L, self.N, self.death_rate);
        other.traits_bound = self.traits_bound.copy();
        other.inviable_index = self.inviable_index;
        return other;
    
    #biomass vector by strain        
    def biomass(self):
        return np.transpose([self.N * self.L]);
    
    #biomass as scalar
    def biomassSum(self):
        return np.sum(self.biomass());
    
    #weighted mean of each trait
    def meanTraits(self):
        bm = self.biomass();
        if sum(bm) > 0:
            return np.sum(self.traits * bm, axis=0) / sum(bm);
        elif np.shape(self.traits)[0] > 0:
            return np.sum(self.traits * bm, axis=0);
        elif np.ndim(self.traits) == 2:
            return np.zeros(np.shape(self.traits)[1]);
        else:
            return 0;
        
    def grow(self, growth_rate):
        self.L *= growth_rate;
     
    def divide(self):
        lgt2 = (self.L > 2);
        self.N[lgt2] *= 2;
        self.L[lgt2] /= 2;
        
        return np.logical_and(np.transpose([lgt2]), self.traits > 0);
    
    def death(self):
        self.N -= commSelect.randFcns.fastbinorv(self.N, self.death_rate);
    
    #mutate dividing cells    
    def mutateTraits(self, pot_mut_index, p_mut, frac_null, sp0, sn0):
        
        nt = np.shape(self.traits)[1]; # number of traits
        for i in range(nt):  # one trait at a time
            if sum(self.N[pot_mut_index[:,i]]) > 0:
                N_mut = commSelect.randFcns.fastbinorv(self.N[pot_mut_index[:,i]], p_mut);
                if sum(N_mut) > 0:
                    L_mut = self.L[pot_mut_index[:,i]];
                    traits_mut = self.traits[pot_mut_index[:,i],:];
                    
                    self.N[pot_mut_index[:,i]] -= N_mut; # remove mutants from orig categories
                    traits_mut = traits_mut[np.nonzero(N_mut),:][0];
                    L_mut = L_mut[np.nonzero(N_mut)];
                    N_mut = N_mut[np.nonzero(N_mut)];
                    
                    #null mutations
                    N_null = commSelect.randFcns.fastbinorv(N_mut, frac_null);
                    N_am = N_mut - N_null;
                    traits_null = traits_mut[np.nonzero(N_null),:][0];
                    L_null = L_mut[np.nonzero(N_null)];
                    N_null = N_null[np.nonzero(N_null)];
                    
                    #add categories to struct
                    if np.size(traits_null)  > 0:
                        traits_null[:,i] = 0;  
                        self.traits = np.append(self.traits, traits_null, axis=0);
                        self.N = np.append(self.N, N_null);
                        self.L = np.append(self.L, L_null);
                        self.n_genos += len(N_null);
                        pot_mut_index = np.append(pot_mut_index, traits_null > 0, axis=0);
                    
                    #active mutations
                    traits_am = traits_mut[np.nonzero(N_am),:][0];
                    L_am = L_mut[np.nonzero(N_am)];
                    N_am = N_am[np.nonzero(N_am)];
                                    
                    traits_split = np.empty((1, nt));
                    L_split = np.array([]);
                    
                    for j in range(len(N_am)):
                        traits_split = np.append(traits_split, \
                            np.repeat([traits_am[j,:]], N_am[j], axis=0), axis=0);
                        L_split = np.append(L_split, np.repeat(L_am[j], N_am[j]));
                    
                    traits_split = traits_split[1:,:];
                    N_split = np.ones(np.shape(L_split)).astype(int);
                    
                    #generate active muts, add categories to struct
                    if np.size(traits_split) > 0:
                        if (traits_split[:,i] == 0).any():
                            print("Error: something is wrong");
                            print(traits_split);
                            return;
                        mult = commSelect.randFcns.mutrndDunham(sp0, sn0, np.shape(traits_split)[0]) + 1;
                        traits_split[:,i] *= mult;
                        
                        self.traits = np.append(self.traits, traits_split, axis=0);
                        self.L = np.append(self.L, L_split);
                        self.N = np.append(self.N, N_split);
                        self.n_genos += len(N_split);
                        pot_mut_index = np.append(pot_mut_index, traits_split > 0, axis=0);
        
        return;
        
    #constrain traits to upper bounds
    def boundTraits(self):
        for i in range(self.n_genos):
            self.traits[i,self.traits[i,:] > self.traits_bound] = \
                self.traits_bound[self.traits[i,:] > self.traits_bound];
      
    #remove strains that cannot grow (zero growth rate or rsrc affinity, etc)
    def clearInviable(self):
        viable = (self.traits[:,self.inviable_index] > 0).all(axis=1);
        n_inviable = len(self.L) - len(np.nonzero(viable)[0]);
        self.traits = self.traits[viable,:];
        self.L = self.L[viable];
        self.N = self.N[viable];
        self.n_genos -= n_inviable;
        
