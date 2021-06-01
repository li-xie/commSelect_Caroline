#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 13:01:17 2020

@author: carolc24
"""


import numpy as np

def normcond(N,p):
    return np.logical_and((N > 9 * p / (1-p)), (N > 9 * (1-p) / p));

def fastbinorv(N,p):
    if (np.size(N) <= 0):
        #print("error: N must not be empty");
        return;
        
    if (np.ndim(N) > 1 and np.shape(N)[1] > 1):
        print("error: N must be a one dimensional array");
        return;
        
    thresh_b2p = 30; # threshold for binomial to poisson dist
    thresh_p2n = 800; # threshold for poisson to normal dist
    
    result = np.zeros(np.shape(N));
    
    if (np.isscalar(p)):
        p = np.ones(np.shape(N)) * p;
        
    if (N.min() > thresh_p2n and normcond(N,p).all()):
        result = np.round(np.random.normal(N * p, np.sqrt(N * p * (1 - p)))).astype(int);
    else:
        bino_ind = np.nonzero(np.logical_and(N > 0, N < thresh_b2p));
        pois_ind = np.nonzero(np.logical_and(N >= thresh_b2p, 
                np.logical_or(N < thresh_p2n, np.logical_not(normcond(N,p)))));
        norm_ind = np.nonzero(np.logical_and(N >= thresh_p2n, normcond(N,p)));
        
        if (np.size(bino_ind) > 0):
            result[bino_ind] = np.random.binomial(N[bino_ind],p[bino_ind]);
            
        result[pois_ind] = np.random.poisson(N[pois_ind] * p[pois_ind]);
        result[norm_ind] = np.round(np.random.normal(N[norm_ind] * p[norm_ind],
              np.sqrt(N[norm_ind] * p[norm_ind] * (1 - p[norm_ind])))).astype(int);
        
    result[result < 0] = result[result < 0] * 0;
    result[result > N] = result[result > N] * 0 + N[result > N];
    
    return result.astype(int);

def mutrndDunham(sp, sn, n):
    u = np.random.rand(n);
    mu = np.zeros(len(u));
    
    idx = np.nonzero(u <= sn / (sp + sn));
    mu[idx] = sn * np.log((sp + sn) * u[idx] / sn );
    idx = np.nonzero(u > sn / (sp + sn));
    mu[idx] = -sp * np.log(1 - (u[idx] - sn/(sp+sn))*(sp+sn)/sp);
    mu[mu<-1] = -1;
    return mu;

# contains information on all variants of this type
class CellType:
    
    #create new cell type
    def __init__(self, traits, L, N):
        if len(traits) == len(N) and len(N) == len(L):
            self.traits = traits;
            self.N = N;
            self.L = L;
            self.n_genos = len(traits);
        elif len(N) == 1 and len(L) == 1:
            self.traits = traits;
            self.N = np.repeat(N, len(traits));
            self.L = np.repeat(L, len(traits));
            self.n_genos = len(traits);
        else:
            print("Error: invalid input (check sizes of vectors)");
            return;
    
    #biomass vector by strain        
    def biomass(self):
        return np.transpose([self.N * self.L]);
    
    #weighted mean of each trait
    def meanTraits(self):
        bm = self.biomass();
        if bm > 0:
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
    
    def death(self, d):
        self.N -= fastbinorv(self.N, d);
    
    #mutate dividing cells    
    def mutateTraits(self, pot_mut_index, p_mut, frac_null, sp0, sn0):
        
        for i in range(self.n_genos):  # one trait at a time
            if sum(self.N[pot_mut_index[:,i]]) > 0:
                N_mut = fastbinorv(self.N[pot_mut_index[:,i]], p_mut);
                if sum(N_mut) > 0:
                    L_mut = self.L[pot_mut_index[:,i]];
                    traits_mut = self.traits[pot_mut_index[:,i],:];
                    
                    self.N[pot_mut_index[:,i]] -= N_mut; # remove mutants from orig categories
                    traits_mut = traits_mut[np.nonzero(N_mut),:][0];
                    L_mut = L_mut[np.nonzero(N_mut)];
                    N_mut = N_mut[np.nonzero(N_mut)];
                    
                    #null mutations
                    N_null = fastbinorv(N_mut, frac_null);
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
                                    
                    traits_split = np.empty((1, self.n_genos));
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
                        mult = mutrndDunham(sp0, sn0, np.shape(traits_split)[0]) + 1;
                        traits_split[:,i] *= mult;
                        
                        self.traits = np.append(self.traits, traits_split, axis=0);
                        self.L = np.append(self.L, L_split);
                        self.N = np.append(self.N, N_split);
                        self.n_genos += len(N_split);
                        pot_mut_index = np.append(pot_mut_index, traits_split > 0, axis=0);
        
        return;
        
    #constrain traits to upper bounds
    def boundTraits(self, traits_bound):
        for i in range(self.n_genos):
            self.traits[i,self.traits[i,:] > traits_bound] = \
                traits_bound[self.traits[i,:] > traits_bound];
      
    #remove strains that cannot grow (zero growth rate or rsrc affinity, etc)
    def clearInviable(self, inviable_index):
        viable = (self.traits[:,inviable_index] > 0).all(axis=1);
        n_inviable = len(self.L) - len(np.nonzero(viable)[0]);
        self.traits = self.traits[viable,:];
        self.L = self.L[viable];
        self.N = self.N[viable];
        self.n_genos -= n_inviable;
        
