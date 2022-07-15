#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:33:56 2020

@author: carolc24
"""

import numpy as np

def normcond(N,p):
    p[p == 0] = 1e-9;  # avoid division errors
    return np.logical_and((N > 9 * p / (1-p)), (N > 9 * (1-p) / p));

def fastbinorv(N,p,rng):
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
        result = np.round(rng.normal(N * p, np.sqrt(N * p * (1 - p)))).astype(int);
    else:
        bino_ind = np.nonzero(np.logical_and(N > 0, N < thresh_b2p));
        pois_ind = np.nonzero(np.logical_and(N >= thresh_b2p, 
                np.logical_or(N < thresh_p2n, np.logical_not(normcond(N,p)))));
        norm_ind = np.nonzero(np.logical_and(N >= thresh_p2n, normcond(N,p)));
        
        if (np.size(bino_ind) > 0):
            result[bino_ind] = rng.binomial(N[bino_ind],p[bino_ind]);
            
        result[pois_ind] = rng.poisson(N[pois_ind] * p[pois_ind]);
        result[norm_ind] = np.round(rng.normal(N[norm_ind] * p[norm_ind],
              np.sqrt(N[norm_ind] * p[norm_ind] * (1 - p[norm_ind])))).astype(int);
        
    result[result < 0] = result[result < 0] * 0;
    result[result > N] = result[result > N] * 0 + N[result > N];
    
    return result.astype(int);

def mutrndDunham(sp, sn, n, rng):
    u = rng.random(n);
    mu = np.zeros(len(u));
    
    idx = np.nonzero(u <= sn / (sp + sn));
    mu[idx] = sn * np.log((sp + sn) * u[idx] / sn );
    idx = np.nonzero(u > sn / (sp + sn));
    mu[idx] = -sp * np.log(1 - (u[idx] - sn/(sp+sn))*(sp+sn)/sp);
    mu[mu<-1] = -1;
    return mu;