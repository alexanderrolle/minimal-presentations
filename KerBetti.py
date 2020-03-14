#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 16:58:20 2020

@author: arolle
"""

import numpy as np

def grid(R):
    labels = list(set(R.labels))
    grid_with_repetition = [(labels[i][0], labels[j][1]) for i in range(len(labels)) for j in range(len(labels))]
    grid = list(set(grid_with_repetition))
    grid.sort()
    return(grid)

def BiRedSub(R, z, pivs) :
        
    cols_reduced_to_zero = 0
    
    Indices = [j for j in range(R.num_cols()) if R.labels[j][0] <= z[0] and R.labels[j][1] == z[1]]
    
    for j in Indices:
        l = R.get_piv(j)
        while l != -1 and pivs[l] != -1 and pivs[l] < j:
            k = pivs[l]
            R.add_column(k,j)
            if R.get_piv(j) == -1:
                cols_reduced_to_zero = cols_reduced_to_zero + 1
            l = R.get_piv(j)
        if l != -1 and (pivs[l] == -1 or pivs[l] > j):
            pivs[l] = j
    
    return(R, pivs, cols_reduced_to_zero)
    
def KerBetti(R) :
    
    betti_0 = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    gr = grid(R)
    
    for z in gr:
        R, pivs, cols_reduced_to_zero = BiRedSub(R, z, pivs)
        if cols_reduced_to_zero > 0:
            betti_0.append((z, cols_reduced_to_zero))
    
    return(R, betti_0)