#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 16:53:11 2020

@author: arolle
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 16:58:20 2020

@author: arolle
"""

import numpy as np
import scipy.sparse
from BiGradedMatrix import *
import time

def grid(R):
    labels = list(set(R.labels))
    grid_with_repetition = [(labels[i][0], labels[j][1]) for i in range(len(labels)) for j in range(len(labels))]
    grid = list(set(grid_with_repetition))
    grid.sort()
    return(grid)

def BiRedSub(R, z, pivs):
        
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
    
def BiRedSub_MinGens(R, z, pivs):
        
    cols_not_reduced_to_zero = []
    
    Indices = [j for j in range(R.num_cols()) if R.labels[j][0] <= z[0] and R.labels[j][1] == z[1]]
    
    for j in Indices:
        l = R.get_piv(j)
        while l != -1 and pivs[l] != -1 and pivs[l] < j:
            k = pivs[l]
            R.add_column(k,j)
            l = R.get_piv(j)
        if l != -1 and R.labels[j] == z:
            cols_not_reduced_to_zero.append(j)
        if l != -1 and (pivs[l] == -1 or pivs[l] > j):
            pivs[l] = j
    
    return(R, pivs, cols_not_reduced_to_zero)
    
def BiRedSubSlave(R, V, z, pivs):
        
    cols_reduced_to_zero = []
    
    Indices = [j for j in range(R.num_cols()) if R.labels[j][0] <= z[0] and R.labels[j][1] == z[1]]
    
    for j in Indices:
        l = R.get_piv(j)
        while l != -1 and pivs[l] != -1 and pivs[l] < j:
            k = pivs[l]
            R.add_column(k,j)
            V.add_column(k,j)
            if R.get_piv(j) == -1:
                cols_reduced_to_zero.append(j)
            l = R.get_piv(j)
        if l != -1 and (pivs[l] == -1 or pivs[l] > j):
            pivs[l] = j
    
    return(R, V, pivs, cols_reduced_to_zero)
    
def KerBetti(R):
    
    betti_0 = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    gr = grid(R)
    
    for z in gr:
        R, pivs, cols_reduced_to_zero = BiRedSub(R, z, pivs)
        if cols_reduced_to_zero > 0:
            betti_0.append((z, cols_reduced_to_zero))
    
    return(R, betti_0)
    
# KerBasis uses scipy.sparse to turn the output into a bigraded matrix
    
def KerBasis(R):
    
    V = R.identity_matrix()
    B_ker = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    gr = grid(R)
    
    for z in gr:
        R, V, pivs, cols_reduced_to_zero = BiRedSubSlave(R, V, z, pivs)
        for j in cols_reduced_to_zero:
            B_ker.append((V.get_col(j), z))
     
    n = R.num_cols()
    m = len(B_ker)
    labels = [B_ker[j][1] for j in range(m)]
    matrix = scipy.sparse.csc_matrix((n,m), dtype='int')
    for j in range(m):
        matrix[:, j] = B_ker[j][0].matrix
        
    B = BiGradedMatrix(labels, matrix)
    
    return(B)
    
def MinGens(R):
    
    S = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    gr = grid(R)
    
    for z in gr:
        R, pivs, cols_not_reduced_to_zero = BiRedSub_MinGens(R, z, pivs)
        for j in cols_not_reduced_to_zero:
            S.append((R.get_col(j), z))
    
    return(S)
    
def pivot_array_of_B_ker(A):
    
    pivots = np.full(shape=A.num_rows(), fill_value=-1, dtype='int', order='C')
    
    for j in range(A.num_cols()):
        pivots[A.get_piv(j)] = j
        
    return(pivots)
    
def Solve(A, pivots, B):
    
    X = A.zero_vector()
    
    for i in reversed(range(A.num_rows())):
        if B.entry_is_nonzero(i, 0):
            if pivots[i] == -1:
                print('I could not find the pivot!')
            else:
                j = pivots[i]
                X.insert_value(j, 0, 1)
                B.matrix = B.matrix - A.get_col(j).matrix
                B.reduce()
                
    return(X)
    
def MinimizePres(P, row_labels):
    
    n = P.num_cols()
    cols_to_remove = []
    rows_to_remove = []
    
    for j in range(n):
        p = P.get_piv(j)
        if P.labels[j] == row_labels[p]:
            cols_to_remove.append(j)
            rows_to_remove.append(p)
            for k in range(j, n):
                P.add_multiple_of_column(P.matrix[p, k], j, k)
                
    cols = [i for i in range(P.num_cols())]
    rows = [i for i in range(P.num_rows())]
    for i in cols_to_remove:
        cols.remove(i)
    for i in rows_to_remove:
        rows.remove(i)
        
    P.labels = [P.labels[i] for i in cols]
    row_labels = [row_labels[i] for i in rows]
    
    P.keep_columns(cols)
    P.keep_rows(rows)
    
    return(P, row_labels)
        
# (delta2, delta1) is an FI-Rep; so, delta2, delta1 are bigraded matrices,
# where delta2 represents a map F^2 -> F^1, 
# and delta1 represents a map F^1 -> F^0
    
def MinimalPres(delta2, delta1):
    
    start = time.time()
    
    S = MinGens(delta2)
    
    MinGens_end = time.time()
    MinGens_time = MinGens_end - start
    
    print("Running MinGens took " + str(MinGens_time))
    
    B = KerBasis(delta1)
    
    KerBasis_end = time.time()
    KerBasis_time = KerBasis_end - MinGens_end
    
    print("Running KerBasis took " + str(KerBasis_time))
    
    
    pivots = pivot_array_of_B_ker(B)
    n = B.num_cols()
    m = len(S)
    labels = [S[j][1] for j in range(m)]
    matrix = scipy.sparse.csc_matrix((n,m), dtype='int')
    P = BiGradedMatrix(labels, matrix)
    
    for j in range(m):
        P.matrix[:, j] = Solve(B, pivots, S[j][0]).matrix
        
    P_end = time.time()
    P_time = P_end - KerBasis_end
    
    print("Computing P took " + str(P_time))
    
    row_labels = B.labels
    P, row_labels = MinimizePres(P, row_labels)
    
    MinimizePres_end = time.time()
    MinimizePres_time = MinimizePres_end - P_end
    Total_time = MinimizePres_end - start
    
    print("Running MinimizePres took " + str(MinimizePres_time))
    print("Total running time was " + str(Total_time))

    return(P, row_labels)