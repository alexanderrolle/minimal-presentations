#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 19:41:37 2020

@author: arolle
"""

import numpy as np
import scipy.sparse
from BiGradedMatrix_v2 import *
import time

def grid_v2(R):
    labels = list(set(R.labels))
    grid_with_repetition = [(labels[i][0], labels[j][1]) for i in range(len(labels)) for j in range(len(labels))]
    grid = list(set(grid_with_repetition))
    grid.sort()
    return(grid)

def BiRedSub_v2(R, z, pivs):
        
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
    
def BiRedSub_MinGens_v2(R, z, pivs):
        
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
    
def BiRedSubSlave_v2(R, V, z, pivs):
        
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
    
def KerBetti_v2(R):
    
    betti_0 = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    gr = grid_v2(R)
    
    for z in gr:
        R, pivs, cols_reduced_to_zero = BiRedSub_v2(R, z, pivs)
        if cols_reduced_to_zero > 0:
            betti_0.append((z, cols_reduced_to_zero))
    
    return(R, betti_0)
    
# uses scipy.sparse to turn the output into a bigraded matrix
def KerBasis_v2(R):
    
    V = R.identity_matrix()
    B_ker = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    gr = grid_v2(R)
    
    for z in gr:
        R, V, pivs, cols_reduced_to_zero = BiRedSubSlave_v2(R, V, z, pivs)
        for j in cols_reduced_to_zero:
            B_ker.append((V.get_col(j), z))
     
    n = R.num_cols()
    m = len(B_ker)
    labels = [B_ker[j][1] for j in range(m)]
    matrix = scipy.sparse.lil_matrix((m,n), dtype='int')
    for j in range(m):
        indices = B_ker[j][0].matrix.rows[0]
        matrix[j, indices] = [1 for i in range(len(indices))]
        
    B = BiGradedMatrix_lil(labels, matrix)
    
    return(B)
    
def MinGens_v2(R):
    
    S = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    gr = grid_v2(R)
    
    for z in gr:
        R, pivs, cols_not_reduced_to_zero = BiRedSub_MinGens_v2(R, z, pivs)
        for j in cols_not_reduced_to_zero:
            S.append((R.get_col(j), z))
    
    return(S)
    
def pivot_array_of_B_ker(A):
    
    pivots = np.full(shape=A.num_rows(), fill_value=-1, dtype='int', order='C')
    
    for j in range(A.num_cols()):
        pivots[A.get_piv(j)] = j
        
    return(pivots)
    
def Solve_v2(A, pivots, B):
    
    X = A.zero_vector()
    
    for i in reversed(range(A.num_rows())):
        if B.entry_is_nonzero(i, 0):
            if pivots[i] == -1:
                print('I could not find the pivot!')
            else:
                j = pivots[i]
                X.insert_value(j, 0, 1)
                col = A.matrix.rows[j]
                B.add_external_column(col, 0)
                
    return(X)
    
def MinimizePres_v2(P, row_labels):
    
    n = P.num_cols()
    cols_to_remove = []
    rows_to_remove = []
    
    for j in range(n):
        p = P.get_piv(j)
        if p == -1:
            return("P contains a zero column!")
        else:
            if P.labels[j] == row_labels[p]:
                cols_to_remove.append(j)
                rows_to_remove.append(p)
                for k in range(j+1, n):
                    if P.entry_is_nonzero(p, k):
                        P.add_column(j, k)
                
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
    
#def MinimizePres(P, row_labels):
#    
#    n = P.num_cols()
#    cols_to_remove = []
#    rows_to_remove = []
#    
#    for j in range(n):
#        p = P.get_piv(j)
#        if P.labels[j] == row_labels[p]:
#            cols_to_remove.append(j)
#            rows_to_remove.append(p)
#            for k in range(j, n):
#                P.add_multiple_of_column(P.matrix[p, k], j, k)
#                
#    cols = [i for i in range(P.num_cols())]
#    rows = [i for i in range(P.num_rows())]
#    for i in cols_to_remove:
#        cols.remove(i)
#    for i in rows_to_remove:
#        rows.remove(i)
#        
#    P.labels = [P.labels[i] for i in cols]
#    row_labels = [row_labels[i] for i in rows]
#    
#    P.keep_columns(cols)
#    P.keep_rows(rows)
#    
#    return(P, row_labels)
#    
#def MinimizePres_test(P, row_labels):
#    
#    start = time.time()
#    
#    n = P.num_cols()
#    cols_to_remove = []
#    rows_to_remove = []
#    
#    for j in range(n):
#        p = P.get_piv(j)
#        if P.labels[j] == row_labels[p]:
#            cols_to_remove.append(j)
#            rows_to_remove.append(p)
#            for k in range(j, n):
#                P.add_multiple_of_column(P.matrix[p, k], j, k)
#                
#    RowOps_end = time.time()
#    RowOps_time = RowOps_end - start
#    print("Row operations took " + str(RowOps_time))
#                
#    cols = [i for i in range(P.num_cols())]
#    rows = [i for i in range(P.num_rows())]
#    for i in cols_to_remove:
#        cols.remove(i)
#    for i in rows_to_remove:
#        rows.remove(i)
#        
#    P.labels = [P.labels[i] for i in cols]
#    row_labels = [row_labels[i] for i in rows]
#    
#    P.keep_columns(cols)
#    P.keep_rows(rows)
#    
#    end = time.time()
#    Total_time = end - start
#    print("Total runtime was " + str(Total_time))
#    
#    return(P, row_labels)
#    
#def MinimizePres_test_v2(P, row_labels):
#    
#    start = time.time()
#    
#    n = P.num_cols()
#    cols_to_remove = []
#    rows_to_remove = []
#    
#    for j in range(n):
#        p = P.get_piv(j)
#        if p == -1:
#            return("P contains a zero column!")
#        else:
#            if P.labels[j] == row_labels[p]:
#                cols_to_remove.append(j)
#                rows_to_remove.append(p)
#                for k in range(j+1, n):
#                    if P.entry_is_nonzero(p, k):
#                        P.add_column(j, k)
#                
#    RowOps_end = time.time()
#    RowOps_time = RowOps_end - start
#    print("Row operations took " + str(RowOps_time))
#                
#    cols = [i for i in range(P.num_cols())]
#    rows = [i for i in range(P.num_rows())]
#    for i in cols_to_remove:
#        cols.remove(i)
#    for i in rows_to_remove:
#        rows.remove(i)
#        
#    P.labels = [P.labels[i] for i in cols]
#    row_labels = [row_labels[i] for i in rows]
#    
#    P.keep_columns(cols)
#    P.keep_rows(rows)
#    
#    end = time.time()
#    Total_time = end - start
#    print("Total runtime was " + str(Total_time))
#    
#    return(P, row_labels)
        
# (delta2, delta1) is an FI-Rep; so, delta2, delta1 are bigraded matrices,
# where delta2 represents a map F^2 -> F^1, 
# and delta1 represents a map F^1 -> F^0
    
def MinimalPres_v2(delta2, delta1):
    
    start = time.time()
    
    S = MinGens_v2(delta2)
    
    MinGens_end = time.time()
    MinGens_time = MinGens_end - start
    
    print("Running MinGens took " + str(MinGens_time))
    
    B = KerBasis_v2(delta1)
    
    KerBasis_end = time.time()
    KerBasis_time = KerBasis_end - MinGens_end
    
    print("Running KerBasis took " + str(KerBasis_time))
    
    pivots = pivot_array_of_B_ker(B)
    n = B.num_cols()
    m = len(S)
    labels = [S[j][1] for j in range(m)]
    matrix = scipy.sparse.lil_matrix((m,n), dtype='int')
    P = BiGradedMatrix_lil(labels, matrix)
    
    for j in range(m):
        P.matrix[j, :] = Solve_v2(B, pivots, S[j][0]).matrix
        
    P_end = time.time()
    P_time = P_end - KerBasis_end
    
    print("Computing P took " + str(P_time))
    
    row_labels = B.labels
    P, row_labels = MinimizePres_v2(P, row_labels)
    
    MinimizePres_end = time.time()
    MinimizePres_time = MinimizePres_end - P_end
    Total_time = MinimizePres_end - start
    
    print("Running MinimizePres took " + str(MinimizePres_time))
    print("Total running time was " + str(Total_time))

    return(P, row_labels)