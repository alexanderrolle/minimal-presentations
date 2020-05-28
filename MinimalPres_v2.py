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
import heapq

def grid_v2(R):
    labels = list(set(R.labels))
    grid_with_repetition = [(labels[i][0], labels[j][1]) for i in range(len(labels)) for j in range(len(labels))]
    grid = list(set(grid_with_repetition))
    grid.sort()
    return(grid)
    
def grid_v3(R):
    
    labels = list(set(R.labels))
    
    x_coords_w_rep = [label[0] for label in labels]
    x_coords = list(set(x_coords_w_rep))
    y_coords_w_rep = [label[1] for label in labels]
    y_coords = list(set(y_coords_w_rep))
    
    grid = [(x, y) for x in x_coords for y in y_coords]
    
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
    
#def BiRedSub_MinGens_v2(R, z, pivs):
#        
#    cols_not_reduced_to_zero = []
#    
#    Indices = [j for j in range(R.num_cols()) if R.labels[j][0] <= z[0] and R.labels[j][1] == z[1]]
#    
#    for j in Indices:
#        l = R.get_piv(j)
#        while l != -1 and pivs[l] != -1 and pivs[l] < j:
#            k = pivs[l]
#            R.add_column(k,j)
#            l = R.get_piv(j)
#        if l != -1 and R.labels[j] == z:
#            cols_not_reduced_to_zero.append(j)
#        if l != -1 and (pivs[l] == -1 or pivs[l] > j):
#            pivs[l] = j
#    
#    return(R, pivs, cols_not_reduced_to_zero)
    
def BiRedSub_MinGens_v2(R, z, pivs, Indices):
        
    cols_not_reduced_to_zero = []
    
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
    
def BiRedSub_MinGens_v3(R, z, pivs, current_columns, X):
        
    cols_not_reduced_to_zero = []
    
    for j in current_columns:
        heapq.heappush(X[z[1]], j)
        
    while len(X[z[1]]) != 0:
        j = heapq.heappop(X[z[1]])
        l = R.get_piv(j)
        while l != -1 and pivs[l] != -1 and pivs[l] < j:
            k = pivs[l]
            R.add_column(k, j)
            l = R.get_piv(j)
        if l != -1 and R.labels[j] == z:
            cols_not_reduced_to_zero.append(j)
        if l != -1 and pivs[l] == -1:
            pivs[l] = j
        if l != -1 and pivs[l] > j:
            p = pivs[l]
            label = R.labels[p]
            heapq.heappush(X[label[1]], p)
            pivs[l] = j
    
    return(R, pivs, X, cols_not_reduced_to_zero)
    
#def BiRedSubSlave_v2(R, V, z, pivs):
#        
#    cols_reduced_to_zero = []
#    
#    Indices = [j for j in range(R.num_cols()) if R.labels[j][0] <= z[0] and R.labels[j][1] == z[1]]
#    
#    for j in Indices:
#        l = R.get_piv(j)
#        while l != -1 and pivs[l] != -1 and pivs[l] < j:
#            k = pivs[l]
#            R.add_column(k,j)
#            V.add_column(k,j)
#            if R.get_piv(j) == -1:
#                cols_reduced_to_zero.append(j)
#            l = R.get_piv(j)
#        if l != -1 and (pivs[l] == -1 or pivs[l] > j):
#            pivs[l] = j
#    
#    return(R, V, pivs, cols_reduced_to_zero)
    
#def BiRedSubSlave_v2(R, V, z, pivs):
#        
#    cols_reduced_to_zero = []
#    
#    Indices = [j for j in range(R.num_cols()) if R.labels[j][0] <= z[0] and R.labels[j][1] == z[1]]
#    
#    for j in Indices:
#        l = R.get_piv(j)
#        while l != -1 and pivs[l] != -1 and pivs[l] < j:
#            k = pivs[l]
#            R.add_column(k,j)
#            V.add_column(k,j)
#            l = R.get_piv(j)
#            if l == -1:
#                cols_reduced_to_zero.append(j)
#        if l != -1 and (pivs[l] == -1 or pivs[l] > j):
#            pivs[l] = j
#    
#    return(R, V, pivs, cols_reduced_to_zero)
    
def BiRedSubSlave_v2(R, V, z, pivs, Indices):
        
    cols_reduced_to_zero = []
    
    for j in Indices:
        l = R.get_piv(j)
        while l != -1 and pivs[l] != -1 and pivs[l] < j:
            k = pivs[l]
            R.add_column(k,j)
            V.add_column(k,j)
            l = R.get_piv(j)
            if l == -1:
                cols_reduced_to_zero.append(j)
        if l != -1 and (pivs[l] == -1 or pivs[l] > j):
            pivs[l] = j
    
    return(R, V, pivs, cols_reduced_to_zero)
    
def BiRedSubSlave_v3(R, V, z, pivs, current_columns, X):
        
    cols_reduced_to_zero = []
    
    for j in current_columns:
        heapq.heappush(X[z[1]], j)
    
    while len(X[z[1]]) != 0:
        j = heapq.heappop(X[z[1]])
        l = R.get_piv(j)
        while l != -1 and pivs[l] != -1 and pivs[l] < j:
            k = pivs[l]
            R.add_column(k, j)
            V.add_column(k, j)
            l = R.get_piv(j)
            if l == -1:
                cols_reduced_to_zero.append(j)
        if l != -1 and pivs[l] == -1:
            pivs[l] = j
        if l != -1 and pivs[l] > j:
            p = pivs[l]
            label = R.labels[p]
            heapq.heappush(X[label[1]], p)
            pivs[l] = j
    
    return(R, V, pivs, X, cols_reduced_to_zero)
    
def KerBetti_v2(R):
    
    betti_0 = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    gr = grid_v2(R)
    
    for z in gr:
        R, pivs, cols_reduced_to_zero = BiRedSub_v2(R, z, pivs)
        if cols_reduced_to_zero > 0:
            betti_0.append((z, cols_reduced_to_zero))
    
    return(R, betti_0)
    
def sort_labels(R):
    
    labels = R.labels
    
    sorted_labels = {}
    current_label = -1
    
    for j, label in enumerate(labels):
        if label[1] == current_label:
            sorted_labels[label[1]].append((j, label[0]))
        else:
            sorted_labels[label[1]] = [(j, label[0])]
            current_label = label[1]
        
    return(sorted_labels)
    
# uses scipy.sparse to turn the output into a bigraded matrix
#def KerBasis_v2(R):
#    
#    V = R.identity_matrix()
#    B_ker = []
#    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
#    gr = grid_v2(R)
#    
#    for z in gr:
#        R, V, pivs, cols_reduced_to_zero = BiRedSubSlave_v2(R, V, z, pivs)
#        for j in cols_reduced_to_zero:
#            B_ker.append((V.get_col(j), z))
#     
#    n = R.num_cols()
#    m = len(B_ker)
#    labels = [B_ker[j][1] for j in range(m)]
#    matrix = scipy.sparse.lil_matrix((m,n), dtype='int')
#    for j in range(m):
#        indices = B_ker[j][0].matrix.rows[0]
#        matrix[j, indices] = [1 for i in range(len(indices))]
#        
#    B = BiGradedMatrix_lil(labels, matrix)
#    
#    return(B)
    
# uses scipy.sparse to turn the output into a bigraded matrix
def KerBasis_v2(R):
    
    V = R.identity_matrix()
    B_ker = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    
    start = time.time()
    gr = grid_v2(R)
    end = time.time()
    
    print("Computing grid took " + str(end - start))
    
    start = time.time()
    sorted_labels = sort_labels(R)
    end = time.time()
    
    print("Sorting labels took " + str(end - start))
    
    for z in gr:
        Indices = []
        for item in sorted_labels[z[1]]:
            if item[1] <= z[0]:
                Indices.append(item[0])
            else:
                break
        R, V, pivs, cols_reduced_to_zero = BiRedSubSlave_v2(R, V, z, pivs, Indices)
        for j in cols_reduced_to_zero:
            B_ker.append((V.get_col(j), z))
    
    start = time.time()
    n = R.num_cols()
    m = len(B_ker)
    labels = [B_ker[j][1] for j in range(m)]
    matrix = scipy.sparse.lil_matrix((m,n), dtype='int')
    for j in range(m):
        indices = B_ker[j][0].matrix.rows[0]
        matrix[j, indices] = [1 for i in range(len(indices))]
        
    B = BiGradedMatrix_lil(labels, matrix)
    end = time.time()
    
    print("Assembling B took " + str(end-start))
    
    return(B)
    
# uses scipy.sparse to turn the output into a bigraded matrix
def KerBasis_v3(R):
    
    V = R.identity_matrix()
    B_ker = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    
#    start = time.time()
    gr = grid_v3(R)
#    end = time.time()
    
#    print("Computing grid took " + str(end - start))
    
#    start = time.time()
    sorted_labels = sort_labels(R)
#    end = time.time()
    
#    print("Sorting labels took " + str(end - start))
    
    for z in gr:
        Indices = []
        for item in sorted_labels[z[1]]:
            if item[1] <= z[0]:
                Indices.append(item[0])
            else:
                break
        R, V, pivs, cols_reduced_to_zero = BiRedSubSlave_v2(R, V, z, pivs, Indices)
        for j in cols_reduced_to_zero:
            B_ker.append((V.matrix.rows[j], z))
    
#    start = time.time()
    n = R.num_cols()
    m = len(B_ker)
    labels = [B_ker[j][1] for j in range(m)]
    matrix = scipy.sparse.lil_matrix((m,n), dtype='int')
    
    for j in range(m):
        matrix.rows[j] = B_ker[j][0]
        matrix.data[j] = [1 for i in range(len(B_ker[j][0]))]
        
    B = BiGradedMatrix_lil(labels, matrix)
#    end = time.time()
    
#    print("Assembling B took " + str(end-start))
    return(B)
    
def columns_by_label(R):
    
    labels = R.labels
    
    columns_by_label = {}
    current_label = (-1, -1)
    
    for j, label in enumerate(labels):
        if label == current_label:
            columns_by_label[label].append(j)
        else:
            columns_by_label[label] = [j]
            current_label = label
            
    return(columns_by_label)
    
# uses scipy.sparse to turn the output into a bigraded matrix
def KerBasis_v4(R):
    
    V = R.identity_matrix()
    B_ker = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    columns = columns_by_label(R)
    
#    start = time.time()
    gr = grid_v3(R)
#    end = time.time()
    
#    print("Computing grid took " + str(end - start))
    
#    start = time.time()
    y_coords_w_rep = []
    for z in gr:
        y_coords_w_rep.append(z[1])
    y_coords = list(set(y_coords_w_rep))
    X = [[] for i in range(len(y_coords))]
#    end = time.time()
    
#    print("Initializing X took " + str(end-start))
    
    for z in gr:
        if z in columns:
            current_columns = columns[z]
        else:
            current_columns = []
        R, V, pivs, X, cols_reduced_to_zero = BiRedSubSlave_v3(R, V, z, pivs, current_columns, X)
        for j in cols_reduced_to_zero:
            B_ker.append((V.matrix.rows[j], z))
     
#    start = time.time()
    n = R.num_cols()
    m = len(B_ker)
    labels = [B_ker[j][1] for j in range(m)]
    matrix = scipy.sparse.lil_matrix((m,n), dtype='int')
    
    for j in range(m):
        matrix.rows[j] = B_ker[j][0]
        matrix.data[j] = [1 for i in range(len(B_ker[j][0]))]
        
    B = BiGradedMatrix_lil(labels, matrix)
#    end = time.time()
    
#    print("Assembling B took " + str(end-start))
    return(B)
    
#def MinGens_v2(R):
#    
#    S = []
#    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
#    gr = grid_v2(R)
#    
#    for z in gr:
#        R, pivs, cols_not_reduced_to_zero = BiRedSub_MinGens_v2(R, z, pivs)
#        for j in cols_not_reduced_to_zero:
#            S.append((R.get_col(j), z))
#    
#    return(S)
    
def MinGens_v2(R):
    
    S = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    gr = grid_v3(R)
    
    sorted_labels = sort_labels(R)
    
    for z in gr:
        Indices = []
        for item in sorted_labels[z[1]]:
            if item[1] <= z[0]:
                Indices.append(item[0])
            else:
                break
        R, pivs, cols_not_reduced_to_zero = BiRedSub_MinGens_v2(R, z, pivs, Indices)
        for j in cols_not_reduced_to_zero:
            S.append((R.get_col(j), z))
    
    return(S)
    
def MinGens_v2a(R):
    
    S = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    gr = grid_v3(R)
    
    sorted_labels = sort_labels(R)
    
    for z in gr:
        Indices = []
        for item in sorted_labels[z[1]]:
            if item[1] <= z[0]:
                Indices.append(item[0])
            else:
                break
        R, pivs, cols_not_reduced_to_zero = BiRedSub_MinGens_v2(R, z, pivs, Indices)
        for j in cols_not_reduced_to_zero:
            S.append((R.matrix.rows[j], z))
    
    return(S)
    
def MinGens_v3(R):
    
    S = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    columns = columns_by_label(R)
    gr = grid_v3(R)
    
    y_coords_w_rep = []
    for z in gr:
        y_coords_w_rep.append(z[1])
    y_coords = list(set(y_coords_w_rep))
    X = [[] for i in range(len(y_coords))]
    
    for z in gr:
        if z in columns:
            current_columns = columns[z]
        else:
            current_columns = []
        R, pivs, X, cols_not_reduced_to_zero = BiRedSub_MinGens_v3(R, z, pivs, current_columns, X)
        for j in cols_not_reduced_to_zero:
            S.append((R.get_col(j), z))
    
    return(S)
    
def MinGens_v4(R):
    
    S = []
    pivs = np.full(shape=R.num_rows(), fill_value=-1, dtype='int', order='C')
    columns = columns_by_label(R)
    gr = grid_v3(R)
    
    y_coords_w_rep = []
    for z in gr:
        y_coords_w_rep.append(z[1])
    y_coords = list(set(y_coords_w_rep))
    X = [[] for i in range(len(y_coords))]
    
    for z in gr:
        if z in columns:
            current_columns = columns[z]
        else:
            current_columns = []
        R, pivs, X, cols_not_reduced_to_zero = BiRedSub_MinGens_v3(R, z, pivs, current_columns, X)
        for j in cols_not_reduced_to_zero:
            S.append((R.matrix.rows[j], z))
    
    return(S)
    
def pivot_array_of_B_ker(A):
    
    pivots = np.full(shape=A.num_rows(), fill_value=-1, dtype='int', order='C')
    
    for j in range(A.num_cols()):
        pivots[A.get_piv(j)] = j
        
    return(pivots)
    
#def Solve_v2(A, pivots, B):
#    
#    X = A.zero_vector()
#    
#    for i in reversed(range(A.num_rows())):
#        if B.entry_is_nonzero(i, 0):
#            if pivots[i] == -1:
#                print('I could not find the pivot!')
#            else:
#                j = pivots[i]
#                X.insert_value(j, 0, 1)
#                col = A.matrix.rows[j]
#                B.add_external_column(col, 0)
#                
#    return(X)
    
def Solve_v2(A, pivots, B):
    
    X = A.zero_vector()
    
    p = B.get_piv(0)
    
    while p != -1:
        k = pivots[p]
        X.insert_value(k, 0, 1)
        col = A.matrix.rows[k]
        B.add_external_column(col, 0)
        p = B.get_piv(0)
        
    return(X)
    
def Solve_v3(A, pivots, B):
    
    X = []
    
    while len(B) != 0:
        p = B[-1]
        k = pivots[p]
        X.append(k)
        col = A.matrix.rows[k]
        B = add_columns(B, col)
        
    X.sort()
        
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
    
def add_columns(col_k, col_j):

    col_sum = []

    length_k = len(col_k)
    length_j = len(col_j)

    position_k = 0
    position_j = 0

    while position_k < length_k or position_j < length_j:
        if position_k >= length_k:
            for i in range(position_j, length_j):
                col_sum.append(col_j[i])
            break
        if position_j >= length_j:
            for i in range(position_k, length_k):
                col_sum.append(col_k[i])
            break
        if col_k[position_k] == col_j[position_j]:
            position_k = position_k + 1
            position_j = position_j + 1
            continue
        if col_k[position_k] < col_j[position_j]:
            col_sum.append(col_k[position_k])
            position_k = position_k + 1
            continue
        if col_k[position_k] > col_j[position_j]:
            col_sum.append(col_j[position_j])
            position_j = position_j + 1
            continue

    return(col_sum)
    
def MinimizePres_v3(P, row_labels):
    
    n = P.num_cols()
    pivs = {}
    
    for j in range(n):
        p = P.get_piv(j)
        if P.labels[j] == row_labels[p]:
            if p in pivs:
                pivs[p].append(j)
            else:
                pivs[p] = [j]
    
    new_cols = []

    for j in range(n):
        
        p = P.get_piv(j)
        
        while P.labels[j] == row_labels[p] and len(pivs[p]) > 1 and pivs[p][0] < j:
            k = pivs[p][0]
            P.add_column(k, j)
            p = P.get_piv(j)
            if pivs.get(p) == None:
                if P.labels[j] == row_labels[p]:
                    pivs[p] = [j]
                else:
                    break
        
        if P.labels[j] != row_labels[p]:
            old_col = P.matrix.rows[j]
            new_col = []
            while len(old_col) != 0:
                p = old_col[-1]
                if p in pivs and pivs[p][0] < j:
                    k = pivs[p][0]
                    col_k = P.matrix.rows[k]
                    old_col = add_columns(old_col, col_k)
#                    if len(old_col) != 0:
#                        if old_col[-1] >= p:
#                            print(j, k, p, old_col, col_k)
#                            return(None)
                else:
                    new_col.append(p)
                    del old_col[-1]
            new_col.sort()
            new_cols.append((j, new_col))
            
    m = len(new_cols)
    k = P.num_rows()
    X = scipy.sparse.lil_matrix((m, k), dtype='int')
    
    for j in range(m):
        X.rows[j] = new_cols[j][1]
        X.data[j] = [1 for i in range(len(new_cols[j][1]))]
        
    labels = []
    
    for col in new_cols:
        labels.append(P.labels[col[0]])
        
    B = BiGradedMatrix_lil(labels, X)

    rows = [i for i in range(k)]
    for i in pivs.keys():
        rows.remove(i)
        
    new_row_labels = [row_labels[i] for i in rows]
    B.keep_rows(rows)
    
    return(B, new_row_labels)
    
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
    
def MinimalPres(delta2, delta1, use_opt1=False, use_opt2=False):
    
    beginning = time.time()
    
    if use_opt1 == True:
        S = MinGens_v3(delta2)
    else:
        S = MinGens_v2(delta2)
        
    end = time.time()    
    print("Running MinGens took " + str(end-beginning))
    
    start = time.time()
    
    if use_opt1 == True:
        B = KerBasis_v4(delta1)
    else:
        B = KerBasis_v3(delta1)
        
    end = time.time()    
    print("Running KerBasis took " + str(end-start))
    
    start = time.time()
    
    pivots = pivot_array_of_B_ker(B)
    n = B.num_cols()
    m = len(S)
    labels = [S[j][1] for j in range(m)]
    matrix = scipy.sparse.lil_matrix((m,n), dtype='int')
    P = BiGradedMatrix_lil(labels, matrix)
    
    for j in range(m):
        P.matrix[j, :] = Solve_v2(B, pivots, S[j][0]).matrix
        
    end = time.time()
    print("Computing P took " + str(end-start))
    
    start = time.time()
    
    if use_opt2 == True:
        row_labels = B.labels
        P, row_labels = MinimizePres_v3(P, row_labels)
    else:
        row_labels = B.labels
        P, row_labels = MinimizePres_v2(P, row_labels)
    
    end = time.time()
    print("Running MinimizePres took " + str(end-start))
    print("Total running time was " + str(end-beginning))

    return(P, row_labels)
    
# (delta2, delta1) is an FI-Rep; so, delta2, delta1 are bigraded matrices,
# where delta2 represents a map F^2 -> F^1, 
# and delta1 represents a map F^1 -> F^0
    
def MinimalPres_v2(delta2, delta1, use_opt1=False, use_opt2=False):
    
    beginning = time.time()
    
    if use_opt1 == True:
        S = MinGens_v4(delta2)
    else:
        S = MinGens_v2a(delta2)
        
    end = time.time()    
    print("Running MinGens took " + str(end-beginning))
    
    start = time.time()
    
    if use_opt1 == True:
        B = KerBasis_v4(delta1)
    else:
        B = KerBasis_v3(delta1)
        
    end = time.time()    
    print("Running KerBasis took " + str(end-start))
    
    start = time.time()
    
    pivots = pivot_array_of_B_ker(B)
    n = B.num_cols()
    m = len(S)
    labels = [S[j][1] for j in range(m)]
    M = scipy.sparse.lil_matrix((m,n), dtype='int')
    
    for j in range(m):
        col = Solve_v3(B, pivots, S[j][0])
        M.rows[j] = col
        M.data[j] = [1 for i in range(len(col))]
        
    P = BiGradedMatrix_lil(labels, M)
        
    end = time.time()
    print("Computing P took " + str(end-start))
    
    start = time.time()
    
    if use_opt2 == True:
        row_labels = B.labels
        P, row_labels = MinimizePres_v3(P, row_labels)
    else:
        row_labels = B.labels
        P, row_labels = MinimizePres_v2(P, row_labels)
    
    end = time.time()
    print("Running MinimizePres took " + str(end-start))
    print("Total running time was " + str(end-beginning))

    return(P, row_labels)