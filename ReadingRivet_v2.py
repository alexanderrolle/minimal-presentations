#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 20:46:31 2020

@author: arolle

This version of Reading Rivet is compatible with BiGradedMatrix_v2
"""

import scipy.sparse
from BiGradedMatrix_v2 import *

def get_first_bigrade(bigrades):
    
    grade_0 = bigrades[0]
    for i in range(1, len(bigrades)):
        grade_0 += bigrades[i]
        if bigrades[i] == ',':
            break
            
    for c in grade_0:
        bigrades = bigrades.replace(c, '', 1)
    
    grade_0 = grade_0.replace('(', '')
    grade_0 = grade_0.replace(',', '')
    grade_0 = int(grade_0)
    
    grade_1 = bigrades[0]
    for i in range(1, len(bigrades)):
        grade_1 += bigrades[i]
        if bigrades[i] == ')':
            break
            
    for c in grade_1:
        bigrades = bigrades.replace(c, '', 1)
        
    grade_1 = grade_1.replace(')', '')
    grade_1 = int(grade_1)

    return((grade_0, grade_1), bigrades)
    
def get_bigrades(bigrades):
    
    bigrades = bigrades.replace(' ', '')
    bigrades = bigrades.replace('|', '')
    bigrades = bigrades.replace('\n', '')
    
    list_of_bigrades = []
    while len(bigrades) > 0:
        bigrade, bigrades = get_first_bigrade(bigrades)
        list_of_bigrades.append(bigrade)
        
    return(list_of_bigrades)

def get_labels_from_rivet(min_pres):
                
    with open(min_pres) as f:
        for num, line in enumerate(f):
            if 'Row bigrades' in line:
                row_bigrades_num = num + 1
                
    with open(min_pres) as f:
        for num, line in enumerate(f):
            if num == row_bigrades_num:
                row_bigrades = line
                
    row_bigrades_list = get_bigrades(row_bigrades)
                
    with open(min_pres) as f:
        for num, line in enumerate(f):
            if 'Column bigrades' in line:
                col_bigrades_num = num + 1
                
    with open(min_pres) as f:
        for num, line in enumerate(f):
            if num == col_bigrades_num:
                col_bigrades = line
                
    col_bigrades_list = get_bigrades(col_bigrades)
    
    return(row_bigrades_list, col_bigrades_list)
    
def permute_matrix(P, row_labels):
    
    enumerated_col_labels = []

    for i in range(len(P.labels)):
        enumerated_col_labels.append((i, P.labels[i]))

    enumerated_col_labels.sort(key=lambda x : (x[1][1], x[1][0]))
    
    resorted_col_labels = sorted(P.labels, key=lambda x : (x[1], x[0]))
    matrix = scipy.sparse.lil_matrix((P.num_cols(), P.num_rows()), dtype='int')
    T = BiGradedMatrix_lil(resorted_col_labels, matrix)
    
    for i in range(P.num_cols()):
        T.matrix[i, :] = P.matrix[enumerated_col_labels[i][0], :]
        
    enumerated_row_labels = []
    
    for i in range(len(row_labels)):
        enumerated_row_labels.append((i, row_labels[i]))
        
    enumerated_row_labels.sort(key=lambda x : (x[1][1], x[1][0]))
    
    resorted_row_labels = sorted(row_labels, key=lambda x : (x[1], x[0]))
    matrix = scipy.sparse.lil_matrix((P.num_cols(), P.num_rows()), dtype='int')
    U = BiGradedMatrix_lil(resorted_col_labels, matrix)
    
    for i in range(P.num_rows()):
        U.matrix[:, i] = T.matrix[:, enumerated_row_labels[i][0]]
        
    U.matrix = U.matrix.transpose()
        
    return(U, resorted_row_labels)
    
def get_first_entry(col):
    entry = col[0]
    for i in range(1, len(col)):
        entry += col[i]
        if col[i] == ' ':
            break
            
    col = col.replace(entry, '', 1)
    entry.replace(' ', '')
    entry = int(entry)
    
    return(entry, col)
    
def get_matrix_from_rivet(min_pres):
    
    with open(min_pres) as f:
        for line in f:
            if 'Number of rows' in line:
                num_rows = line.replace('Number of rows:', '')
                num_rows = int(num_rows)
                
    with open(min_pres) as f:
        for line in f:
            if 'Number of columns' in line:
                num_cols = line.replace('Number of columns:', '')
                num_cols = int(num_cols)
    
    with open(min_pres) as f:
        for num, line in enumerate(f):
            if 'Column bigrades' in line:
                matrix_start = num + 2
            
    cols = []

    with open(min_pres) as f:
        for num, line in enumerate(f):
            if num >= matrix_start:
                col = line.replace('\n', '')
                cols.append(col)
                
    cols_list = []
    for col in cols:
        col_list = []
        while len(col) > 0:
            entry, col = get_first_entry(col)
            col_list.append(entry)

        cols_list.append(col_list)
        
    data = []
    for col in cols_list:
        data.append([1 for i in range(len(col))])
        
    lil_matrix = scipy.sparse.lil_matrix((num_cols, num_rows), dtype = 'int')
    for i in range(num_cols):
        lil_matrix[i, cols_list[i]] = data[i]

    csc_matrix = scipy.sparse.csc_matrix(lil_matrix)
    csc_matrix = csc_matrix.transpose()
    return(csc_matrix)