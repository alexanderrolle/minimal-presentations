#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 12:36:37 2020

@author: arolle

This version of the BiGradedMatrix class implements sparse matrices
using scipy.sparse.lil_matrix

"""

#import numpy as np
import scipy.sparse
#import heapq

class BiGradedMatrix_lil:
    
    def __init__(self, labels, matrix):
        self.labels = labels
        self.matrix = matrix
        
    def num_rows(self):
        return(self.matrix.shape[1])
    
    def num_cols(self):
        return(self.matrix.shape[0])
        
    def identity_matrix(self):
        n = self.num_cols()
        I = scipy.sparse.identity(n=n, dtype='int', format='lil')
        Id = BiGradedMatrix_lil(labels=self.labels, matrix=I)
        return(Id)
        
    def zero_vector(self):
        m = self.num_cols()
        x = scipy.sparse.lil_matrix((1, m), dtype='int')
        X = BiGradedMatrix_lil(labels=[0], matrix=x)
        return(X)
        
    def get_col(self, j):
        col = self.matrix.getrow(j)
        label = self.labels[j]
        R = BiGradedMatrix_lil(label, col)
        return(R)
        
    def entry_is_nonzero(self, i, j):
        return(self.matrix[(j, i)] != 0)
        
#    def entry_is_nonzero(self, i, j):
#        row_entries = self.matrix.getrow(j).rows[0]
#        for row in row_entries:
#            if row == i:
#                return(True)
#            
#        return(False)
        
    # insert one at position (n, m)
#    def insert_one(self, n, m):
#        col = self.matrix.rows[m]
#        
#        if len(col) == 0:
#            self.matrix.rows[m].append(n)
#            return
#        else:
#            for i in range(len(col)):
#                if col[i] == n:
#                    return
#                if col[i] > n:
#                    self.matrix.rows[m].insert(i, n)
#                    return
#            
#        self.matrix.rows[m].append(n)
        
    #insert n at position (i, j)           
    def insert_value(self, i, j, n):
        self.matrix[(j, i)] = n
        
    # add column k to column j
#    def add_column(self, k, j):
#        row_k = self.matrix.getrow(k).rows[0]
#        row_j = self.matrix.getrow(j).rows[0]
#        
#        x = heapq.merge(row_k, row_j)
#        
#        y = []
#        last_index = -1
#        for i in x:
#            if i == last_index:
#                y.remove(i)
#            else:
#                last_index = i
#                y.append(i)
#                
#        num_rows = self.matrix.shape[1]
#        self.matrix[j, range(num_rows)] = [0 for i in range(num_rows)]
#        self.matrix[j, y] = [1 for i in range(len(y))]
    
    # add column k to column j
    def add_column(self, k, j):
        
        col_k = self.matrix.rows[k]
        col_j = self.matrix.rows[j]
        
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
                
        self.matrix.rows[j] = col_sum
        self.matrix.data[j] = [1 for i in range(len(col_sum))]
        
    # add col to column j
    def add_external_column(self, col, j):
        
        col_k = col
        col_j = self.matrix.rows[j]
        
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
                
        self.matrix.rows[j] = col_sum
        self.matrix.data[j] = [1 for i in range(len(col_sum))]
        
    def get_piv(self, j):
        row = self.matrix.rows[j]
        if len(row) == 0:
            return(-1)
        else:
            return(row[-1])
            
    def keep_columns(self, cols):
        self.matrix = self.matrix[cols, :]
        
    def keep_rows(self, rows):
        self.matrix = self.matrix[:, rows]