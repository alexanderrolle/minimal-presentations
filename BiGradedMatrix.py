#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 16:31:33 2020

@author: arolle
"""

import numpy as np
import scipy.sparse

class BiGradedMatrix:
    
    def __init__(self, labels, matrix):
        self.labels = labels
        self.matrix = matrix
        
    def num_rows(self):
        return(self.matrix.shape[0])
    
    def num_cols(self):
        return(self.matrix.shape[1])
        
    # reduce the entries in a matrix mod 2
    def reduce(self):
        M = self.matrix
        N = scipy.sparse.csc_matrix((M.data % 2, M.indices, M.indptr), shape = M.shape)
        scipy.sparse.csc_matrix.eliminate_zeros(N)
        self.matrix = N
        
    # add column k to column j, 
    # by multiplying on the right with an elementary matrix
    def add_column(self, k, j):
        
        m = self.num_cols()
        data = np.array([1])
        row = np.array([k])
        col = np.array([j])
        J = scipy.sparse.csc_matrix((data, (row, col)), shape=(m, m), dtype='int')
        I = scipy.sparse.identity(n=m, dtype='int', format='csc')
        
        self.matrix = self.matrix * (I + J)
        self.reduce()
        
    def get_piv(self, j):
        M = self.matrix
        col = scipy.sparse.csc_matrix.getcol(M, j)
        if len(col.indices) == 0:
            return(-1)
        else:
            return(max(col.indices))