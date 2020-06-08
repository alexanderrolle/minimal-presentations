#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 11:40:17 2020

@author: arolle
"""

from FIRep_v3 import *
import numpy as np
import scipy.sparse
from BiGradedMatrix_v2 import *
import time
from MinimalPres_v2 import *

R, Q = FIRep_v3("code/minimal_presentations/k_fold_filtrations/bifiltrations/bifiltration_random100_1.txt", 1)

S = MinGens_v4(R)

c_1 = 3136
c_2 = 3137

print("Before reparameterizing, column c_1 is " + str(S[c_1]))
print("Before reparameterizing, column c_2 is " + str(S[c_2]))
print("Label of row 3010 is " + str(Q.labels[3010]))
print("Label of row 3035 is " + str(Q.labels[3035]))

B = KerBasis_v4(Q)

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

row_labels = B.labels

print("After reparameterizing, column c_1 is " + str(P.matrix.rows[c_1]))
print("After reparameterizing, column c_2 is " + str(P.matrix.rows[c_2]))
print("Label of c_1 is " + str(P.labels[c_1]))
print("Label of c_2 is " + str(P.labels[c_2]))
print("Label of row 2848 is " + str(row_labels[2848]))