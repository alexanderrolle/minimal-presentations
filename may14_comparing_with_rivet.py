#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 13:49:42 2020

@author: arolle
"""

import time
from FIRep_v3 import *
from MinimalPres_v2 import *
from ReadingRivet_v2 import *

start = time.time()
R, Q = FIRep_v3("code/minimal_presentations/k_fold_filtrations/bifiltrations/bifiltration_random100_1.txt", 1)
end = time.time()
print("Producing FIRep took " + str(end-start))

P, row_labels = MinimalPres_v2(R, Q)

U, resorted_row_labels = permute_matrix(P, row_labels)

rivet_row_labels, rivet_col_labels = get_labels_from_rivet("code/minimal_presentations/rivet_output/rivet_random100_1.txt")

if resorted_row_labels == rivet_row_labels:
    print("Row labels agree with rivet")
else:
    print("Row labels do not agree with rivet")
    
if U.labels == rivet_col_labels:
    print("Column labels agree with rivet")
else:
    print("Column labels do not agree with rivet")
    
rivet_matrix = get_matrix_from_rivet("code/minimal_presentations/rivet_output/rivet_random100_1.txt")

if U.matrix.toarray().tolist() == rivet_matrix.toarray().tolist():
    print("Matrix agrees with rivet")
else:
    print("Matrix does not agree with rivet")