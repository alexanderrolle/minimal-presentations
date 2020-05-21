#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 11:40:03 2020

@author: arolle
"""

import time
from FIRep_v3 import *
from MinimalPres_v2 import *
from ReadingRivet_v2 import *

bifiltration = "code/minimal_presentations/k_fold_filtrations/bifiltrations/bifiltration_random100_1.txt"

start = time.time()
R, Q = FIRep_v3(bifiltration, 1)
end = time.time()
print("Producing FIRep took " + str(end-start))

P, row_labels = MinimalPres(R, Q, use_opt=True)

U, resorted_row_labels = permute_matrix(P, row_labels)

rivet_output = "code/minimal_presentations/rivet_output/rivet_random100_1.txt"

rivet_row_labels, rivet_col_labels = get_labels_from_rivet(rivet_output)

if resorted_row_labels == rivet_row_labels:
    print("Row labels agree with rivet")
else:
    print("Row labels do not agree with rivet")
    
if U.labels == rivet_col_labels:
    print("Column labels agree with rivet")
else:
    print("Column labels do not agree with rivet")
    
rivet_matrix = get_matrix_from_rivet(rivet_output)

if U.matrix.toarray().tolist() == rivet_matrix.toarray().tolist():
    print("Matrix agrees with rivet")
else:
    print("Matrix does not agree with rivet")