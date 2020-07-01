#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 19:41:48 2020

@author: arolle
"""

from function_Rips import *
import numpy as np
from sklearn.neighbors import KernelDensity
import time

# Choose number of points and ambient dimension for point cloud
num_points = 100
ambient_dimension = 2

# Choose dimension for persistent homology
homology_dimension = 1

# Choose the bandwidth for the kernel density estimate
band = 0.2

# Choose name of file containing the bifiltration
#file_name_1 = "july1_fR_2.txt"
#file_name_2 = "july1_fR_2_int.txt"
file_name_3 = "july1_FIRep_2.txt"

start = time.time()

points = np.random.rand(num_points, ambient_dimension)
kde = KernelDensity(kernel='gaussian', bandwidth=band).fit(points)
function_vals = np.exp(kde.score_samples(points))

#number_of_simplices = function_Rips(points, function_vals, homology_dimension, file_name_1)
#number_of_simplices_2 = function_Rips_v2(points, function_vals, homology_dimension, file_name_2)
function_Rips_FIRep(points, function_vals, homology_dimension, file_name_3)

end = time.time()

#print("Number of simplices is " + str(number_of_simplices))
print("Runtime is " + str(end-start))