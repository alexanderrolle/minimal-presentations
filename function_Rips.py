#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 19:33:50 2020

@author: arolle
"""

from itertools import combinations
from sklearn.metrics.pairwise import euclidean_distances

# Computes the p-1, p, and p+1 simplices of the function-Rips bifiltration,
# and prints them to a file.
#
# Input: 
# (1) a numpy array, or list, of points in euclidean space;
# (2) a list of function values for each point;
# (3) the homology dimension p;
# (4) the name of the file to print to.
#
# The function returns the number of simplices it computes.
    
def function_Rips(points, function_vals, p, file_name):
    
    distances = euclidean_distances(points)
    
    simplices_to_print = []
    
    p_minus_one_simplices = combinations(range(len(points)), p)

    for simplex in p_minus_one_simplices:
        simplex_to_print = ""
        bidegree = bidegree_of_simplex(simplex, function_vals, distances)
        for vertex in simplex:
            simplex_to_print += str(vertex)
            simplex_to_print += " "
        simplex_to_print += "; "
        simplex_to_print += str(bidegree[0])
        simplex_to_print += " "
        simplex_to_print += str(bidegree[1])
        simplices_to_print.append(simplex_to_print)
        
    p_simplices = combinations(range(len(points)), p+1)

    for simplex in p_simplices:
        simplex_to_print = ""
        bidegree = bidegree_of_simplex(simplex, function_vals, distances)
        for vertex in simplex:
            simplex_to_print += str(vertex)
            simplex_to_print += " "
        simplex_to_print += "; "
        simplex_to_print += str(bidegree[0])
        simplex_to_print += " "
        simplex_to_print += str(bidegree[1])
        simplices_to_print.append(simplex_to_print)

    p_plus_one_simplices = combinations(range(len(points)), p+2)

    for simplex in p_plus_one_simplices:
        simplex_to_print = ""
        bidegree = bidegree_of_simplex(simplex, function_vals, distances)
        for vertex in simplex:
            simplex_to_print += str(vertex)
            simplex_to_print += " "
        simplex_to_print += "; "
        simplex_to_print += str(bidegree[0])
        simplex_to_print += " "
        simplex_to_print += str(bidegree[1])
        simplices_to_print.append(simplex_to_print)

    f = open(file_name, "w")
    for simplex in simplices_to_print:
        print(simplex, file=f)
    f.close()
    
    return(len(simplices_to_print))
    
def bidegree_of_simplex(simplex, function_vals, distances):
    
    simplex_function_vals = []
    for vertex in simplex:
        simplex_function_vals.append(function_vals[vertex])
    
    k = min(simplex_function_vals)
    
    
    if len(simplex) == 1:
        s = 0
    else:
        simplex_distances = []
        for pair in combinations(simplex, 2):
            simplex_distances.append(distances[pair])
        
        s = max(simplex_distances)
    
    return((s, k))