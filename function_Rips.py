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
        
    function_vals_negative = []
    for k in function_vals:
        function_vals_negative.append(-k)
    
    distances = euclidean_distances(points)
    
    simplices_to_print = []
    
    p_minus_one_simplices = combinations(range(len(points)), p)

    for simplex in p_minus_one_simplices:
        simplex_to_print = ""
        bidegree = bidegree_of_simplex(simplex, function_vals_negative, distances)
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
        bidegree = bidegree_of_simplex(simplex, function_vals_negative, distances)
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
        bidegree = bidegree_of_simplex(simplex, function_vals_negative, distances)
        for vertex in simplex:
            simplex_to_print += str(vertex)
            simplex_to_print += " "
        simplex_to_print += "; "
        simplex_to_print += str(bidegree[0])
        simplex_to_print += " "
        simplex_to_print += str(bidegree[1])
        simplices_to_print.append(simplex_to_print)

    f = open(file_name, "w")
    print("bifiltration", file=f)
    print("s", file=f)
    print("k", file=f)
    for simplex in simplices_to_print:
        print(simplex, file=f)
    f.close()
    
    return(len(simplices_to_print))
    
# The difference between function_Rips and function_Rips_v2
# is that v2 reindexes the k values so they are integers.    
    
def function_Rips_v2(points, function_vals, p, file_name):
    
    function_vals_negative = []
    for k in function_vals:
        function_vals_negative.append(-k)
    
    function_vals_distinct = list(set(function_vals_negative))
    function_vals_distinct.sort()
    function_vals_dict = {function_vals_distinct[i] : i for i in range(len(function_vals_distinct))}
    
    distances = euclidean_distances(points)
    
    simplices_to_print = []
    
    p_minus_one_simplices = combinations(range(len(points)), p)

    for simplex in p_minus_one_simplices:
        simplex_to_print = ""
        bidegree = bidegree_of_simplex(simplex, function_vals_negative, distances)
        for vertex in simplex:
            simplex_to_print += str(vertex)
            simplex_to_print += " "
        simplex_to_print += "; "
        simplex_to_print += str(bidegree[0])
        simplex_to_print += " "
        simplex_to_print += str(function_vals_dict[bidegree[1]])
        simplices_to_print.append(simplex_to_print)
        
    p_simplices = combinations(range(len(points)), p+1)

    for simplex in p_simplices:
        simplex_to_print = ""
        bidegree = bidegree_of_simplex(simplex, function_vals_negative, distances)
        for vertex in simplex:
            simplex_to_print += str(vertex)
            simplex_to_print += " "
        simplex_to_print += "; "
        simplex_to_print += str(bidegree[0])
        simplex_to_print += " "
        simplex_to_print += str(function_vals_dict[bidegree[1]])
        simplices_to_print.append(simplex_to_print)

    p_plus_one_simplices = combinations(range(len(points)), p+2)

    for simplex in p_plus_one_simplices:
        simplex_to_print = ""
        bidegree = bidegree_of_simplex(simplex, function_vals_negative, distances)
        for vertex in simplex:
            simplex_to_print += str(vertex)
            simplex_to_print += " "
        simplex_to_print += "; "
        simplex_to_print += str(bidegree[0])
        simplex_to_print += " "
        simplex_to_print += str(function_vals_dict[bidegree[1]])
        simplices_to_print.append(simplex_to_print)

    f = open(file_name, "w")
    print("bifiltration", file=f)
    print("s", file=f)
    print("k", file=f)
    for simplex in simplices_to_print:
        print(simplex, file=f)
    f.close()
    
    return(len(simplices_to_print))
    
def bidegree_of_simplex(simplex, function_vals_negative, distances):
    
    simplex_function_vals = []
    for vertex in simplex:
        simplex_function_vals.append(function_vals_negative[vertex])
    
    k = max(simplex_function_vals)
    
    
    if len(simplex) == 1:
        s = 0
    else:
        simplex_distances = []
        for pair in combinations(simplex, 2):
            simplex_distances.append(distances[pair])
        
        s = max(simplex_distances)
    
    return((s, k))
    
# Computes an FIRep of H_p of the function-Rips bifiltration,
# and prints it to a file (in the Rivet format).
#
# Input: 
# (1) a numpy array, or list, of points in euclidean space;
# (2) a list of function values for each point;
# (3) the homology dimension p;
# (4) the name of the file to print to.
    
def function_Rips_FIRep(points, function_vals, p, file_name):
    
    function_vals_negative = []
    for k in function_vals:
        function_vals_negative.append(-k)
    
    distances = euclidean_distances(points)
    
    p_minus_one_simplices = combinations(range(len(points)), p)
    
    p_minus_one_simplices_list = []
    for simplex in p_minus_one_simplices:
        p_minus_one_simplices_list.append((simplex, bidegree_of_simplex(simplex, function_vals_negative, distances)))
    p_minus_one_simplices_list.sort(key=lambda x : (x[1][1], x[1][0]))
    
    p_minus_one_simplices_dict = {}
    for counter, item in enumerate(p_minus_one_simplices_list):
        p_minus_one_simplices_dict[item[0]] = counter
    
    
    p_simplices = combinations(range(len(points)), p+1)
    
    p_simplices_list = []
    for simplex in p_simplices:
        p_simplices_list.append((simplex, bidegree_of_simplex(simplex, function_vals_negative, distances)))
    p_simplices_list.sort(key=lambda x : (x[1][1], x[1][0]))
    
    p_simplices_dict = {}
    for counter, item in enumerate(p_simplices_list):
        p_simplices_dict[item[0]] = counter
        
    
    p_plus_one_simplices = combinations(range(len(points)), p+2)
    
    p_plus_one_simplices_list = []
    for simplex in p_plus_one_simplices:
        p_plus_one_simplices_list.append((simplex, bidegree_of_simplex(simplex, function_vals_negative, distances)))
    p_plus_one_simplices_list.sort(key=lambda x : (x[1][1], x[1][0]))
    
#    p_plus_one_simplices_dict = {}
#    for counter, item in enumerate(p_plus_one_simplices_list):
#        p_plus_one_simplices_dict[item[0]] = counter
        
    
    f = open(file_name, "w")
    print("firep", file=f)
    print("s", file=f)
    print("k", file=f)
    print(str(len(p_plus_one_simplices_list)) + ' ' + str(len(p_simplices_list)) + ' ' + str(len(p_minus_one_simplices_list)), file=f)
    
    for item in p_plus_one_simplices_list:
        
        line_to_print = ''
        line_to_print += str(item[1][0])
        line_to_print += ' '
        line_to_print += str(item[1][1])
        line_to_print += ' ;'
        
        col = []
        vertices = item[0]
    
        for j in range(len(vertices)):
            face = [vertices[k] for k in range(len(vertices)) if k != j]
            row_index = p_simplices_dict[tuple(face)]
            col.append(row_index)
        col.sort()
        
        for row_index in col:
            line_to_print += ' '
            line_to_print += str(row_index)
            
        print(line_to_print, file=f)
        
    for item in p_simplices_list:
        
        line_to_print = ''
        line_to_print += str(item[1][0])
        line_to_print += ' '
        line_to_print += str(item[1][1])
        line_to_print += ' ;'
        
        col = []
        vertices = item[0]
    
        for j in range(len(vertices)):
            face = [vertices[k] for k in range(len(vertices)) if k != j]
            row_index = p_minus_one_simplices_dict[tuple(face)]
            col.append(row_index)
        col.sort()
        
        for row_index in col:
            line_to_print += ' '
            line_to_print += str(row_index)
            
        print(line_to_print, file=f)
    
    f.close()
