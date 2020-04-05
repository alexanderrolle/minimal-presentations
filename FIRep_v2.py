#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 17:12:13 2020

@author: arolle
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 20:13:07 2020

@author: arolle
"""

import scipy.sparse
from BiGradedMatrix import *

def simplex_dimension(simplex):
    
    if simplex[0] == 'b' or simplex[0] == 'C' or simplex[0] == 'O':
        return(-1)
    
    num_spaces = 0
    for c in simplex:
        if c == " ":
            num_spaces = num_spaces + 1
        if c == ";":
            break
    return(num_spaces-1)

def get_first_vertex(simplex):
    
    vertex = []
        
    for c in simplex:
        vertex.append(c)
        if c == " ":
            break
    
    for c in vertex:
        simplex.remove(c)
        
    return(vertex)
    
def get_vertex_chars(simplex):
    
    vertices = []
    
    for c in simplex:
        vertices.append(c)
        if c == ";":
            break
    return(vertices)

def get_vertices(simplex):
    
    vertex_chars = get_vertex_chars(simplex)
    
    n = simplex_dimension(vertex_chars)
    vertices = []
    
    for i in range(n+1):
        vertex = get_first_vertex(vertex_chars)
        vertex.remove(" ")
        vertex_as_string = vertex[0]
        for i in range(1, len(vertex)):
            vertex_as_string += vertex[i]
        vertex = int(vertex_as_string)
        vertices.append(vertex)
        
    return(vertices)
    
def get_label(simplex):
    
    label_characters = []
    
    for c in simplex:
        label_characters.append(c)
    
    vertices = get_vertex_chars(simplex)
    
    for c in vertices:
        label_characters.remove(c)
    
    label_characters.remove(" ")

    r_as_list = []
    for c in label_characters:
        r_as_list.append(c)
        if c == " ":
            break
        
    r_as_list.remove(" ")

    r_as_string = r_as_list[0]
    for i in range(1, len(r_as_list)):
        r_as_string += r_as_list[i]

    r = float(r_as_string)
    
    for c in r_as_list:
        label_characters.remove(c)
        
    label_characters.remove(" ")
    label_characters.remove("\n")
    
    k_as_string = label_characters[0]
    for i in range(1, len(label_characters)):
        k_as_string += label_characters[i]
        
    k = int(k_as_string)
    
    return((r,k))
    
def get_simplices_in_colex_order(bifiltration, p):
    
    p_simplices = []

    with open(bifiltration) as f:
        for line in f:
            if simplex_dimension(line) == p:
                p_simplices.append((get_vertices(line), get_label(line)))
                
    p_simplices.sort(key=lambda x : (x[1][1], x[1][0]))
    
    return(p_simplices)
    
# This uses scipy.sparse to construct sparse matrices
def MatrixRep(bifiltration, p):
    
    p_plus_one_simplices = get_simplices_in_colex_order(bifiltration, p+1)
    p_simplices = get_simplices_in_colex_order(bifiltration, p)
    
    n = len(p_simplices)
    m = len(p_plus_one_simplices)

    labels = []
    matrix = scipy.sparse.csc_matrix((n,m), dtype='int')

    for i in range(m):
        simplex = p_plus_one_simplices[i]
    
        label = simplex[1]
        labels.append(label)
    
        col = scipy.sparse.csc_matrix((n,1), dtype='int')
        vertices = simplex[0]
        boundary = []
    
        for j in range(p+2):
            boundary_simplex = [vertices[k] for k in range(p+2) if k != j]
            boundary.append(boundary_simplex)
        
        for sigma in boundary:
            rows = [index for index in range(n) if set(p_simplices[index][0]) == set(sigma)]
        
            if len(rows) != 1:
                print('multiple row indices!')
            else:
                row = rows[0]
                col[(row, 0)] = 1
    
        matrix[:, i] = col
    
    R = BiGradedMatrix(labels, matrix)  
    
    return(R)

def FIRep(bifiltration, p):
    
    if p == 0:
        R = MatrixRep(bifiltration, p)
        return(R)
        
    else:
        R = MatrixRep(bifiltration, p)
        Q = MatrixRep(bifiltration, p - 1)
        return(R, Q)