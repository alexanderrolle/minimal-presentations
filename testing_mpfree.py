#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 15:26:20 2020

@author: arolle
"""

import math

def compare_bigrades(col_bigrades_1, row_bigrades_1, col_bigrades_2, row_bigrades_2, tol):
    
    if len(col_bigrades_1) != len(col_bigrades_2):
        return('number of columns is different!')
        
    if len(row_bigrades_1) != len(row_bigrades_2):
        return('number of rows is different!')

    col_comparisons = []
    
    for i in range(len(col_bigrades_1)):
        col_comparisons.append(math.isclose(col_bigrades_1[i][0], col_bigrades_2[i][0], rel_tol=tol))
        col_comparisons.append(math.isclose(col_bigrades_1[i][1], col_bigrades_2[i][1], rel_tol=tol))
        
    col_match = all(col_comparisons)
    
    row_comparisons = []
    
    for i in range(len(row_bigrades_1)):
        row_comparisons.append(math.isclose(row_bigrades_1[i][0], row_bigrades_2[i][0], rel_tol=tol))
        row_comparisons.append(math.isclose(row_bigrades_1[i][1], row_bigrades_2[i][1], rel_tol=tol))
        
    row_match = all(row_comparisons)
    
    match = all([col_match, row_match])
    
    return(match)

#########################################################################################
### reading mpfree's output #############################################################
#########################################################################################

def mpfree_grades(firep):
    
    with open(firep) as f:
        for num, line in enumerate(f):
            if num == 3:
                num_cols, num_rows = get_mpfree_dimensions(line)
    
    mpfree_col_bigrades = []
    mpfree_row_bigrades = []
    
    with open(firep) as f:
        for num, line in enumerate(f):
            if 3 < num and num < 4 + num_cols:
                mpfree_col_bigrades.append(get_mpfree_bigrade(line))
                
    with open(firep) as f:
        for num, line in enumerate(f):
            if 3 + num_cols < num and num < 4 + num_cols + num_rows:
                mpfree_row_bigrades.append(get_mpfree_bigrade(line))
                
    return(mpfree_col_bigrades, mpfree_row_bigrades)

def get_mpfree_bigrade(line):
    
    x_grade = ''
    y_grade = ''
    
    for i in range(len(line)):
        if line[i] == ' ':
            break
        else:
            x_grade += line[i]
            
    for i in range(len(x_grade) + 1, len(line)):
        if line[i] == ' ':
            break
        else:
            y_grade += line[i]
            
    x_grade = float(x_grade)
    y_grade = float(y_grade)
    
    return((x_grade, y_grade))           
                
                
def get_mpfree_dimensions(line):
    
    num_cols = ''
    num_rows = ''
    
    for i in range(len(line)):
        if line[i] == ' ':
            break
        else:
            num_cols += line[i]
            
    for i in range(len(num_cols) + 1, len(line)):
        if line[i] == ' ':
            break
        else:
            num_rows += line[i]
            
    num_cols = int(num_cols)
    num_rows = int(num_rows)
    
    return(num_cols, num_rows)

#########################################################################################
### reading rivet's output ##############################################################
#########################################################################################

def rivet_grades(rivet_output):
    
    # get bigrades from rivet (as integers)
    
    rivet_int_row_bigrades, rivet_int_col_bigrades = get_labels_from_rivet(rivet_output)
    
    # put them in lexicographical order
    
    rivet_int_row_bigrades.sort()
    rivet_int_col_bigrades.sort()
    
    # convert integer grades to real grades
    
    with open(rivet_output) as f:
        for num, line in enumerate(f):
            if 'y-grades' in line:
                y_start = num
                
    rivet_row_bigrades = []
    rivet_col_bigrades = []
    
    for bigrade in rivet_int_row_bigrades:
        rivet_row_bigrades.append(bigrades_real_value(rivet_output, bigrade[0], bigrade[1], y_start))
    
    for bigrade in rivet_int_col_bigrades:
        rivet_col_bigrades.append(bigrades_real_value(rivet_output, bigrade[0], bigrade[1], y_start))
    
    return(rivet_col_bigrades, rivet_row_bigrades)

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

def get_labels_from_rivet(rivet_output):
                
    with open(rivet_output) as f:
        for num, line in enumerate(f):
            if 'Row bigrades' in line:
                row_bigrades_num = num + 1
                
    with open(rivet_output) as f:
        for num, line in enumerate(f):
            if num == row_bigrades_num:
                row_bigrades = line
                
    row_bigrades_list = get_bigrades(row_bigrades)
                
    with open(rivet_output) as f:
        for num, line in enumerate(f):
            if 'Column bigrades' in line:
                col_bigrades_num = num + 1
                
    with open(rivet_output) as f:
        for num, line in enumerate(f):
            if num == col_bigrades_num:
                col_bigrades = line
                
    col_bigrades_list = get_bigrades(col_bigrades)
    
    return(row_bigrades_list, col_bigrades_list)
    
def frac_to_float(line):
    
    negative = False
    
    if line[0] == '-':
        negative = True
        line = line.replace('-', '')
    
    numerator = ''
    denominator = ''
    
    for i in range(len(line)):
        if line[i] == '/':
            break
        else:
            numerator += line[i]
            
    denominator_exists = False
    
    for i in range(len(line)):
        if line[i] == '/':
            denominator_exists = True
            break
        
    if denominator_exists:
    
        denominator_start = False        
                
        for i in range(len(line)-1):
            if line[i] == '/':
                denominator_start = True
                continue
            if denominator_start:
                denominator += line[i]
                
    else:
        
        denominator = '1'
            
    numerator = float(numerator)
    denominator = float(denominator)
    value = numerator / denominator
    
    if negative:
        
        value = -1 * value
            
    return(value)
    
def bigrades_real_value(rivet_output, x_grade, y_grade, y_start):
    
    with open(rivet_output) as f:
        for num, line in enumerate(f):
            if num == 1 + x_grade:
                x_value = frac_to_float(line)
                
    with open(rivet_output) as f:
        for num, line in enumerate(f):
            if num == 1 + y_start + y_grade:
                y_value = frac_to_float(line)
                
    return((x_value, y_value))