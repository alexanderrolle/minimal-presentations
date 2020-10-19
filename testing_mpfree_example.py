#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 13:15:07 2020

@author: arolle
"""

from testing_mpfree import *

data_name = 'function_rips_sample'

# mpfree output with default options
mpfree_output = data_name + '_output.firep'

# lists of column and row bigrades from mpfree with default options
mpfree_col_bigrades, mpfree_row_bigrades = mpfree_grades(mpfree_output)    


# compare with output from other options
# 0 = --no-chunk
# 1 = --clearing
# 2 = --bit_tree_pivot_column
for option in ['0', '1', '2', '01', '02', '12', '012']:
    mpfree_output_with_option = data_name + '_output_' + option + '.firep'
    
    # check if output files are identical
    if open(mpfree_output,'r').read() == open(mpfree_output_with_option,'r').read():
        print('output matches with option ' + option)
    
    # check if lists of bigrades match
    else:
        print('output does not match with option ' + option)
        mpfree_with_option_col_bigrades, mpfree_with_option_row_bigrades = mpfree_grades(mpfree_output_with_option)
        match = compare_bigrades(mpfree_col_bigrades, mpfree_row_bigrades, 
                                 mpfree_with_option_col_bigrades, mpfree_with_option_row_bigrades, 1e-9)
        print('lists of bigrades match option ' + option +': ' + str(match))


# compare lists of bigrades with output from rivet
rivet_output = 'rivet_' + data_name + '_output.txt'
rivet_col_bigrades, rivet_row_bigrades = rivet_grades(rivet_output)
match = compare_bigrades(mpfree_col_bigrades, mpfree_row_bigrades, 
                       rivet_col_bigrades, rivet_row_bigrades, 1e-9)
print('lists of bigrades match rivet: ' + str(match))