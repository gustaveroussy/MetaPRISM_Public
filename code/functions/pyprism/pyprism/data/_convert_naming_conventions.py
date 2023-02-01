# -*- coding: utf-8 -*-
"""
@created: Nov 25 2020
@modified: Nov 25 2020
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Naming conventions
    CapitalizedWords (CamelCase)
    Capitalized_Under_Scores (Camel_Snake_Case)
    mixedCase (camelCase)
    under_scores (snake_case)

User functions
    any_to_cw
    any_to_cu
    any_to_mc
    any_to_us

References
----------
https://www.oreilly.com/library/view/python-cookbook/0596001673/ch03s16.html
"""

import re

# User-level ===========================================================================================================

def any_to_cw(x):
    return _lw_to_cw(_any_to_lw(x))

def any_to_cu(x):
    return _lw_to_cu(_any_to_lw(x))

def any_to_mc(x):
    return _lw_to_mc(_any_to_lw(x))

def any_to_us(x):
    return _lw_to_us(_any_to_lw(x))

# Internal =============================================================================================================

def _upper2cw(x):
    return ''.join(map(str.capitalize, re.sub(r'(?<=[a-z])[A-Z]|(?<!^)[A-Z](?=[a-z])', r" \g<0>", x).split()))

def _lw_to_us(x):
    return '_'.join(x)

def _lw_to_cw(x):
    return ''.join(map(str.capitalize,x))

def _lw_to_cu(x):
    return '_'.join(map(str.capitalize,x))

def _lw_to_mc(x):
    x = list(x)
    return x[0]+''.join(map(str.capitalize,x[1:]))

def _any_to_lw(x):  # any format of identifier to list of lowercased words
    # First, see if there are dashes:
    lw = str.split(x, '-')
    if len(lw)>1: x = x.replace('-', '_')

    # Second, see if there are dots:
    lw = str.split(x, '.')
    if len(lw)>1: x = x.replace('.', '_')

    # See if there are underscores:
    lw = str.split(x, '_')
    if len(lw)>1: return map(str.lower, lw)

    # No. Then uppercase letters are the splitters:
    # Note: cases with whole words in upper cases like ABCYou are first
    # transformed into AbcYou.
    pieces = re.split('([A-Z])', _upper2cw(x))

    # Ensure first word follows the same rules as the others:
    if pieces[0]: pieces = [''] + pieces
    else: pieces = pieces[1:]

    # Join two by two, lowercasing the splitters as you go
    return [pieces[i].lower(  )+pieces[i+1] for i in range(0,len(pieces),2)]
