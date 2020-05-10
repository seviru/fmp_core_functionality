#!/usr/bin/env python3
"""Different functions to calculate the node scores for a given
tree. Score would measure how similar the sequences in a node are,
given different criteria.
"""

from collections import Counter
from math import factorial as fact
import itertools
import src.utils as utils
import sys

def get_scorer(calculus_algorithm):
    """ Returns the desired calculus function for node score calculus
    given the calculus algorithm
    """
    try:
        if calculus_algorithm == "simple":
            scorer = simple_calculus
        
        elif calculus_algorithm == "all_vs_all":
            scorer = all_vs_all_calculus
        
        elif calculus_algorithm == "whole_annotation_simple":
            scorer = whole_annotation_simple_calculus
        
        elif calculus_algorithm == "whole_annotation_all_vs_all":
            scorer = whole_annotation_all_vs_all_calculus

        elif calculus_algorithm == "all_vs_all_means":
            scorer = all_vs_all_calculus_means
            
        elif calculus_algorithm == "whole_annotation_all_vs_all_means":
            scorer = whole_annotation_all_vs_all_calculus_means
    except:
        sys.stderr.write("Error at choosing score (scoring_functions.get_scorer).\n")
        sys.exit(1)

    return scorer


def simple_calculus (first_branch_matrix, second_branch_matrix):
    """Compares the 2 branches between them, taking into account only
    differences between branches multiple times. score>=0 (0 being equal nodes).
    """
    try:
        score = 0
        for position, item in enumerate(first_branch_matrix):
            leaf_number = len(first_branch_matrix[position]) + len(second_branch_matrix[position])
            differences = utils.dict_diff(Counter(item), Counter(second_branch_matrix[position]))
            if (len(first_branch_matrix[position]) == 0 or len(second_branch_matrix[position]) == 0):   # If one of the compared positions is always a gap
                continue                                                                                # we skip taking into account this leaf for score calcule.
            score += sum( [differences[different_aminoacid] / leaf_number for different_aminoacid in differences ] )
        score = round(score, 2)
    except:
        sys.stderr.write("Error at execution of calculus function (scoring_functions.simple_calculus).\n")
        sys.exit(1)

    return score


def whole_annotation_simple_calculus (first_branch_matrix, second_branch_matrix):
    """Compares the 2 branches between them, taking into account only
    differences between branches multiple times, taking all the positions
    from the annotation as an only one. score>=0 (0 being equal nodes).
    """
    try:
        merged_matrix = []
        merged_branch = []
        leaf_number = len(first_branch_matrix[0]) + len(second_branch_matrix[0])
        score = 0
        for position, item in enumerate(first_branch_matrix[0]):
            merged_matrix.append(item + first_branch_matrix[1][position])
        for position, item in enumerate(second_branch_matrix[0]):
            merged_branch.append(item + second_branch_matrix[1][position])
        merged_matrix = [merged_matrix] + [merged_branch]

        differences = utils.dict_diff(Counter(merged_matrix[0]), Counter(merged_matrix[1]))
        score += sum( [differences[different_aminoacid] / leaf_number for different_aminoacid in differences ] )
        score = round(score, 2)
    except:
        sys.stderr.write("Error at execution of calculus function (scoring_functions.whoile_annotation.simple_calculus).\n")
        sys.exit(1)

    return score


def all_vs_all_calculus (first_branch_matrix, second_branch_matrix):
    """Compares a given position with all the position inside and outside Its
    branch. Sensitive to node size. score=[0-1] (0 being equal nodes).
    """
    try:
        aminoacid_matrix = []
        score = 0
        leaf_number = len(first_branch_matrix[0]) + len(second_branch_matrix[0])
        for position in range(len(first_branch_matrix)):
            aminoacid_matrix.append(first_branch_matrix[position] + second_branch_matrix[position])
        for position_aminoacids in aminoacid_matrix:
            for aa1, aa2 in itertools.combinations(position_aminoacids, 2):
                if aa1 != aa2:
                    # score += 1/leaf_number
                    score += 1/(fact(leaf_number)/(fact(2)*fact(leaf_number-2)))
        score = round(score, 2)
    except:
        sys.stderr.write("Error at execution of calculus function (scoring_functions.all_vs_all_calculus).\n")
        sys.exit(1)

    return score


def whole_annotation_all_vs_all_calculus (first_branch_matrix, second_branch_matrix):
    """Compares a given position with all the position inside and outside Its
    branch, taking all the positions from the annotation as an only one. 
    Sensitive to node size. score=[0-1] (0 being equal nodes).
    """
    try:
        merged_matrix = []
        leaf_number = len(first_branch_matrix[0]) + len(second_branch_matrix[0])
        score = 0
        for position, item in enumerate(first_branch_matrix[0]):
            merged_matrix.append(item + first_branch_matrix[1][position])
        for position, item in enumerate(second_branch_matrix[0]):
            merged_matrix.append(item + second_branch_matrix[1][position])
        for aa1, aa2 in itertools.combinations(merged_matrix, 2):
            if aa1 != aa2:
                score += 1/(fact(leaf_number)/(fact(2)*fact(leaf_number-2))) 
        score = round(score, 2)
    except:
        sys.stderr.write("Error at execution of calculus function (scoring_functions.whole_annotation_all_vs_all_calculus).\n")
        sys.exit(1)
    
    return score


def all_vs_all_calculus_means (first_branch_matrix, second_branch_matrix):
    """Compares a given position with all the position inside and outside Its
    branch, giving a mean value for diversity. Sensitive to node size. score=[0-1] (0 being equal nodes).
    """
    try:
        position_means = []
        for position in range(len(first_branch_matrix)):
            position_aminoacids = first_branch_matrix[position] + second_branch_matrix[position]        
            position_comparison = []
            for aa1, aa2 in itertools.combinations(position_aminoacids, 2):
                if aa1 == aa2:
                    position_comparison.append(0)
                else:
                    position_comparison.append(1)
            position_means.append(sum(position_comparison)/len(position_comparison))
        score = sum(position_means)/len(position_means)
    except:
        sys.stderr.write("Error at execution of calculus function (scoring_functions.all_vs_all_calculus_means).\n")
        sys.exit(1)

    return score


def whole_annotation_all_vs_all_calculus_means (first_branch_matrix, second_branch_matrix):
    """Compares a given position with all the position inside and outside Its
    branch, taking all the positions from the annotation as an only one, and 
    giving a mean value for diversity. Sensitive to node size. score=[0-1] (0 being equal nodes).
    """
    try:
        merged_matrix = []
        score = 0
        for position, item in enumerate(first_branch_matrix[0]):
            merged_matrix.append(item + first_branch_matrix[1][position])
        for position, item in enumerate(second_branch_matrix[0]):
            merged_matrix.append(item + second_branch_matrix[1][position])
        position_comparison = [] 
        for aa1, aa2 in itertools.combinations(merged_matrix, 2):
            if aa1 == aa2:
                position_comparison.append(0)
            else:
                position_comparison.append(1)
        score = sum(position_comparison)/len(position_comparison)
        score = round(score, 2)
    except:
        sys.stderr.write("Error at execution of calculus function (scoring_functions.whole_annotation_all_vs_all_calculus_means).\n")
        sys.exit(1)

    return score


## END