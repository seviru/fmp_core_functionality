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
import src.feature_processing as fp

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
            leaf_number = 0
            for aminoacid in (first_branch_matrix[position] + second_branch_matrix[position]):
                if aminoacid != "":
                    leaf_number += 1
            differences = utils.dict_diff(Counter(item), Counter(second_branch_matrix[position]))
            if "" in differences:
                del differences[""]
            if (all(aminoacid == "" for aminoacid in first_branch_matrix[position]) or 
                all(aminoacid == "" for aminoacid in second_branch_matrix[position])) == True:  # If one of the compared positions is always a gap
                continue
            score += sum( [differences[different_aminoacid] / leaf_number for different_aminoacid in differences ] )
            print(f"score:::{score}")
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
        leaf_number = len(first_branch_matrix[0]) + len(second_branch_matrix[0])
        score = 0
        first_branch_merged = fp.merge_annotation(first_branch_matrix)
        second_branch_merged = fp.merge_annotation(second_branch_matrix)
        merged_matrix = [first_branch_merged] + [second_branch_merged]
        differences = utils.dict_diff(Counter(merged_matrix[0]), Counter(merged_matrix[1]))
        score += sum( [differences[different_aminoacid] / leaf_number for different_aminoacid in differences ] )
        score = round(score, 2)
    except:
        sys.stderr.write("Error at execution of calculus function (scoring_functions.whole_annotation.simple_calculus).\n")
        sys.exit(1)

    return score


def all_vs_all_calculus (first_branch_matrix, second_branch_matrix):
    """Compares a given position with all the position inside and outside Its
    branch. Sensitive to node size. score>0 (0 being equal nodes).
    """
    try:
        aminoacid_matrix = []
        score = 0
        for position in range(len(first_branch_matrix)):
            if (len(first_branch_matrix[position]) == 0 or len(second_branch_matrix[position]) == 0): # Don't need it right now, but saving it just in case
                continue
            else:
                aminoacid_matrix.append(first_branch_matrix[position] + second_branch_matrix[position])
        for position_aminoacids in aminoacid_matrix:
            leaf_number = len(position_aminoacids)
            for aa1, aa2 in itertools.combinations(position_aminoacids, 2):
                if aa1 != aa2:
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
        first_branch_merged = fp.merge_annotation(first_branch_matrix)
        second_branch_merged = fp.merge_annotation(second_branch_matrix)
        merged_matrix = first_branch_merged + second_branch_merged
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
            if (len(first_branch_matrix[position]) == 0 or len(second_branch_matrix[position]) == 0):
                continue
            else:
                position_aminoacids = first_branch_matrix[position] + second_branch_matrix[position]        
                position_comparison = []
                for aa1, aa2 in itertools.combinations(position_aminoacids, 2):
                    if aa1 == aa2:
                        position_comparison.append(0)
                    else:
                        position_comparison.append(1)
                position_means.append(sum(position_comparison)/len(position_comparison))
        if len(position_means) != 0:
            score = sum(position_means)/len(position_means)
        else:
            score = 0
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
        first_branch_merged = fp.merge_annotation(first_branch_matrix)
        second_branch_merged = fp.merge_annotation(second_branch_matrix)
        merged_matrix = first_branch_merged + second_branch_merged
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