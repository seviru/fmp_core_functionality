#!/usr/bin/env python3
"""Functions to retrieve information from our data cluster.
"""

import src.scoring_functions as sf
import sys

def check_features (table_file, evalue_threshold, uniprot_info):
    """Returns a set of the present annotated features for a given table
    and within the evalue threshold.
    """
    try:
        uniprot_hit_set = {}
        feature_set = {}
        uniprot_hit_set = set(uniprot_hit_set)
        feature_set = set(feature_set)
        
        for line in table_file:
            hit_type = line.split("\t")[2]
            hit_name = line.split("\t")[3]
            hit_evalue = line.split("\t")[4]
            if hit_evalue == "-":
                hit_evalue = 0.0
            else:
                hit_evalue = float(hit_evalue)
            if (hit_type == ("spb" or "spe") and
                hit_name not in uniprot_hit_set and
                (evalue_threshold is None or hit_evalue <= evalue_threshold)):
                for feature in uniprot_info[hit_name]["FT"]: # TO CHECK FEATURES WHEN WE PICK A NEW FILE MANUALLY
                    if feature["ft"] not in feature_set:
                        feature_set.add(feature["ft"])
                uniprot_hit_set.update(hit_name)
    except:
        sys.stderr.write("Error at checking data features (feature_processing.check_features).\n")
        sys.exit(1)
        
    return feature_set


def retrieve_features (feature_list, table_file, evalue_threshold, uniprot_info):
    """Returns a dictionary containing all the specified feature information, and
    a deletion list for all the annotated leaves which don't pass the cut of the evalue_threshold
    """
    try:
        uniprot_hit_hash = {}
        leaf_deleting_list = []
        leaf_saving_list = []

        for line in table_file:
            hit_type = line.split("\t")[2]
            hit_name = line.split("\t")[3]
            hit_evalue = line.split("\t")[4]
            if hit_evalue == "-":
                hit_evalue = 0.0
            else:
                hit_evalue = float(hit_evalue)

            if (hit_type == ("spb" or "spe") and 
                hit_name not in uniprot_hit_hash and 
                hit_evalue <= evalue_threshold):
                features_newlist = []
                for feature in uniprot_info[hit_name]["FT"]:
                    for feature_tag in feature_list:
                        if feature["ft"] == feature_tag:
                            features_newlist.append(feature)
                if len(features_newlist) > 0:
                    uniprot_hit_hash[hit_name] = features_newlist
                leaf_saving_list.append(hit_name)
            elif hit_evalue > evalue_threshold:
                leaf_deleting_list.append(hit_name)         
        leaf_deleting_list = set(leaf_deleting_list) - set(leaf_saving_list)
    except:
        sys.stderr.write("Error at retrieving data features (feature_processing.retrieve_features).\n")
        sys.exit(1)

    return uniprot_hit_hash, leaf_deleting_list


def get_alignment_position (sequence_position, sequence):
    """Transforms a given sequence position into the real alignment
    position for our specific case.
    """
    try:
        alignment_position = 0
        aminoacid_counted = 0
        for aminoacid in sequence:
            if aminoacid.isalpha():
                aminoacid_counted += 1
            if not (aminoacid_counted <= sequence_position or aminoacid_counted == 0):
                break
            alignment_position += 1
    except:
        sys.stderr.write("Error at retrieving alignment position (feature_processing.get_alignment_position).\n")
        sys.exit(1)

    return alignment_position


def get_positions_matrix (feature_hash, tree):
    """Returns a matrix of all the feature positions in our alignment.
    """
    try:
        position_matrix = []
        for unigene in feature_hash:
            unigene_sequence = (tree&unigene).sequence
            for feature in feature_hash[unigene]:
                feature_start = int(feature["s"])-1
                feature_end = int(feature["e"])-1
                alignment_feature_start = get_alignment_position(feature_start, unigene_sequence)
                alignment_feature_end = get_alignment_position(feature_end, unigene_sequence)
                for position in range (alignment_feature_start, alignment_feature_end+1):
                    position_matrix.append(position)
        position_matrix = sorted(list(set(position_matrix)))
    except:
        sys.stderr.write("Error at retrieving position matrix (feature_processing.get_position_matrix).\n")
        sys.exit(1)

    return position_matrix


def calculate_node_score (node, position_matrix, calculus_algorithm, differentiate_gaps):
    """Calculates the node score for a given tree branch, with different
    algorithm variations. If user wants to ignore positions that are gaps,
    the only algorithm that does this at the moment (06/05/2020) is "simple_calculus"
    (calculus_algorithm = "simple").
    """
    try:
        branch_matrix = []
        for branch in node.get_children():
            aminoacid_matrix = []
            for position in position_matrix:
                position_aminoacids = []
                for leaf in branch.iter_leaves():
                    if (leaf.sequence[position] == "-" and differentiate_gaps == "Y"):
                        continue
                    else:
                        position_aminoacids.append(leaf.sequence[position])

                aminoacid_matrix.append(position_aminoacids)
            branch_matrix.append(aminoacid_matrix)
        first_branch_matrix = branch_matrix[0]
        second_branch_matrix = branch_matrix[1]
        scorer = sf.get_scorer(calculus_algorithm)
        score = scorer(first_branch_matrix, second_branch_matrix)
    except:
        sys.stderr.write("Error at calculating node score (feature_processing.calculate_node_score).\n")
        sys.exit(1)

    return score