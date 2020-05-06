#!/usr/bin/env python3
"""Functions to retrieve information from our data cluster.
"""

def check_features (table_file, evalue_threshold, uniprot_info):
    """Returns a set of the present annotated features for a given table
    and within the evalue threshold.
    """
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
            hit_evalue <= evalue_threshold):
            for feature in uniprot_info[hit_name]["FT"]: # TO CHECK FEATURES WHEN WE PICK A NEW FILE MANUALLY
                if feature["ft"] not in feature_set:
                    print(feature["ft"])
                    feature_set.add(feature["ft"])
                    print(feature_set)
            uniprot_hit_set.update(hit_name)
        
    return feature_set


def retrieve_features (feature_tag, table_file, evalue_threshold, uniprot_info):
    """Returns a dictionary containing all the specified feature information, and
    a deletion list for all the annotated leaves which don't pass the cut of the evalue_threshold
    """
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
                if feature["ft"] == feature_tag:
                    features_newlist.append(feature)
            if len(features_newlist) > 0:
                uniprot_hit_hash[hit_name] = features_newlist
            leaf_saving_list.append(hit_name)
        elif hit_evalue > evalue_threshold:
            leaf_deleting_list.append(hit_name)
            
    leaf_deleting_list = set(leaf_deleting_list) - set(leaf_saving_list)

    return uniprot_hit_hash, leaf_deleting_list

    
## END