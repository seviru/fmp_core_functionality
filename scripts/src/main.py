#!/usr/bin/env python3
"""Core of the script, which gives us a node score dictionary and an output tree given different 
input values by the user, as well as different intermediate useful information.
"""

import src.scoring_functions as sf
import json
import re
from ete3 import PhyloTree, TreeStyle, SeqMotifFace, TextFace

def get_alignment_position (sequence_position, sequence):
    """Transforms a given sequence position into the real alignment
    position for our specific case.
    """
    alignment_position = 0
    aminoacid_counted = 0
    for aminoacid in sequence:
        if aminoacid.isalpha():
            aminoacid_counted += 1
        if not (aminoacid_counted <= sequence_position or aminoacid_counted == 0):
            break
        alignment_position += 1

    return alignment_position


def get_positions_matrix (feature_hash, tree):
    """Returns a matrix of all the feature positions in our alignment.
    """
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

    return position_matrix


def calculate_node_score (node, position_matrix, calculus_algorithm, ignore_gaps):
    """Calculates the node score for a given tree branch, with different
    algorithm variations. If user wants to ignore positions that are gaps,
    the only algorithm that does this at the moment (06/05/2020) is "simple_calculus"
    (calculus_algorithm = "simple").
    """
    branch_matrix = []
    for branch in node.get_children():
        aminoacid_matrix = []
        for position in position_matrix:
            position_aminoacids = []
            for leaf in branch.iter_leaves():
                if (leaf.sequence[position] == "-" and ignore_gaps == True and calculus_algorithm == "simple"):
                    continue
                else:
                    position_aminoacids.append(leaf.sequence[position])

            aminoacid_matrix.append(position_aminoacids)
        branch_matrix.append(aminoacid_matrix)
    
    first_branch_matrix = branch_matrix[0]
    second_branch_matrix = branch_matrix[1]
    
    scorer = sf.get_scorer(calculus_algorithm)
    score = scorer(first_branch_matrix, second_branch_matrix)

    return score



### MAIN 

tree_path = #IMPORTAR UN ARGUMENTO
alignment_path = #IMPORTAR OTRO ARGUMENTO
table_path = #IMPORTAR OTRO ARGUMENTO
uniprot_path = #IMPORTAR OTRO ARGUMENTO

try:
    with open(table_path, "r") as table_file:
        table_info = table_file.readlines()
    table_file.close()
except:
    print ("No table path.")
        
try:
    with open(alignment_path, "r") as alignment_file:
        alignment_info = alignment_file.read()
    alignment_file.close()
except:     
    print ("No alignment path.")

uniprot_info = {}
with open(uniprot_path, "r") as uniprot_file:
    for line in uniprot_file:
        uniprot_entry = json.loads(line)
        uniprot_info.update(uniprot_entry)
    uniprot_file.close


### AQUI ABAJO, "BINDING" Y EL VALOR DE EVALUE DEBERIAN SER ARGUMENTOS TRAIDOS
### IF FEATURE = ALL, PENSAR LO QUE QUEREMOS HACER
uniprot_hit_hash, leaf_delete_list = retrieve_features("BINDING", table_info, 1e-100)
tree = PhyloTree(tree_path, alignment=alignment_info, alg_format="fasta")
md = tree.get_midpoint_outgroup()
tree.set_outgroup(md)
position_matrix = get_positions_matrix(uniprot_hit_hash, tree)
ts = TreeStyle()
ts.layout_fn = lambda x: True

node_number = 0
node_scores = {}

for leaf in tree.iter_leaves():
    if leaf.name in leaf_delete_list:
        leaf.delete()
        
for node in tree.traverse():
    if node.is_leaf() == True:
        draw_position = 0
        for position in position_matrix:
            seqFace = SeqMotifFace(node.sequence[position], seq_format="seq")
            (tree&node.name).add_face(seqFace, draw_position, "aligned")
            draw_position += 1  
    else:
        node_score = round(calculate_node_score(node, position_matrix, "simple"), 2) ###AQUI EL ALGORITMO DEBERIA VENIR DE UN ARGUMENTO?
        node.add_feature("node_score", node_score)
        node_scores[node_number] = node_score
        node_number += 1
        if node_score > 0:
            score_face = TextFace(node.node_score)
            node.add_face(score_face, 0, "branch-top")
        
### AQUI ABAJO EL ARCHIVODE SALIDA DEBERIA VENIR DE LA LINEA DE ARGUMENTOS
tree.render("mytree_deleted.png", w=800, h=2000, units="px", tree_style=ts)




## END