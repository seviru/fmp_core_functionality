#!/usr/bin/env python3
"""Script to handle the arguments that the user inputs and that
will launch different functionalities in the main script
"""
import argparse

parser = argparse.ArgumentParser(description="Handles the desired output information by the user and asks for the required input files.")
parser.add_argument("--tree_infile", type=str, required=True,
                    help="Tree file in Newick format.")
parser.add_argument("--alignment_infile", type=str, required=True,
                    help="Alignment corresponding to the tree file, with the tree leaf names as headers.")
parser.add_argument("--table_infile", type=str, required=True,
                    help="Table file with the required format.")
parser.add_argument("--uniprot_infile", type=str, required=True,
                    help="File containing the uniprot annotations for our data cluster.")
parser.add_argument("--min_evalue", type=float, default=1e-10,
                    help="Minimum evalue to take into account an uniprot hit.")
parser.add_argument("--annotation_feature", type=str, required=True,
                    help="Feature for which we want to represent our tree or get the node scores.")
parser.add_argument("--tree_outfile", type=str, default="tree_image",
                    help="File where we will save the image of our phylogenetic tree with the secuence in the specified annotation\
                        and the node scores.")
parser.add_argument("--node_score_table", type=str, default="node_score.tsv",
                    help="File where we will store the node scores for our tree.")

args = parser.parse_args()
## END