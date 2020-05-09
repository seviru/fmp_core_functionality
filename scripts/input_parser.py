#!/usr/bin/env python3
"""Script to handle the arguments that the user inputs and that
will launch different functionalities in the main script
"""
import argparse
import src.main_class as main

parser = argparse.ArgumentParser(description="Handles the desired output information by the user and asks for the required input files.")
infile_args = parser.add_argument_group("INFILES", "Files of input data for our script to process.")
outfile_args = parser.add_argument_group("OUTFILES", "Files where our script will write the output.")
parameter_args = parser.add_argument_group("SCRIPT PARAMETERS", "Different parameters which will modify how the script works.")


infile_args.add_argument("-i", "--tree_infile", type=str, required=True,
                    help="Tree file in Newick format.")
infile_args.add_argument("-a", "--alignment_infile", type=str, required=True,
                    help="Alignment corresponding to the tree file, with the tree leaf names as headers.")
infile_args.add_argument("-t", "--table_infile", type=str, required=True,
                    help="Table file with the required format.")
infile_args.add_argument("-u", "--uniprot_infile", type=str, required=True,
                    help="File containing the uniprot annotations for our data cluster.")
outfile_args.add_argument("-o", "--tree_outfile", type=str, default="tree_image.png",
                    help="File where we will save the image of our phylogenetic tree with the secuence in the specified annotation\
                          and the node scores.")
outfile_args.add_argument("-n", "--node_score_table", type=str, default="node_score.tsv",
                    help="File where we will store the node scores for our tree.")
parameter_args.add_argument("-c", "--calculus_algorithm", type=str, default="simple", 
                    choices=["simple", "all_vs_all", "whole_annotation_simple", "whole_annotation_all_vs_all", "all_vs_all_means", "whole_annotation_all_vs_all_means"]
                    help="Select which calculation algorithm you want to utilize to calculate the node scores in your tree.")
parameter_args.add_argument("-g", "--ignore_gap_positions", type=str, default="N", choices=["Y", "N"],
                    help="Decide if you want your calculus algorithm to ignore positions with gaps. At the moment, only 'simple'\
                          calculus algorithm is capable of doing this.")
parameter_args.add_argument("-e", "--min_evalue", type=float, default=1e-10,
                    help="Minimum evalue to take into account an uniprot hit.")
parameter_args.add_argument("-f", "--annotation_feature", type=str, default="ALL",
                    help="Feature for which we want to represent our tree or get the node scores.")



args = parser.parse_args()

###AQUI CREO LA INSTANCIA DE LA CLASE Y LLAMO A SUS METODOS
case_study = main.FeatureStudy(args.tree_infile, args.alignment_infile, args.table_infile,
                           args.uniprot_infile, args.annotation_feature, args.min_evalue,
                           args.tree_outfile, args.node_score_table, args.calculus_algorithm,
                           args.ignore_gap_positions)

# print(case_study.study_features)
# case_study.calculate_nodes()

# case_study.design_tree()

# case_study.processed_tree.render(case_study.tree_out, w=800, h=2000, units="px", tree_style=ts)
# ts = TreeStyle()
# ts.layout_fn = lambda x: True

case_study.plot_tree()

## END