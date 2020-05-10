#!/usr/bin/env python3
"""Script to handle the arguments that the user inputs and that
will launch different functionalities in the main script
"""
import argparse
import src.main_class as main
import sys

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
outfile_args.add_argument("-o", "--tree_outfile", type=str,
                    help="File where we will save the image of our phylogenetic tree with the secuence in the specified annotation\
                          and the node scores.")
outfile_args.add_argument("-n", "--node_score_table", type=str,
                    help="File where we will store the node scores for our tree.")
parameter_args.add_argument("-c", "--calculus_algorithm", type=str, default="simple", 
                    choices=["simple", "all_vs_all", "whole_annotation_simple", "whole_annotation_all_vs_all", "all_vs_all_means", "whole_annotation_all_vs_all_means"],
                    help="Select which calculation algorithm you want to utilize to calculate the node scores in your tree.")
parameter_args.add_argument("-g", "--ignore_gap_positions", type=str, default="N", choices=["Y", "N"],
                    help="Decide if you want your calculus algorithm to ignore positions with gaps. At the moment, only 'simple'\
                          calculus algorithm is capable of doing this.")
parameter_args.add_argument("-e", "--min_evalue", type=float, default=1e-10,
                    help="Minimum evalue to take into account an uniprot hit.")
parameter_args.add_argument("-f", "--annotation_feature", type=str, default="ALL",
                    help="Feature for which we want to represent our tree or get the node scores.")

args = parser.parse_args()

if args.tree_outfile is None and args.node_score_table is None:
      sys.stderr.write("Not specified any outfile. Specify at least tree_outfile or node_score_table.\n")
      sys.exit(1)

print("Creating study instance...")
case_study = main.FeatureStudy(args.tree_infile, args.alignment_infile, args.table_infile,
                           args.uniprot_infile, args.annotation_feature, args.min_evalue,
                           args.calculus_algorithm, args.ignore_gap_positions)

print("Designing tree.")
case_study.design_tree()

if args.tree_outfile is not None:
      print("Plotting tree.")
      case_study.plot_tree(outfile=f"{args.tree_outfile}.png", tree=case_study.processed_tree)
if args.node_score_table is not None:
      print("Writing nodes to a file.")
      case_study.write_node_file(outfile=args.node_score_table)

print("Study finished.")

## END