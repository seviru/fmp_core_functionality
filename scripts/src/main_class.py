### TRYING TO PUT EVERYTHING INSIDE A CLASS, STILL WORK IN PROGRESS ###
"""Class which will allow us to create FeatureStudy objects, that we will be able
to operate via different methods in order to analyze them.
"""

import json
import sys
import src.feature_processing as fp
from ete3 import PhyloTree, TreeStyle, SeqMotifFace, TextFace

class FeatureStudy:
    def __init__(self, tree_path, alignment_path, table_path,
                uniprot_path, annotation_features, min_evalue, 
                tree_outfile, node_score_table, node_score_algorithm, 
                ignore_gap_positions):
        try:
            self.tree_in = tree_path
            self.min_eval = min_evalue
            self.tree_out = tree_outfile
            self.score_table = node_score_table
            self.calc_alg = node_score_algorithm
            self.ignore_gaps = ignore_gap_positions
            self.align_in, self.table_in, self.uniprot_in, self.study_features = self.__setup__(alignment_path, 
                                                                                table_path, 
                                                                                uniprot_path, 
                                                                                annotation_features)
        except:
            print("Error at instance initialization.")
            sys.stderr.write("Error at instance initialization.")
            sys.exit(1)


    def __setup__(self, alignment_path, table_path,
                uniprot_path, annotation_features):
        """Method to check if the infiles exist, to process the features 
        the user want to study, and if the parameter "ignore_gaps" 
        corresponds with an algorithm that can actually handle it.
        """
        try:
            try:
                with open(alignment_path, "r") as alignment_file:
                    alignment_info = alignment_file.read()
                alignment_file.close()
            except:
                print("Alignment file not found.")
                sys.stderr.write("Alignment file not found.")
                sys.exit(1)
            
            try:
                with open(table_path, "r") as table_file:
                    table_info = table_file.readlines()
                table_file.close()
            except:
                print("Table file not found.")
                sys.stderr.write("Table file not found.")
                sys.exit(1)

            try:
                uniprot_info = {}
                with open(uniprot_path, "r") as uniprot_file:
                    for line in uniprot_file:
                        uniprot_entry = json.loads(line)
                        uniprot_info.update(uniprot_entry)
                uniprot_file.close
            except:
                print("Uniprot file not found.")
                sys.stderr.write("Uniprot file not found.")
                sys.exit(1)

            try:
                if annotation_features == "ALL":
                    annotation_features = fp.check_features(table_info, self.min_eval, uniprot_info)
                    if "CHAIN" in annotation_features:
                        annotation_features.remove("CHAIN")
                elif "," in annotation_features:
                    annotation_features = set([feature.upper() for feature in annotation_features.split(",")])
                else:
                    pass
            except:
                print("Feature unpacking gone wrong.")
                sys.stderr.write("Feature unpacking gone wrong.")
                sys.exit(1)

            if self.calc_alg != "simple" and self.ignore_gaps == "Y":
                print("Only calculus algorithm supporting gap differentiation is 'simple'.")
                sys.stderr.write("Only calculus algorithm supporting gap differentiation is 'simple'.")
                sys.exit(1)
        
        except:
            print("Error at instance setup.")
            sys.stderr.write("Error at instance setup.")
            sys.exit(1)
                
        return alignment_info, table_info, uniprot_info, annotation_features


    def calculate_nodes(self):
        """Method to calculate the different internal node scores
        for a given calculus method, and store those values both in
        a dictionary (if the user wants to) and in an instance
        of a processed tree.
        """
        try:
            uniprot_hit_hash, leaf_deleting_list = fp.retrieve_features(self.study_features, self.table_in, self.min_eval, self.uniprot_in)
            tree = PhyloTree(self.tree_in, alignment=self.align_in, alg_format="fasta")
            md = tree.get_midpoint_outgroup()
            tree.set_outgroup(md)
            position_matrix = fp.get_positions_matrix(uniprot_hit_hash, tree)
            node_number = 0
            node_scores = {}
            for leaf in tree.iter_leaves():
                if leaf.name in leaf_deleting_list:
                    leaf.delete()
            for node in tree.traverse():
                if node.is_leaf() == False:
                    node_score = round(fp.calculate_node_score(node, position_matrix, self.calc_alg, self.ignore_gaps), 2) ###AQUI EL ALGORITMO DEBERIA VENIR DE UN ARGUMENTO?
                    node.add_feature("node_score", node_score)
                    node_scores[node_number] = node_score
                    node_number += 1
            self.processed_tree = tree
            self.processed_tree__position_matrix = position_matrix
        except:
            print("Error at calculating nodes.")
            sys.stderr.write("Error at calculating nodes.")
            sys.exit(1)

        return node_scores


    def design_tree(self, node_scores="CALCULATE", plot_threshold=0):
        """Method that allows us to design a tree with It's node scores.
        """
        try:
            if node_scores == "CALCULATE":
                self.calculate_nodes()
            else: # DUDAS EN ESTA OPCION
                tree = PhyloTree(self.tree_in, alignment=self.align_in, alg_format="fasta")
                node_number = 0
                for node in tree.traverse():
                    node.add_feature("node_number", node_number)
                    node.add_feature("node_score", node_scores[node_number])
                    node_number += 1
                self.processed_tree = tree
            
            for node in self.processed_tree.traverse():
                if node.is_leaf() == True:
                    draw_position = 0
                    for position in self.processed_tree__position_matrix:
                        seqFace = SeqMotifFace(node.sequence[position], seq_format="seq")
                        (self.processed_tree&node.name).add_face(seqFace, draw_position, "aligned")
                        draw_position += 1  
                else:
                    if node.node_score > plot_threshold:
                        score_face = TextFace(node.node_score)
                        node.add_face(score_face, 0, "branch-top")
        except:
            print("Error at designing tree.")
            sys.stderr.write("Error at designing tree.")
            sys.exit(1)

        return


    def plot_tree(self, tree="DESIGN", width=800, heigth=2000, units="px", plot_threshold=0): # Units may be “px”: pixels, “mm”: millimeters, “in”: inches
        """Method to plot a tree to a file given a set of
        parameters.
        """
        try:
            ts = TreeStyle()
            ts.layout_fn = lambda x: True
            if tree == "DESIGN":
                self.design_tree(plot_threshold=plot_threshold)
                tree = self.processed_tree
            tree.render(self.tree_out, w=width, h=heigth, units=units, tree_style=ts)
        except:
            print("Error at plotting tree.")
            sys.stderr.write("Error at plotting tree.")
            sys.exit(1)

## END