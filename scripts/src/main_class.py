### TRYING TO PUT EVERYTHING INSIDE A CLASS, STILL WORK IN PROGRESS ###
"""Class which will allow us to create FeatureStudy objects, that we will be able
to operate via different methods in order to analyze them.
"""

import json
import sys
import os.path
import src.feature_processing as fp
from src import config
from ete3 import PhyloTree, TreeStyle, SeqMotifFace, TextFace


class FeatureStudy:
    def __init__(self, tree_path, alignment_path, 
                node_score_algorithm, differentiate_gap_positions,
                table_path=None, uniprot_path=None, 
                annotation_features=None, min_evalue=None,
                position_matrix=None):
        try:
            self.align_in = alignment_path
            self.tree_in = tree_path
            self.calc_alg = node_score_algorithm
            self.differentiate_gaps = differentiate_gap_positions
            self.position_matrix = position_matrix
            self.__setup_basic__
            if self.position_matrix == None:
                self.min_eval = min_evalue
                self.table_info, self.uniprot_info, self.study_features, self.all_features = self.__setup_features__(table_path, 
                                                                                                            uniprot_path, 
                                                                                                            annotation_features)
        except:
            sys.stderr.write("Error at instance initialization.\n")
            sys.exit(1)


# SEPARAR SETUPS

    def __setup_basic__(self):
        """Method to check if the different files exist, as well
        as the calculus algorithm.
        """
        try:
            try:
                os.path.isfile(self.align_in)
            except:
                print("Alignment file not found.")
                sys.stderr.write("Alignment file not found.")
                sys.exit(1)
            
            try:
                os.path.isfile(self.tree_in)
            except:
                print("Tree file not found.")
                sys.stderr.write("Tree file not found.")
                sys.exit(1)

            if self.calc_alg not in {"simple", "all_vs_all", "all_vs_all_means"} and self.differentiate_gaps == "Y":
                sys.stderr.write("Only calculus algorithm supporting gap differentiation are 'simple', 'all_vs_all', 'all_vs_all_means'.\n")
                sys.exit(1)
        except:
            sys.stderr.write("Error at basic setup.")
            sys.exit(1)
        return

    def __setup_features__(self, table_path,
                uniprot_path, annotation_features):
        """Method to process the features 
        the user want to study.
        """
        try:
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
                uniprot_hit_set = set()
                for line in table_info:
                    hit_name = line.split("\t")[3]
                    if hit_name not in uniprot_hit_set:
                        uniprot_hit_set.update([hit_name])               
                with open(uniprot_path, "r") as uniprot_file:
                    for line in uniprot_file:
                        uniprot_entry = json.loads(line)
                        for unigene in uniprot_entry.keys():
                            if unigene in uniprot_hit_set:
                                uniprot_info.update(uniprot_entry)
                uniprot_file.close
            except:
                print("Uniprot file not found.")
                sys.stderr.write("Uniprot file not found.")
                sys.exit(1)

            try:
                feature_collection = fp.check_features(table_info, None, uniprot_info)
                all_features = fp.check_features(table_info, self.min_eval, uniprot_info)
                temporal_features = set()
                if "CHAIN" in all_features:
                    all_features.remove("CHAIN")
                if annotation_features == "ALL":
                    annotation_features = all_features
                else:
                    if "," in annotation_features:
                        annotation_features = set([feature.upper() for feature in annotation_features.split(",")])
                    else:
                        annotation_features = set([annotation_features.upper()])                    
                    temporal_features.update(annotation_features & all_features)
                    not_found_annotations = annotation_features - temporal_features
                    annotation_features = temporal_features
                    if len(not_found_annotations) > 0:
                        sys.stderr.write(f"The features {not_found_annotations} were not found for the given parameters.\n")
                if len(annotation_features) == 0:
                    sys.stderr.write("No features found for the given parameters.\n")
                    sys.exit(1)
            except:
                sys.stderr.write("Feature unpacking gone wrong.\n")
                sys.exit(1)

        except:
            sys.stderr.write("Error at feature setup.\n")
            sys.exit(1)
                
        return table_info, uniprot_info, annotation_features, feature_collection

    def update_features(self, update_parameters):
        """Method to update our case study when
        new parameters arrive.
        """
        if update_parameters["calc_alg"][0] is not "": # IF PARAMETER HAS BEEN MODIFIED
            self.calc_alg = update_parameters["calc_alg"][0]
        if "features" in update_parameters:
            self.study_features = set(update_parameters["features"])
        if update_parameters["evalue"][0] is not "":
            self.min_eval = float(update_parameters["evalue"][0])
        if (config.calculus_algorithms[self.calc_alg]["differentiate_gaps"]) == "N": # TO ENSURE THEY DONT CHANGE GAP PARAMETER IN A NOT ALLOWED ALGOORITHM
            self.differentiate_gaps = "N"
        else:
            if "diff_gaps" in update_parameters and update_parameters["diff_gaps"][0] is not "":
                self.differentiate_gaps = update_parameters["diff_gaps"][0]
        self.position_matrix = None

        return


    def calculate_nodes(self):
        """Method to calculate the different internal node scores
        for a given calculus method, and store those values both in
        a dictionary (if the user wants to) and in an instance
        of a processed tree.
        """
        try:
            # print(f"""Computing tree with the following parameters: 
            # - STUDY FEATURES: {self.study_features}
            # - EVALUE THRESHOLD: {self.min_eval}
            # - CALCULUS ALGORITHM: {self.calc_alg}
            # - DIFFERENTIATE GAPS: {self.differentiate_gaps}""")

            tree = PhyloTree(self.tree_in, alignment=self.align_in, alg_format="fasta")
            md = tree.get_midpoint_outgroup()
            tree.set_outgroup(md)
            leaf_deleting_list = set()
            if self.position_matrix == None:
                uniprot_hit_hash, leaf_deleting_list = fp.retrieve_features(self.study_features, self.table_info, self.min_eval, self.uniprot_info)
                self.position_matrix = fp.get_positions_matrix(uniprot_hit_hash, tree)
            node_number = 0
            node_scores = {}
            for leaf in tree.iter_leaves():
                if leaf.name in leaf_deleting_list:
                    leaf.delete()
            for index, node in enumerate(tree.traverse("preorder")):
                node._nid = index
                if node.is_leaf() == False:
                    node_score = round(fp.calculate_node_score(node, self.position_matrix, self.calc_alg, self.differentiate_gaps), 2) ###AQUI EL ALGORITMO DEBERIA VENIR DE UN ARGUMENTO?
                    node.add_feature("node_score", node_score)
                    node_scores[node_number] = node_score
                    node_number += 1
            self.processed_tree = tree
            self.node_scores = node_scores
        except:
            sys.stderr.write("Error at calculating nodes.\n")
            sys.exit(1)

        return 


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
                    for position in self.position_matrix:
                        seqFace = SeqMotifFace(node.sequence[position], seq_format="seq")
                        (self.processed_tree&node.name).add_face(seqFace, draw_position, "aligned")
                        draw_position += 1  
                else:
                    if node.node_score > plot_threshold:
                        score_face = TextFace(node.node_score)
                        node.add_face(score_face, 0, "branch-top")
        except:
            sys.stderr.write("Error at designing tree.\n")
            sys.exit(1)

        return


    def plot_tree(self, outfile, tree="DESIGN", width=800, heigth=2000, units="px", plot_threshold=0):
        """Method to plot a tree to a file given a set of
        parameters. Units may be “px”: pixels, “mm”: millimeters, “in”: inches.
        Plot threshold referes to the minimum score a node must have to plot
        it in the tree image.
        """
        try:
            ts = TreeStyle()
            ts.layout_fn = lambda x: True
            if tree == "DESIGN":
                self.design_tree(plot_threshold=plot_threshold)
                tree = self.processed_tree
            tree.render(outfile, w=width, h=heigth, units=units, tree_style=ts)
        except:
            sys.stderr.write("Error at plotting tree.\n")
            sys.exit(1)


    def write_node_file(self, outfile):
        """Method to write our node score list to a file
        """
        try:
            with open(outfile, 'w') as file:
                for node in self.node_scores:
                    file.write(f"{node}\t{self.node_scores[node]}\n")
                file.close()
        except:
            sys.stderr.write("Error at writing nodes to file.\n")
            sys.exit(1)


## END