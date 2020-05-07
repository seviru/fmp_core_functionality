### TRYING TO PUT EVERYTHING INSIDE A CLASS, STILL WORK IN PROGRESS ###



class FeatureStudy:
    def __init__(self, tree_infile, alignment_infile, table_infile,
                uniprot_infile, min_evalue, annotation_feature,
                tree_outfile, node_score_table)
        self.tree_in = tree_infile
        self.align_in = alignment_infile
        self.table_in = table_infile
        self.uniprot_in = uniprot_infile
        self.min_eval = min_evalue
        self.study_feature = annotation_feature
        self.tree_out = tree_outfile
        self.score_table = node_score_table
        self.calc_alg = node_score_algorithm
        self.ignore_gaps = ignore_gap_positions


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