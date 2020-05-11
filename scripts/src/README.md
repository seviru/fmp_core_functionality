# Folder containing the main algorithm to be called by the user, as well as a src folder where all functionalities are developed.

## __init__.py
Empty file needed fot Python to recognise this folder as a module.

## feature_processing.py
File containing different functions that carry on the feature processing for a given case study, such as:
    - check_features: Check what features are available for our dataset.
    - retrieve_features: Retrieve the feature information for out dataset.
    - get_alignment_position: Transform a given sequence position in the equivalent alignment position.
    - get_positions_matrix: Get a matrix of the equivalent positions in our alignmetn for all the required features.
    - calculate_node_score: Calculate the node score for the given annotations.

## scoring_functions.py
File containing different functions to computate the node score:
    - get_scorer: For a given calculus algorithm, retrieves the needed calculus function.
    - simple_calculus: Compares the different aminoacids between the two branches.
    - whole_annotation_simple_calculus: Same as simple_calculus but takes all the positions as if it was just one.
    - all_vs_all_calculus: Compares the different aminoacids between the branches and inside the branch itself.
    - whole_annotation_all_vs_all_calculus: Same as all_vs_all_calculus but takes all the positions as if it was just one.
    - all_vs_all_calculus_means: Compares the different aminoacids between the branches and inside the branch itself, giving it a mean value.
    - whole_annotation_all_vs_all_calculus_means: Same as all_vs_all_calculus_means but takes all the positions as if it was just one.

## utils.py
File containing functions with flexible utilities, such as:
    - dict_diff: Returns a dictionary result of the elements only present in one
    of the inputs. The key point here is that returns the key and values pairs without operating with them.

## main_class.py
File containing the class FeatureStudy definition, as well as the different methods that call the functions to carry on the analysis. It has a __setup__ method as well that checks that all the needed inputs are correctly given, and that the analysis can actually be carried out.