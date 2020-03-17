###############################################################
### Script to load the trees using ete3, cutting the        ###
### functional sites following the swissprot annotation.    ###
###############################################################

# Manually do: module load ETE/3.1.1-foss-2019a-Python-3.7.2
# Manually do: module load scikit-learn/0.20.3-foss-2019a


### USING THE CLUSTER NAME ###
import sys 
cluster_name = sys.argv[1]
partition_number = sys.argv[2]
BASE_PATH = "../data/partitions"

### FILE HANDLING THE INPUT AND OUTPUT FILES ###
from pathlib import Path
tree_infile = f"{BASE_PATH}/{partition_number}/trees/{cluster_name}.tree"
alignment_infile = f"{BASE_PATH}/{partition_number}/alignments/{cluster_name}.fas.alg"    # Set the file output path
table_infile = f"{BASE_PATH}/{partition_number}/tables/{cluster_name}.tsv"

### TREE WORKOUT ###
from ete3 import PhyloTree
tree = PhyloTree(tree_infile, alignment=alignment_infile, alg_format="fasta")

# print (tree)
