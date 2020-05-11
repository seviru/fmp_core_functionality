# Folder containing different example input files.

## alignment.example
Alignment file in .fasta format. The headers will be the same as the leaf names for our tree.example file. For more information regarding .fasta format, check: https://en.wikipedia.org/wiki/FASTA_format

## table.example
Tabulated file containing the different hits that our individual sequences have to a known uniprot protein. Format example:
<SEQUENCE_NAME><\t><SEQUENCE_SOURCE><\t><HIT_TYPE><\t><HIT_NAME><\t><HIT_EVALUE><\t><HIT_SCORE><\t><HIT_IDENTITY><\t><QUERY_COVERY><\t><TARGET_COVERY>
"""
cat table.example | head -n 2
294_737_599 HEMA    spb	A3CL77	2.32e-48	169	0.966	0.8476190476190476  0.21342925659472423
147_661_409	HEMA	spb	A3CL77	3.59e-127	401	0.956	1.0	0.49640287769784175
"""

## tree.example
File containing our tree in Newick format. For more information regarding this format, check https://en.wikipedia.org/wiki/Newick_format

## uniprot.example
File containing all the information for our annotated proteins. Format example:
{<UNIPROT_IDENTIFIER>:{"ID":<PROTEIN_NAME>, "FT":[<FEATURES_LIST>], "SQ": <SECUENCE>}}
FEATURES_LIST format:
[{"ft": <FEATURE_TYPE>, "s": <FEATURE_START>, "e": <FEATURE_END>, "ann": <FEATURE_ANNOTATION>}]
"""
cat uniprot.example | head -n 1
{"A3CL77": {"ID": "HEM1_STRSV", "EC": "1.2.1.70", "FT": [{"ft": "CHAIN", "s": "1", "e": "1", "ann": "Glutamyl-tRNA reductase./FTId=PRO_1000004706."}, {"ft": "NP_BIND", "s": "189", "e": "189", "ann": "NADP. {ECO:0000255|HAMAP-Rule:MF_00087}."}, {"ft": "REGION", "s": "49", "e": "49", "ann": "Substrate binding. {ECO:0000255|HAMAP-Rule:MF_00087}."}, {"ft": "REGION", "s": "114", "e": "114", "ann": "Substrate binding. {ECO:0000255|HAMAP-Rule:MF_00087}."}, {"ft": "ACT_SITE", "s": "50", "e": "50", "ann": "Nucleophile. {ECO:0000255|HAMAP-Rule:MF_00087}."}, {"ft": "BINDING", "s": "109", "e": "109", "ann": "Substrate. {ECO:0000255|HAMAP-Rule:MF_00087}."}, {"ft": "BINDING", "s": "120", "e": "120", "ann": "Substrate. {ECO:0000255|HAMAP-Rule:MF_00087}."}, {"ft": "SITE", "s": "99", "e": "99", "ann": "Important for activity.{ECO:0000255|HAMAP-Rule:MF_00087}."}], "SQ": "MHLLYVGLTHRETPLTILEKAHFSDQEGLKALKLLKREKSILENIILSTCNRTELYLVVDQLHTGRYYSKHFLADWFQIPVKELEEYLVFREGDEALRHLLRVSIGLESKIVGESQVLGQLKQAFLTAQDAGTTGIVLNQAFKQALTFAKRMHDTYRINDRPISIGLTAIQELDRMGLDYSTKKIAVIGLGEIGQLVTKYALQRPFESVMLLNRTVSKAQAFLTEDRVSAHGWDELEEVLADADVVFSAVKTEEYIIFPSMLKEGAIVFDLCLPRSCHPSSSLKLYNIENLTNQLEQYKAERQEIAGRIALEIDEELVKFADWRQQLGIIPLIQEIRDKALEAQASAMESLNRKIPDLTEREQKQISKHMKSIINQVLKEPILQLKELSVGEHSDYDIALIAKIFGLHRERGKDEGH"}}
"""