# DiffMut-W
DiffMut-W was designed to detect oncogenes by comparing the mutaitonal profile between cancer population and natural population.

The method is still on progress, please send comments and feedbacks to zhanghaoyang0@hotmail.com

# Requirements
Operation system: Unix/Linux

python 3 (with os, sys, numpy, pandas, math packagegs installed)

# Data preparation
inputfile is a VCF file with the INFO annotated by refGene and dbnsfp33a database through ANNOVAR (http://annovar.openbioinformatics.org/en/latest/)

INFO fields are encoded as a semicolon-separated series of short keys with optional values.
Please add gene name(at first field), mutatioin type and patient id before anotating on ANNOVAR:

1	152760163	.	G	A	.	.	KPRP;Missense_Mutation;TCGA-BH-A0HO-01A-11W-A050-09

After annotation, the vcf will be like:

1	152760163	.	G	A	.	.	KPRP;Missense_Mutation;TCGA-BH-A0HO-01A-11W-A050-09;ANNOVAR_DATE=2018-04-16;...;M-CAP_score=0.005;...

... indicate information omitted here.
M-CAP_score is the prior mutation weight we used in analysis. 
Please put yuur inputfile file in data folder.
# Usage
python calculate.py data/inputfile top
where inpute file is your file name, top is the number of output

Example:
python calculate.py data/BRAC.vcf 50

Result files are:
1. a csv with genes names and coresponding uEMDs. 
2. a csv with genes have not results (because they have not corresponding 1000G gene with missense_mutation).
