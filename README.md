# DGAT-onco
DGAT-onco was designed to detect oncogenes by comparing the functional impacts profile between cancer population and natural population.  
Comments and feedbacks: zhanghaoyang0@hotmail.com  
# Requirements
Operation system: Unix/Linux  
Python 3  
Modules: os, sys, numpy, pandas, math, stats  
# Example files
* ```BRCA.maf``` is a MAF file like:  
Hugo_Symbol	Entrez_Gene_Id	Center	NCBI_Build	Chromosome	...  
FOXD3	27022	WUGSC	GRCh38	chr1	...  
...  

* ```BRCA.vcf```  is a VCF file like:  
1	152760163	.	G	A	.	.	KPRP;TCGA-BH-A0HO-01A-11W-A050-09;Missense_Mutation  
3	172133474	.	A	G	.	.	FNDC3B;TCGA-BH-A0HO-01A-11W-A050-09;Missense_Mutation  
...  

* ```BRCA.hg38_multianno.vcf```  is an annotated VCF file like:  
1	152760163	.	G	A	.	.	KPRP;TCGA-BH-A0HO-01A-11W-A050-09;Missense_Mutation;ANNOVAR_DATE=2018-04-16;SIFT_score=.;SIFT_converted_rankscore=.;SIFT_pred=.;Polyphen2_HDIV_score=0.002 ...  
3	172133474	.	A	G	.	.	FNDC3B;TCGA-BH-A0HO-01A-11W-A050-09;Missense_Mutation;ANNOVAR_DATE=2018-04-16;SIFT_score=0.006;SIFT_converted_rankscore=0.614;SIFT_pred=D;Polyphen2_HDIV_score=0.619;Polyphen2_HDIV_rankscore=0.380;Polyphen2_HDIV_pred=P;Polyphen2_HVAR_score=0.311 ...  
...  

* ```BRCA_mutmatrix.csv``` is a mutation functional impacts matrix like:  
gene,TCGA-3C-AAAU-01A-11D-A41F-09,TCGA-3C-AALI-01A-11D-A41F-09 ...  
A1CF,7.833920877399138e-05,7.833920877399138e-05 ...  
...  

* ```BRCA_uemd``` is the calculation result like:  
gene,uemd  
A1CF,1.8041124150158794e-16  
A2M,0.0  
...  

# Pre-calculated functional impact files from 1000Genomics data(natural population)
Totaly 25 files of functional impacts depending to different scoring function: including:  
```'SIFT_converted_rankscore'```, ```'Polyphen2_HDIV_rankscore'```, ```'Polyphen2_HVAR_rankscore'```, ```'LRT_converted_rankscore'```,
 ```'MutationTaster_converted_rankscore'```, ```'MutationAssessor_score_rankscore'```, ```'FATHMM_converted_rankscore'```,
 ```'PROVEAN_converted_rankscore'```, ```'VEST3_rankscore'```, ```'MetaSVM_rankscore'```, ```'MetaLR_rankscore'```, ```'M-CAP_rankscore'```,
 ```'CADD_raw_rankscore'```, ```'CADD_phred'```, ```'DANN_rankscore'```, ```'fathmm-MKL_coding_rankscore'```, ```'GenoCanyon_score_rankscore'```, ```'integrated_fitCons_score_rankscore'```,
 ```'GERP++_RS_rankscore'```, ```'phyloP100way_vertebrate_rankscore'```, ```'phyloP20way_mammalian_rankscore'```, ```'phastCons100way_vertebrate_rankscore'```,
 ```'phastCons20way_mammalian_rankscore'```, ```'SiPhy_29way_logOdds_rankscore'```, ```'constant'```  
 Note that ```'constant'```  mean the weight of each mutation are equal.
# Script files
* ```maf_to_vcf.py```　　
This script generate a VCF of missense mutations from a MAF file  

Usage: python3 maf_to_vcf.py input_file output_file  
```input_file```: passing a address of inpute file  
```output_file```: passing a address of output file  
Example: python3 maf_to_vcf.py data/BRCA.maf ./BRCA.vcf  

* ```vcf_to_mutmatrix.py```

This script generate a mutation matrix based on a pathogenic score from VCF file annotated by dbnsfp33a dataset. 

usage: python3 vcf_to_mutmatrix.py input_file output_file score impute_missing  
```input_file```: passing a address of inpute file  
```output_file```: passing a address of output file  
```score```: passing a scoring function to score mutations  
"impute mutations": passing a value from 0-1 to score mutations without score annotation  
example: python3 vcf_to_mutmatrix.py data/BRCA.hg38_multianno.vcf data/BRCA_mutmatrix.csv M-CAP_rankscore 0  

Note that 25 scoring functions were provided (see above).

* ```calculate_uemd.py``` 
This script generate uemd and its significance of each gene from a mutation matrix file.  

Usage: python3 calculate_uemd.py input_file output_file score sig_test  
```input_file```: passing a address of inpute file  
```output_file```: passing a address of output file  　
```score```: passing a scoring function to selected scored mutations in 1000G data  
```sig_test```: passing 1 for significance test, other 0  
Example: python3 calculate_uemd.py data/BRCA_mutmatrix.csv data/BRCA_uemd  M-CAP_rankscore 0  

After calculation, the result will be like:  
## Annotation
The VCF generated from ```maf_to_vcf.py```　should be annotated by  dbnsfp33a dataset.  

Depending your MAF build, the annotation can be conducted by ANNOVAR('http://annovar.openbioinformatics.org/en/latest/user-guide/download/') by following commands:  

table_annovar.pl your_input humandb/ -buildver your_build -out your_output -remove -protocol dbnsfp33a -operation f -nastring . -vcfinput  

