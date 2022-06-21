# DGAT-onco
DGAT-onco was designed to detect oncogenes by comparing the functional impact-based profiles between cancer population and natural population.  

It is developped based on the framework of DiffMut (https://github.com/Singh-Lab/Differential-Mutation-Analysis), a frequency-based method.  
DiffMut detects oncogenes by comparing the frequency-based distribution between cancer and natural population. But it treats each mutation equally.  
Thus, we uses a scoring function to consider the differnce of their functional impact and make this method.  

Reference: https://ieeexplore.ieee.org/document/9669388

Comments and feedbacks: zhanghaoyang0@hotmail.com  
  
# Overview  
![image](https://github.com/zhanghaoyang0/DGAT-onco/blob/master/fig/overview.png)
# Requirements
* R (≥4.0.3)  
Packages: dplyr (≥1.0.7), tidyr (≥1.1.3) 
* ANNOVAR  
Databases: refGene, dbnsfp35c
* Linux (to perform ANNOVAR)
# Example files (hg19-based mutation data, in the "example" folder)
* ```BRCA.maf```, a Mutation Annotation Format (MAF) file, is the input of our method, .  
* We only need 5 columns, as below. It looks like:  
Chromosome	Start_Position	Tumor_Seq_Allele1	Tumor_Seq_Allele2	Tumor_Sample_Barcode  
chr1	100535194	G	C	TCGA-A7-A0DA-01A-31D-A10Y-09  
chr1	100547556	C	G	TCGA-AR-A2LE-01A-11D-A17W-09  
...  
* ```BRCA_for_annovar.vcf```, a Variant Call Format (VCF) file prepared for ANNOVAR, also can be the input of our method.
* It can be converted from the MAF file, with the function ```make_vcf_for_annovar```.  
* It looks like:  
chr1	100535194	.	G	C	.	.	TCGA-A7-A0DA-01A-31D-A10Y-09  
chr1	100547556	.	C	G	.	.	TCGA-AR-A2LE-01A-11D-A17W-09  
...  
* ```BRCA.hg19_multianno.vcf```  is a VCF file annotated with ANNOVAR. 
* It looks like:  
Chr	Start	End	Ref	Alt	Func.refGene	Gene.refGene	GeneDetail.refGene	ExonicFunc.refGene	AAChange.refGene	SIFT_score ...  
chr1	100535194	100535194	G	C	exonic	MFSD14A	.	nonsynonymous SNV	MFSD14A:NM_033055:exon8:c.G844C:p.V282L	0.266 ...  
chr1	100547556	100547556	C	G	intronic	MFSD14A	. .	.	. ...  
...  
* ```BRCA_result.txt``` is the uEMD calculated for each gene:  
* It looks like:  
gene	uemd  
A1CF	0  
A2M	0  
...  
# Pre-calculated functional impact-based distributions of natural (1000 genomics [1000g]) population (in the "data" folder)
* ```bins_1000g.rdata```, an Rdata of impact-based distributions of 1000g population.  
* 20 scoring functions are provided, they are: 
```"SIFT_converted_rankscore"```, ```"LRT_converted_rankscore"```, ```"MutationTaster_converted_rankscore"```, ```"MutationAssessor_score_rankscore"```, ```"FATHMM_converted_rankscore"```, ```"PROVEAN_converted_rankscore"```, ```"MetaSVM_rankscore"```, ```"MetaLR_rankscore"```, ```"M.CAP_rankscore"```, ```"MutPred_rankscore"```, ```"fathmm.MKL_coding_rankscore"```, ```"GenoCanyon_score_rankscore"```, ```"integrated_fitCons_score_rankscore"```, ```"GERP++_RS_rankscore"``` , ```"phyloP100way_vertebrate_rankscore"```, ```"phyloP20way_mammalian_rankscore"```, ```"phastCons100way_vertebrate_rankscore"```, ```"phastCons20way_mammalian_rankscore"```, ```"SiPhy_29way_logOdds_rankscore"```, ```"constant_rankscore"```  
* ```"constant_rankscore"``` means each mutation has equal (constant) weight. 
# Example
## Load functions and packages  
source('./code/dgat-onco.r')  
names(bins_1000g) # list scoring functions
## Make vcf for ANNOVAR  
make_vcf_for_annovar(chr_col = 'Chromosome', pos_col = 'Start_Position', ref_col = 'Tumor_Seq_Allele1',
                     alt_col = 'Tumor_Seq_Allele2', id_col = 'Tumor_Sample_Barcode',
                     file_in = './example/BRCA.maf', file_out = './example/BRCA_for_annovar.vcf')
* This function generate a VCF file from a MAF input.  
## ANNOVAR (on Linux)
file_in="./example/BRCA_for_annovar.vcf"  
file_out='./example/BRCA'  
path_annovar='path_of_your_annovar' # path of your annovar  

$path_annovar\table_annovar.pl \  
$file_in \  
$path_annovar\humandb/ \  
-buildver hg19 \  
-out $file_out \  
-remove -protocol refGene,dbnsfp35c -operation g,f -nastring . -vcfinput -polish
## Calculate mutation summation score for each cancer individual
summ = calculate_summation_score(file_in = './example/BRCA.hg19_multianno.txt', id_col = 'Otherinfo11',
                                 gene_col = 'Gene.refGene', func_col= 'ExonicFunc.refGene')  
colnames(summ)[grepl('score', colnames(summ))] # list scoring functions
## Calculate uEMD based on a scoring function
calculate_uemd(summ, file_out = './example/BRCA_result.txt', score = 'MutPred_rankscore')  
# If your input is hg38-based  
Just set the ```-buildver``` flag in ANNOVAR to hg38 and the rest is the same.
