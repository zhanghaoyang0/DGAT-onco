#!/usr/bin/env python3
# this script generate a mutation matrix based on a pathogenic score from a maf file
# usage: python3 vcf_to_mutmatrix.py input_file output_file score impute_missing
# "input_file": passing a address of inpute file
# "output_file": passing a address of output file
# "score": passing a scoring function to score mutations
# "impute mutations": passing a value from 0-1 to score mutations without score annotation
# example: python3 vcf_to_mutmatrix.py data/BRCA.hg38_multianno.vcf data/BRCA_mutmatrix.csv M-CAP_rankscore 0
import sys
import pandas as pd
import os

all_score = ['SIFT_converted_rankscore', 'Polyphen2_HDIV_rankscore', 'Polyphen2_HVAR_rankscore', 'LRT_converted_rankscore',
 'MutationTaster_converted_rankscore', 'MutationAssessor_score_rankscore', 'FATHMM_converted_rankscore',
 'PROVEAN_converted_rankscore', 'VEST3_rankscore', 'MetaSVM_rankscore', 'MetaLR_rankscore', 'M-CAP_rankscore',
 'CADD_raw_rankscore', 'CADD_phred', 'DANN_rankscore', 'fathmm-MKL_coding_rankscore', 'GenoCanyon_score_rankscore', 'integrated_fitCons_score_rankscore',
 'GERP++_RS_rankscore', 'phyloP100way_vertebrate_rankscore', 'phyloP20way_mammalian_rankscore', 'phastCons100way_vertebrate_rankscore',
 'phastCons20way_mammalian_rankscore', 'SiPhy_29way_logOdds_rankscore', 'constant']

print('NOTICE: this script generate a mutation matrix based on a pathogenic score from a maf file\n'
      'USAGE: python3 vcf_to_mutmatrix.py input_file output_path output_name score impute_missing\n'
      '"input_file": passing a address of inpute file\n'
      '"output_file": passing a address of output file\n'
      '"score": passing a scoring function to score mutations\n'
      '"impute mutations": passing a value from 0-1 to score mutations without score annotation\n'
      'NOTICE: score is a scoring function to generate weight of each mutation, including:\n'+str(all_score)+
      'EXAMPLE: python3 vcf_to_mutmatrix.py data/BRCA.hg38_multianno.vcf data/BRCA_mutmatrix.csv M-CAP_rankscore 0'
)

if len(sys.argv) != 5:
    print('ERROR: please check your number of argument!')
    sys.exit()

input_file = sys.argv[1]
output_file = sys.argv[2]
score = sys.argv[3]
impute = float(sys.argv[4])

if impute > 1 or impute < 0:
    print('ERROR: impute_missing should be in [0, 1]!')
    sys.exit()
if score not in all_score:
    print('ERROR: score should be in')
    print(str(all_score)+'!')
    sys.exit()

print('loading your annotated vcf...')
vcf = pd.read_csv(os.getcwd()+'/'+input_file, sep='\t', header=None)
if ('ANNOVAR' not in vcf.iloc[0,7]):
    ('ERROR: please check if your vcf have been annotated by ANNOVAR!')
    sys.exit()
print('done!')
gene = pd.DataFrame([x.split(';')[0] for x in vcf.iloc[:, -1]], columns=['gene'])
id = pd.DataFrame([x.split(';')[1] for x in vcf.iloc[:, -1]], columns=['id'])
print('n of id: '+str(len(set([x[0] for x in id.values])))+', n of gene: '+str(len(set([x[0] for x in gene.values])))+', n of missense mutation: '+str(vcf.shape[0]))
print('summarizing...')

# separate and filter score columns
item = vcf.iloc[:, -1].str.split(';', expand=True).iloc[:, 4:-1] # annotation part
colname = []
for col in range(0, item.shape[1]):
    colname += [item.iloc[:, col].str.split('=', expand=True)[0][1]]
    item.iloc[:, col] = item.iloc[:, col].str.split('=', expand=True)[1]
item.columns = colname
item = item.replace({'.': None})  # replace '.' with NA
item['constant'] = 1
mut_score = pd.concat([id, gene, item[score]], axis=1).fillna(impute).rename(columns={score: "score"})
mut_score['score'] = pd.to_numeric(mut_score['score'])
mut_matrix = pd.pivot_table(mut_score, values='score', index='gene', columns='id', aggfunc='sum')
mut_matrix.fillna(0, inplace=True)
mut_matrix = mut_matrix.apply(lambda x: x.rank(method='min') / mut_matrix.shape[0])  # rank normalize
print('done!\nsaving...')
mut_matrix.to_csv(os.getcwd()+'/'+output_file)
print('outfile '+output_file+' have been writed!')