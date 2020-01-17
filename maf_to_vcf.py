#!/usr/bin/env python3
# this script generate a vcf of missense mutations from a maf file
# usage: python3 maf_to_vcf.py input_file output_file
# "input_file": passing a address of inpute file
# "output_file": passing a address of output file
# example: python3 maf_to_vcf.py data/BRCA.maf ./BRCA.vcf
import sys
import pandas as pd
import os

print('NOTICE: this script generate a vcf of missense mutations from a maf file\n'
      'USAGE: python3 maf_to_vcf.py input_file output_path output_name\n'
      '"input_file": passing a address of inpute file\n'
      '"output_file": passing a address of output file\n'
      'EXAMPLE: python3 maf_to_vcf.py data/BRCA.maf ./BRCA.vcf')

if len(sys.argv) != 3:
    print('ERROR: please check your number of argument!')
    sys.exit()

input_file = sys.argv[1]
output_file = sys.argv[2]

print('loading your maf...')
maf = pd.read_csv(os.getcwd()+'/'+input_file, sep='\t', low_memory=False)
print('done!')

if 'Variant_Classification' not in maf.columns:
    print('ERROR: no col names Variant_Classification found in your maf!\n')
    sys.exit()
if 'Missense_Mutation' not in set(maf['Variant_Classification']):
    print('ERROR: no Missense_Mutation found in your Variant_Classification column!')
    sys.exit()

maf = maf.loc[maf['Variant_Classification'] == 'Missense_Mutation', :].reset_index(drop=True)
if 'chr' in maf['Chromosome'][0]:
    chr = maf['Chromosome'].str.split('chr', expand=True)[1]
else:
    chr = maf['Chromosome']

print('making a vcf...')
dot = pd.DataFrame(['.']*maf.shape[0])
info = pd.DataFrame([gene+';'+sample_id+';'+mut_type for (gene, sample_id, mut_type) in zip(maf['Hugo_Symbol'], maf['Tumor_Sample_Barcode'], maf['Variant_Classification'])])
vcf = pd.concat([chr, maf['Start_Position'], dot, maf['Reference_Allele'], maf['Tumor_Seq_Allele2'], dot, dot, info], axis=1)
vcf.to_csv(os.getcwd()+'/'+output_file, sep='\t', header=None, index=False)

n_id = str(len(set(maf['Tumor_Sample_Barcode'])))
n_gene = str(len(set(maf['Hugo_Symbol'])))
n_mismut = str(maf.shape[0])
print('done!\n'
      'n of id: '+n_id+', n of gene: '+n_gene+', n of missense mutation: '+n_mismut+
      '\noutfile ' + output_file + ' have been writed!')
