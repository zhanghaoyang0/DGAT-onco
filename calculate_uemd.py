#!/usr/bin/env python3
# this script generate uemd and its significance of each gene from a mutation matrix file
# usage: python3 calculate_uemd.py input_file output_file score sig_test
# "input_file": passing a address of inpute file
# "output_file": passing a address of output file
# "score": passing a scoring function to selected scored mutations in 1000G data
# "sig_test": passing 1 for significance test, other 0
# example: python3 calculate_uemd.py data/BRCA_mutmatrix.csv data/BRCA_uemd  M-CAP_rankscore 0

import sys
import pandas as pd
import os
import math
import numpy as np
from scipy import stats

def getuemd(health_bin, cancer_vector): # input two list get uemd
    cancer_bin = [0] * 100
    for k in cancer_vector:
        cancer_bin[math.ceil(k * 100) - 1] += 1
    cancer_bin = [x / len(cancer_vector) for x in cancer_bin]
    move = 0
    score = 0
    for k in reversed(range(100)):
        move += cancer_bin[k] - health_bin[k]
        score += max(0, move)
    return float(score)

all_score = ['SIFT_converted_rankscore', 'Polyphen2_HDIV_rankscore', 'Polyphen2_HVAR_rankscore', 'LRT_converted_rankscore',
 'MutationTaster_converted_rankscore', 'MutationAssessor_score_rankscore', 'FATHMM_converted_rankscore',
 'PROVEAN_converted_rankscore', 'VEST3_rankscore', 'MetaSVM_rankscore', 'MetaLR_rankscore', 'M-CAP_rankscore',
 'CADD_raw_rankscore', 'CADD_phred', 'DANN_rankscore', 'fathmm-MKL_coding_rankscore', 'GenoCanyon_score_rankscore', 'integrated_fitCons_score_rankscore',
 'GERP++_RS_rankscore', 'phyloP100way_vertebrate_rankscore', 'phyloP20way_mammalian_rankscore', 'phastCons100way_vertebrate_rankscore',
 'phastCons20way_mammalian_rankscore', 'SiPhy_29way_logOdds_rankscore', 'constant']

print('NOTICE: this script generate uemd and its significance of each gene from a mutation matrix file\n'
      'USAGE: python3 calculate_uemd.py input_file output_file sig_test\n'
      '"input_file": passing a address of inpute file\n'
      '"output_file": passing a address of output file\n'
      '"score": passing a scoring function to selected scored mutations in 1000G data\n'
      '"sig_test": passing 1 for significance test, other 0\n'
      'NOTICE: score is a scoring function to generate weight of each mutation, including:\n'+str(all_score)+
      '\nEXAMPLE: python3 calculate_uemd.py data/BRCA_mutmatrix.csv data/BRCA_uemd  M-CAP_rankscore 0\n')

input_file = sys.argv[1]
output_file = sys.argv[2]
score = sys.argv[3]
sig_test = sys.argv[4]

print('loading mutation matrix...')
mutmatrix = pd.read_csv(os.getcwd()+'/'+input_file, index_col=0)
print('done!\nloading mutaiton profile from 1000g data...')
natbin = pd.read_csv('data/1000g/' + score + '_bin', sep='\t', header=0)
print('done!\ncalculating uemd...')

gene_list = list(sorted(set(natbin.columns).intersection(set(mutmatrix.index))))
cancer_vector_list = [mutmatrix.loc[mutmatrix.index==x, ].values.tolist()[0] for x in gene_list]
health_bin_list = [natbin[x].values.tolist() for x in gene_list]
uemd_list = [getuemd(x, y) for x,y in zip(health_bin_list, cancer_vector_list)]
uemd = pd.DataFrame(list(zip(gene_list, uemd_list)), columns=['gene', 'uemd'])
print('done!')

if sig_test == '1':
    print('doing significance test, it will take some time...')
    ranuemd = pd.DataFrame()
    for shift in range(2):
        print('shifting mutation rank and calculate ranuemd... '+str(shift + 1) + '/5')
        np.random.seed(shift)
        mutmatrix_shift = mutmatrix.copy()
        mutmatrix_shift.apply(np.random.shuffle, axis=0)
        cancer_vector_list = [mutmatrix_shift.loc[mutmatrix_shift.index == x,].values.tolist()[0] for x in gene_list]
        uemd_list = [getuemd(x, y) for x, y in zip(health_bin_list, cancer_vector_list)]
        ranuemd = pd.concat([ranuemd, pd.DataFrame([gene_list, [shift]*len(gene_list), uemd_list]).T])
    ranuemd.columns = ['gene', 'shift', 'ranuemd']
    ranuemd.columns = ['gene', 'shift', 'ranuemd']
    print('done!\ndoing significance test...')
    # one sample one sided t-test
    t = [float(stats.ttest_1samp(ranuemd.loc[ranuemd['gene']==x, 'ranuemd'].values.tolist(),
                           uemd.loc[uemd['gene']==x, 'uemd'].values.tolist())[0]) for x in gene_list]
    p = [float(stats.ttest_1samp(ranuemd.loc[ranuemd['gene']==x, 'ranuemd'].values.tolist(),
                           uemd.loc[uemd['gene']==x, 'uemd'].values.tolist())[1]/2) for x in gene_list]
    p = [0 if x==-math.inf else 1 if (x==math.inf or math.isnan(x)) else x for x in p]
    p = [1-y if x > 0 else y for x, y in zip(t, p)]
    uemd['p'] = p
    n_test = len(p) - sum(1 if math.isnan(x) else 0 for x in p)  # n test = total test - number of nan result
    uemd = uemd.sort_values(by=['p'])
    uemd['rank'] = range(1, uemd.shape[0] + 1)
    uemd['q'] = uemd['p'] * n_test / uemd['rank']
    uemd = uemd[['gene', 'uemd', 'q']].sort_values(by=['gene'])
    print('done!')
uemd.to_csv(os.getcwd()+'/'+output_file)
print(output_file + ' have been writed!')
