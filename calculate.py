#!/usr/bin/env python3

import pandas as pd
import numpy as np
import math
import os
import sys

# load data
if len(sys.argv) == 1:
    print('No inputfile has been specified!')
    os._exit()
if len(sys.argv) == 2:
    top = 50
else:
    top = sys.argv[2]
inputfile_path = os.getcwd() + '/' + sys.argv[1]


def getuemd(inputfile_path, top=50):
    # load data
    print('loading data ...')
    path = inputfile_path[:-len(inputfile_path.split('/')[-1])]
    filename = inputfile_path.split('/')[-1]
    os.chdir(path)
    d_1000G = np.load('1000G_M-CAP_weighted.npy').item()  # load 1000G mutation data(pre calculated for uemd measurement)
    df = pd.read_csv(inputfile_path, sep='\t', header=None)  # load data
    weight = 'M-CAP_rankscore'  # use M-CAP as prior weight
    if 'Hugo_Symbol' not in df.iloc[0,7]:
        print('break because not genename in your data! you should present it like "Hugo_Symbol = genename" in your INFO fields')
        return
    if df.shape[0] == 0:
        print('break because not Missense_Mutation in your data!')
        return
    if weight not in df.iloc[0, 7]:
        print('break because' + weight + 'not in your annotation!')
        return
    print('done!')
    # make mutation score matrix
    print('making mutation matrix ...')
    df = df[df.iloc[:, 7].str.contains('Missense_Mutation')]  # select Missense_Mutation
    item = df.iloc[:, 7]
    score = [x.split(weight)[1].split(';')[0].split('=')[1] for x in item]  # get M-CAP score
    score = [float(x) if x != '.' else 0 for x in score]  # impute missing weight by 0

    gene = [x.split('Hugo_Symbol=')[1].split(';')[0] for x in item]  # get gene
    id = [x.split('TCGA-')[1].split(';')[0] for x in item]  # get tcga id
    df = pd.DataFrame([id, gene, score]).T
    df.columns = ['id', 'gene', 'score']
    df = pd.pivot_table(df, values='score', index='gene', columns='id', aggfunc='sum')
    df.fillna(0, inplace=True)
    df = df.apply(lambda x: x.rank() / df.shape[0])  # rank normalized
    print('done!')
    # calculate uemd
    print('calculating uemd ...')
    d_uemd = dict()  # uemd dict
    print(str(len(set(df.index) - set(
        d_1000G.keys()))) + ' genes have not results because they are not found in 1000G gene with missense_mutation')
    missing_gene = pd.DataFrame([set(df.index) - set(d_1000G.keys())]).T
    missing_gene.to_csv('../result/' + filename + '_genewithoutresult.csv', header=False, index=False)  # write the gene without result
    for i in set(df.index) & set(d_1000G.keys()):
        health_bin = d_1000G[i]
        cancer_vector = df.loc[df.index == i].values.T
        cancer_bin = [0] * 100
        for j in cancer_vector:
            cancer_bin[math.ceil(j * 100) - 1] += 1
        cancer_bin = [x / len(cancer_vector) for x in cancer_bin]
        move = 0
        uemd = 0
        for k in reversed(range(100)):
            move += cancer_bin[k] - health_bin[k]
            uemd += max(0, move)
        d_uemd[i] = uemd
    print('done!')
    out = pd.DataFrame([d_uemd]).T
    out.columns = ['uEMD']
    out.sort_values('uEMD', ascending=0, inplace=True)
    if int(top) >= out.shape[0]:
        print('top gene number higher than measured genes, output all records')
        out.to_csv('../result/' + filename + '_result_all.csv')
    else:
        print('output ' + top + ' genes')
        out = out.head(int(top))
        out.to_csv('../result/' + filename + '_result_' + top +'.csv')

getuemd(inputfile_path, top)


