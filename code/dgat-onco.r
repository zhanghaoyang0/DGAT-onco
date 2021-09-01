#=====================================================================================
# annovar
#=====================================================================================
## annovar
file_in="./example/BRCA_for_annovar.vcf"
file_out='./example/BRCA'
path_annovar='/home/yanglab_data/user/zhanghy/soft_slurm/annovar/'  

$path_annovar\table_annovar.pl \
$file_in \
$path_annovar\humandb/ \
-buildver hg19 \
-out $file_out \
-remove -protocol refGene,dbnsfp35c -operation g,f -nastring . -vcfinput -polish
#=====================================================================================
# load functions and packages
#=====================================================================================
library(dplyr)
library(tidyr)

make_vcf_for_annovar = function(chr_col = 'Chromosome', pos_col = 'Start_Position', ref_col = 'Tumor_Seq_Allele1', 
                                alt_col = 'Tumor_Seq_Allele1', id_col = 'Tumor_Sample_Barcode', 
                                file_in = './example/BRCA.maf', file_out = './example/BRCA_for_annovar.vcf'){
  cat(paste0('loading input ', file_in, ' ...\n'))
  df = read.table(file_in, header=1)
  cat('making vcf from maf ...\n')
  vcf = df
  vcf$dot = '.'
  vcf = vcf[,c(chr_col, pos_col, 'dot', ref_col, alt_col, 'dot', 'dot', id_col)]
  write.table(vcf, file_out, sep='\t', row.names=F, col.names=F, quote=F)
  cat('done!\n')
  cat('saving output ...\n')
  cat(paste0('done!\noutput saved in ', file_out, '\n'))
}

calculate_summation_score = function(file_in = './example/BRCA.hg19_multianno.txt', id_col = 'Otherinfo11', 
                              gene_col = 'Gene.refGene', func_col= 'ExonicFunc.refGene'){
  cat(paste0('loading input ', file_in, ' ...\n'))
  df = read.table(file_in, header=1, sep='\t', quote='', check.names=F)
  df = df[df[,func_col]=='nonsynonymous SNV',]
  cat('done!\n')
  
  if (sum(colnames(df)%in%all_of(c(id_col, gene_col, func_col)))!=3){
    stop(paste0('at least one  of ', id_col, ', ', gene_col, ', or', func_col, ' columns were not in your input file!'))
  }
  
  if (nrow(df)==0){
    stop('not nonsynonymous SNV in your input file!')
  }
  
  df = df%>%rename(id=all_of(id_col), gene=all_of(gene_col))
  df = df[!grepl(';', df$gene),] # drop inter-genic variant
  df$constant_rankscore = 1
  score_col = colnames(df)[grepl('rank', colnames(df))]
  df = df[,c('id', 'gene', score_col)]
  df[,score_col][df[,score_col]=='.'] = '0'
  df[,score_col] = sapply(df[,score_col], as.numeric)
  
  cat('scoring function in your input:\n\n')
  cat(paste0(score_col, collapse='\n'))
  cat('\n\nnote that "constant_rankscore" means every mutation has a equal weight\n')
  
  cat('\ncalculating mutation scores for each individual ...\n
      ')
  summ = aggregate(.~id+gene, df, FUN=sum)
  cat('\ndone!\n')
  return(summ)
}

calculate_uemd = function(summ, file_out = './example/BRCA_result.txt',
                          score = 'MutPred_rankscore'){
  cat(paste0('calculating uemd using ', score, ' scoring function ... \n'))
  mat = summ[,c('id', 'gene', score)]
  mat = spread(mat, gene, score)
  mat$id = NULL # drop id
  mat[is.na(mat)] = 0
  bins_cancer = apply(t(apply(mat, 1, fastRank)), 2, bins)
  colnames(bins_cancer) = colnames(mat)
  
  append_genes = setdiff(colnames(bins_cancer), colnames(bins_1000g[[score]])) # if cancer gene not in 1000g, assume a distribution as c(1,rep(0,99))
  bins_append = matrix(rep(c(1,rep(0,99)), length(append_genes)), nrow=100, byrow=F)
  colnames(bins_append) = append_genes
  bins_1000g_append = cbind(bins_1000g[[score]], bins_append)
  uemd = sapply(colnames(bins_cancer), function(x) uEMD(bins_cancer[,x], bins_1000g_append[,x]))
  
  res = data.frame(gene=colnames(bins_cancer), uemd=uemd)
  row.names(res) = NULL
  write.table(res, file_out, row.names=F, sep='\t', quote=F)
  cat(paste0('done!\nresult saved in ', file_out, '\n'))
}

# fastRank, bins, uEMD functions were from diffmut
fastRank = function(x) { 
  x[x!=0] = rank(x[x!=0], ties.method="min")+length(which(x==0)) 
  as.numeric(x)/length(x) 
}

bins = function(v, p=100) {
  l = length(v)
  nBins = rep(0,p)
  for(val in v){
    nBins[ceiling(val*(p-1))+1] = nBins[ceiling(val*(p-1))+1]+1
  }
  nBins = nBins/l
  nBins
}

uEMD = function(tBins, nBins, p=100) {
  sum = 0
  move = 0
  for(i in p:2){
    move = move+tBins[i]-nBins[i]
    sum = sum+max(move,0)
  }
  sum
}

#=====================================================================================
# example
#=====================================================================================
## load functional impact-based distribution of natural population (1000 genomics[1000g]) 
file_bin_1000g = './data/bins_1000g.rdata'
load(file_bin_1000g)
cat('scoring function in the pre-calculated 1000g bins:\n\n')
cat(paste0(names(bins_1000g), collapse='\n'))
cat('\n\nnote that "constant_rankscore" means every mutation has a equal weight\n')

## make vcf for annovar
make_vcf_for_annovar(chr_col = 'Chromosome', pos_col = 'Start_Position', ref_col = 'Tumor_Seq_Allele1', 
                     alt_col = 'Tumor_Seq_Allele2', id_col = 'Tumor_Sample_Barcode', 
                     file_in = './example/BRCA.maf', file_out = './example/BRCA_for_annovar.vcf')

## after annovar, calculate mutation summation score for each cancer individual
summ = calculate_summation_score(file_in = './example/BRCA.hg19_multianno.txt', id_col = 'Otherinfo11', 
                                 gene_col = 'Gene.refGene', func_col= 'ExonicFunc.refGene')

## calculate uemd
calculate_uemd(summ, file_out = './example/BRCA_result.txt', score = 'MutPred_rankscore')
