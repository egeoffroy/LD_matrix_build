library(dplyr)
library(data.table)
source('/home/rschubert1/scratch/test_coloc_for_elyse/05_coloc.R')


args <- commandArgs(trailingOnly=TRUE)
print(args)
pop <- args[1]
phenotype <- args[2]
eqtl <- args[3]
gwas<- args[4]
LD_dir<- args[5]
gene_id<-args[6]
pop_size<-as.numeric(args[7])
ld_matrix<-args[8]
gwas_size<-as.numeric(args[9])


genes<-readRDS(gene_id)

main(eqtl=eqtl, 
     gwas=gwas, 
     directory=LD_dir,  
     mode = 'bse', 
     gene_list=genes, 
     eqtlGeneCol='gene_id', 
     eqtlSNPCol='variant_id', 
     eqtlMAFCol='maf', 
     eqtlSeCol=6, 
     eqtlBetaCol=5, 
     eqtlSampleSize=pop_size, 
     gwasSNPCol=1, 
     LD=ld_matrix, 
     method="cond", 
     gwasBetaCol=2, 
     gwasSeCol=3, 
     gwasSampleSize=gwas_size, 
     outFile=paste('/home/rschubert1/scratch/test_coloc_for_elyse/', pop, '_', phenotype, '.txt', sep = ''), ld_header = 'T')

