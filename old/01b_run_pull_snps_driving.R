setwd("/home/egeoffroy/LD_matrix/")
source("01_pull_snps_driving.R")
"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
chrom <- argv[1]
pops <- c('AFA', 'HIS', 'ALL', 'CAU', 'CHN')
pip <- '0.001'

for(pop1 in pops){
if ( pop1 == "ALL"){
        expr="/data/rschubert1/TOPMED_Proteome/ALL/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_ALL_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids_sorted.txt.gz"
} else {
        expr="/data/rschubert1/TOPMED_Proteome/" %&% pop1 %&% "/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_" %&% pop1 %&% "_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids.txt.gz"
}
snp_annot_file <- "/home/rschubert1/data/TOPMED_Proteome/" %&% pop1 %&% "/PCAIR_modeling/01genotypes/preddb_input/uniq_pred_db_hg38.chr" %&% chrom %&% ".maf0.01.R20.8.anno.txt.gz"
gene_annot_file <- "/data/rschubert1/TOPMED_Proteome/annotation_all_aptamers_ENSG.txt"
genotype_file <- "/home/rschubert1/data/TOPMED_Proteome/" %&% pop1 %&% "/PCAIR_modeling/01genotypes/preddb_input/uniq_pred_db_hg38.chr" %&% chrom %&% ".maf0.01.R20.8.geno.txt.gz"
expression_file <- expr
protein_file <- "/data/rschubert1/TOPMED_Proteome/" %&% pop1 %&% "/PCAIR_modeling/05dapg/summary_dapg_out/summary_snps.txt"
prefix <- "/home/rschubert1/data/TOPMED_Proteome/" %&% pop1 %&% "/PCAIR_modeling/06Elastic_net/dapg_" %&% pip %&% "_T_out/" %&% pop1 %&% "_chr" %&% chrom %&% "_" %&% "PIP_" %&% pip %&% "_clus_filt_T"
#"/home/ryan/topmed/multiomic_modeling/output/CAU_PBMC_chr" %&% chrom %&% "_multiomic_models.txt"
#prefix <- "/home/ryan/topmed/proteome/dapg_net/redo_AFA/PIP_" %&% pip %&% "_" %&% suffix %&% "/" %&% pop1 %&% "/" %&% pop1 %&% "_chr" %&% chrom %&% "_" %&% "PIP_" %&% pip %&% "_" %&% suffix

main(snp_annot_file=snp_annot_file,
        gene_annot_file,
        genotype_file=genotype_file,
        expression_file=expression_file,
        pop=pop1,
        chrom=as.numeric(chrom),
        prefix=prefix,
        null_testing=FALSE,
        protein_file=protein_file)
}
