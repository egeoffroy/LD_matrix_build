#We need to figure out how to make LD file for GWAS dataset

if(!require("remotes"))
   install.packages("remotes") # if necessary
library(remotes)
install_github("chr1swallace/coloc")

GWAS_dataset=list(pvalues=GWAS_data$pval,N=36000,s=0.2,MAF=GWAS_data$minor_allele_frequency,beta=GWAS_data$beta,varbeta=GWAS_data$se,type="cc",snp=GWAS_data$rsid,LD=NULL)
eQTL_dataset=list(eQTL_df$pval_nominal,N=1000,MAF=eQTL_df$maf,beta=eQTL_df$slope,varbeta=eQTL_df$slope_se,type="quant",snp=eQTL_df$rs_id_dbSNP151_GRCh38p7,LD=NULL)
col_analysis=coloc.signals(dataset1 = GWAS_dataset,dataset2 = eQTL_dataset,method ="mask",mode="iterative")
