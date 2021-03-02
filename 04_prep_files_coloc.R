# This is where the pQTL and GWAS SS files will be formatted for coloc v2 
# Author: Elyse Geoffroy

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(R.utils))

"%&%" = function(a,b) paste(a,b,sep="") #'WBC', 'Platelet'
#phenos <- c('Total_cholesterol', 'HDL_cholesterol', 'LDL_cholesterol', 'BMI', 'Waist-hip50', 'Waist-hip51', 'Waist-hip52', 'Chronic_kidney', 'End_renal', 'Diastolic_blood', 'Systolic_blood', 'Diabetes', 'Triglyceride', 'Mean_corpuscular_hemoglobin', 'Hemoglobin', 'Height', 'Glomerular', 'Fasting_insulin', 'Fasting_glucose', 'PR_interval', 'QRS_duration', 'QT_interval', 'Coffee', 'Smoking') # make this into a parameter for pipeline
#phenos <- c('PR_interval', 'QRS_duration', 'QT_interval', 'Coffee', 'Smoking')
phenos <- c('BMI')
chrs <- c(1:22)
#pops <- c('AFA', 'ALL', 'CAU', 'CHN', 'HIS')
#pop_sample_size <- c(183, 971, 416, 71, 301)
pops <- c( 'CAU')
pop_sample_size <-  c( 416)
for(pop in 1:length(pops)){ #read in pop's .frq file for MAF
  # some file locations may change
  frq <- fread(paste("/home/egeoffroy/LD_matrix/", pops[pop], "_prot_hg38.frq", sep = ''))
  frq <- frq %>% dplyr::select(SNP, MAF)

  for(pheno in phenos){ #read in GWAS output file
    if(!file.exists("/home/egeoffroy/LD_matrix/coloc/GWAS_TOPMED" %&% pops[pop] %&% "_" %&% pheno %&% ".txt") & !file.exists("/home/egeoffroy/LD_matrix/coloc/pQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt")){
    GWAS_result <- fread(paste("/home/egeoffroy/Wojcik/Wojcik_build38/gzipped_versions/WojcikG_", pheno, ".txt.gz", sep = '') , header = T)
    GWAS_result$chr_pos <- paste(gsub("chr", "", GWAS_result$chromosome), GWAS_result$base_pair_location, sep = ":")
    GWAS_for_COLOC <- GWAS_result %>% dplyr::select(SNP_hg38, Beta, SE, `Effect-allele-frequency`, `Sample-size`) #subset to COLOC input

    colnames(GWAS_for_COLOC) <- c("panel_variant_id", "effect_size", "standard_error", "frequency", "sample_size")
    GWAS_for_COLOC <- GWAS_for_COLOC[complete.cases(GWAS_for_COLOC),] #COLOC does not like missing values
    #GWAS_write <- data.frame(panel_variant_id = character(), effect_size = numeric(), standard_error = numeric(), frequency = numeric(), sample_size = numeric(), stringsAsFactors = F)
    GWAS_write <- GWAS_for_COLOC

    pQTL_write <- data.frame(gene_id = character(), variant_id = character(), maf = numeric(), pval_nominal = numeric(), slope = numeric(), slope_se = numeric(), stringsAsFactors = F)

    for(chr in chrs){ # it may be more useful to also have chromosome as a parameter for the overall pipeline
      mpqtl <- fread("/home/egeoffroy/LD_matrix/cis_eQTLs_" %&% pops[pop] %&% "_WG_all_cis.txt.gz", nThread = 40) #read in matrix eQTL results
      mpqtl$se <- mpqtl$beta / mpqtl$statistic #make your own standard error since it's not in the meQTL output
      mpqtl$n_samples <- pop_sample_size[pop]
      mpQTL_for_COLOC <- left_join(mpqtl, frq, by = c("snps" = "SNP")) #add freq to COLOC input
      mpQTL_for_COLOC <- mpQTL_for_COLOC %>% dplyr::select(gene, snps, MAF, pvalue, beta, se) #subset to COLOC input
      colnames(mpQTL_for_COLOC) <- c("gene_id", "variant_id", "maf", "pval_nominal", "slope", "slope_se")
      mpQTL_for_COLOC <- mpQTL_for_COLOC[complete.cases(mpQTL_for_COLOC),]

      #GWAS_write <- rbind(GWAS_write, GWAS_for_COLOC_chr)
      pQTL_write <- rbind(pQTL_write, mpQTL_for_COLOC)
    }

    snps_in_both <- intersect(GWAS_write$panel_variant_id, pQTL_write$variant_id) #is there a better way to do this? Probably. 
    GWAS_write <- subset(GWAS_write, panel_variant_id %in% snps_in_both)
    print(head(GWAS_write))
    pQTL_write <- subset(pQTL_write, variant_id %in% snps_in_both)
    pQTL_write <- pQTL_write[order(pQTL_write$gene_id),]

#    if(!file.exists("/home/egeoffroy/LD_matrix/coloc/pQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt")){
    # write out the files that will be used to run coloc
    	fwrite(pQTL_write, "/home/egeoffroy/LD_matrix/coloc/pQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", quote = F, sep = "\t", na = "NA", row.names = F, col.names = T)
    	gzip("/home/egeoffroy/LD_matrix/coloc/pQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", destname = "/home/egeoffroy/LD_matrix/coloc/pQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt.gz") #script may only take .gz values so can't hurt to be too careful
    	fwrite(GWAS_write, "/home/egeoffroy/LD_matrix/coloc/GWAS_TOPMED" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
    	gzip("/home/egeoffroy/LD_matrix/coloc/GWAS_TOPMED" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", "/home/egeoffroy/LD_matrix/coloc/GWAS_TOPMED" %&% pops[pop] %&% "_" %&% pheno %&% ".txt.gz")
    	print("Completed with " %&% pops[pop] %&% ", for " %&% pheno %&% ".")
    } else{
	next
    }
  }
}
