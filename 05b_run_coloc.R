library(dplyr)
library(data.table)
source('/home/egeoffroy/LD_matrix/05_coloc.R')

# should have an input list of the genes to test in which populations and which chrom they are on?
phenotype <- 'C-reactive'
gene_id <- 'SL004783_ENSG00000163221.4' 
gene <- gsub("\\..*","",gene_id)
print(gene)
#pops <- c('AFA', 'ALL', 'CAU', 'CHN', 'HIS')
pops <- c('CHN')
pop_sample_sizes <- c(71)
#pop_sample_sizes <- c(183, 971, 416, 71, 301)

for(pop in pops){
	for(pop_size in pop_sample_sizes){
		eqtl <- paste('/home/egeoffroy/LD_matrix/coloc/pQTL_', pop, '_', phenotype, '.txt.gz', sep = '')
		gwas <- paste('/home/egeoffroy/LD_matrix/coloc/GWAS_TOPMED', pop, '_', phenotype, '.txt.gz', sep = '')
		ld_matrix <- list.files(paste('/home/egeoffroy/LD_matrix/', pop, '_10kb_LDMatrix/', sep = ''), pattern = gene, full.names = T)[[1]][1]
		#ld_matrix <- paste('/home/egeoffroy/LD_matrix/', pop, '_10kb_LDMatrix/', pop, '_chr_', chrom, '_', gene, '_10kb_LD.ld.gz', sep = '')
		print(ld_matrix)
		F_gwas <- fread(gwas, header = T, stringsAsFactors=F)
		gwas_size <- F_gwas$`sample_size`[1]
		print(gwas_size)
		main(eqtl=eqtl, gwas=gwas, mode = 'bse', gene=gene_id, eqtlGeneCol='gene_id', eqtlSNPCol='variant_id', eqtlMAFCol='maf', eqtlSeCol=6, eqtlBetaCol=5, eqtlSampleSize= pop_size, gwasSNPCol=1, LD=ld_matrix, gwasBetaCol=2, gwasSeCol=3, gwasSampleSize=gwas_size, outFile=paste(pop, '_', phenotype, '_', gene, '.txt', sep = ''))
	}
}
