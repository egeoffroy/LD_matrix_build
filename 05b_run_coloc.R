library(dplyr)
library(data.table)
source('/home/egeoffroy/LD_matrix/05_coloc.R')

# should have an input list of the genes to test in which populations and which chrom they are on?

coloc_analysis <- function(gene_id=NULL, pop=NULL, pop_size=NULL, phenotype=NULL){
	gene <- gsub("\\..*","",gene_id)
	print(gene)

	eqtl <- paste('/home/egeoffroy/LD_matrix/coloc/pQTL_', pop, '_', phenotype, '.txt.gz', sep = '')
	gwas <- paste('/home/egeoffroy/LD_matrix/coloc/GWAS_TOPMED', pop, '_', phenotype, '.txt.gz', sep = '')
	
        ld_matrix <- 'T'
		if(!is.null(ld_matrix)){
					F_gwas <- fread(gwas, header = T, stringsAsFactors=F)
					print(F_gwas)
					gwas_size <- F_gwas$`sample_size`[1]
					print(gwas_size)
					if(pop == 'CAU'){ # for some reason CAU uses pvalues instead of bse values
#						main(eqtl=eqtl, gwas=gwas, directory=paste('/home/egeoffroy/LD_matrix/', pop, '_1Mb_coords_LDMatrix', sep =''),  mode = 'p', gene_list=gene_id, eqtlGeneCol='gene_id', eqtlSNPCol='variant_id', eqtlMAFCol='maf', eqtlSeCol=6, eqtlBetaCol=5, eqtlSampleSize= pop_size, gwasSNPCol=1, LD=ld_matrix, method="cond", gwasBetaCol=2, gwasSeCol=3, gwasSampleSize=gwas_size, outFile=paste('/home/egeoffroy/LD_matrix/coloc_output/', pop, '/', pop, '_', phenotype, '.txt', sep = ''), ld_header = 'T')
						eqtl1 <- fread(eqtl, header = T, stringsAsFactors=F)
						gwas1 <- fread(gwas, header = T, stringsAsFactors=F)
						eqtl1$variant_id <- str_replace_all(eqtl1$variant_id, 'chr', '')
						gwas1$panel_variant_id <- str_replace_all(gwas1$panel_variant_id, 'chr', '')
						write.table(eqtl1, eqtl, quote= F, row.names=F)
						write.table(gwas1, gwas, quote=F, row.names=F)
					} 
						main(eqtl=eqtl, gwas=gwas, directory=paste('/home/egeoffroy/LD_matrix/', pop, '_1Mb_coords_LDMatrix', sep =''),  mode = 'bse', gene_list=gene_id, eqtlGeneCol='gene_id', eqtlSNPCol='variant_id', eqtlMAFCol='maf', eqtlSeCol=6, eqtlBetaCol=5, eqtlSampleSize= pop_size, gwasSNPCol=1, LD=ld_matrix, method="cond", gwasBetaCol=2, gwasSeCol=3, gwasSampleSize=gwas_size, outFile=paste('/home/egeoffroy/LD_matrix/coloc_output/', pop, '/', pop, '_', phenotype, '.txt', sep = ''), ld_header = 'T')
#					}
				}
}
