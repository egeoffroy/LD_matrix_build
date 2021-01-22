suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
source('/home/egeoffroy/LD_matrix/05_coloc.R')

phenotype <- 'Total_cholesterol'
ld_matrix <- '/home/egeoffroy/LD_matrix/5kb_of_gene/AFA_chr_19_SL000276_ENSG00000130203_5kb_LD.ld.gz'
pops <- c('AFA', 'ALL', 'CAU', 'CHN', 'HIS')
pop_sample_sizes <- c(183, 971, 416, 71, 301)

for(pop in pops){
        for(pop_size in pop_sample_sizes){
                eqtl <- paste('/home/egeoffroy/LD_matrix/coloc/pQTL_', pop, '_Total_cholesterol.txt.gz', sep = '')
                gwas <- paste('/home/egeoffroy/LD_matrix/coloc/GWAS_TOPMED', pop, '_Total_cholesterol.txt.gz', sep = '')
                F_gwas <- fread(gwas, header = T, stringsAsFactors=F)
                gwas_size <- F_gwas$`sample_size`[1]
                print(gwas_size)
                main(eqtl, gwas, gene='SL000276_ENSG00000130203', eqtlGeneCol= 'gene_id', eqtlMAFCol='maf', eqtlSeCol='se', eqtlBetaCol='slope', eqtlSampleSize= pop_size,
                        gwasSNPCol='panel_variant_id',gwasMAFCol='frequency', LD=ld_matrix,
                        gwasBetaCol='effect_size', gwasSeCol='standard_error', gwasSampleSize=gwas_size)
        }
}
