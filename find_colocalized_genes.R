library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)

files <- list.files(path = '/home/rschubert1/scratch/test_coloc_for_elyse', pattern = '.txt', full.names=T)
print(files)

coloc_combined_output <- data.frame()

for(file in files){
  if(!str_detect(file, 'QTL_effects') & !str_detect(file, 'GWAS_effects')){
	if(str_detect(file, 'AFA')){pop <- 'AFA'}
        if(str_detect(file, 'ALL')){pop<-'ALL'}
        if(str_detect(file, 'CAU')){pop<-'CAU'}
        if(str_detect(file, 'CHN')){pop<-'CHN'}
        if(str_detect(file, 'HIS')){pop<-'HIS'}
        print(pop)

	pheno <- str_remove(file, '/home/rschubert1/scratch/test_coloc_for_elyse')
	pheno <- str_remove(pheno, '.txt')
	pheno <- str_remove(pheno, paste('/', pop, '_', sep =''))
	print(pheno)

	file1 <- fread(file, header = T, stringsAsFactors=F)
	print(head(file1))
	file1$Population <- rep(pop, nrow(file1))
	file1$Phenotype <- rep(pheno, nrow(file1))

	coloc_combined_output <- rbind(coloc_combined_output, file1)
   }
}

print(coloc_combined_output)
write.table(coloc_combined_output, '/home/egeoffroy/LD_matrix/coloc_combined_output.txt', row.names=F, quote=F)



# Make p4 distribution figures
for(pheno in unique(coloc_combined_output$Phenotype)){
	sub <- coloc_combined_output %>% filter(Phenotype == pheno)
	p <- ggplot(sub, aes(Population, PP.H4.abf, fill = Population)) + geom_violin() + ggtitle(str_replace_all(pheno, '_', ' ')) + xlab("Topmed Population") + ylab("COLOC P4 Distribution") + ylim(0,1)+ ggsave(paste('TOPMED_', pheno, '_violin_p4.png', sep = ''))
}



#Colocalized genes 
P4 <- coloc_combined_output %>% filter(PP.H4.abf > 0.5)
print(P4)
write.table(P4, '/home/egeoffroy/LD_matrix/colocalized_topmed.txt', row.names=F, quote=F)
P4_pairs <- P4 %>% select(gene,Phenotype)
print(unique(P4_pairs))


P3 <- coloc_combined_output %>% filter(PP.H3.abf > 0.5)
write.table(P3, '/home/egeoffroy/LD_matrix/independent_topmed.txt', row.names=F, quote=F)
