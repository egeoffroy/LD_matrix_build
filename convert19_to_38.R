# This script removes the alleles and the chr from hg38 to hg19 for CAU
# Author: Elyse Geoffroy
library(data.table)
library(stringr)

all_files <- list.files('/home/egeoffroy/LD_matrix/1Mb_of_gene_coords', pattern = 'CAU', full.names=T)
print(head(all_files))

for(file in all_files){
	file1 <- fread(file, header = F, stringsAsFactors=F, sep = ':')
	file1$V1 <- str_replace_all(file1$V1, 'chr', '')
	file1 <- paste(file1$V1, file1$V2, sep = ':')
	print(head(file1))
	write.table(file1, file, row.names=F, quote=F, col.names=F)
}

