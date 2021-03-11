library(dplyr)
library(data.table)
library(stringr)

files <- list.files('/home/egeoffroy/LD_matrix/1Mb_of_gene_coords', pattern = 'ALL', full.names=T)
for(file in files){
	file1 <- fread(file, header = F, stringsAsFactors=F, sep = ':')
	file1 <- paste(file1$V1, file1$V2, sep = ':')
	print(head(file1))
	write.table(file1, file, row.names=F, col.names=F, quote=F)
}
