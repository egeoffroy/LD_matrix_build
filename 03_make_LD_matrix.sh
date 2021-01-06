#!/bin/bash
# This file makes the LD matrices for colocalization using plink with the dosage bed/bim/fam files and a list of sig snps to include inthe matrix for each gene. 
pop=$1 #first argument is the population abbreviation

for chr in {1..22} #for chromosomes 1 through 22 --> may want to break up in chunks
do
	for file in $"/home/egeoffroy/LD_matrix/5kb_of_gene"/${pop}_chr_${chr}*
	do
		filename="${file##*/}"
		filename="${filename%%.*}"
		echo ${filename}
		./plink --bfile ${pop}_chr${chr}_dose --r bin yes-really --extract ${file} --out ${filename}_LD
	done
done
