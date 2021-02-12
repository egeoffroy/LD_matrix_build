#!/bin/bash
# This file makes the LD matrices for colocalization using plink with the dosage bed/bim/fam files and a list of sig snps to include inthe matrix for each gene. 
pop=$1 #first argument is the population abbreviation

for chr in {1..22} #{1..22} #for chromosomes 1 through 22 --> may want to break up in chunks
do
	for file in $"/home/egeoffroy/LD_matrix/10kb_of_gene"/${pop}_chr_${chr}_*
	do
		filename="${file##*/}"
		filename="${filename%%.*}"
		echo ${filename}
		./plink --bfile /home/egeoffroy/LD_matrix/${pop}/${pop}_chr${chr}_dose --r square gz yes-really --extract ${file} --out ${pop}_10kb_LDMatrix/${filename}_10kb_LD
	done
done
