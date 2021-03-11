#!/bin/bash
# Author: Elyse Geoffroy
# This file makes the LD matrices for colocalization using plink with the dosage bed/bim/fam files and a list of sig snps to include inthe matrix for each gene. 
pop=$1 #first argument is the population abbreviation

for chr in {1..22} #{1..22} #for chromosomes 1 through 22 --> may want to break up in chunks
do
	for file in $"/home/egeoffroy/LD_matrix/1Mb_of_gene_coords"/${pop}_chr_${chr}_*
	do
		filename="${file##*/}"
		#filename="${filename%%.*}"
		substring=('_1Mb_of_gene.txt')
		filename=${filename%"${substring}"}
		echo ${filename}
		./plink --bfile /home/rschubert1/data/TOPMED_Proteome/${pop}/PCAIR_modeling/01genotypes/dosages/hg38chr${chr}.dosage --r square gz yes-really --extract ${file} --write-snplist --out ${pop}_1Mb_coords_LDMatrix/${filename}_1Mb_LD
	done
done
