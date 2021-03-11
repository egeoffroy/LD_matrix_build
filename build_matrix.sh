#!/bin/bash
pop=$1
chr=$2
./plink --vcf /home/egeoffroy/LD_matrix/${pop}_chr${chr}.dose.vcf --make-bed --out /home/egeoffroy/LD_matrix/${pop}_chr${chr}_dose
./plink --bfile ${pop}_chr${chr}_dose --r bin yes-really --extract /home/egeoffroy/LD_matrix/AFA_SL000437_ENSG00000257017.9_sig_gene_snps.txt --out ${pop}_chr${chr}_LD  
