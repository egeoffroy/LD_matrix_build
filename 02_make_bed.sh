#!/bin/bash

pop=$1
for chr in {7..9}
do
	./plink --vcf /home/egeoffroy/MESA_TOPMED_Imputation/${pop}/chr${chr}.dose.vcf --make-bed --out /home/egeoffroy/LD_matrix/${pop}_chr${chr}_dose
done
