#!/bin/bash

pop=$1
for chr in {1..22}
do
	./plink --vcf /home/egeoffroy/MESA_TOPMED_Imputation/${pop}/chr${chr}.dose.vcf.gz --make-bed --out /home/egeoffroy/LD_matrix/${pop}/${pop}_chr${chr}_dose
done
