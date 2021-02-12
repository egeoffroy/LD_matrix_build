#!/bin/bash
pop=$1

for chr in {1..22}
do
	./plink --vcf /home/egeoffroy/LD_matrix/${pop}_chr${chr}.dose.vcf --make-bed --out /home/egeoffroy/LD_matrix/${pop}_chr${chr}_dose
done
