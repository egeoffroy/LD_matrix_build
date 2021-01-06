#!/bin/bash

for i in {1..22}
do 
	nohup Rscript run_pull_snps_driving.R $i > nohup_10kb_chrom${i}.out &
done
