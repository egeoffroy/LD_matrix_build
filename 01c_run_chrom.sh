#!/bin/bash

for i in {1..22}
do 
	nohup Rscript 01b_run_pull_snps_driving.R $i > nohup_1Mb_chrom${i}_2.out &
done
