#!/bin/bash
while read -r sample
do
	echo $sample
	sbatch ./code/fastqc.job $sample
	echo ${sample}_submitted
done
