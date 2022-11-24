#!/bin/bash
while read -r sample
do
	echo $sample
	sbatch ./code/trim.job $sample
	echo ${sample}_submitted
done
