#!/bin/bash
while read -r sample
do
	echo $sample
	sbatch ./code/subread_align.job $sample
	echo ${sample}_submitted
done
