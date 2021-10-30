#!/bin/bash

#loop through vcf-files counting SNPs
VCF_PATH=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS/
TARGET_FILE=SNP_count_6.1.txt

touch $TARGET_FILE
> $TARGET_FILE

for vcfFile in $(ls -tr $VCF_PATH/*.vcf); do
	name=$(basename $vcfFile)
	number=$(grep -cv '#' $vcfFile)
	echo -e $name '\t' $number '\t' >> $TARGET_FILE
	max=$(sed -r 's/.* ([0-9]+\.*[0-9]*).*?/\1/' $TARGET_FILE | head -1)
	percentage=$(echo "$number/$max*100" | bc -l | head -c 4)
	echo -e $percentage >> $TARGET_FILE

done
