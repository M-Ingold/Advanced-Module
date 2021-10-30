#!/bin/bash

#loop through vcf-files counting SNPs
VCF_PATH=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_solyntus
TARGET_FILE=SNP_count_solyntus.txt

touch $TARGET_FILE
> $TARGET_FILE

max=$(grep -cv '#' $VCF_PATH/freebayes_192_samples_09_solyntus.vcf)

for vcfFile in $(ls -tr $VCF_PATH/*.vcf); do
	name=$(basename $vcfFile)
	number=$(grep -cv '#' $vcfFile)
	percentage=$(echo "$number/$max*100" | bc -l | head -c 4)
	echo -e $name '\t' $number '\t' $percentage >> $TARGET_FILE
	#max=$(sed -r 's/.* ([0-9]+\.*[0-9]*).*?/\1/' $TARGET_FILE | head -1)

done
