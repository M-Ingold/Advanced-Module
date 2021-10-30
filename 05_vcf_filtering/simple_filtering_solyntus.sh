#!/bin/bash

PATH_IN=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_solyntus

# Filtering out non-SNP variants
vcffilter -s -f "TYPE = snp" $PATH_IN/freebayes_192_samples_09_solyntus.vcf > $PATH_IN/freebayes_192_samples_09_SNPs.vcf

# Filtering out samples that weren't phenotyped
bcftools view --samples-file samples.txt $PATH_IN/freebayes_192_samples_09_SNPs.vcf > $PATH_IN/freebayes_135_samples_09_SNPs.vcf

# Filtering for biallelic SNPs
bcftools view -m2 -M2 $PATH_IN/freebayes_135_samples_09_SNPs.vcf > $PATH_IN/freebayes_135_samples_09_biallelic_SNPs.vcf

# Filtering out sites with variant call quality < 30
vcftools --vcf $PATH_IN/freebayes_135_samples_09_biallelic_SNPs.vcf --minQ 30 --recode --recode-INFO-all \
--out $PATH_IN/freebayes_135_samples_09_SNPs_QUAL_30

# Filtering by read depth of all samples per SNP: 30 ≤ Depth ≤ 10 × d
# First determine mean sequencing depth d
vcftools --vcf $PATH_IN/freebayes_135_samples_09_SNPs_QUAL_30.recode.vcf --site-mean-depth
d=$(grep -v 'nan' out.ldepth.mean | awk '{ total += $3 } END { print total/NR }')
d10=$(echo "$d *10" | bc -l)
# Actual filtering by depth
vcftools --vcf $PATH_IN/freebayes_135_samples_09_SNPs_QUAL_30.recode.vcf --min-meanDP 30 --max-meanDP $d10 --recode --recode-INFO-all \
--out $PATH_IN/freebayes_135_samples_09_SNPs_QUAL_30_depth

# Filtering out samples with depth < 62
gatk VariantFiltration -V $PATH_IN/freebayes_135_samples_09_SNPs_QUAL_30_depth.recode.vcf -filter "DP > 60" --filter-name "MINDEPTH" \
-O $PATH_IN/freebayes_135_samples_09_SNPs_QUAL_30_depth_depth_gt_61.vcf

# Filtering out SNPs <90% genotyped
vcftools --vcf $PATH_IN/freebayes_135_samples_09_SNPs_QUAL_30_depth_depth_gt_61.vcf --max-missing 0.9 --recode --recode-INFO-all \
--out $PATH_IN/freebayes_135_samples_09_SNPs_QUAL_30_depth_depth_gt_61_missing_0.9.vcf

# At least one read per strand
vcffilter -s -f "SAF > 0 & SAR > 0" $PATH_IN/freebayes_135_samples_09_SNPs_QUAL_30_depth_depth_gt_61_missing_0.9.vcf.recode.vcf \
> $PATH_IN/freebayes_135_samples_09_SNPs_QUAL_30_depth_depth_gt_61_missing_0.9_1_read.vcf

# Minor Allele called at least once
bcftools view --min-ac=1:minor $PATH_IN/freebayes_135_samples_09_SNPs_QUAL_30_depth_depth_gt_61_missing_0.9_1_read.vcf > $PATH_IN/solyntus_filtered.vcf
