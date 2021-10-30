#!/bin/bash

################################################################
#                                                              #
# Master script for applying all selected filters sequencially #
#                                                              #
################################################################

# In and out for SNP filtering
IN_SNP=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_tetraploid/freebayes_192_samples_chr09_tetraploid.vcf
OUT_SNP=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_tetraploid/freebayes_192_samples_chr09_tetraploid_SNPs.vcf

# In and out for phenotype filtering
IN_PHENO=$OUT_SNP
OUT_PHENO=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_tetraploid/freebayes_135_samples_chr09_tetraploid_SNPs.vcf

# In and ouput for biallelic SNP filtering
IN_BIALLEL=$OUT_PHENO
OUT_BIALLEL=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_tetraploid/freebayes_135_samples_chr09_tetraploid_biallelic_SNPs.vcf

# In and output for QUAL filtering
IN_QUAL=$OUT_BIALLEL
OUT_QUAL=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_tetraploid/freebayes_135_samples_chr09_tetraploid_biallelic_SNPs_QUAL_30.vcf

IN_DEPTH=$OUT_QUAL
OUT_DEPTH=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_tetraploid/freebayes_135_samples_chr09_tetraploid_biallelic_SNPs_QUAL_30_depth.vcf

IN_BLANK=$OUT_DEPTH
OUT_BLANK=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_tetraploid/freebayes_135_samples_chr09_tetraploid_biallelic_SNPs_QUAL_30_depth_0.9_blanked.vcf

IN_BIAS=$OUT_BLANK
OUT_BIAS=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_tetraploid/freebayes_135_samples_chr09_tetraploid_biallelic_SNPs_QUAL_30_depth_0.9_blanked_1_read.vcf

IN_HET=$OUT_BIAS
OUT_HET=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_tetraploid/tetraploid_filtered.vcf

# Filtering for SNPs
#vcffilter -s -f "TYPE = snp" $IN_SNP > $OUT_SNP

# Only phenotyped data, samples.txt generated using R
#bcftools view --samples-file samples.txt $IN_PHENO > $OUT_PHENO

# Only biallelic SNPs
touch $OUT_BIALLEL
python3 filter_VCF_GWAS_biallelic_SNPs.py $IN_BIALLEL $OUT_BIALLEL

# QUAL >= 30
sh filter_VCF_GWAS_QUAL.sh $IN_QUAL $OUT_QUAL

# Create statistics for average depth for subsequent filtering of max depth
#IN_PATH_SUBSET=../../data/vcf_GWAS_input
#IN_VCF_FILE=$IN_PATH_SUBSET/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30.vcf
#OUT_STAT_PATH=../../analysis/variant_calling_GWAS_SNPS
#mkdir -p $OUT_STAT_PATH
#OUT_STAT_NAME=$OUT_STAT_PATH/SNPS_135-samples_Biall_SNPs_QUAL-30
# Extracting information from the subsetted SNP VCF file
vcftools --vcf $OUT_QUAL --site-mean-depth
#vcftools --vcf $IN_VCF_FILE --depth --out $OUT_STAT_NAME
#vcftools --vcf $IN_VCF_FILE --site-quality --out $OUT_STAT_NAME


# Depth filter
touch $OUT_DEPTH
python3 filter_VCF_GWAS_depth_per_site.py $IN_DEPTH $OUT_DEPTH


# Blanking "bad" genotypes and filtering variants with too few genotyped samples 
touch $OUT_BLANK
python3 filter_VCF_GWAS_Missing.py $IN_BLANK $OUT_BLANK


# Min. 1 read per strand per variant
sh filter_VCF_GWAS_StrandBias.sh $IN_BIAS $OUT_BIAS

# Min. 1 ALT and 1 REF allele in the population at this site
touch $OUT_HET
python3 filter_VCF_GWAS_HET.py $IN_HET $OUT_HET