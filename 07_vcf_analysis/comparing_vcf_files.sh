#!/bin/bash

PATH_404=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_4.4/404_lifted_over.vcf
PATH_404_unlifted=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_4.4/404_filtered.vcf
PATH_61=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS/61_filtered.vcf
PATH_SOLYNTUS=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_solyntus/solyntus_lifted_over.vcf
PATH_TETRAPLOID=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_tetraploid/tetraploid_filtered.vcf

:'
# bgzip files for processing
bgzip $PATH_404 > $PATH_404.gz
#bgzip $PATH_61 > $PATH_61.gz
#bgzip $PATH_SOLYNTUS > $PATH_SOLYNTUS.gz
#bgzip $PATH_TETRAPLOID > $PATH_TETRAPLOID.gz
#bgzip $PATH_404_unlifted > $PATH_404_unlifted.gz

bcftools index -f $PATH_404.gz
#bcftools index -f $PATH_61.gz #(cant be indexed?)
#bcftools index -f $PATH_SOLYNTUS.gz

#bcftools
bcftools isec $PATH_404.gz $PATH_61.gz -p isec_404_to_61
#bcftools isec $PATH_404.gz ../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_4.4/404_lifted_over_transanno.vcf.gz -p isec404_404
# -> only 115 SNPs mapped differently between VcfLiftover and Transanno!
#bcftools isec $PATH_SOLYNTUS.gz $PATH_61.gz -p isec_solyntus_to_61
bcftools isec $PATH_SOLYNTUS.gz $PATH_404.gz -p isec_solyntus_404
bcftools isec -n=3 $PATH_SOLYNTUS.gz $PATH_404.gz $PATH_61.gz -p isec_SNPs_in_all

#bcftools index -f $PATH_TETRAPLOID.gz
#bcftools index -f $PATH_404_unlifted.gz
#bcftools isec $PATH_TETRAPLOID.gz $PATH_404_unlifted.gz -p isec_tetraploid_to_404

# F#bcftools isec $PATH_SOLYNTUS.gz $PATH_61.gz -p isec_solyntus_to_61
#touch rename-chrs.txt
#echo "StSOLv1.1ch09 chr09" > rename-chrs.txt
#bcftools annotate --rename-chrs rename-chrs.txt $PATH_SOLYNTUS/solyntus_filtered.vcf > $PATH_SOLYNTUS/solyntus_filtered_chr_renamed.vcf

# SnpSift
#SnpSift concordance -v $PATH_404 $PATH_61 > concordance404_61.txt
#SnpSift concordance -v $PATH_SOLYNTUS/solyntus_filtered_chr_renamed.vcf $PATH_61 > concordancesolyntus_61.txt
#SnpSift concordance -v $PATH_SOLYNTUS/solyntus_filtered_chr_renamed.vcf $PATH_404 > concordancesolyntus_404.txt
'

# extract positions from isec vcfs to txt to load in R
bcftools query -f '%POS\t%CHROM\t%POS\n' isec_solyntus_to_61/solyntus.vcf > isec_solyntus_to_61/snps_solyntus.txt
bcftools query -f '%POS\t%CHROM\t%POS\n' isec_solyntus_to_61/61.vcf > isec_solyntus_to_61/snps_61.txt
bcftools query -f '%POS\t%CHROM\t%POS\n' isec_solyntus_to_61/both.vcf > isec_solyntus_to_61/snps_solyntus_61.txt

bcftools query -f '%POS\t%CHROM\t%POS\n' isec_404_to_61/404.vcf > isec_404_to_61/snps_404.txt
bcftools query -f '%POS\t%CHROM\t%POS\n' isec_404_to_61/61.vcf > isec_404_to_61/snps_61.txt
bcftools query -f '%POS\t%CHROM\t%POS\n' isec_404_to_61/both.vcf > isec_404_to_61/snps_404_61.txt

bcftools query -f '%POS\t%CHROM\t%POS\n' isec_solyntus_404/404.vcf > isec_solyntus_404/snps_404.txt
bcftools query -f '%POS\t%CHROM\t%POS\n' isec_solyntus_404/solyntus.vcf > isec_solyntus_404/snps_solyntus.txt
bcftools query -f '%POS\t%CHROM\t%POS\n' isec_solyntus_404/both.vcf > isec_solyntus_404/snps_solyntus_404.txt