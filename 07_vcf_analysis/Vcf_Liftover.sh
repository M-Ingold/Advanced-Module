#!/bin/bash

# Lift over SNPs from 4.0.4 and Solyntus to 6.1
# Should allow accurate concordance analysis and may show differences in alignment regions
# For generation of chain files see https://github.com/wurmlab/flo

# Generate dict file needed for Liftover
#java -jar /home/rna/genetools/picard.jar CreateSequenceDictionary -R ../../data/References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fasta

: '
java -jar /home/rna/genetools/picard.jar LiftoverVcf \
     I=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_4.4/404_filtered.vcf \
     O=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_4.4/404_lifted_over.vcf \
     CHAIN=flo_404_61_chains/run/combined.chn.sorted \
     REJECT=rejected_variants.vcf \
     R=../../data/References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fasta \
     LIFTOVER_MIN_MATCH=0.1 \
     RECOVER_SWAPPED_REF_ALT=true
'

: '
java -jar /home/rna/genetools/picard.jar LiftoverVcf \
     I=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_solyntus/solyntus_filtered.vcf \
     O=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_solyntus/solyntus_lifted_over.vcf \
     CHAIN=flo_solyntus_61_chains/run/combined.chn.sorted \
     REJECT=rejected_variants_solyntus.vcf \
     R=../../data/References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fasta \
     LIFTOVER_MIN_MATCH=0.1 \
     RECOVER_SWAPPED_REF_ALT=true \
     VERBOSITY=DEBUG
# Not working well, trying transanno instead!
'

: '
# Chain File generation using transanno

#minimap2 -cx asm5 --cs ../../data/References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fasta ../../data/References/SolyntusV1.1.fasta > Sol_to_61.paf
#/home/rna/genetools/transanno-v0.2.4-x86_64-unknown-linux-musl/transanno minimap2chain Sol_to_61.paf --output Sol_to_61.chain

minimap2 -cx asm5 --cs ../../data/References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fasta ../../data/References/potato_dm_v404_all_pm_un.fasta > 404_to_61.paf
/home/rna/genetools/transanno-v0.2.4-x86_64-unknown-linux-musl/transanno minimap2chain 404_to_61.paf --output 404_to_61.chain




# Liftover using transanno
/home/rna/genetools/transanno-v0.2.4-x86_64-unknown-linux-musl/transanno liftvcf -m --chain Sol_to_61.chain \
-o ../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_solyntus/solyntus_lifted_over_transanno.vcf.gz \
--query ../../data/References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fasta \
--reference ../../data/References/SolyntusV1.1.fasta \
--vcf ../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_solyntus/solyntus_filtered.vcf.gz \
--fail FAILED_liftover.vcf.gz
'

#Chain 6.1 to Solyntus
minimap2 -cx asm5 --cs ../../data/References/SolyntusV1.1.fasta ../../data/References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fasta > 61_to_Sol.paf
/home/rna/genetools/transanno-v0.2.4-x86_64-unknown-linux-musl/transanno minimap2chain 61_to_Sol.paf --output 61_to_Sol.chain

:'
# LiftoverVcf using transanno chain file
java -jar /home/rna/genetools/picard.jar LiftoverVcf \
     I=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_4.4/404_filtered.vcf \
     O=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_4.4/404_lifted_over.vcf \
     CHAIN=404_to_61.chain \
     REJECT=rejected_variants_404.vcf \
     R=../../data/References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fasta \
     LIFTOVER_MIN_MATCH=0.1 \
     RECOVER_SWAPPED_REF_ALT=true
#5807 / 5837 (99.4860%)
'
:'
java -jar /home/rna/genetools/picard.jar LiftoverVcf \
     I=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_solyntus/solyntus_filtered.vcf \
     O=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_solyntus/solyntus_lifted_over.vcf \
     CHAIN=Sol_to_61.chain \
     REJECT=rejected_variants_solyntus.vcf \
     R=../../data/References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fasta \
     LIFTOVER_MIN_MATCH=0.1 \
     RECOVER_SWAPPED_REF_ALT=true 
#5274 / 6001 (87.8854%)
'

java -jar /home/rna/genetools/picard.jar LiftoverVcf \
     I=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS/61_filtered.vcf \
     O=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS/61_lifted_over.vcf \
     CHAIN=61_to_Sol.chain \
     REJECT=rejected_variants_61.vcf \
     R=../../data/References/SolyntusV1.1.fasta \
     LIFTOVER_MIN_MATCH=0.1 \
     RECOVER_SWAPPED_REF_ALT=true 
#5207 / 5766 (90.3052%)