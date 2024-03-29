#!/bin/bash

# This program calls variants in BAM files generated by bwa mem and sorted and indexed by samtools.

# Input and output for bamaddrg
PATH_IN_BAMS=../../data/GWAS_data/Final_BIG_GBS/alignment_data
PATH_BAM_LIST=bam_list_for_bamaddrg.txt

# Input and output for Freebayes
PATH_REF_FILE=../../data/References/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa
PATH_VCF_OUT=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS
mkdir -p $PATH_VCF_OUT
PATH_STAT_OUT=../../analysis/variant_calls_freebayes
mkdir -p $PATH_STAT_OUT

# Creation of input file for bamadrrg
# Creating a space separated list of input BAM files of the form "-b $BAM_FILE_PATH -s $SAMPLE_NAME

for folder in ${PATH_IN_BAMS}/Sample*; do
	sample=$(basename $folder)
	bamFile=$(realpath ${folder}/${sample}_sorted.bam)
	echo -n "-b ${bamFile} -s ${sample} " >> $PATH_BAM_LIST
done

# Calculating number of total input samples
TOTAL_SAMPLES=$(ls -l $PATH_IN_BAMS | grep -c ^d)

# Defining color variables for printing progress output!
LGREEN='\033[1;32m'
NC='\033[0m'

# Loop containing the execution of bamadrrg and freebayes #
# for chromosome 9                                        #

# Output name for the concatenated BAM file
BIG_BAM_FILE=$PATH_IN_BAMS/all_samples_09.bam
# Output path and name for the VCF file for each chromosome
VCF_OUT_FILE=$PATH_VCF_OUT/freebayes_${TOTAL_SAMPLES}_samples_09.vcf
# Specifying file name for compuation times for each chromosome
COMP=${PATH_STAT_OUT}/computation_times_09.txt

# Running bamaddrg for creating sample names for the 192 alignment files in PATH_BAM_LIST
# Test for existence of output file
if [ -f "$BIG_BAM_FILE" ]; then
	echo "$BIG_BAM_FILE already exists for Chr9, skipping bamaddrg."
else 
	BAMADDRG_INPUT=$(head -n 1 ${PATH_BAM_LIST}) # This is necessary for storing the file content in a variable, as expected from bamadrrg
	bdd_start=`date +%s`
	bamaddrg --clear --region 09 $BAMADDRG_INPUT > $BIG_BAM_FILE
	bdd_end=`date +%s`
	bdd_time=$(($bdd_end-bdd_start))
	echo "bamaddrg sample concatenation and naming time in s ($TOTAL_SAMPLES samples and whole 09):	$bdd_time" >> $COMP
	echo -e "${LGREEN}BAM concatenation of $TOTAL_SAMPLES samples and for chromosome 9 finished.${NC}"
fi

# Indexing the stacked BAM file for variant calling using Freebayes
idx_start=`date +%s`
samtools index $BIG_BAM_FILE
idx_end=`date +%s`
idx_time=$(($idx_end-idx_start))
echo "indexing time of the stacked BAM file in s ($TOTAL_SAMPLES samples and whole chromosome 09):	$idx_time" >> $COMP
echo -e "${LGREEN}Indexing of BAM file containing $TOTAL_SAMPLES samples and for chromosome 09 finished.${NC}"

# Running freebayes on chr09
fb_start=`date +%s`
freebayes \
	--fasta-reference $PATH_REF_FILE \
        --genotype-qualities \
    	--ploidy 2 \
    	--region chr09 \
    	--use-duplicate-reads \
    	--vcf $VCF_OUT_FILE \
    	${BIG_BAM_FILE}
fb_end=`date +%s`
fb_time=$(($fb_end-fb_start))
echo "freebayes variant calling time in s ($TOTAL_SAMPLES samples and 9):	$fb_time" >> $COMP
echo -e "${LGREEN}Variant call of $TOTAL_SAMPLES samples and for chromosome 9 finished.${NC}"


# Recording the total computation time
echo "total time in s:	$((bdd_time+idx_time+fb_time))" >> $COMP
echo -e "${LGREEN}Variant calling pipeline for $TOTAL_SAMPLES samples and for chromosome 9 completely done!${NC}"


