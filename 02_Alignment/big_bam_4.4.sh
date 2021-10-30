PATH_IN_BAMS=../../data/GWAS_data/Final_BIG_GBS/alignment_data_4.4
PATH_BAM_LIST=bam_list_for_bamaddrg_4.4.txt

# Input and output for Freebayes
PATH_REF_FILE=../../data/References/potato_dm_v404_all_pm_un.fasta
PATH_VCF_OUT=../../data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_4.4

# Creating a space separated list of input BAM files of the form "-b $BAM_FILE_PATH -s $SAMPLE_NAME
for folder in ${PATH_IN_BAMS}/Sample*; do
	sample=$(basename $folder)
	bamFile=$(realpath ${folder}/${sample}_sorted.bam)
	echo -n "-b ${bamFile} -s ${sample} " >> $PATH_BAM_LIST
done


TOTAL_SAMPLES=$(ls -l $PATH_IN_BAMS | grep -c ^d)

BIG_BAM_FILE=$PATH_IN_BAMS/all_samples.bam

BAMADDRG_INPUT=$(head -n 1 ${PATH_BAM_LIST}) # This is necessary for storing the file content in a variable, as expected from bamadrrg

bamaddrg --clear $BAMADDRG_INPUT > $BIG_BAM_FILE

samtools index $BIG_BAM_FILE
