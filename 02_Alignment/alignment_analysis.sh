#!/bin/bash

# interactive mode needed!
# analyse aligned files, input parent folder of folders with bam files

#mkdir ../../analysis

BAMFOLDER61=/media/rna/CYSTOPHORA/GWAS/data/GWAS_data/Final_BIG_GBS/alignment_data
BAMFOLDER404=/media/rna/CYSTOPHORA/GWAS/data/GWAS_data/Final_BIG_GBS/alignment_data_4.4
BAMFOLDERSOL=/media/rna/CYSTOPHORA/GWAS/data/GWAS_data/Final_BIG_GBS/alignment_data_solyntus

# For allSamples File
# by default the genome is split into 400 windows in which coverage etc is averaged

#qualimap bamqc -bam $BAMFOLDER61/all_samples_09.bam -outdir /media/rna/CYSTOPHORA/GWAS/analysis/$(basename $BAMFOLDER61) \
#--java-mem-size=50G -gff /media/rna/CYSTOPHORA/GWAS/data/References/DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3 \
#--output-genome-coverage /media/rna/CYSTOPHORA/GWAS/analysis/$(basename $BAMFOLDER61)/coverage --paint-chromosome-limits

qualimap bamqc -bam $BAMFOLDER404/all_samples.bam -outdir /media/rna/CYSTOPHORA/GWAS/analysis/$(basename $BAMFOLDER404) \
--java-mem-size=50G -gff /media/rna/CYSTOPHORA/GWAS/data/References/4.0.4_annotation.gff \
--output-genome-coverage /media/rna/CYSTOPHORA/GWAS/analysis/$(basename $BAMFOLDER404)/coverage --paint-chromosome-limits

#qualimap bamqc -bam $BAMFOLDERSOL/all_samples_09.bam -outdir /media/rna/CYSTOPHORA/GWAS/analysis/$(basename $BAMFOLDERSOL) \
#--java-mem-size=50G -gff /media/rna/CYSTOPHORA/GWAS/data/References/GeMoMaAnnotationsDM.gff \
#--output-genome-coverage /media/rna/CYSTOPHORA/GWAS/analysis/$(basename $BAMFOLDERSOL)/coverage --paint-chromosome-limits


# For each file
# better way by using multi-bamqc, config file with paths to files: http://qualimap.conesalab.org/doc_html/command_line.html#cmdline-multibamqc