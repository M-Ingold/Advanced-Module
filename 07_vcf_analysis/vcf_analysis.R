library(vcfR)
library(ape)
library(here)

set_here(path = "/media/rna/CYSTOPHORA/GWAS")

# Analysis of 6.1
vcf_path <- file.path('data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS/61_filtered.vcf')
dna_path <- file.path('data/References/6.1_chr09.fa')

vcf_file_6.1 <- vcfR::read.vcfR(vcf_path)
dna_file_6.1 <- ape::read.dna(dna_path, format = "fasta")

chrom_6.1 <- create.chromR(name = 'Chr09_6.1', vcf = vcf_file_6.1, seq = dna_file_6.1)
chrom_6.1 <- proc.chromR(chrom_6.1, verbose=TRUE)

mqm_61 <- extract.info(chrom_6.1, "MQM")
head(mqm_61)

#plot(chrom_6.1)
png(filename = "Chr9_Coverage_61.png", 1000, 1200, res = 100)
chromoqc(chrom_6.1, dp.alpha = 90)
dev.off()

# Analysis of 4.0.4
vcf_path <- file.path('data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_4.4/404_filtered.vcf')
dna_path <- file.path('data/References/4.0.4_chr09.fasta')

vcf_file_4.0.4 <- vcfR::read.vcfR(vcf_path)
dna_file_4.0.4 <- ape::read.dna(dna_path, format = "fasta")

chrom_4.0.4 <- create.chromR(name = 'Chr09_4.0.4', vcf = vcf_file_4.0.4, seq = dna_file_4.0.4)
chrom_4.0.4 <- proc.chromR(chrom_4.0.4, verbose=TRUE)
#plot(chrom_4.0.4)

png(filename = "Chr9_Coverage_404.png", 1000, 1200, res = 100)
chromoqc(chrom_4.0.4, dp.alpha = 90)
dev.off()


# Analysis of Solyntus
vcf_path <- file.path('data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_solyntus/solyntus_filtered.vcf')
dna_path <- file.path('data/References/solyntus_chr09.fasta')

vcf_file_solyntus <- vcfR::read.vcfR(vcf_path)
dna_file_solyntus <- ape::read.dna(dna_path, format = "fasta")

chrom_solyntus <- create.chromR(name = 'Chr09_Solyntus', vcf = vcf_file_solyntus, seq = dna_file_solyntus)
chrom_solyntus <- proc.chromR(chrom_solyntus, verbose=TRUE)
#plot(chrom_solyntus)
png(filename = "Chr9_Coverage_Solyntus.png", 1000, 1200, res = 100)
par(mar = c(6, 10, 6, 6))
chromoqc(chrom_solyntus, dp.alpha = 90)
dev.off()

