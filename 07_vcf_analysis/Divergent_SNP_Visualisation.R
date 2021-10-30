library(vcfR)
library(ape)
library(chromoMap)
library(karyoploteR)
library(VariantAnnotation)
library(VennDiagram)
library(png)

vcf_path <- file.path('/media/rna/CYSTOPHORA/GWAS/scripts/vcf_analysis/isec_404_to_601/404_SNPs.vcf')
dna_path <- file.path('data/References/6.1_chr09.fa')

read.delim("/media/rna/CYSTOPHORA/GWAS/data/References/DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3", header=F, comment.char="#") -> gff
gff <- subset(gff, V1 == "chr09")


#vcf_file_6.1 <- vcfR::read.vcfR(vcf_path)
#dna_file_6.1 <- ape::read.dna(dna_path, format = "fasta")

#chrom_6.1 <- create.chromR(name = 'Chr09_6.1', vcf = vcf_file_6.1, seq = dna_file_6.1)
#chrom_6.1 <- proc.chromR(chrom_6.1, verbose=TRUE, win.size = 100000)

#chromoqc(chrom_6.1)

#generiert mit bcftools query
#chr_file = "/media/rna/CYSTOPHORA/GWAS/scripts/vcf_analysis/isec_404_to_601/chr.txt"
#anno_file = "/media/rna/CYSTOPHORA/GWAS/scripts/vcf_analysis/isec_404_to_601/404_pos.txt"

#chromoMap(chr_file,anno_file)
#          data_type = "numeric",
#          plots = "bar")

# VCF as GRanges
#vcf <- readVcf(vcf_path, )
#vcfGRanges <- vcfWhich(vcf)

pos_404 <- read.delim("/media/rna/CYSTOPHORA/GWAS/scripts/vcf_analysis/isec_404_to_61/snps_404.txt", header=F)
pos_404$V4 <- pos_404$V3 + 1
colnames(pos_404) <- c("n", "chr", "start", "end")
grange404 = makeGRangesFromDataFrame(pos_404)

pos_61 <- read.delim("/media/rna/CYSTOPHORA/GWAS/scripts/vcf_analysis/isec_404_to_61/snps_61.txt", header=FALSE)
pos_61$V4 <- pos_61$V3 + 1
colnames(pos_61) <- c("n", "chr", "start", "end")
grange61 = makeGRangesFromDataFrame(pos_61)

pos_both <- read.delim("/media/rna/CYSTOPHORA/GWAS/scripts/vcf_analysis/isec_404_to_61/snps_404_61.txt", header=FALSE)
pos_both$V4 <- pos_both$V3 + 1
colnames(pos_both) <- c("n", "chr", "start", "end")
grange_both = makeGRangesFromDataFrame(pos_both)

custom.genome <- toGRanges(data.frame(chr="chr09", start=1, end=67600300))

png(filename = "SNP_density_61_404.png", 1000, 800)

kp <- plotKaryotype(genome = custom.genome)
kpPlotDensity(kp, data=grange404, window.size = 100000, r0=0, r1=0.3, col="blue")
kpPlotDensity(kp, data=grange61, window.size = 100000, r0=0.33, r1=0.63, col="red")
kpPlotDensity(kp, data=grange_both, window.size = 100000, r0=0.66, r1=1, col="orchid")
kpAddLabels(kp, labels = " 4.0.4", r0=-0.1, r1=0.3, cex = 1.5)
kpAddLabels(kp, labels = "6.1", r0=0.23, r1=0.63, cex = 1.5)
kpAddLabels(kp, labels = "Both", r0=0.56, r1=1, cex = 1.5)
kpAddMainTitle(kp, "SNP Density of DM 4.0.4 vs. DM 6.1", cex = 2)

dev.off()

#solyntus 6.1 comparison
pos_61_sol <- read.delim("/media/rna/CYSTOPHORA/GWAS/scripts/vcf_analysis/isec_solyntus_to_61/snps_61.txt", header=FALSE)
pos_61_sol$V4 <- pos_61_sol$V3 + 1
colnames(pos_61_sol) <- c("n", "chr", "start", "end")
grange61_sol = makeGRangesFromDataFrame(pos_61_sol)

pos_sol <- read.delim("/media/rna/CYSTOPHORA/GWAS/scripts/vcf_analysis/isec_solyntus_to_61/snps_solyntus.txt", header=FALSE)
pos_sol$V4 <- pos_sol$V3 + 1
colnames(pos_sol) <- c("n", "chr", "start", "end")
grangesol = makeGRangesFromDataFrame(pos_sol)

pos_61_and_sol <- read.delim("/media/rna/CYSTOPHORA/GWAS/scripts/vcf_analysis/isec_solyntus_to_61/snps_solyntus_61.txt", header=FALSE)
pos_61_and_sol$V4 <- pos_61_and_sol$V3 + 1
colnames(pos_61_and_sol) <- c("n", "chr", "start", "end")
grange61_and_sol = makeGRangesFromDataFrame(pos_61_and_sol)

png(filename = "SNP_density_61_solyntus.png", 1000, 800)

kp <- plotKaryotype(genome = custom.genome)
kpPlotDensity(kp, data=grangesol, window.size = 100000, r0=0, r1=0.3, col="blue")
kpPlotDensity(kp, data=grange61_sol, window.size = 100000, r0=0.33, r1=0.63, col="red")
kpPlotDensity(kp, data=grange61_and_sol, window.size = 100000, r0=0.66, r1=1, col="orchid")
kpAddLabels(kp, labels = " Sol", r0=-0.1, r1=0.3, cex = 1.5)
kpAddLabels(kp, labels = "6.1", r0=0.23, r1=0.63, cex = 1.5)
kpAddLabels(kp, labels = "Both", r0=0.56, r1=1, cex = 1.5)
kpAddMainTitle(kp, "SNP Density of Solyntus vs. 6.1", cex = 2)

dev.off()
# Venn Diagram
#all_SNPs_61 <- read.delim("/media/rna/CYSTOPHORA/GWAS/data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS/61_filtered.txt")
#all_SNPs_404 <- read.delim("/media/rna/CYSTOPHORA/GWAS/data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_4.4/404_lifted_over_transanno.txt")
#all_SNPs_sol <- read.delim("/media/rna/CYSTOPHORA/GWAS/data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_solyntus/solyntus_lifted_over.txt")
#venn.diagram(x = list(all_SNPs_404$X10093298.1, all_SNPs_61$X119023.1, all_SNPs_sol$X30933981.1), filename = 'venn.png')


png(filename = "common_snps.png", 1000, 1000)
draw.triple.venn(area1 = 5766 , area2 = 5807, area3 = 5274, 
                 n12 = 5412, n13 = 4230, n23 = 4178, n123 = 4124, 
                 category = c('6.1', '4.0.4', 'Solyntus'),
                 main = "SNPs", main.pos = c(0.5, 1.05),
                 output = TRUE ,
                 imagetype="png" ,
                 #height = 480 , 
                 #width = 480 , 
                 #resolution = 300,
                 #compression = "lzw",
                 lwd = 3,
                 col=c("#440154ff", '#21908dff', '#fde725ff'),
                 fill = c("#440154ff",'#21908dff','#fde725ff'),
                 alpha = c(0.3,0.3,0.3),
                 cex = 2.5,
                 cat.cex = 3,
                 fontfamily = "sans",
                 cat.default.pos = "outer",
                 #cat.pos = c(-27, 27, 135),
                 #cat.dist = c(0.055, 0.055, 0.085),
                 #cat.fontfamily = "sans",
                 #cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
                 rotation = 1)
dev.off()

png(filename = "common_snps_ploidy.png", 1000, 1000)
draw.pairwise.venn(area1 = 5837, area2 = 2996, 
                   cross.area = 1863,
                   cex = 2.2,
                   cat.cex = 2.2,
                   category = c("diploid", "tetraploid"),
                   fill = c("#440154ff",'#21908dff'),
                   main = "SNPs", main.pos = c(0.5, 1.05),
                   alpha = c(0.3,0.3))
dev.off()