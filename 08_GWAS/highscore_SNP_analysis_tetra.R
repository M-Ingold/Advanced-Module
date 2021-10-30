library(vcfR)

phenofile <- "./data/GWAS_data/Additional_data/GWAS_input.csv"
phenoData <- read_csv(phenofile)


vcfTetra <- read.vcfR("./data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_tetraploid/tetraploid_filtered.vcf")
dosagetetra <- extract.gt(vcfTetra, as.numeric = F)

dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/0/0", 0)
dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/0/1", 1)
dosagetetra <- replace(dosagetetra, dosagetetra == "0/0/1/1", 2)
dosagetetra <- replace(dosagetetra, dosagetetra == "0/1/1/1", 3)
dosagetetra <- replace(dosagetetra, dosagetetra == "1/1/1/1", 4)


snpvisualisation <- function(snpbp) {
    snp <- paste("chr09_",snpbp, sep = "")
    plot(jitter(as.numeric(dosagetetra[snp, ]), factor = 0.2), phenoData$Trait, 
         xlab = "Dosage", ylab = "Phenotype Score", xlim=c(-0.1,4.1), cex.lab=1.5, 
         cex.axis=1.5, cex.main=1.5, cex.sub=1.5, main = snp)
    abline(lm(phenoData$Trait ~ as.numeric(dosagetetra[snp, ])))
}

# layout(matrix(c(1,2,3), byrow = T))

png(filename = "snps_tetraploid_additive.png", width = 1000, height = 2000)

par(mfrow=c(5,2))

for (i in 1:length(addscoretetra$BP)){
  snpvisualisation(addscoretetra$BP[i])
}

dev.off()


png(filename = "snps_tetraploid_general.png", width = 1000, height = 2000)

par(mfrow=c(5,2))

for (i in 1:length(genscoretetra$BP)){
  snpvisualisation(genscoretetra$BP[i])
}

dev.off()

png(filename = "snps_tetraploid_duplex.png", width = 1000, height = 2000)

par(mfrow=c(5,2))

for (i in 1:length(twoscoretetra$BP)){
  snpvisualisation(twoscoretetra$BP[i])
}

dev.off()


#snpvisualisation(5596081)

#snpvisualisation(51373565)

png(filename = "5596081_tetra.png", width = 1000, height = 700)
plot1 <- snpvisualisation(5596081)
dev.off()

png(filename = "6058273_tetra.png", width = 1000, height = 700)

plot2 <- snpvisualisation(6058273)
dev.off()
