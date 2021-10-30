library(vcfR)

vcf <- read.vcfR("./data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_4.4/404_filtered.vcf")
phenofile <- "./data/GWAS_data/Additional_data/GWAS_input.csv"
phenoData <- read_csv(phenofile)

dosage <- extract.gt(vcf, as.numeric = F)

dosage <- replace(dosage, dosage == "0/0", 0)
dosage <- replace(dosage, dosage == "0/1", 1)
dosage <- replace(dosage, dosage == "1/1", 2)



snpvisualisation <- function(snpbp){
snp <- paste("chr09_",snpbp, sep = "")
plot(jitter(as.numeric(dosage[snp, ]), factor = 0.2), 
     phenoData$Trait, xlab = "Dosage", ylab = "Phenotype Score", xlim=c(-0.1,2.1), 
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, main = snp)
abline(lm(phenoData$Trait ~ as.numeric(dosage[snp, ])))
}



# layout(matrix(c(1,2,3), byrow = T))

png(filename = "snps_diploid_additive.png", width = 1000, height = 2000)

par(mfrow=c(5,2))

for (i in 1:length(addscore$BP)){
  snpvisualisation(addscore$BP[i])
}

dev.off()

png(filename = "snps_diploid_simplex.png", width = 1000, height = 2000)

par(mfrow=c(5,2))

for (i in 1:length(onescore$BP)){
  snpvisualisation(onescore$BP[i])
}

dev.off()


#snpvisualisation(5596081)

#snpvisualisation(51373565)
png(filename = "5596081_dip.png", width = 1000, height = 700)

plot1 <- snpvisualisation(5596081)

dev.off()

png(filename = "6058273_dip.png", width = 1000, height = 700)

plot2 <- snpvisualisation(6058273)
dev.off()
