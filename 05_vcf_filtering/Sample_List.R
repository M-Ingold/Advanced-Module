library(readxl)
genotypes <- read_excel("/media/rna/CYSTOPHORA/GWAS/data/GWAS_data/Additional_data/22-01-2021_Zuordnung_genotypes_SB.aktualisiert.xlsx", range = "E1:G190")
genotypesdf <- data.frame(genotypes)
relevantIDs <- genotypesdf[genotypesdf$Stärke...verfügbar == 'ja',]
t <- "Sample_"
newlist <- paste(t, relevantIDs$Identifier, sep = "")
write.table(newlist, file = "samples.txt", sep = "\n",
            row.names = FALSE, col.names = F, quote = F)