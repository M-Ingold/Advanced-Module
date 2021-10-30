library(GWASpoly)
library(gemma2)
library(tidyverse)
library(heritability)
library(cowplot)
#library(ldsep)
library(here)

#########################################
#                                       #
# Reading the geno- and phenotype files #
#                                       #
#########################################

# Files generated using MultiGWAS

genofile <- "/media/rna/CYSTOPHORA/GWAS/data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS/61_filtered_GWASPoly.csv"
phenofile <- "/media/rna/CYSTOPHORA/GWAS/data/GWAS_data/Additional_data/GWAS_input.csv"
phenoData <- read_csv(phenofile)
genoData <- read.csv(genofile)

# setdiff(colnames(genoData), phenoData$Individual)
# setdiff(phenoData$Individual,colnames(genoData))

data <- read.GWASpoly(ploidy=2,pheno.file=phenofile,geno.file=genofile,
                      format="ACGT",n.traits=1,delim=",")

###############################################
#                                             #
# Creating kinship matrix and performing GWAS #
#                                             #
###############################################

# Set the kinship matrix for polygenic effects
data <- set.K(data) #LOCO=F weil nur ein Chromosom?

# Further incorporation of population structure
params <- set.params(MAF = 0.01)

#Performing the GWAS multiple times
dataAdd <- GWASpoly(data = data, models = c("additive"),
                    traits = "Trait", params = params)
dataOne <- GWASpoly(data = data, models = c("1-dom"),
                    traits = "Trait", params = params)
#dataTwo <- GWASpoly(data = data, models = c("2-dom"),
#                    traits = "Trait", params = params)
# dataDipGen <- GWASpoly(data = data, models = c("diplo-general"),
#                        traits = "Trait", params = params)
# dataDipAdd <- GWASpoly(data = data, models = c("diplo-additive"),
#                        traits = "Trait", params = params)
dataGen <- GWASpoly(data = data, models = c("general"),
                    traits = "Trait", params = params)

####################################################################################
#                                                                                  #
# Heritability estimation, lambda estimation and GWAS sig. threshold determination #
#                                                                                  #
####################################################################################

# Estimating narrow-sense heritability for all models
# Additive model
#h2AddEst <- marker_h2_means(data.vector = dataAdd@pheno$Trait, geno.vector = dataAdd@pheno$Individual,
#                            K = dataAdd@K, max.iter = 1000); h2Add <- h2AddEst$h2; h2AddInt <- h2AddEst$conf.int1
# # Simplex dominant
# h2OneEst <- marker_h2(data.vector = dataOne@pheno$Trait, geno.vector = dataOne@pheno$Individual,
#                       K = dataOne@K, max.iter = 1000); h2One <- h2OneEst$h2; h2OneInt <- h2OneEst$conf.int1
# # Duplex dominant
# h2TwoEst <- marker_h2(data.vector = dataTwo@pheno$Trait, geno.vector = dataTwo@pheno$Individual,
#                       K = dataTwo@K, max.iter = 1000); h2Two <- h2TwoEst$h2; h2TwoInt <- h2TwoEst$conf.int1
# # Diploid General
# h2DipGenEst <- marker_h2(data.vector = dataDipGen@pheno$Trait, geno.vector = dataDipGen@pheno$Individual,
#                          K = dataDipGen@K, max.iter = 1000); h2DipGen <- h2DipGenEst$h2; h2DipGenInt <- h2DipGenEst$conf.int1
# # Diploid Additive
# h2DipAddEst <- marker_h2(data.vector = dataDipAdd@pheno$Trait, geno.vector = dataDipAdd@pheno$Individual,
#                          K = dataDipAdd@K, max.iter = 1000); h2DipAdd <- h2DipAddEst$h2; h2DipAddInt <- h2DipAddEst$conf.int1
# # General
# h2GenEst <- marker_h2(data.vector = dataGen@pheno$Trait, geno.vector = dataGen@pheno$Individual,
#                       K = dataGen@K, max.iter = 1000); h2Gen <- h2GenEst$h2; h2GenInt <- h2GenEst$conf.int1

# Setting significance thresholds
dataAdd <- set.threshold(data = dataAdd, method = "M.eff", level = 0.05)
#dataAdd <- set.threshold(data = dataAdd, method = "Bonferroni", level = 0.05)
#dataAdd <- set.threshold(data = dataAdd, method = "FDR", level = 0.05)
#dataAdd <- set.threshold(data = dataAdd, method = "permute", level = 0.05, n.permute = 100, n.core = 11)

dataOne <- set.threshold(data = dataOne, method = "M.eff", level = 0.05)
#dataOne <- set.threshold(data = dataOne, method = "Bonferroni", level = 0.05)
#dataOne <- set.threshold(data = dataOne, method = "FDR", level = 0.05)
#dataOne <- set.threshold(data = dataOne, method = "permute", level = 0.05, n.permute = 100, n.core = 11)

# dataTwo <- set.threshold(data = dataTwo, method = "M.eff", level = 0.05)
# dataTwo <- set.threshold(data = dataTwo, method = "Bonferroni", level = 0.05)
# dataTwo <- set.threshold(data = dataTwo, method = "FDR", level = 0.05)
# dataTwo <- set.threshold(data = dataTwo, method = "permute", level = 0.05, n.permute = 100)

# dataDipGen <- set.threshold(data = dataDipGen, method = "M.eff", level = 0.05)
# dataDipGen <- set.threshold(data = dataDipGen, method = "Bonferroni", level = 0.05)
# dataDipGen <- set.threshold(data = dataDipGen, method = "FDR", level = 0.05)
# #dataDipGen <- set.threshold(data = dataDipGen, method = "permute", level = 0.05, n.permute = 100, n.core = 11)
# 
# dataDipAdd <- set.threshold(data = dataDipAdd, method = "M.eff", level = 0.05)
# dataDipAdd <- set.threshold(data = dataDipAdd, method = "Bonferroni", level = 0.05)
# dataDipAdd <- set.threshold(data = dataDipAdd, method = "FDR", level = 0.05)
# #dataDipAdd <- set.threshold(data = dataDipAdd, method = "permute", level = 0.05, n.permute = 100, n.core = 11)

dataGen <- set.threshold(data = dataGen, method = "M.eff", level = 0.05)
#dataGen <- set.threshold(data = dataGen, method = "Bonferroni", level = 0.05)
#dataGen <- set.threshold(data = dataGen, method = "FDR", level = 0.05)
#dataGen <- set.threshold(data = dataGen, method = "permute", level = 0.05, n.permute = 100, n.core = 11)

# write results to the output files
# Additive model
# outFileScoresAdd <- here('data', 'GWAS_results', 'scoresAdd.csv')
# write.GWASpoly(data = dataAdd, trait = 'Trait', filename = outFileScoresAdd, what = 'scores', delim = ',')
# outFileEffectsAdd <- here('data', 'GWAS_results', 'effectsAdd.csv')
# write.GWASpoly(data = dataAdd, trait = 'Trait', filename = outFileEffectsAdd, what = 'effects', delim = ',')
# # Simplex dominant model
# outFileScoresOne <- here('data', 'GWAS_results', 'scoresOne.csv')
# write.GWASpoly(data = dataOne, trait = 'Trait', filename = outFileScoresOne, what = 'scores', delim = ',')
# outFileEffectsOne <- here('data', 'GWAS_results', 'effectsOne.csv')
# write.GWASpoly(data = dataOne, trait = 'Trait', filename = outFileEffectsOne, what = 'effects', delim = ',')
# # # Diplo-general model
# # outFileScoresDipGen <- here('data', 'GWAS_results', 'scoresDipGen.csv')
# # write.GWASpoly(data = dataDipGen, trait = 'Trait', filename = outFileScoresDipGen, what = 'scores', delim = ',')
# # # Diplo-additive model
# # outFileScoresDipAdd <- here('data', 'GWAS_results', 'scoresDipAdd.csv')
# # write.GWASpoly(data = dataDipAdd, trait = 'Trait', filename = outFileScoresDipAdd, what = 'scores', delim = ',')
# # outFileEffectsDipAdd <- here('data', 'GWAS_results', 'effectsDipAdd.csv')
# # write.GWASpoly(data = dataDipAdd, trait = 'Trait', filename = outFileEffectsDipAdd, what = 'effects', delim = ',')
# # General model
# outFileScoresGen <- here('data', 'GWAS_results', 'scoresGen.csv')
# write.GWASpoly(data = dataGen, trait = 'Trait', filename = outFileScoresGen, what = 'scores', delim = ',')

library(ldsep)
library(here)
library(qqplotr)
library(qqman)

# AddSNPs <- read.csv(here('data', 'GWAS_results', 'scoresAdd.csv'))
# OneSNPs <- as.character(unlist(read.csv(here('data', 'GWAS_results', 'scoresOne.csv'), header = F), use.names = F))
# DipGenSNPs <- as.character(unlist(read.csv(here('data', 'GWAS_results', 'scoresDipGen.csv'), header = F), use.names = F))
# DipAddSNPs <- as.character(unlist(read.csv(here('data', 'GWAS_results', 'scoresDipAdd.csv'), header = F), use.names = F))
# GenSNPs <- as.character(unlist(read.csv(here('data', 'GWAS_results', 'scoresGen.csv'), header = F), use.names = F))

#######################################
#                                     #
# Creating manhattan plots with qqman #
#                                     #
#######################################

# Create dataframes from the output of GWASpoly for each gene action model
AddFrame <- cbind(dataAdd@map$Marker, dataAdd@map$Chrom, dataAdd@map$Position, dataAdd@scores$Trait)
colnames(AddFrame) <- c("SNP", "CHR", "BP", "score")
AddFrame$CHR <- substr(AddFrame$CHR, 4,5)
AddFrame$SNP <- as.character(AddFrame$SNP)
AddFrame$CHR <- as.numeric(AddFrame$CHR)

OneFrame <- cbind(dataOne@map$Marker, dataOne@map$Chrom, dataOne@map$Position, dataOne@scores$Trait)
colnames(OneFrame) <- c("SNP", "CHR", "BP", "score1", "score2")
OneFrame$CHR <- substr(OneFrame$CHR, 4,5)
OneFrame$SNP <- as.character(OneFrame$SNP)
OneFrame$CHR <- as.numeric(OneFrame$CHR)
OneFrame$score <- rep(NA, nrow(OneFrame))
for (i in 1:nrow(OneFrame)) {
  if (!is.na(OneFrame$score1)[i] & is.na(OneFrame$score2)[i]) {
    OneFrame$score[i] = OneFrame$score1[i]
  }
  else if (is.na(OneFrame$score1)[i] & !is.na(OneFrame$score2)[i]) {
    OneFrame$score[i] = OneFrame$score2[i]
  }
  else if (!is.na(OneFrame$score1)[i] & !is.na(OneFrame$score2)[i]) {
    if (OneFrame$score1[i] >= OneFrame$score2[i]) {
      OneFrame$score[i] = OneFrame$score1[i]
    }
    else {
      OneFrame$score[i] = OneFrame$score2[i]
    }
  }
  else {
    OneFrame$score[i] = NA
  }
}
# 
# DipGenFrame <- cbind(dataDipGen@map$Marker, dataDipGen@map$Chrom, dataDipGen@map$Position, dataDipGen@scores$Trait)
# colnames(DipGenFrame) <- c("SNP", "CHR", "BP", "score")
# DipGenFrame$CHR <- substr(DipGenFrame$CHR, 4,5)
# DipGenFrame$SNP <- as.character(DipGenFrame$SNP)
# DipGenFrame$CHR <- as.numeric(DipGenFrame$CHR)
# 
# DipAddFrame <- cbind(dataDipAdd@map$Marker, dataDipAdd@map$Chrom, dataDipAdd@map$Position, dataDipAdd@scores$Trait)
# colnames(DipAddFrame) <- c("SNP", "CHR", "BP", "score")
# DipAddFrame$CHR <- substr(DipAddFrame$CHR, 4,5)
# DipAddFrame$SNP <- as.character(DipAddFrame$SNP)
# DipAddFrame$CHR <- as.numeric(DipAddFrame$CHR)

GenFrame <- cbind(dataGen@map$Marker, dataGen@map$Chrom, dataGen@map$Position, dataGen@scores$Trait)
colnames(GenFrame) <- c("SNP", "CHR", "BP", "score")
GenFrame$CHR <- substr(GenFrame$CHR, 4,5)
GenFrame$SNP <- as.character(GenFrame$SNP)
GenFrame$CHR <- as.numeric(GenFrame$CHR)


png(filename='manhattan_plot_61.png', res=130, width = 1000, height = 700)
# Defining layout for big plot
layout(matrix(c(1,2,3), byrow = F),
       heights =c(6,6,7))

# Creating the manhattan plots
par(mar = c(2, 4.2, 0.5, 2))
manAdd <- manhattan(AddFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 7), 
                    genomewideline = 4.5, suggestiveline = F)
title("Additive", line = -1.5)

manOne <- manhattan(OneFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 7), 
                    genomewideline = 4.43, suggestiveline = F)
title("Simplex", line = -1.5)

# manDipGen <- manhattan(DipGenFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 7), 
#                        genomewideline = 4.5, suggestiveline = F)
# title("Diploid-general", line = -1.5)
# 
# manDipAdd <- manhattan(DipAddFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 7), 
#                        genomewideline = 4.5, suggestiveline = F)
# title("Diploid-additive", line = -1.5)

#par(mar = c(4, 4.2, 0.5, 2))
manGen <- manhattan(GenFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 7), 
                    genomewideline = 4.5, suggestiveline = F)
title("General", line = -1.5)

dev.off()

# Calculating and plotting inflation rates for all models
pAdd <- 10^-na.omit(unlist(dataAdd@scores))
lambdaAdd <- median(qchisq(pAdd, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
pOne <- 10^-na.omit(unlist(dataOne@scores))
lambdaOne <- median(qchisq(pOne, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
#pTwo <- 10^-na.omit(unlist(dataTwo@scores))
#lambdaTwo <- median(qchisq(pTwo, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
# pDipGen <- 10^-na.omit(unlist(dataDipGen@scores))
# lambdaDipGen <- median(qchisq(pDipGen, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
# pDipAdd <- 10^-na.omit(unlist(dataDipAdd@scores))
# lambdaDipAdd <- median(qchisq(pDipAdd, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
pGen <- 10^-na.omit(unlist(dataGen@scores))
lambdaGen <- median(qchisq(pGen, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)

####################
#                  #
# qq-plot plotting #
#                  #
####################

qqAdd <- qq.plot(dataAdd, trait = "Trait") + 
  ggtitle("Additive") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "none",
        strip.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  annotate("text", label = "lambda[GC]^A == '1.11'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
  geom_hline(aes(yintercept=4.5), color = "#1B9E77") + 
  xlim(0, 5) + 
  ylim(0, 6.5)

qqOne <- qq.plot(dataOne, trait = "Trait") + 
  ggtitle("Simplex") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = c(0.5,0.9),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  annotate("text", label = "lambda[GC]^S == '0.99'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
  geom_hline(aes(yintercept=4.43), color = "#1B9E77") + 
  xlim(0, 5) + 
  ylim(0, 6.5)

# qqTwo <- qq.plot(dataTwo, trait = "Trait") +
#   ggtitle("Duplex") +
#   theme(plot.title = element_text(hjust = 0.5, size = 15),
#         legend.position = c(0.5,0.9),
#         legend.title = element_blank(),
#         legend.background = element_blank(),
#         legend.key = element_blank(),
#         axis.title.x = element_text(size = 10),
#         axis.title.y = element_text(size = 10),
#         strip.text.x = element_blank()) +
#   annotate("text", label = "lambda[GC]^D == '1.01'", parse = TRUE, x = 0.7, y = 6, size = 4.5) +
#   geom_hline(aes(yintercept=4.79), color = "#1B9E77") +
#   xlim(0, 5) +
#   ylim(0, 6.5)

# qqDipGen <- qq.plot(dataDipGen, trait = "Trait") + 
#   ggtitle("Diploid-general") + 
#   theme(plot.title = element_text(hjust = 0.5, size = 15),
#         legend.position = "none",
#         axis.title.x = element_text(size = 10),
#         axis.title.y = element_text(size = 10),
#         strip.text.x = element_blank()) + 
#   annotate("text", label = "lambda[GC]^DG == '0.98'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
#   geom_hline(aes(yintercept=4.5), color = "#1B9E77") + 
#   xlim(0, 5) + 
#   ylim(0, 6.5)
# 
# qqDipAdd <- qq.plot(dataDipAdd, trait = "Trait") + 
#   ggtitle("Diploid-additive") + 
#   theme(plot.title = element_text(hjust = 0.5, size = 15),
#         legend.position = "none",
#         axis.title.x = element_text(size = 10),
#         axis.title.y = element_text(size = 10),
#         strip.text.x = element_blank()) + 
#   annotate("text", label = "lambda[GC]^DA == '1.11'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
#   geom_hline(aes(yintercept=4.5), color = "#1B9E77") + 
#   xlim(0, 5) + 
#   ylim(0, 6.5)

qqGen <- qq.plot(dataGen, trait = "Trait") + 
  ggtitle("General") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "none",
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text.x = element_blank()) + 
  annotate("text", label = "lambda[GC]^G == '0.98'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
  geom_hline(aes(yintercept=4.5), color = "#1B9E77") + 
  xlim(0, 5) + 
  ylim(0, 6.5)

png(filename='qq_plots_diploid.png', res=130, width = 1000, height = 1000)
# Defining layout for big plot
#par(bg = 'white', mfrow = c(6,1))
#layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = T),
#       heights =c(6,6,6,6,6,7))
plot_grid(qqAdd, qqOne, qqGen,
          nrow = 2, ncol = 2, labels = "auto", label_size = 20)
dev.off()