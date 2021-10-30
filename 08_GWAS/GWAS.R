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

genofile <- "/media/rna/CYSTOPHORA/GWAS/data/GWAS_data/Final_BIG_GBS/variant_calls_GWAS_using_4.4/404_filtered_GWASPoly.csv"
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
#                      K = dataAdd@K, max.iter = 1000); h2Add <- h2AddEst$h2; h2AddInt <- h2AddEst$conf.int1
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
dataAdd <- set.threshold(data = dataAdd, method = "Bonferroni", level = 0.05)
dataAdd <- set.threshold(data = dataAdd, method = "FDR", level = 0.05)
#dataAdd <- set.threshold(data = dataAdd, method = "permute", level = 0.05, n.permute = 100, n.core = 11)

dataOne <- set.threshold(data = dataOne, method = "M.eff", level = 0.05)
dataOne <- set.threshold(data = dataOne, method = "Bonferroni", level = 0.05)
dataOne <- set.threshold(data = dataOne, method = "FDR", level = 0.05)
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
dataGen <- set.threshold(data = dataGen, method = "Bonferroni", level = 0.05)
dataGen <- set.threshold(data = dataGen, method = "FDR", level = 0.05)
#dataGen <- set.threshold(data = dataGen, method = "permute", level = 0.05, n.permute = 100, n.core = 11)

# write results to the output files
# Additive model
outFileScoresAdd <- here('data', 'GWAS_results', 'scoresAdd.csv')
write.GWASpoly(data = dataAdd, trait = 'Trait', filename = outFileScoresAdd, what = 'scores', delim = ',')
outFileEffectsAdd <- here('data', 'GWAS_results', 'effectsAdd.csv')
write.GWASpoly(data = dataAdd, trait = 'Trait', filename = outFileEffectsAdd, what = 'effects', delim = ',')
# Simplex dominant model
outFileScoresOne <- here('data', 'GWAS_results', 'scoresOne.csv')
write.GWASpoly(data = dataOne, trait = 'Trait', filename = outFileScoresOne, what = 'scores', delim = ',')
outFileEffectsOne <- here('data', 'GWAS_results', 'effectsOne.csv')
write.GWASpoly(data = dataOne, trait = 'Trait', filename = outFileEffectsOne, what = 'effects', delim = ',')
# # Diplo-general model
# outFileScoresDipGen <- here('data', 'GWAS_results', 'scoresDipGen.csv')
# write.GWASpoly(data = dataDipGen, trait = 'Trait', filename = outFileScoresDipGen, what = 'scores', delim = ',')
# # Diplo-additive model
# outFileScoresDipAdd <- here('data', 'GWAS_results', 'scoresDipAdd.csv')
# write.GWASpoly(data = dataDipAdd, trait = 'Trait', filename = outFileScoresDipAdd, what = 'scores', delim = ',')
# outFileEffectsDipAdd <- here('data', 'GWAS_results', 'effectsDipAdd.csv')
# write.GWASpoly(data = dataDipAdd, trait = 'Trait', filename = outFileEffectsDipAdd, what = 'effects', delim = ',')
# General model
outFileScoresGen <- here('data', 'GWAS_results', 'scoresGen.csv')
write.GWASpoly(data = dataGen, trait = 'Trait', filename = outFileScoresGen, what = 'scores', delim = ',')
