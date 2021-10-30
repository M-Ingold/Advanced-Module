library(readxl)
library(reshape2)
library(multcompView)
library(tidyverse)
library(compare)

data <- read_excel("/media/rna/CYSTOPHORA/GWAS/data/GWAS_data/Additional_data/Screeningdaten_Set 1-6.xlsx", skip = 1)
sampleNames <- read_excel("/media/rna/CYSTOPHORA/GWAS/data/GWAS_data/Additional_data/22-01-2021_Zuordnung_genotypes_SB_aktualisiert_verbessert.xlsx")

# Fixing Nr column
data$Nr <- as.numeric(data$Nr.)
colnames(data)[1] <- "Nr"

# Only taking the data for which we have unique identifier of genotyped samples in the 'sampleNames'
data <- data[1:191,]
sampleNames <- sampleNames[1:189,c(1,5)]

#setdiff(data$`3`, sampleNames$VARIETY)

# Deleting rows from the data with unmeasured phenotype, for matching the sample identifier
data <- data[-c(139,156), ]

# Specifying different data frames to concatenate when desired
commonData <- data[1:12] # Data common to all measured traits
dataStarch <- data[34:41]
newColNamesStarch <- c('CR1', 'CR2', 'CR3', 'CR4', 'HR1', 'HR2', 'HR3', 'HR4')
colnames(dataStarch) <- newColNamesStarch

# Appending the common data by the identifier for the GWAS
sampleNames <- sampleNames[order(sampleNames$Nr), ]
commonData <- merge(sampleNames, commonData, by = "Nr")

# Separating the starch data into heat and control and calculate average of all observed data points for each category
dataStarchControl <- dataStarch %>%
  select(CR1, CR2, CR3, CR4)
dataStarchControl <- dataStarchControl %>% 
  mutate(CR = rowMeans(select(., starts_with("CR")), na.rm = TRUE))

dataStarchHeat <- dataStarch %>%
  select(HR1, HR2, HR3, HR4)
dataStarchHeat <- dataStarchHeat %>% 
  mutate(HR = rowMeans(select(., starts_with("HR")), na.rm = TRUE))

# Creating a new data frame containing only the necessary columns for the starch content
dataStarchNew <- cbind(commonData, dataStarchHeat, dataStarchControl)
# Creating percentage ratio of starch concentration for Heat/Control
dataStarchNew$Perc <- 100 * (dataStarchNew$HR / dataStarchNew$CR)

# Selecting only experiments 1 to 4 for GWAS
expStarchGwas <- dataStarchNew %>%
  filter(Exp. == "Exp_01" | Exp. == "Exp_02" | Exp. == "Exp_03" | Exp. == "Exp_04" | Identifier == "P2-A12-191-LI" | Identifier == "P2-B02-106-AMS")
# Removing rows with NAs in the 'Perc" column'
expStarchGwas <- subset(expStarchGwas, (!is.na(expStarchGwas$Perc)))

phenoGWAS <- expStarchGwas %>%
  select(Identifier, Perc)
colnames(phenoGWAS) <- c("Individual", "Trait")

# Adopting individual names to fit to the VCF files
phenoGWAS$Individual <- paste0("Sample_", phenoGWAS$Individual)

phenoGWAS$Trait <- round(phenoGWAS$Trait,2)

#phenoGWAS$Individual <- str_replace_all(phenoGWAS$Individual, "-", ".")

# Remove 045-AG, as it was wrongly marked as not phenotyped, so it was filtered out
phenoGWAS <- subset(phenoGWAS, Individual!="Sample_P1-E06-045-AG")

write.csv(phenoGWAS, "/media/rna/CYSTOPHORA/GWAS/data/GWAS_data/Additional_data/GWAS_input.csv", row.names=FALSE, quote=FALSE)