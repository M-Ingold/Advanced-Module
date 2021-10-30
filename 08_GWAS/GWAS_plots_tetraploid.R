library(ldsep)
library(here)
library(qqplotr)
library(qqman)


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

TwoFrame <- cbind(dataTwo@map$Marker, dataTwo@map$Chrom, dataTwo@map$Position, dataTwo@scores$Trait)
colnames(TwoFrame) <- c("SNP", "CHR", "BP", "score1", "score2")
TwoFrame$CHR <- substr(TwoFrame$CHR, 4,5)
TwoFrame$SNP <- as.character(TwoFrame$SNP)
TwoFrame$CHR <- as.numeric(TwoFrame$CHR)
TwoFrame$score <- rep(NA, nrow(TwoFrame))
for (i in 1:nrow(TwoFrame)) {
  if (!is.na(TwoFrame$score1)[i] & is.na(TwoFrame$score2)[i]) {
    TwoFrame$score[i] = TwoFrame$score1[i]
  }
  else if (is.na(TwoFrame$score1)[i] & !is.na(TwoFrame$score2)[i]) {
    TwoFrame$score[i] = TwoFrame$score2[i]
  }
  else if (!is.na(TwoFrame$score1)[i] & !is.na(TwoFrame$score2)[i]) {
    if (TwoFrame$score1[i] >= TwoFrame$score2[i]) {
      TwoFrame$score[i] = TwoFrame$score1[i]
    }
    else {
      TwoFrame$score[i] = TwoFrame$score2[i]
    }
  }
  else {
    TwoFrame$score[i] = NA
  }
}

DipGenFrame <- cbind(dataDipGen@map$Marker, dataDipGen@map$Chrom, dataDipGen@map$Position, dataDipGen@scores$Trait)
colnames(DipGenFrame) <- c("SNP", "CHR", "BP", "score")
DipGenFrame$CHR <- substr(DipGenFrame$CHR, 4,5)
DipGenFrame$SNP <- as.character(DipGenFrame$SNP)
DipGenFrame$CHR <- as.numeric(DipGenFrame$CHR)

DipAddFrame <- cbind(dataDipAdd@map$Marker, dataDipAdd@map$Chrom, dataDipAdd@map$Position, dataDipAdd@scores$Trait)
colnames(DipAddFrame) <- c("SNP", "CHR", "BP", "score")
DipAddFrame$CHR <- substr(DipAddFrame$CHR, 4,5)
DipAddFrame$SNP <- as.character(DipAddFrame$SNP)
DipAddFrame$CHR <- as.numeric(DipAddFrame$CHR)

GenFrame <- cbind(dataGen@map$Marker, dataGen@map$Chrom, dataGen@map$Position, dataGen@scores$Trait)
colnames(GenFrame) <- c("SNP", "CHR", "BP", "score")
GenFrame$CHR <- substr(GenFrame$CHR, 4,5)
GenFrame$SNP <- as.character(GenFrame$SNP)
GenFrame$CHR <- as.numeric(GenFrame$CHR)

# SNPs to be marked
# addSNPs = c(5552536,5596081,5596242,6489191,7096264,7096341,7110360,55797680)
# oneSNPs = c(5552536,5596081,5596242,7096264,7096341,7110360,55797680)
# twoSNPs = c(55569227)
# dipgenSNPs = c(5552536,5596081,5596242,6489191,7096264,7096341,7110360)
# dipaddSNPs = c(5552536,5596081,5596242,6489191,7096264,7096341,7110360,54075466)
# genSNPs = c(5552536, 5596242)
# 
# addRanges <- list()
# oneRanges <- list()
# twoRanges <- list()
# dipgenRanges <- list()
# dipaddRanges <- list()
# genRanges <- list()
# 
# for(i in 1:length(addSNPs)) {
#   n = c(addSNPs[i]-670000, addSNPs[i]+670000)
#   addRanges[[i]] <- n
# }
# 
# AddFrameMarker <- data.frame(row.names = row.names(AddFrame))
# 
# for(i in 1:length(addRanges)) {
#   AddFrameMarker <- rbind(AddFrameMarker, subset(AddFrame, BP > addRanges[[i]][1] & BP < addRanges[[i]][2]))
# }
# 
# 
# for(i in 1:length(oneSNPs)) {
#   n = c(oneSNPs[i]-670000, oneSNPs[i]+670000)
#   oneRanges[[i]] <- n
# }
# 
# OneFrameMarker <- data.frame(row.names = row.names(OneFrame))
# 
# for(i in 1:length(oneRanges)) {
#   OneFrameMarker <- rbind(OneFrameMarker, subset(OneFrame, BP > oneRanges[[i]][1] & BP < oneRanges[[i]][2]))
# }
# 
# for(i in 1:length(twoSNPs)) {
#   n = c(twoSNPs[i]-670000, twoSNPs[i]+670000)
#   twoRanges[[i]] <- n
# }
# 
# TwoFrameMarker <- data.frame(row.names = row.names(TwoFrame))
# 
# for(i in 1:length(twoRanges)) {
#   TwoFrameMarker <- rbind(TwoFrameMarker, subset(TwoFrame, BP > twoRanges[[i]][1] & BP < twoRanges[[i]][2]))
# }
# 
# for(i in 1:length(dipgenSNPs)) {
#   n = c(dipgenSNPs[i]-670000, dipgenSNPs[i]+670000)
#   dipgenRanges[[i]] <- n
# }
# 
# DipgenFrameMarker <- data.frame(row.names = row.names(DipGenFrame))
# 
# for(i in 1:length(dipgenRanges)) {
#   DipgenFrameMarker <- rbind(DipgenFrameMarker, subset(DipGenFrame, BP > dipgenRanges[[i]][1] & BP < dipgenRanges[[i]][2]))
# }
# 
# for(i in 1:length(dipaddSNPs)) {
#   n = c(dipaddSNPs[i]-670000, dipaddSNPs[i]+670000)
#   dipaddRanges[[i]] <- n
# }
# 
# dipaddFrameMarker <- data.frame(row.names = row.names(DipAddFrame))
# 
# for(i in 1:length(dipaddRanges)) {
#   dipaddFrameMarker <- rbind(dipaddFrameMarker, subset(DipAddFrame, BP > dipaddRanges[[i]][1] & BP < dipaddRanges[[i]][2]))
# }
# 
# for(i in 1:length(genSNPs)) {
#   n = c(genSNPs[i]-670000, genSNPs[i]+670000)
#   genRanges[[i]] <- n
# }
# 
# GenFrameMarker <- data.frame(row.names = row.names(GenFrame))
# 
# for(i in 1:length(genRanges)) {
#   GenFrameMarker <- rbind(GenFrameMarker, subset(GenFrame, BP > genRanges[[i]][1] & BP < genRanges[[i]][2]))
# }



# find 5 highest scoring SNPs per model
highscores <- function(Frame){
  tempscores <- sort(Frame$score, decreasing = T)[1:10]
  return(subset(Frame, score >= min(tempscores)))
}

snparea <- function(topscorers, targetFrame){
  Ranges <- list()
  for(i in 1:nrow(topscorers)) {
    n = c(topscorers$BP[i]-670000, topscorers$BP[i]+670000)
    Ranges[[i]] <- n
  }
  
  FrameMarker <- data.frame(row.names = row.names(targetFrame))
  
  for(i in 1:nrow(topscorers)) {
    FrameMarker <- rbind(FrameMarker, subset(targetFrame, BP > Ranges[[i]][1] & BP < Ranges[[i]][2]))
  }
  return(FrameMarker)
}

addscoretetra <- highscores(AddFrame)
addarea <- snparea(addscoretetra, AddFrame)

onescoretetra <- highscores(OneFrame)
onearea <- snparea(onescoretetra, OneFrame)

twoscoretetra <- highscores(TwoFrame)
twoarea <- snparea(twoscoretetra, TwoFrame)

dipaddscoretetra <- highscores(DipAddFrame)
dipaddarea <- snparea(dipaddscoretetra, DipAddFrame)

dipgenscoretetra <- highscores(DipGenFrame)
dipgenarea <- snparea(dipgenscoretetra, DipGenFrame)

genscoretetra <- highscores(GenFrame)
genarea <- snparea(genscoretetra, GenFrame)









png(filename='manhattan_plot_tetraploid_short.png', res=130, width = 1000, height = 700)
# Defining layout for big plot
layout(matrix(c(1,2,3), byrow = F),
       heights =c(6,6,7))

# Creating the manhattan plots
par(mar = c(2, 4.2, 0.5, 2))
manAdd <- manhattan(AddFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 6), 
                    genomewideline = dataAdd@threshold, suggestiveline = min(addscoretetra$score), highlight = addarea$SNP)
title("Additive", line = -1.5)

manOne <- manhattan(OneFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 6), 
                    genomewideline = dataOne@threshold, suggestiveline = min(onescoretetra$score), highlight = onearea$SNP)
title("Simplex", line = -1.5)

# manTwo <- manhattan(TwoFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 6), 
#                     genomewideline = dataTwo@threshold, suggestiveline = min(twoscoretetra$score), highlight = twoarea$SNP)
# title("Duplex", line = -1.5)
# 
# manDipGen <- manhattan(DipGenFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 6), 
#                        genomewideline = dataDipGen@threshold, suggestiveline = min(dipgenscoretetra$score), highlight = dipgenarea$SNP)
# title("Diploid-general", line = -1.5)
# 
# manDipAdd <- manhattan(DipAddFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 6), 
#                        genomewideline = dataDipAdd@threshold, suggestiveline = min(dipaddscoretetra$score), highlight = dipaddarea$SNP)
# title("Diploid-additive", line = -1.5)

#par(mar = c(4, 4.2, 0.5, 2))
manGen <- manhattan(GenFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 6), 
                    genomewideline = dataGen@threshold, suggestiveline = min(genscoretetra$score), highlight = genarea$SNP)
title("General", line = -1.5)

dev.off()

# Calculating and plotting inflation rates for all models
pAdd <- 10^-na.omit(unlist(dataAdd@scores))
lambdaAdd <- median(qchisq(pAdd, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
pOne <- 10^-na.omit(unlist(dataOne@scores))
lambdaOne <- median(qchisq(pOne, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
pTwo <- 10^-na.omit(unlist(dataTwo@scores))
lambdaTwo <- median(qchisq(pTwo, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
pDipGen <- 10^-na.omit(unlist(dataDipGen@scores))
lambdaDipGen <- median(qchisq(pDipGen, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
pDipAdd <- 10^-na.omit(unlist(dataDipAdd@scores))
lambdaDipAdd <- median(qchisq(pDipAdd, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
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
  annotate("text", label = "lambda[GC]^A == '0.98'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
  geom_hline(aes(yintercept=4.08), color = "#1B9E77") + 
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
  annotate("text", label = "lambda[GC]^S == '0.92'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
  geom_hline(aes(yintercept=3.97), color = "#1B9E77") + 
  xlim(0, 5) + 
  ylim(0, 6.5)

qqTwo <- qq.plot(dataTwo, trait = "Trait") +
  ggtitle("Duplex") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = c(0.5,0.9),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text.x = element_blank()) +
  annotate("text", label = "lambda[GC]^D == '1.00'", parse = TRUE, x = 0.7, y = 6, size = 4.5) +
  geom_hline(aes(yintercept=3.99), color = "#1B9E77") +
  xlim(0, 5) +
  ylim(0, 6.5)

qqDipGen <- qq.plot(dataDipGen, trait = "Trait") + 
  ggtitle("Diploid-general") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "none",
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text.x = element_blank()) + 
  annotate("text", label = "lambda[GC]^DG == '0.93'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
  geom_hline(aes(yintercept=4.08), color = "#1B9E77") + 
  xlim(0, 5) + 
  ylim(0, 6.5)

qqDipAdd <- qq.plot(dataDipAdd, trait = "Trait") + 
  ggtitle("Diploid-additive") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "none",
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text.x = element_blank()) + 
  annotate("text", label = "lambda[GC]^DA == '1.08'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
  geom_hline(aes(yintercept=4.08), color = "#1B9E77") + 
  xlim(0, 5) + 
  ylim(0, 6.5)

qqGen <- qq.plot(dataGen, trait = "Trait") + 
  ggtitle("General") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "none",
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text.x = element_blank()) + 
  annotate("text", label = "lambda[GC]^G == '0.84'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
  geom_hline(aes(yintercept=4.08), color = "#1B9E77") + 
  xlim(0, 5) + 
  ylim(0, 6.5)

png(filename='qq_plots_tetraploid.png', res=130, width = 1000, height = 1300)
# Defining layout for big plot
#par(bg = 'white', mfrow = c(6,1))
#layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = T),
#       heights =c(6,6,6,6,6,7))
plot_grid(qqAdd, qqOne, qqTwo, qqDipGen, qqDipAdd, qqGen,
          nrow = 3, ncol = 2, labels = "auto", label_size = 20)
dev.off()

png(filename='qq_plots_tetraploid_short.png', res=130, width = 1000, height = 1000)

plot_grid(qqAdd, qqOne, qqGen,
          nrow = 2, ncol = 2, labels = "auto", label_size = 20)
dev.off()