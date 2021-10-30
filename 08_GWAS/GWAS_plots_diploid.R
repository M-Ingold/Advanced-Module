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
# 
# # SNPs to be marked
# addSNPs = c(5552536,5596081,5596242,6489191,7096264,7096341,7110360,55797680)
# oneSNPs = c(5552536,5596081,5596242,7096264,7096341,7110360,55797680)
# genSNPs = c(5552536, 5596242)
# 
# addRanges <- list()
# oneRanges <- list()
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
# for(i in 1:length(addRanges)) {
#   OneFrameMarker <- rbind(OneFrameMarker, subset(OneFrame, BP > oneRanges[[i]][1] & BP < oneRanges[[i]][2]))
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

addscore <- highscores(AddFrame)
addarea <- snparea(addscore, AddFrame)

onescore <- highscores(OneFrame)
onearea <- snparea(onescore, OneFrame)

genscore <- highscores(GenFrame)
genarea <- snparea(genscore, GenFrame)

# png(filename='manhattan_plot_diploid_correct.png', res=130, width = 1000, height = 700)
# # Defining layout for big plot
# layout(matrix(c(1,2,3), byrow = F),
#        heights =c(6,6,7))
# 
# # Creating the manhattan plots
# par(mar = c(2, 4.2, 0.5, 2))
# manAdd <- manhattan(AddFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 6), 
#                     genomewideline = dataAdd@threshold, suggestiveline = min(addscore$score), highlight = addarea$SNP)
# title("Additive", line = -1.5)
# 
# manOne <- manhattan(OneFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 6), 
#                     genomewideline = dataOne@threshold, suggestiveline = min(onescore$score), highlight = onearea$SNP)
# title("Simplex", line = -1.5)
# 
# # manDipGen <- manhattan(DipGenFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 7), 
# #                        genomewideline = 4.5, suggestiveline = F)
# # title("Diploid-general", line = -1.5)
# # 
# # manDipAdd <- manhattan(DipAddFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 7), 
# #                        genomewideline = 4.5, suggestiveline = F)
# # title("Diploid-additive", line = -1.5)
# 
# #par(mar = c(4, 4.2, 0.5, 2))
# manGen <- manhattan(GenFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), ylim = c(0, 6), 
#                     genomewideline = dataGen@threshold, suggestiveline = min(genscore$score), highlight = genarea$SNP)
# title("General", line = -1.5)
# 
# dev.off()
# 
# # Calculating and plotting inflation rates for all models
# pAdd <- 10^-na.omit(unlist(dataAdd@scores))
# lambdaAdd <- median(qchisq(pAdd, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
# pOne <- 10^-na.omit(unlist(dataOne@scores))
# lambdaOne <- median(qchisq(pOne, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
# #pTwo <- 10^-na.omit(unlist(dataTwo@scores))
# #lambdaTwo <- median(qchisq(pTwo, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
# # pDipGen <- 10^-na.omit(unlist(dataDipGen@scores))
# # lambdaDipGen <- median(qchisq(pDipGen, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
# # pDipAdd <- 10^-na.omit(unlist(dataDipAdd@scores))
# # lambdaDipAdd <- median(qchisq(pDipAdd, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
# pGen <- 10^-na.omit(unlist(dataGen@scores))
# lambdaGen <- median(qchisq(pGen, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)

# ####################
# #                  #
# # qq-plot plotting #
# #                  #
# ####################
# 
# qqAdd <- qq.plot(dataAdd, trait = "Trait") + 
#   ggtitle("Additive") + 
#   theme(plot.title = element_text(hjust = 0.5, size = 15),
#         legend.position = "none",
#         strip.text.x = element_blank(),
#         axis.title.x = element_text(size = 10),
#         axis.title.y = element_text(size = 10)) + 
#   annotate("text", label = "lambda[GC]^A == lambdaAdd", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
#   geom_hline(aes(yintercept=4.5), color = "#1B9E77") + 
#   xlim(0, 5) + 
#   ylim(0, 6.5)
# 
# qqOne <- qq.plot(dataOne, trait = "Trait") + 
#   ggtitle("Simplex") + 
#   theme(plot.title = element_text(hjust = 0.5, size = 15),
#         legend.position = c(0.5,0.9),
#         legend.title = element_blank(),
#         legend.background = element_blank(),
#         legend.key = element_blank(),
#         strip.text.x = element_blank(),
#         axis.title.x = element_text(size = 10),
#         axis.title.y = element_text(size = 10)) + 
#   annotate("text", label = "lambda[GC]^S == lambdaOne", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
#   geom_hline(aes(yintercept=4.43), color = "#1B9E77") + 
#   xlim(0, 5) + 
#   ylim(0, 6.5)
# 
# # qqTwo <- qq.plot(dataTwo, trait = "Trait") +
# #   ggtitle("Duplex") +
# #   theme(plot.title = element_text(hjust = 0.5, size = 15),
# #         legend.position = c(0.5,0.9),
# #         legend.title = element_blank(),
# #         legend.background = element_blank(),
# #         legend.key = element_blank(),
# #         axis.title.x = element_text(size = 10),
# #         axis.title.y = element_text(size = 10),
# #         strip.text.x = element_blank()) +
# #   annotate("text", label = "lambda[GC]^D == '1.01'", parse = TRUE, x = 0.7, y = 6, size = 4.5) +
# #   geom_hline(aes(yintercept=4.79), color = "#1B9E77") +
# #   xlim(0, 5) +
# #   ylim(0, 6.5)
# 
# # qqDipGen <- qq.plot(dataDipGen, trait = "Trait") + 
# #   ggtitle("Diploid-general") + 
# #   theme(plot.title = element_text(hjust = 0.5, size = 15),
# #         legend.position = "none",
# #         axis.title.x = element_text(size = 10),
# #         axis.title.y = element_text(size = 10),
# #         strip.text.x = element_blank()) + 
# #   annotate("text", label = "lambda[GC]^DG == '0.98'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
# #   geom_hline(aes(yintercept=4.5), color = "#1B9E77") + 
# #   xlim(0, 5) + 
# #   ylim(0, 6.5)
# # 
# # qqDipAdd <- qq.plot(dataDipAdd, trait = "Trait") + 
# #   ggtitle("Diploid-additive") + 
# #   theme(plot.title = element_text(hjust = 0.5, size = 15),
# #         legend.position = "none",
# #         axis.title.x = element_text(size = 10),
# #         axis.title.y = element_text(size = 10),
# #         strip.text.x = element_blank()) + 
# #   annotate("text", label = "lambda[GC]^DA == '1.11'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
# #   geom_hline(aes(yintercept=4.5), color = "#1B9E77") + 
# #   xlim(0, 5) + 
# #   ylim(0, 6.5)
# 
# qqGen <- qq.plot(dataGen, trait = "Trait") + 
#   ggtitle("General") + 
#   theme(plot.title = element_text(hjust = 0.5, size = 15),
#         legend.position = "none",
#         axis.title.x = element_text(size = 10),
#         axis.title.y = element_text(size = 10),
#         strip.text.x = element_blank()) + 
#   annotate("text", label = "lambda[GC]^G == lambdaGen", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
#   geom_hline(aes(yintercept=4.5), color = "#1B9E77") + 
#   xlim(0, 5) + 
#   ylim(0, 6.5)
# 
# png(filename='qq_plots_diploid_correct.png', res=130, width = 1000, height = 1000)
# # Defining layout for big plot
# #par(bg = 'white', mfrow = c(6,1))
# #layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = T),
# #       heights =c(6,6,6,6,6,7))
# plot_grid(qqAdd, qqOne, qqGen,
#           nrow = 2, ncol = 2, labels = "auto", label_size = 20)
# dev.off()