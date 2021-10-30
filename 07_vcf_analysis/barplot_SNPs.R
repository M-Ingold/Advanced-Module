library(ggplot2)
library(gridExtra)
library(ggpubr)

df <- data.frame(genome=c('DM 4.0.4', 'DM 6.1', 'Solyntus'),
                 unfiltered_SNPs=c(179444,198709,213405),
                 filtered_SNPs=c(5837,5766,6506))
# df$genome <- factor(df$genome, levels = c('DM 4.0.4', 'DM 6.1', 'Solyntus'))
# 
# bp1 <- ggbarplot(df, x='genome', y='unfiltered_SNPs',
#                  xlab = 'Referene Genome',
#                  ylab = 'Unfiltered Variants',
#                  #color = "black",
#                  palette = "jco",            
#                  sort.val = "asc",
#                  ggtheme = theme_gray())
#                  
# 
# bp2 <- ggbarplot(df, x='genome', y='filtered_SNPs',
#                 xlab = 'Referene Genome',
#                 ylab = 'Filtered SNPs',
#                 #color = "black",
#                 palette = "jco",            
#                 sort.val = "asc",
#                 ggtheme = theme_gray(),
#                 #order = c('DM 4.0.4', 'DM 6.1', 'Solyntus')
#                 sort.by.groups = F)
# 
# bp1 = bp1 + font("xlab", size = 15) +
#   font("ylab", size = 15) +
#   font("xy.text", size = 13)
# 
# bp2 = bp2 + font("xlab", size = 15) +
#   font("ylab", size = 15) +
#   font("xy.text", size = 13)

bp1 <- ggplot(df, aes(x=genome, y=unfiltered_SNPs)) + 
  geom_bar(stat = "identity", width=0.5) +
  labs(x="Reference Genome", y = "Unfiltered Variants") +
  font("ylab", size = 15) +
  font("xlab", size = 15) +
  font("xy.text", size = 13)

bp2 <- ggplot(df, aes(x=genome, y=filtered_SNPs), xlab = 'Referene Genome', ylab = 'Filtered SNPs') + 
  geom_bar(stat = "identity", width=0.5) +
  labs(x="Reference Genome", y = "Filtered SNPs") +
  font("ylab", size = 15) +
  font("xlab", size = 15) +
  font("xy.text", size = 13)


png(filename = 'barplot_SNP_count.png', width = 1000, height = 500)
ggarrange(bp1, bp2)
dev.off()