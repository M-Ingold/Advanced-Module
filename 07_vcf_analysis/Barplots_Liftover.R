library(ggpubr)



df <- data.frame(Status=rep(c("Kept","Lost"), each=2),
                 Genome=rep(c("DM 4.0.4", "Solyntus"),2),
                 SNPs=c(5807,5274,30,727))

df <- df[order(df$Status, decreasing = F),]

png(filename = "SNPs_lifted_over.png", width=500, height = 500)
ggbarplot(df, "Genome", "SNPs",
          fill = "Status", color = "Status", palette = "grey",
          label = TRUE, lab.col = "white", lab.pos = "in", legend = "right") +
  ggtitle("SNPs lifted over") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
# ggplot(df[order(df$Status, decreasing = F),],
#        aes(Genome, SNPs, fill=factor(Status))) +
#   geom_bar(position="stack",stat = "identity")