library(ggplot2)
library(ggpubr)
library(bubbleHeatmap)

OG_data = read.csv("PostPandemic_75.csv")
Corr_data = read.csv("Covid19_BM_Heatmap.csv")

########################
#Heatmap
########################
Corr = subset(Corr_data, Corr_data$Group == "Correlation")
rownames(Corr) = Corr[,1]
Corr$X = NULL; Corr$Group = NULL 
Corr[] <- lapply(Corr, function(x) { as.numeric(as.character(x))})
Corr_M = as.matrix(Corr)

P = subset(Corr_data, Corr_data$Group == "P")
P = P[-1,]
rownames(P) = P[,1]
P$X = NULL; P$Group = NULL 
P[] <- lapply(P, function(x) { -log10(as.numeric(as.character(x)))})
P_M = as.matrix(P)

tree <-  bubbleHeatmap(Corr_M, P_M, 
                       legendTitles = c("-log10(P)", "ρ"))

################################################
#Correlation Scatterplots 
################################################
A_Pseudo = ggplot(OG_data, aes(y = log10(Pseudo.IC50), x = log10(MSD.IgA.Spike)))+ 
  geom_point(shape = 21, size = 3, fill = "goldenrod") + 
  theme_bw() + 
  coord_cartesian(ylim = c(0,2.75), xlim = c(0.5,3.7)) + 
  ylab(expression("Pseudovirus: Log10 IC"[50]))+ 
  xlab(expression("Spike IgA: Log10 BAU/ML"))+
  geom_text(data = OG_data[1,], aes(y = log10(220), x = log10(18), 
             label = c("ρ=0.740\np<0.001"))); A_Pseudo 

A_FRNT = ggplot(OG_data, aes(y = log10(FRNT.50), x = log10(MSD.IgA.Spike)))+ 
  geom_point(shape = 21, size = 3, fill = "mediumseagreen") + 
  theme_bw() + 
  coord_cartesian(ylim = c(0,2.75), xlim = c(0.5,3.7)) + 
  ylab(expression("FRNT-AF647: Log10 FRNT"[50]))+ 
  xlab(expression("Spike IgA: Log10 BAU/ML"))+
  geom_text(data = OG_data[1,], aes(y = log10(220), x = log10(18), 
             label = c("ρ=0.840\np<0.001")))


GFP_FRNT = ggplot(OG_data, aes(y = log10(FRNT.eGFP.50), x = log10(MSD.IgA.Spike)))+ 
  geom_point(shape = 21, size = 3, fill = "darkseagreen4") + 
  theme_bw() + 
  coord_cartesian(ylim = c(0,2.75), xlim = c(0.5,3.7)) + 
  ylab(expression("FRNT-eGFP: Log10 FRNT"[50]))+ 
  xlab(expression("Spike IgA: Log10 BAU/ML"))+
  geom_text(data = OG_data[1,], aes(y = log10(220), x = log10(18), 
             label = c("ρ=0.804\np<0.001"))); GFP_FRNT

###################################################
#Save Images 
###################################################
jpeg("Heatmap.jpeg", res = 400, height = 1500, width = 2000)
grid.newpage()
grid.draw(tree)
dev.off()

jpeg("IgA_Corr.jpeg", res = 400, height = 1500, width = 3000)
ggarrange(A_Pseudo, A_FRNT,GFP_FRNT, ncol = 3)
dev.off()
