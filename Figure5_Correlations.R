library(ggplot2)
library(ggpubr)

setwd("C:/Users/PMilletich/OneDrive - University of Maryland School of Medicine/Desktop/Assorted/BM_Covid19/Final_Data.Scripts/")
OG_data = read.csv("PostPandemic_75.csv")

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


jpeg("IgA_Corr.jpeg", res = 400, height = 1500, width = 3000)
ggarrange(A_Pseudo, A_FRNT,GFP_FRNT, ncol = 3)
dev.off()
