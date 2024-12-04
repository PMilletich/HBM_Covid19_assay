library(tidyr)
library(ggplot2)
library(ggpubr)

OG_data = read.csv("PostPandemic_75.csv")

#Covert Wide to Long
OG_data_L <- gather(OG_data, Type, Measurement, FRNT.eGFP.50:MSD.IgA.Spike, factor_key=TRUE)


########################################################
#Spike Distribution 
########################################################
#Data Subset
OG_data_L_Spike = subset(OG_data_L, OG_data_L$Type %in% c("MSD.IgG.Spike","MSD.IgA.Spike"))
#Data Formatting 
OG_data_L_Spike$Ig = ifelse(grepl("IgA", OG_data_L_Spike$Type), "IgA", "IgG")
OG_data_L_Spike$Ig = ifelse(grepl("IgA", OG_data_L_Spike$Type), "IgA", "IgG")
OG_data_L_Spike$Group = gsub("MSD.IgG.Spike", "Spike IgG", OG_data_L_Spike$Type)
OG_data_L_Spike$Group = gsub("MSD.IgA.Spike", "Spike IgA", OG_data_L_Spike$Group)
OG_data_L_Spike$Spike = "Spike"

OG_data_L_Spike$Group  = factor(OG_data_L_Spike$Group, levels = c("Spike IgG", "Spike IgA"))
IG = ggplot(OG_data_L_Spike, aes(x = Group, y = log10(Measurement), fill = Ig, group = Group)) + 
  #Add Horizontal lines for cutoff from Figure 3
  geom_text(data = OG_data_L_Spike[OG_data_L_Spike$Spike == "Spike",], 
            aes(x = "Spike IgG", y = log10(0.03), label = "_ _ _ _ _ _ _ _ _ _"), 
            color = "coral") + 
  geom_text(data = OG_data_L_Spike[OG_data_L_Spike$Spike == "Spike",], 
            aes(x = "Spike IgA", y = log10(3.7),label = "_ _ _ _ _ _ _ _ _ _"), 
            color = "cornflowerblue") + 
  #Add points for each individual 
  geom_jitter(shape = 21, width = 0.3, size = 4) + 

  #Add Error bars and median line
  geom_point(stat = "summary", fun = "median", size = 10,
             shape = "-", color = "black", fill = "black") + 
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                linewidth = 1, 
                fun.args = list(mult = 1),
                position =  position_dodge(width = 0.9), color = "black") + 
  #Facet wrap to match Panel B
  facet_grid(~Spike, scale = "free_x", space = "free_x") + 
  #Aesthetics 
  ylab("Log10 BAU/mL") + 
  theme_bw() + 
  scale_fill_manual(breaks = c("IgG", "IgA"), 
                    values = c("darksalmon", "lightskyblue")) + 
  theme(axis.text.y = element_text(size = 11), axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        legend.position = "none",
        strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm"), 
                                    size = 14)); IG



########################################################
#Nucleocapsid cutoff
########################################################
OG_data_L_N = subset(OG_data_L, OG_data_L$Type %in% c("MSD.IgA.N", "MSD.IgG.N"))
OG_data_L_N$SelfReported= ifelse(OG_data_L_N$SR_Covid == 1, 
                                 "Infected", "Uninfected")
OG_data_L_N = subset(OG_data_L_N, is.na(OG_data_L_N$SelfReported) == F)
#Cutoffs from Figure 3
A_Cutoff = 12.7528; G_Cutoff = 0.0113

OG_data_L_N$Trend = ifelse(OG_data_L_N$Type == "MSD.IgA.N" & OG_data_L_N$Measurement <= A_Cutoff, "Negative: IgA", 
                           ifelse(OG_data_L_N$Type == "MSD.IgA.N" & OG_data_L_N$Measurement > A_Cutoff, "Positive: IgA", 
                                  ifelse(OG_data_L_N$Type == "MSD.IgG.N" & OG_data_L_N$Measurement <= G_Cutoff, "Negative: IgG", 
                                         ifelse(OG_data_L_N$Type == "MSD.IgG.N" & OG_data_L_N$Measurement > G_Cutoff, "Positive: IgG", NA ))))

#Summarize Data
# OG_data_L_N[OG_data_L_N$isotype == "N IgA",] %>%
#   group_by(SelfReported) %>%
#   summarise(
#     n = n(),
#     median = median(Measurement),
#     Q1 = quantile(Measurement, probs = c(0.25)), 
#     Q2 = quantile(Measurement, probs = c(0.75))
#   )

OG_data_L_N$isotype = gsub("MSD.IgG.N", "N IgG", OG_data_L_N$Type)
OG_data_L_N$isotype = gsub("MSD.IgG.A", "N IgA", OG_data_L_N$isotype)

OG_data_L_N$isotype = factor(OG_data_L_N$isotype, levels = c("N IgG", "N IgA"))
set.seed(1)
N_cutoff = ggplot(OG_data_L_N, aes(x = SelfReported,
                                   y = log10(Measurement), 
                                   fill = Trend, shape = Trend,
                                   group = SelfReported))  + 
  #Add cutoffs from Figure 3 
  geom_hline(data = subset(OG_data_L_N, OG_data_L_N$isotype == "N IgG")[1,], 
             aes(yintercept = log10(G_Cutoff)), linetype = "dashed", color = "red3") + 
  geom_hline(data = subset(OG_data_L_N, OG_data_L_N$isotype == "N IgA")[1,], 
             aes(yintercept = log10(A_Cutoff)), linetype = "dashed", color = "darkblue") + 
  #Add points for each indivdiual 
  geom_jitter(width = 0.25, size = 3) +
  #Add Error bars and median 
  geom_point(stat = "summary", fun = "median", size = 10, 
             shape = "-", color= "black", fill = "black", 
             position = position_dodge(0.8)) + 
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                linewidth = 1, color = "black",
                width = 0.4, position = position_dodge(0.8), 
                fun.args = list(mult = 1))+ 
  #Compare between infected and uninfected 
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))), 
                     label.x = 1.35, label.y = 2.4, 
                     size = 5) + 
  #Facet wrap to separate IgG and IgA 
  facet_wrap(~isotype) + 
  #Change shape and color 
  scale_shape_manual("Trend", breaks = c("Positive: IgG", "Positive: IgA", 
                                         "Negative: IgG", "Negative: IgA"), 
                     values = c(21, 21, 24, 24)) + 
  scale_fill_manual("Trend", breaks = c("Positive: IgG", "Positive: IgA", 
                                        "Negative: IgG", "Negative: IgA"), 
                    values = c("darksalmon", "lightskyblue", 
                               "indianred1", "royalblue1")) + 
  #Update Aesthetics 
  theme_bw() + 
  ylab("Log10 BAU/mL") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 15),
        legend.position = "bottom",
        legend.title = element_blank(), 
        legend.text = element_text(size = 10), 
        strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm"), 
                                    size = 14)); N_cutoff


########################################################
#Pseudovirus
########################################################

OG_data_L_Pseudo = subset(OG_data_L, OG_data_L$Type == "Pseudo.IC50")
OG_data_L_Pseudo$Type = "Pseudovirus"
Pseudo_plot = ggplot(OG_data_L_Pseudo, aes(x = Type, y = log10(Measurement))) + 
  #LLOQ 
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dotted") + 
  geom_violin( fill = "goldenrod", alpha = 0.1) + 
  geom_jitter(width = 0.3, size = 4, fill = "goldenrod", shape = 21)+ 
  #Error bars
  geom_point(stat = "summary", fun = "median", size = 10,
             shape = "-", color = "black", fill = "black") + 
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                linewidth = 1, width  = 0.5, 
                fun.args = list(mult = 1), color = "black") + 
  #Aesthetics 
  coord_cartesian(ylim = c(0, 3)) + 
  ylab(expression("Log10 IC"[50])) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15), axis.title.x = element_blank(), 
        legend.position = "none",
        strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm"), 
                                    size = 14)); Pseudo_plot


########################################################
#FRNT Boxplots 
########################################################
OG_data_L_FRNT = subset(OG_data_L, OG_data_L$Type %in% c("FRNT.eGFP.50","FRNT.50"))

OG_data_L_FRNT$Type = gsub("FRNT.eGFP.50", "FRNT-eGFP", OG_data_L_FRNT$Type )
OG_data_L_FRNT$Type = gsub("FRNT.50", "FRNT-AF647", OG_data_L_FRNT$Type )
OG_data_L_FRNT$Type = factor(OG_data_L_FRNT$Type , 
                             levels = c("FRNT-eGFP", "FRNT-AF647"))

OG_data_L_FRNT$Threshold = ifelse(OG_data_L_FRNT$Measurement == (min(OG_data_L_FRNT$Measurement, na.rm = T)), 
                                  "Threshold", "Over")

FRNT_plot = ggplot(OG_data_L_FRNT, aes(x = Type, y = log10(Measurement))) + 
  #LLOQ
  geom_hline(yintercept = log10(min(OG_data_L_FRNT$Measurement, na.rm = T)),
             linewidth = 0.5, linetype = "dotted") +
  #Points and violin
  geom_violin( fill = "mediumseagreen", alpha = 0.1) + 
  geom_jitter(width = 0.3, size = 4, fill = "mediumseagreen", shape = 21)+ 
  #Error bars 
  geom_point(stat = "summary", fun = "median", size = 10,
             shape = "-", color = "black", fill = "black") + 
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                linewidth = 1, width  = 0.5, 
                fun.args = list(mult = 1), color = "black") + 
  #Compare two distributions 
  stat_compare_means(method = "wilcox", 
                     aes(label = paste0("p = ", after_stat(p.format))), 
                     label.x = 1.35, size = 5) + 
  #Aesthetics 
  theme_bw() + 
  ylab(expression("Log10 FRNT"[50])) + 
  coord_cartesian(ylim = c(0, 3)) + 
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 15), axis.title.x = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm"), 
                                    size = 14)); FRNT_plot

########################################################
#FRNT Correlations 
########################################################
OG_data_L_FRNT_corr = data.frame(Y = OG_data_L_FRNT[OG_data_L_FRNT$Type == "FRNT-AF647","Measurement"],
                                 X = OG_data_L_FRNT[OG_data_L_FRNT$Type == "FRNT-eGFP","Measurement"])

FRNT_plot_corr = ggplot(OG_data_L_FRNT_corr, aes(x = log10(X), y = log10(Y))) +

  geom_point(size = 4, fill = "mediumseagreen", shape = 21)+ 

  #Add text from correlation:... stat_cor()
  geom_text(data = OG_data_L_FRNT_corr[1,],
            aes(x = 0.75, y = 2.5, label = "ρ=0.95    \np<0.0001"),
            size = 5) + 
  #Aesthetics 
  ylab(expression("Log10 FRNT"[50]*"-AF647")) +
  xlab(expression("Log10 FRNT"[50]*"-eGFP")) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.title = element_blank()); FRNT_plot_corr


########################################################
#Save Image 
########################################################
jpeg("./Figure4.jpeg", res = 400, height = 4000, width = 5000)
ggarrange(ggarrange(IG, N_cutoff, ncol = 2, widths = c(1,1.5), align = "h"),
          ggarrange(Pseudo_plot, FRNT_plot, FRNT_plot_corr, align = "h",
                    ncol = 3),
          nrow = 2)
dev.off()
