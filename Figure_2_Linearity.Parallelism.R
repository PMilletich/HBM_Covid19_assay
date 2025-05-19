library(ggplot2); library(ggpubr); library(dplyr)

# ---- Linearity ---- 
Linearity_data = read.csv("Linearity.csv")

IgA = subset(Linearity_data, Linearity_data$Antibody == "IgA")
IgG = subset(Linearity_data, Linearity_data$Antibody == "IgG")




IgA_plot = ggplot(IgA, aes(x = log10(BAU.mL), y = log10(ECL), 
                           shape = Antigen, group = Antigen, 
                           fill = Antibody)) + 
  theme_bw() + theme(legend.position = "none") + 
  geom_line() + 
  geom_point(size = 3) + 
  stat_cor(aes(label = after_stat(rr.label)), r.digits = 3, color = "black",
           label.x  = -2.5, method= "pearson") + 
  scale_fill_manual(breaks = c("IgA", "IgG"), 
                    values = c("lightskyblue", "darksalmon")) + 
  scale_shape_manual(breaks = c("N", "RBD", "Spike"),
                     values = c(22, 24, 21)) + 
  xlab("Log10 BAU/mL") + 
  ylab("Log10 ECL");IgA_plot

IgG_plot = ggplot(IgG, aes(x = log10(BAU.mL), y = log10(ECL), 
                           shape = Antigen, group = Antigen, 
                           fill = Antibody)) + 
  theme_bw() + theme(legend.position = "none") + 
  geom_line() + 
  geom_point(size = 3) + 
  stat_cor(aes(label = after_stat(rr.label)), r.digits = 3,
           label.x = -3.5) + 
  scale_fill_manual(breaks = c("IgA", "IgG"), 
                    values = c("lightskyblue", "darksalmon")) + 
  scale_shape_manual(breaks = c("N", "RBD", "Spike"),
                     values = c(22, 24, 21)) + 
  xlab("Log10 BAU/mL") + 
  ylab("Log10 ECL")


Linearity = ggarrange(IgG_plot, IgA_plot, ncol = 2); Linearity


####################################### -
# ---- Parallelism ---- 
####################################### -

Ab = "IgA"



for (Ab in c("IgA", "IgG")) {
  if (Ab == "IgA") {
    current_data = read.csv("./Parallelism/Parallelism_IgA_2.csv")

  } else {
    current_data = read.csv("./Parallelism/Parallelism_IgG_2.csv")
  }
  
  
  #Show only data in detection range or the RS1 samples 
  current_data_1 = subset(current_data, current_data$Detection.Range == "In.Detection.Range" | 
                            current_data$Sample == "RS1")
  #Convert to log scale for dilution and BAU/mL
  current_data_1$log.BAU = log10(round(as.numeric(current_data_1$BAU.mL),5))
  current_data_1$log.Dilution = log10(round(as.numeric(current_data_1$Inverse.of.dilution),5))
  
  
  #Remove those with NAs without BAU/mL [unusedRS1s]
  current_data_1 = subset(current_data_1, is.na(current_data_1$log.BAU) == F)
  current_data_1$ID_Date = paste(current_data_1$Date, current_data_1$Sample)
  current_data_1$Ag = factor(current_data_1$Ag, levels = c("Spike", "RBD", "N"))
  
  #########################################
  current_antigen = "RBD"
  for (current_antigen in c("RBD", "Spike", "N")) {
    antigen_subset_1 = subset(current_data_1, current_data_1$Ag == current_antigen )
    table(antigen_subset_1$Date)
    #Remove most concentrated dilution from each sample which tend to be the most variable
     if (Ab == "IgA") {
       antigen_subset_1$log.Dilution = ifelse(antigen_subset_1$Sample == "RS1", antigen_subset_1$log.Dilution,
                                              ifelse(antigen_subset_1$log.Dilution < -3.7, NA, 
                                                     antigen_subset_1$log.Dilution))
       antigen_subset_1 = antigen_subset_1[is.na(antigen_subset_1$log.Dilution) == F,]
      for (current_sample in unique(antigen_subset_1$ID_Date)) {
        current_sample_df = subset(antigen_subset_1, antigen_subset_1$ID_Date == current_sample)
        
        if (current_sample_df[1,]$Sample != "RS1") {
          #Max
          current_sample_df_max = subset(current_sample_df, current_sample_df$log.Dilution ==
                                           max(current_sample_df$log.Dilution))
          antigen_subset_1 = antigen_subset_1[!(antigen_subset_1$ID_Date == current_sample &
                                                  antigen_subset_1$log.BAU == current_sample_df_max$log.BAU),]

        }
      }
     } else if (current_antigen == "RBD") {
       antigen_subset_1$log.Dilution = ifelse(antigen_subset_1$Sample == "RS1", antigen_subset_1$log.Dilution,
                                                ifelse(antigen_subset_1$log.Dilution > -2, NA, 
                                                       antigen_subset_1$log.Dilution))
    }

    Count_Samples = data.frame(table(antigen_subset_1$ID_Date))
    Count_Samples = subset(Count_Samples, Count_Samples$Freq >= 3)
    
    antigen_subset_1.1 = antigen_subset_1[antigen_subset_1$ID_Date  %in% Count_Samples$Var1,]
    antigen_subset_1.1 = antigen_subset_1.1[is.finite(antigen_subset_1.1$log.BAU),]
    
    
    antigen_subset_1.1 = subset(antigen_subset_1.1, is.na(antigen_subset_1.1$log.Dilution) == F)

    if (Ab == 'IgA') {
      antigen_subset_1.1 = subset(antigen_subset_1.1,  antigen_subset_1.1$Sample %in% 
                                  c("Pt.8.Day.38", "217C.post.dose.2",
                                    "UI-002-25MAY2022", "UI-003-25MAY2022", "RS1"  ))
      if (current_antigen == "N") {
        antigen_subset_1.1 = subset(antigen_subset_1.1, antigen_subset_1.1$Sample != "RSRB1507-IAD-COVID19")
        
      }
    } else {
      antigen_subset_1.1 = subset(antigen_subset_1.1, antigen_subset_1.1$ID_Date
                                  %in% c("6.7.2022 RS1",
                                         "6.7.2022 217C.post.dose.2",
                                         "6.7.2022 HP-0094782-211203",
                                         "6.7.2022 HP-0094782-220215",
                                         "6.7.2022 UI-001-25MAY2022" ))
    }
    
    antigen_subset_1.1 = subset(antigen_subset_1.1, is.infinite(antigen_subset_1.1$log.BAU) == F)
    antigen_subset_1.1 = subset(antigen_subset_1.1, is.infinite(antigen_subset_1.1$log.Dilution) == F)
    antigen_subset_ANCOVA_subset = anova(lm(log.BAU ~ log.Dilution * ID_Date, data = antigen_subset_1.1, ))
    antigen_subset_ANCOVA_subset = antigen_subset_ANCOVA_subset$`Pr(>F)`[3]
    antigen_subset_ANCOVA_subset = ifelse(antigen_subset_ANCOVA_subset > 0.001, 
                                          round(antigen_subset_ANCOVA_subset, 3), 
                                          formatC(antigen_subset_ANCOVA_subset, format = "e", digits = 2))
    print(antigen_subset_ANCOVA_subset)
    
    if (Ab == "IgA") {
      color_Ab = "lightskyblue"; ylim = 0.5; xlim = -3
    } else {
      color_Ab = "darksalmon"; ylim = 0; xlim = -4
    }
    
    antigen_subset_1.1$Ab_Shape = ifelse(antigen_subset_1.1$Sample == "RS1", 
                                      paste("RS1", Ab, sep = "\n"), 
                                      paste("HBM", Ab, sep = "\n"))
    
    Subset_antigen_subset_1 = ggplot(antigen_subset_1.1, aes(x = log.Dilution, y = log.BAU, 
                                                             group = ID_Date, shape = Ab_Shape)) + 
      geom_line( ) +
      geom_point(size = 3, fill = color_Ab) + theme_bw() +
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank()) + 
      ggtitle(paste(current_antigen, Ab)) + 
      scale_shape_manual(breaks = c("RS1\nIgG", "RS1\nIgA", "HBM\nIgG", "HBM\nIgA"), 
                         values = c( 24, 24, 21, 21)) + 
      annotate("text", label = paste("p=", antigen_subset_ANCOVA_subset, sep= ""),
               x = xlim, y = ylim, size = 4, colour = "black"); Subset_antigen_subset_1
    
    assign(paste(current_antigen, Ab, "plot", sep = "_"), Subset_antigen_subset_1)
    print(paste(current_antigen, Ab,  "plot", sep = "_"))
    
  }
}

jpeg("Parallelism_Figure2.jpeg", res = 600, height = 5000, width = 6000)
ggarrange(Linearity, ggarrange(Spike_IgG_plot, RBD_IgG_plot, N_IgG_plot, 
          Spike_IgA_plot, RBD_IgA_plot, N_IgA_plot,
          ncol = 3, nrow = 2, common.legend = T, legend = "none"), 
          nrow = 2, heights = c(1, 2))
dev.off()          
