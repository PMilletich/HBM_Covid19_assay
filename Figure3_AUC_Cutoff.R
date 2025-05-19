library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(cutpointr) #install.packages("cutpointr")
library(pROC)

#################
#Data upload
N_Data_OG = read.csv("N_Pre.Post.csv")
PostPandemic = read.csv("PostPandemic_236.csv")
S_data_OG = read.csv("S_Pre.Post.csv")


####################################################################################
#Figure 3, Generation of Cutoffs 
####################################################################################

###############################################################
#Step 1. AUC 
###############################################################
#Step 1.1, Nucleocapsid: Data Formatting 
Infection = subset(PostPandemic, PostPandemic$antigen == "Nucleocapsid")
Infection = subset(Infection, Infection$Cov == 1)
Infection = data.frame( "PATID" = Infection$guspec, 
                          "Time" = "Self-Reported", 
                          "N_IgG" = Infection[Infection$isotype == "IgG","bau_ml_conversion"], 
                          "N_IgA" = Infection[Infection$isotype == "IgA","bau_ml_conversion"])


Pre_Pandemic = N_Data_OG[N_Data_OG$Time == "Pre-Pandemic",]

AUC_N_Data = rbind(Pre_Pandemic, Infection)
AUC_N_Data$log_G = log10(AUC_N_Data$N_IgG)
AUC_N_Data$log_A = log10(AUC_N_Data$N_IgA)

AUC_N_Data$Vacc = ifelse(AUC_N_Data$Time == "Pre-Pandemic", 0, 1)

# Print summary data Example
# AUC_N_Data %>%
#   group_by(Time) %>%
#   summarise(
#     n = n(),
#     mean = mean(N_IgA), sd = sd(N_IgA),
#     min = min(N_IgA),median = median(N_IgA), max = max(N_IgA)   )


#####################
#Step 1.2, Nucleocapsid: ROC 
#IgG
N_IgG_Cutoff = summary(cutpointr(AUC_N_Data, N_IgG, Vacc))
N_IgG_Cutoff = as.numeric(N_IgG_Cutoff$cutpointr[[1]][2])

#IgA
N_IgA_Cutoff = summary(cutpointr(AUC_N_Data, N_IgA, Vacc))
N_IgA_Cutoff = as.numeric(N_IgA_Cutoff$cutpointr[[1]][2])

ROC_N_G = pROC::roc(AUC_N_Data$Vacc, AUC_N_Data$log_G, direction = "<")
ROC_N_A = pROC::roc(AUC_N_Data$Vacc, AUC_N_Data$log_A, direction = "<")


plot.roc(AUC_N_Data$Vacc, AUC_N_Data$N_IgG, print.thres = (N_IgG_Cutoff), print.auc = T)
#         xlim = c(0.5, 0.25))



#####################
#Step 1.3; Spike Data Formatting  
S_data_OG$log_G = log10(S_data_OG$S_IgG)
S_data_OG$log_A = log10(S_data_OG$S_IgA)

#####################
#Step 1.4; Spike ROC 
#IgA
ROC_S_G = pROC::roc(S_data_OG$Vaccinated, S_data_OG$log_G, direction = "<")
ROC_S_A = pROC::roc(S_data_OG$Vaccinated, S_data_OG$log_A, direction = "<")

#like.coordinates <- coords(ROC_S_G, c(-Inf, sort(S_IgG_Cutoff), Inf), input="threshold", ret=c("specificity", "sensitivity"))

plot.roc(S_data_OG$Vaccinated, S_data_OG$log_G, print.thres = 0.03)


###############################################################
#Step 2. Distribution of N
###############################################################

#####################
#Step 2.1: Data Formatting 
N_Data_Pre = subset(N_Data_OG, N_Data_OG$Time == "Pre-Pandemic")
N_Data_Pre$Time = "Controls"
head(N_Data_Pre)
#Wide to long
N_Data_Pre = gather(N_Data_Pre, "Type", "Measure", N_IgG:N_IgA, factor_key=TRUE)

#Post Pandemic formatting 
N_Data_Post_OG = subset(PostPandemic, PostPandemic$antigen == "Nucleocapsid")
N_Data_Post_OG$Time = ifelse(N_Data_Post_OG$Cov == 1, "Vaccinated\nwith Infection", 
                             "Vaccinated\nwithout Infection")
N_Data_Post_OG$Type = paste("N_", N_Data_Post_OG$isotype, sep = "")
N_Data_Post_OG = N_Data_Post_OG[, c("guspec", "Time", "Type", "bau_ml_conversion")]
colnames(N_Data_Post_OG) = c("PATID", "Time", "Type", "Measure")

N_Data_Combined = rbind(N_Data_Pre, N_Data_Post_OG)
table(N_Data_Combined$Type)
table(N_Data_Combined$Time)

#Using cutoffs established above 
N_Data_Combined$Trend = ifelse(N_Data_Combined$Type == "N_IgG" & N_Data_Combined$Measure > 0.011, 
                               "Positive: IgG", 
                               ifelse(N_Data_Combined$Type == "N_IgG" & N_Data_Combined$Measure <= 0.011, 
                                      "Negative: IgG", 
                                      ifelse(N_Data_Combined$Type == "N_IgA" & N_Data_Combined$Measure > 12.75, 
                                             "Positive: IgA", "Negative: IgA")))
table(N_Data_Combined$Trend, N_Data_Combined$Time)
#Clean up variables 
N_Data_Combined$Type = gsub("N_IgG", "Nucleocapsid IgG", N_Data_Combined$Type)
N_Data_Combined$Type = gsub("N_IgA", "Nucleocapsid IgA", N_Data_Combined$Type)

#Re-level factors for graphing 
N_Data_Combined$Type = factor(N_Data_Combined$Type, 
                              levels = c("Nucleocapsid IgG", "Nucleocapsid IgA"))
N_Data_Combined$Time_1 = factor(N_Data_Combined$Time, levels = 
                              c("Controls",
                                "Vaccinated\nwith Infection",
                                "Vaccinated\nwithout Infection"))
#Create Comparison list for stat_compare_means
compare_list = list(c("Controls", "Vaccinated\nwith Infection"), 
                    c("Vaccinated\nwithout Infection","Vaccinated\nwith Infection"),
                    c("Controls","Vaccinated\nwithout Infection"))


N = ggplot(N_Data_Combined, aes(x = Time_1, y = log10(Measure), 
                            shape = Trend, fill = Trend, group = Time_1))  + 
  #Create horizontal lines for cutoffs 
  geom_hline(data = N_Data_Combined[N_Data_Combined$Type == "Nucleocapsid IgA",],
             aes(yintercept = log10(0.053)),
             color = "darkblue", linetype = "dotted")+
  geom_hline(data = N_Data_Combined[N_Data_Combined$Type == "Nucleocapsid IgG",],
             aes(yintercept = log10(0.00043)),
             color = "red4", linetype = "dotted")+
  geom_hline(data = subset(N_Data_Combined, N_Data_Combined$Type == "Nucleocapsid IgA")[1,], 
             aes(yintercept = log10(12.75)), linetype = "dashed", color = "darkblue") + 
  geom_hline(data = subset(N_Data_Combined, N_Data_Combined$Type == "Nucleocapsid IgG")[1,], 
             aes(yintercept = log10(0.011)), linetype = "dashed", color = "red3") +  
  #Scatter plot of individual values 
  geom_point(size= 3, position = position_jitterdodge(0.65))+ 
  #Manually add Error bars 
  geom_point(stat = "summary", fun = "median", size = 8, 
             shape = "-", color= "black", fill = "black", 
             position = position_dodge(0.8)) + 
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                linewidth = 1, color = "black",
                width = 0.4, position = position_dodge(0.8), 
                fun.args = list(mult = 1)) + 
  #Set Y axis 
  ylab("Log 10 BAU/mL") + 
  #Compare Groups
  stat_compare_means(comparisons = compare_list, method = "wilcox", 
                     label.y = c(2.1, 2.8, 3.5),
                     size = 4, tip.length = 0, 
                     symnum.args=list(
                       cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                       symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "ns"))) + 
  #Facet by N IgA or IgG
  facet_wrap(~Type, scale = "free_y") + 
  #Standardize y axis
  coord_cartesian(ylim = c(-4, 4.5)) + 
  #Manually set point shape and color 
  scale_shape_manual("Trend", breaks = c("Positive: IgG", "Positive: IgA", 
                                         "Negative: IgG", "Negative: IgA"), 
                     values = c(21, 21, 24, 24)) + 
  scale_fill_manual("Trend", breaks = c("Positive: IgG", "Positive: IgA", 
                                        "Negative: IgG", "Negative: IgA"), 
                    values = c("darksalmon", "lightskyblue", 
                               "indianred1", "royalblue1")) + 
  #Set theme for aesthetics 
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 11),
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm"), 
                                    size = 12)); N

###############################################################
#Step 3. Distribution of S
###############################################################

S_data = gather(S_data_OG, "Type", "Measure", S_IgG:S_IgA, factor_key=TRUE)
S_data$Trend = ifelse(S_data$Type == "S_IgA" & S_data$Measure > 3.7, "Positive: IgA", 
                      ifelse(S_data$Type == "S_IgA" & S_data$Measure <= 3.7, "Negative: IgA", 
                             ifelse(S_data$Type == "S_IgG" & S_data$Measure > 0.03, "Positive: IgG",
                                    "Negative: IgG")))

S_data$Time = ifelse(S_data$Vaccinated == "no", "Controls", "Vaccinated±Infection")
S_data$Time_Ig = paste(S_data$Time, S_data$Ig)
table(S_data$Trend)


S_data$Time = factor(S_data$Time, levels = c("Controls", "Vaccinated±Infection"))
S_data$Time_Ig = factor(S_data$Time_Ig, levels = c("Pre-Pandemic IgG", "Pre-Pandemic IgA", 
                                                   "Post-Pandemic IgG", "Post-Pandemic IgA"))
S_data$Type = gsub("S_IgG", "Spike IgG", S_data$Type)
S_data$Type = gsub("S_IgA", "Spike IgA", S_data$Type)

S_data$Type = factor(S_data$Type , levels = c("Spike IgG", "Spike IgA"))


S = ggplot(S_data, aes(x = Time, y = log10(Measure), 
                       group = Time_Ig, fill = Trend, shape = Trend)) + 
  #Horizontal lines for LLOQ and LOD 
  geom_hline(data = S_data[S_data$Type == "Spike IgA",],
             aes(yintercept = log10(0.032)),
             color = "darkblue", linetype = "dotted")+
  geom_hline(data = S_data[S_data$Type == "Spike IgG",],
             aes(yintercept = log10(0.0018)),
             color = "red4", linetype = "dotted")+
  geom_hline(data = S_data[S_data$Type == "Spike IgA",],
             aes(yintercept = log10(3.7)), 
             color = "darkblue", linetype = "dashed")+
  geom_hline(data = S_data[S_data$Type == "Spike IgG",],
             aes(yintercept = log10(0.03)), 
             color = "red4", linetype = "dashed")+ 
  #Scatterplot of individual's values 
  geom_point(size= 3, position = position_jitterdodge(0.65))+ 
  #Wrap by Antibody Type
  facet_wrap(~Type, scale = "free_y") + 
  #Manual Y axis
  coord_cartesian(ylim = c(-4, 3.75)) + 
  #Manually add error bars 
  geom_point(stat = "summary", fun = "median", size = 8, 
             shape = "-", color= "black", fill = "black", 
             position = position_dodge(0.8))  +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                linewidth = 1, color = "black",
                width = 0.4, position = position_dodge(0.8), 
                fun.args = list(mult = 1)) + 
  #Add pvalues for statistical summary 
  geom_text(data = subset(S_data, S_data$Type == "Spike IgA")[1,], 
            aes(x = 1.5, y = 3.5) , size  = 4,
            label = "p<0.001")+ 
  geom_text(data = subset(S_data, S_data$Type == "Spike IgG")[1,], 
            aes(x = 1.5, y = 3.5) , size = 4,
            label = "p<0.001") +
  #Add Aesthetics and themes 
  ylab("Log 10 BAU/mL") + 
  theme_bw()+ 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 11),
        legend.position = "bottom",
        legend.title = element_blank(), 
        legend.text = element_text(size = 13), 
        strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm"), 
                                    size = 12)) + 
  scale_shape_manual("Trend", breaks = c("Positive: IgG", "Positive: IgA", 
                                         "Negative: IgG", "Negative: IgA"), 
                     values = c(21, 21, 24, 24)) + 
  scale_fill_manual("Trend", breaks = c("Positive: IgG", "Positive: IgA", 
                                        "Negative: IgG", "Negative: IgA"), 
                    values = c("darksalmon", "lightskyblue", 
                               "indianred1", "royalblue1")) ; S

###############################################################
#Step 4. Save Figures 
###############################################################

jpeg("Figure_3_SN_Cutoff.jpeg", res = 500, height = 3500, width = 4000)
ggarrange(S, N, common.legend = T, 
          align = "v",
          nrow = 2, legend = "bottom")
dev.off()

jpeg("AUC_Spike_Nucleo.jpeg", res = 400, height = 3000, width = 3000 )
par(mfrow = c(2, 2))
plot(ROC_S_G, print.auc = T); plot(ROC_S_A, print.auc = T);plot(ROC_N_G, print.auc = T); plot(ROC_N_A, print.auc = T) 
dev.off()

