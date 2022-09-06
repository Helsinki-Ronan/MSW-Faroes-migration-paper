



                    ############################################################
                    #                                                          #
                    #       Genetic stock identification reveals greater       #
                    #                use of an oceanic feeding                 #
                    #             ground around the Faroe Islands              #
                    #         by multi-sea winter Atlantic salmon----          #
                    #                                                          #
                    ############################################################



setwd("C:\\Users\\osuronan\\Dropbox\\Research\\Faroese salmon\\Data\\Data for analysis")

library(readxl)
library(ggplot2)
library(readr)
library(ggridges)
library(dplyr)
library(gridExtra)
library(MCMCglmm)                    
library(scales)
library(ggforce)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)
library(ggsn)
library(stringr)
library(ggOceanMaps)
library(ggOceanMapsData)
library(gtsummary)
                    


                    
                    ############################################################
                    #                                                          #
                    #          1.0 Load data and check structure----           #
                    #                                                          #
                    ############################################################

# 1.1 Fishery data----
fishery_data<- read_csv("faroese_fisheries_data_for_analysis.csv")
summary(fishery_data)

# 1.2 Recode variables----
fishery_data$season<- as.factor(as.character(fishery_data$season))

fishery_data$netting_record<- as.integer(fishery_data$netting_record)

fishery_data$sea_age<- as.factor(as.character(fishery_data$sea_age))

fishery_data$num_prey<- as.integer(fishery_data$num_prey)
fishery_data$num_fish<- as.integer(fishery_data$num_fish)
fishery_data$num_crus<- as.integer(fishery_data$num_crus)

fishery_data$prey_group<- as.factor(as.character(fishery_data$prey_group))

fishery_data$sex<- as.factor(as.character(fishery_data$sex))

fishery_data$location<- as.factor(as.character(fishery_data$location))

fishery_data$best_estimate<-
  as.factor(as.character(fishery_data$best_estimate))

fishery_data$fishing_season<- 
  as.factor(as.character(fishery_data$fishing_season))

summary(fishery_data)

# 1.3 Exclude North American samples----
fishery_data_2<- subset(fishery_data, best_estimate != "North_America")
table(fishery_data_2$best_estimate, useNA = "always")
fishery_data_2$best_estimate<- factor(fishery_data_2$best_estimate)
table(fishery_data_2$best_estimate, useNA = "always")




                    ############################################################
                    #                                                          #
                    #     2.0 Grouping of fishery data by gross sea age----    #
                    #                                                          #
                    ############################################################

table(fishery_data_2$best_estimate, fishery_data_2$sea_age)

# 2.1 Create a gross sea age variable----
fishery_data_2$gross_sea_age<- 
  as.factor(ifelse(fishery_data_2$sea_age == 1, "1SW", "MSW"))

summary(fishery_data_2)

# 2.2 Group fishery data by gross sea age, nested within 
# each reporting group----
fishery_data_grouped_1<- fishery_data_2 %>%
  group_by(best_estimate) %>%
  subset(gross_sea_age == "1SW") %>%
  summarise(n_1SW = length(gross_sea_age))

# 2.3 Add in zero counts of 1SW fish manually----
fishery_data_grouped_1b<- data.frame(matrix(ncol = 2, nrow = 6))
names(fishery_data_grouped_1b)[1]<- "best_estimate"
names(fishery_data_grouped_1b)[2]<- "n_1SW"

fishery_data_grouped_1b[1, 1]<- "Eastern_Barents_White Sea"
fishery_data_grouped_1b[1, 2]<- fishery_data_grouped_1[1, 2]
fishery_data_grouped_1b[2, 1]<- "EasternFinnmark_Kola_P"
fishery_data_grouped_1b[2, 2]<- 0
fishery_data_grouped_1b[3, 1]<- "England_Ireland"
fishery_data_grouped_1b[3, 2]<- fishery_data_grouped_1[2, 2]
fishery_data_grouped_1b[4, 1]<- "Northern_Norway"
fishery_data_grouped_1b[4, 2]<- 0
fishery_data_grouped_1b[5, 1]<- "Southern_Norway"
fishery_data_grouped_1b[5, 2]<- fishery_data_grouped_1[3, 2]
fishery_data_grouped_1b[6, 1]<- "Teno"
fishery_data_grouped_1b[6, 2]<- 0

fishery_data_grouped_2<- fishery_data_2 %>%
  group_by(best_estimate) %>%
  subset(gross_sea_age == "MSW") %>%
  summarise(n_MSW = length(gross_sea_age))

fishery_data_grouped_3<- 
  cbind.data.frame(fishery_data_grouped_1b, fishery_data_grouped_2[, 2])

fishery_data_grouped_3$total_abundance<- 
  fishery_data_grouped_3$n_1SW + fishery_data_grouped_3$n_MSW

fishery_data_grouped_3$prop_MSW<- 
  fishery_data_grouped_3$n_MSW / fishery_data_grouped_3$total_abundance

fishery_data_grouped_3$prop_oneSW<- 1 - fishery_data_grouped_3$prop_MSW




                    ############################################################
                    #                                                          #
                    #             3.0 Estimate uncertainty around              #
                    #             sea ages from fisheries data----             #
                    #                                                          #
                    ############################################################

# 3.1 Create dataframe----
binom<- rbind(c(107, 10, 324, 115, 785, 24), c(3, 0, 130, 0, 32, 0))

d<- data.frame(pMSW = rep(rep(c(1, 0),6), binom),
           group = rep(c("EBWS","EFKP", "EI", "NN", "SN", "T"), colSums(binom)))

# 3.2 Set prior----
n<- 6
prior1<- list(B = list(mu = rep(0,6), V = diag(n) *
                         
                         + (1 + pi^2/3)), R = list(V = 1, fix = 1))

# 3.3 Run model----
mod_1<- MCMCglmm(pMSW~-1+group, data=d, verbose=TRUE, family="categorical",
                
                nitt=2000000, burnin=500000, thin=150, prior=prior1)

# Purpose of second model run is to check convergence of chains----
mod_2<- MCMCglmm(pMSW~-1+group, data=d, verbose=TRUE, family="categorical",
               
                 nitt=2000000, burnin=500000, thin=150, prior=prior1)

# 3.4 Save models----
setwd("C:\\Users\\osuronan\\Dropbox\\Research\\Faroese salmon\\Data\\Data for analysis\\model")
#saveRDS(mod_1, file = "uncertainty_estimates_1.rds")
#saveRDS(mod_2, file = "uncertainty_estimates_2.rds")

# 3.5 Load models and check summary statistics and convergence----
mod_1<- readRDS("uncertainty_estimates_1.rds")
mod_2<- readRDS("uncertainty_estimates_2.rds")

summary(mod_1)
plot(mod_1$Sol)
diag(autocorr(mod_1$Sol)[2,,]) 
autocorr.plot(mod_1$Sol)

summary(mod_2)
plot(mod_2$Sol)
diag(autocorr(mod_2$Sol)[2,,]) 
autocorr.plot(mod_2$Sol)
gelman.diag(mcmc.list(as.mcmc(mod_1$Sol), as.mcmc(mod_2$Sol)))

# 3.6 Build dataframe of point estimates and credible intervals----
summary(mod_1$Sol)

model_estimates<- data.frame(matrix(nrow = 6, ncol = 4))
colnames(model_estimates)<- c("reporting_group", "median", "LCI", "UCI")

model_estimates[1, 1]<- "EBWS"
model_estimates[2, 1]<- "EFKP"
model_estimates[3, 1]<- "EI"
model_estimates[4, 1]<- "NN"
model_estimates[5, 1]<- "SN"
model_estimates[6, 1]<- "Teno"

model_estimates[1, 2]<- 3.824 
model_estimates[2, 2]<- 3.147     
model_estimates[3, 2]<- 1.095     
model_estimates[4, 2]<- 5.260     
model_estimates[5, 2]<- 3.619     
model_estimates[6, 2]<- 3.874     

model_estimates[1, 3]<- 2.8859     
model_estimates[2, 3]<- 1.1199     
model_estimates[3, 3]<- 0.8573     
model_estimates[4, 3]<- 3.7748     
model_estimates[5, 3]<- 3.2661     
model_estimates[6, 3]<- 2.1893     

model_estimates[1, 4]<- 5.009     
model_estimates[2, 4]<- 5.826    
model_estimates[3, 4]<- 1.342    
model_estimates[4, 4]<- 7.504     
model_estimates[5, 4]<- 4.001     
model_estimates[6, 4]<- 6.377     

# 3.7 Get logit of point estimates and 2.5% and 97.5% quantiles----
# Build function that converts data from link to data scale----
logit_to_mean<- function(x) {
  
  1/(1+exp(-x))
  
}

# Check function works----
1/(1+exp(-3.824))
logit_to_mean(3.824)

model_estimates$logit_median<- logit_to_mean(model_estimates$median)
model_estimates$logit_LCI<- logit_to_mean(model_estimates$LCI)
model_estimates$logit_UCI<- logit_to_mean(model_estimates$UCI)




                    ############################################################
                    #                                                          #
                    #                  4.0 Markov chains for                   #
                    #           Faroes MSW proportion estimates----            #
                    #                                                          #
                    ############################################################

# 4.1 Transpose and get full Markov chains of proportion estimates----
# EBWS----
prop_EBWS_MSW<-
  as.data.frame((sort(logit_to_mean(mod_1$Sol[, 1]))))

prop_EBWS_1SW<-
  as.data.frame((sort(1-logit_to_mean(mod_1$Sol[, 1]))))

# EFNK----
prop_EFNK_MSW<-
  as.data.frame((sort(logit_to_mean(mod_1$Sol[, 2]))))

prop_EFNK_1SW<-
  as.data.frame((sort(1-logit_to_mean(mod_1$Sol[, 2]))))

# IUK----
prop_IUK_MSW<-
  as.data.frame((sort(logit_to_mean(mod_1$Sol[, 3]))))

prop_IUK_1SW<-
  as.data.frame((sort(1-logit_to_mean(mod_1$Sol[, 3]))))

# NN----
prop_NN_MSW<-
  as.data.frame((sort(logit_to_mean(mod_1$Sol[, 4]))))

prop_NN_1SW<-
  as.data.frame((sort(1-logit_to_mean(mod_1$Sol[, 4]))))

# SN----
prop_SN_MSW<-
  as.data.frame((sort(logit_to_mean(mod_1$Sol[, 5]))))

prop_SN_1SW<-
  as.data.frame((sort(1-logit_to_mean(mod_1$Sol[, 5]))))

# Teno----
prop_T_MSW<-
  as.data.frame((sort(logit_to_mean(mod_1$Sol[, 6]))))

prop_T_1SW<-
  as.data.frame((sort(1-logit_to_mean(mod_1$Sol[, 6]))))

# 4.2 Build dataframe of Faroes' Markov chains----
faroes_prop_markov<- cbind.data.frame(prop_EBWS_MSW, prop_EBWS_1SW, 
                                      prop_EFNK_MSW, prop_EFNK_1SW,
                                      prop_IUK_MSW, prop_IUK_1SW, 
                                      prop_NN_MSW,  prop_NN_1SW, 
                                      prop_SN_MSW, prop_SN_1SW, 
                                      prop_T_MSW, prop_T_1SW)

colnames(faroes_prop_markov)<- c("prop_EBWS_MSW_f", "prop_EBWS_1SW_f", 
                                 "prop_EFNK_MSW_f", "prop_EFNK_1SW_f",
                                 "prop_IUK_MSW_f", "prop_IUK_1SW_f",
                                 "prop_NN_MSW_f", "prop_NN_1SW_f",
                                 "prop_SN_MSW_f", "prop_SN_1SW_f",
                                 "prop_T_MSW_f", "prop_T_1SW_f")




                    ############################################################
                    #                                                          #
                    #   5.0 ICES deep sea pre-fisheries abundance data----     #
                    #                                                          #
                    ############################################################

load("C:/Users/osuronan/Dropbox/Research/Faroese salmon/Data/Data for analysis/pfa.RData")

# 5.1 Calculate 1SW and MSW proportions for each reporting group---- 

# EBWS----
# 1SW----
totals_EBWS_1SW_t<- t(Russia_East_pfa1SW)

totals_EBWS_1SW_totals<- as.data.frame(totals_EBWS_1SW_t[, 1] + 
                                       totals_EBWS_1SW_t[, 2])

# MSW----
totals_EBWS_MSW_t<- t(Russia_East_pfaMSW)

totals_EBWS_MSW_totals<- as.data.frame(totals_EBWS_MSW_t[, 1] + 
                                       totals_EBWS_MSW_t[, 2])

prop_MSW_EBWS_pfa<- as.data.frame(totals_EBWS_MSW_totals[, 1] / 
                    (totals_EBWS_1SW_totals[, 1] + totals_EBWS_MSW_totals[, 1]))

prop_1SW_EBWS_pfa<- 1-prop_MSW_EBWS_pfa


# EFNK----
# 1SW----
totals_EFNK_1SW_t<- t(Russia_KolaBarents_pfa1SW)

totals_EFNK_1SW_totals<- as.data.frame(totals_EFNK_1SW_t[, 1] + 
                                       totals_EFNK_1SW_t[, 2])

# MSW----
totals_EFNK_MSW_t<- t(Russia_KolaBarents_pfaMSW)

totals_EFNK_MSW_totals<- as.data.frame(totals_EFNK_MSW_t[, 1] + 
                                       totals_EFNK_MSW_t[, 2])

prop_MSW_EFNK_pfa<- as.data.frame(totals_EFNK_MSW_totals[, 1] / 
                    (totals_EFNK_1SW_totals[, 1] + totals_EFNK_MSW_totals[, 1]))

prop_1SW_EFNK_pfa<- 1-prop_MSW_EFNK_pfa


# IUK---
# 1SW----
totals_IUK_1SW_t<- t(UK_Ireland_pfa1SW)

totals_IUK_1SW_totals<- as.data.frame(totals_IUK_1SW_t[, 1] + 
                                      totals_IUK_1SW_t[, 2])

# MSW----
totals_IUK_MSW_t<- t(UK_Ireland_pfaMSW)

totals_IUK_MSW_totals<- as.data.frame(totals_IUK_MSW_t[, 1] +
                                      totals_IUK_MSW_t[, 2])

prop_MSW_IUK_pfa<- as.data.frame(totals_IUK_MSW_totals[, 1] / 
                       (totals_IUK_1SW_totals[, 1] + totals_IUK_MSW_totals[, 1]))

prop_1SW_IUK_pfa<- 1-prop_MSW_IUK_pfa


# NN----
# 1SW----
totals_NN_1SW_t<- t(North_Norway_pfa1SW)

totals_NN_1SW_totals<- as.data.frame(totals_NN_1SW_t[, 1] + 
                                     totals_NN_1SW_t[, 2])

# MSW----
totals_NN_MSW_t<- t(North_Norway_pfaMSW)

totals_NN_MSW_totals<- as.data.frame(totals_NN_MSW_t[, 1] +
                                     totals_NN_MSW_t[, 2])

prop_MSW_NN_pfa<- as.data.frame(totals_NN_MSW_totals[, 1] /
                        (totals_NN_1SW_totals[, 1] + totals_NN_MSW_totals[, 1]))

prop_1SW_NN_pfa<- 1-prop_MSW_NN_pfa


# SN----
# 1SW----
totals_SN_1SW_t<- t(South_Norway_pfa1SW)

totals_SN_1SW_totals<- as.data.frame(totals_SN_1SW_t[, 1] + 
                                     totals_SN_1SW_t[, 2])

# MSW----
totals_SN_MSW_t<- t(South_Norway_pfaMSW)

totals_SN_MSW_totals<- as.data.frame(totals_SN_MSW_t[, 1] +
                                     totals_SN_MSW_t[, 2])

prop_MSW_SN_pfa<- as.data.frame(totals_SN_MSW_totals[, 1] / 
                        (totals_SN_1SW_totals[, 1] + totals_SN_MSW_totals[, 1]))

prop_1SW_SN_pfa<- 1-prop_MSW_SN_pfa

# Teno----
# 1SW----
totals_T_1SW_t<- t(Tana_pfa1SW)
totals_T_1SW_totals<- as.data.frame(totals_T_1SW_t[, 1] + 
                                    totals_T_1SW_t[, 2])

# MSW----
totals_T_MSW_t<- t(Tana_pfaMSW)
totals_T_MSW_totals<- as.data.frame(totals_T_MSW_t[, 1] + 
                                    totals_T_MSW_t[, 2])

prop_MSW_T_pfa<- as.data.frame(totals_T_MSW_totals[, 1] /
                          (totals_T_1SW_totals[, 1] + totals_T_MSW_totals[, 1]))

prop_1SW_T_pfa<- 1-prop_MSW_T_pfa


coastal_prop_markov<- cbind.data.frame(prop_MSW_EBWS_pfa, prop_1SW_EBWS_pfa, 
                                       prop_MSW_EFNK_pfa, prop_1SW_EFNK_pfa,
                                       prop_MSW_IUK_pfa, prop_1SW_IUK_pfa, 
                                       prop_MSW_NN_pfa, prop_1SW_NN_pfa, 
                                       prop_MSW_SN_pfa, prop_1SW_SN_pfa, 
                                       prop_MSW_T_pfa, prop_1SW_T_pfa)

colnames(coastal_prop_markov)<- c("prop_MSW_EBWS_pfa", "prop_1SW_EBWS_pfa", 
                                  "prop_MSW_EFNK_pfa", "prop_1SW_EFNK_pfa",
                                  "prop_MSW_IUK_pfa", "prop_1SW_IUK_pfa",
                                  "prop_MSW_NN_pfa", "prop_1SW_NN_pfa",
                                  "prop_MSW_SN_pfa", "prop_1SW_SN_pfa",
                                  "prop_MSW_T_pfa", "prop_1SW_T_pfa")

# 5.2 Make both dataframes of the same dimensions----
set.seed(123)
faroes_prop_markov_2<- faroes_prop_markov[-sample(1), ]

props_f_c<- cbind.data.frame(faroes_prop_markov_2, coastal_prop_markov)
dim(props_f_c)

head(props_f_c)

# 5.3 Point estimate and 95% credible interval for PFA of each
# reporting group----

# EBWS----
prop_EBWS_MSW_pfa_PE<-
  as.data.frame(t(sort(prop_MSW_EBWS_pfa[, 1])[c(250,5000,9750)]))

# EFNK----
prop_EFNK_MS_pfa_PE<-
  as.data.frame(t(sort(prop_MSW_EFNK_pfa[, 1])[c(250,5000,9750)]))

# IUK----
prop_IUK_MSW_pfa_PE<-
  as.data.frame(t(sort(prop_MSW_IUK_pfa[, 1])[c(250,5000,9750)]))

# NN----
prop_NN_MSW_pfa_PE<-
  as.data.frame(t(sort(prop_MSW_NN_pfa[, 1])[c(250,5000,9750)]))

# SN----
prop_SN_MSW_pfa_PE<-
  as.data.frame(t(sort(prop_MSW_SN_pfa[, 1])[c(250,5000,9750)]))

# Teno----
prop_T_MSW_pfa_PE<-
  as.data.frame(t(sort(prop_MSW_T_pfa[, 1])[c(250,5000,9750)]))

# 5.4 Make dataframe of all deep seas pre-fisheries abundance 
# proportion estimates----
MSW_props_all<- rbind.data.frame(prop_EBWS_MSW_pfa_PE, prop_EFNK_MS_pfa_PE, 
                                 prop_IUK_MSW_pfa_PE, prop_NN_MSW_pfa_PE,
                                 prop_SN_MSW_pfa_PE, prop_T_MSW_pfa_PE)

colnames(MSW_props_all)<- c("LCI", "median", "UCI")


MSW_props_all$region<- as.factor(model_estimates[, 1])
MSW_props_all$habitat<- as.factor(c("coastal"))

model_estimates$habitat<- as.factor(c("faroes"))

props_MSW_plotting<- data.frame(model_estimates[, c(6, 5, 7, 1, 8)])
colnames(props_MSW_plotting)<- c("LCI", "median", "UCI", "region", "habitat")

props_MSW_plotting_2<- rbind.data.frame(props_MSW_plotting, MSW_props_all)

props_MSW_plotting_2$region<- 
  as.factor(as.character(props_MSW_plotting_2$region))
  
# 5.5 Plot Figure 2a----
# 5.5.1 Arrange factor levels in the order in which they are to be plotted----
props_MSW_plotting_2$region<- 
  factor(props_MSW_plotting_2$region,
         levels = c("EI", "SN", "NN", "Teno", "EFKP", "EBWS"))

# 5.5.2 Plot----
fig_2a<- ggplot(props_MSW_plotting_2,
                aes(x=as.factor(region), y=median, group = region))+
  
  geom_point(position= position_dodge(width = .9),
             size = 4,
             aes(shape = habitat, colour = region)) +
  
  geom_errorbar(aes(ymin = LCI, ymax = UCI, colour = region), 
                position="dodge", size = 1, width = 0)+
  
 xlab("")+
 ylab("")+ 
 theme_bw()+
 scale_color_manual(values = c("black", "black", "black", 
                               "black", "black", "black"))+
  
 scale_y_continuous(limits = c(0, 1), labels = label_number(accuracy = 0.01))+
  
 theme(plot.margin = unit(c(1,1,1,1), units = , "cm"),
       text = element_text(size = 20),
       axis.text = element_text( colour = "black"),
      legend.position = "none") +
  scale_x_discrete(labels=c("Ireland & UK",
                            "Southern Norway",
                            "Northern Norway", 
                            "Teno", 
                            "Eastern Finnmark\nNorth Kola", 
                            "Eastern Barents\nWhite Sea")) +
  scale_shape_manual(values = c(19, 1))


setwd("C:\\Users\\osuronan\\Dropbox\\Research\\Faroese salmon\\Faroese population assignment work\\Figures and tables\\Tables")

#write.csv2(props_MSW_plotting_2, "proportion_estimates.csv")




                    ############################################################
                    #                                                          #
                    #          6.0 Consistency of pattern analysis----         #
                    #                                                          #
                    ############################################################

load("C:/Users/osuronan/Dropbox/Research/Faroese salmon/Data/Data for analysis/pfa1SW2SW3SW.RData")

# 6.1 Northern Norway----
# For each sea age, calculate the column totals----
North_Norway_pfa1SW<- as.data.frame(North_Norway_pfa1SW)

# 1SW----
North_Norway_pfa1SW_totals<- North_Norway_pfa1SW %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total")))

# 2SW----
North_Norway_pfa2SW<- as.data.frame(North_Norway_pfa2SW)

North_Norway_pfa2SW_totals<- North_Norway_pfa2SW %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total")))

# 3SW----
North_Norway_pfa3SW<- as.data.frame(North_Norway_pfa3SW)

North_Norway_pfa3SW_totals<- North_Norway_pfa3SW %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total")))

totals_NN<- rbind(North_Norway_pfa1SW_totals[3, ],
                  North_Norway_pfa2SW_totals[3, ], 
                  North_Norway_pfa3SW_totals[3, ])

totals_NN_2<- totals_NN %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total")))

row.names(totals_NN_2)<- 1:nrow(totals_NN_2)

# 6.2 Southern Norway----
South_Norway_pfa1SW<- as.data.frame(South_Norway_pfa1SW)

# 1SW----
South_Norway_pfa1SW_totals<- South_Norway_pfa1SW %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total")))

# 2SW----
South_Norway_pfa2SW<- as.data.frame(South_Norway_pfa2SW)

South_Norway_pfa2SW_totals<- South_Norway_pfa2SW %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total")))

# 3SW----
South_Norway_pfa3SW<- as.data.frame(South_Norway_pfa3SW)

South_Norway_pfa3SW_totals<- South_Norway_pfa3SW %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total")))

totals_SN<- rbind(South_Norway_pfa1SW_totals[3, ],
                  South_Norway_pfa2SW_totals[3, ], 
                  South_Norway_pfa3SW_totals[3, ])

totals_SN_2<- totals_SN %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total")))

row.names(totals_SN_2)<- 1:nrow(totals_SN_2)

# 6.3 Get relevant fishery data from Faroes samples----
fishery_data_3<- subset(fishery_data,
                        best_estimate == "Northern_Norway" | 
                        best_estimate == "Southern_Norway")

fishery_data_3$best_estimate<- factor(fishery_data_3$best_estimate)

table(fishery_data_3$sea_age, fishery_data_3$best_estimate, useNA = "always")


# 6.4 Estimate uncertainty around sea ages from fisheries data----

# 6.4.1 Create dataframe----
binom<- rbind(c(291, 71), c(494, 44))

d<- data.frame(p = rep(rep(c(1, 0),2), binom),
               group = rep(c("SN", "NN"), colSums(binom)))

table(d)
head(d)
table(d$p, d$group)

# 6.4.2 Design prior----
n<- 2
prior_2<- list(B = list(mu = rep(0,2), V = diag(n) *
                         
                         + (1 + pi^2/3)), R = list(V = 1, fix = 1))

# 6.4.3 Run model----
mod_Norway_1<- MCMCglmm(p~-1+group, data=d, verbose=TRUE, family="categorical",
                    nitt=2000000, burnin=500000, thin=150, prior=prior_2)

summary(mod_Norway_1)
plot(mod_Norway_1$Sol)
diag(autocorr(mod_Norway_1$Sol)[2,,]) 
autocorr.plot(mod_Norway_1$Sol)

# NN 3SW vs 2SW proportions in Faroes----
head(mod_Norway_1$Sol[, 1])

# SN 3SW vs 2SW proportions in Faroes----
head(mod_Norway_1$Sol[, 2])

# 6.4.4 Save model----
setwd("C:\\Users\\osuronan\\Dropbox\\Research\\Faroese salmon\\Data\\Data for analysis\\model")
#saveRDS(mod_Norway_1, file = "uncertainty_estimates_Norway.rds")

mod_Norway_1<- readRDS("uncertainty_estimates_Norway.rds")


# 6.5 Calculate likelihoods for Southern Norway----
summary(mod_Norway_1)

# 6.5.1 Deep sea pre-fisheries proportions for 2SW and 3SW salmon---- 
totals_SN_T<- t(totals_SN_2)
head(totals_SN_T)

# 6.5.2 Ratio of 3SW to 2SW in the high seas for Southern Norway----
SN_3SWto2SWprop_highseas_SN<- totals_SN_T[,3] / totals_SN_T[,2]  

# 6.5.3 Ratio of 3SW to 2SW around the Faroes for Southern Norway----
SN_3SWto2SWprop_faroes_SN<- 
 (1-logit_to_mean(mod_Norway_1$Sol[, 2])) / logit_to_mean(mod_Norway_1$Sol[, 2]) 


sort(SN_3SWto2SWprop_faroes_SN)[c(250,5000,9750)]
# 1.9 times less 3SW around Faroes compared to 2SW for Southern Norway----

sort(SN_3SWto2SWprop_highseas_SN)[c(250,5000,9750)]
# 0.23 times less 3SW in the high seas compared to 2SW for Southern Norway----

SN_3SWto2SWprop_faroes_SN_2<- SN_3SWto2SWprop_faroes_SN[-sample(1)]

increase_in_3SW_prop_in_MSW_group_SN<- 
  SN_3SWto2SWprop_faroes_SN_2 / SN_3SWto2SWprop_highseas_SN 

sort(increase_in_3SW_prop_in_MSW_group_SN)[c(250,5000,9750)]
# 3SW fish are 8.31 times more likely around the Faroes compared to in the deep 
# seas for Southern Norway----


# 6.6 Calculate likelihoods for Northern Norway----

# 6.6.1 Deep sea pre-fisheries proportions for 2SW and 3SW salmon---- 
totals_NN_T<- t(totals_NN_2)
head(totals_NN_T)

# 6.6.2 Ratio of 3SW to 2SW in the high seas for Northern Norway----
SN_3SWto2SWprop_highseas_NN<- totals_NN_T[,3] / totals_NN_T[,2]

# 6.6.3 Ratio of 3SW to 2SW around the Faroes for Northern Norway----
NN_3SWto2SWprop_faroes_NN<- 
 (1-logit_to_mean(mod_Norway_1$Sol[, 1])) / logit_to_mean(mod_Norway_1$Sol[, 1]) 

sort(NN_3SWto2SWprop_faroes_NN)[c(250,5000,9750)] 

sort(SN_3SWto2SWprop_highseas_NN)[c(250,5000,9750)] 

NN_3SWto2SWprop_faroes_NN_2<- NN_3SWto2SWprop_faroes_NN[-sample(1)]

increse_in_3SW_prop_in_MSW_group_NN<-
  NN_3SWto2SWprop_faroes_NN_2 / SN_3SWto2SWprop_highseas_NN 

sort(increse_in_3SW_prop_in_MSW_group_NN)[c(250,5000,9750)] 
# 3SW fish are 1.76 times more likely around the Faroes compared to in the deep 
# seas for Northern Norway----


# 6.7 Arrange data for plotting----
prop2SW_FO_NN<-
  sort(1-logit_to_mean(mod_Norway_1$Sol[, 1]))[c(250,5000,9750)]

prop3SW_FO_NN<-
  sort(logit_to_mean(mod_Norway_1$Sol[, 1]))[c(250,5000,9750)]

prop2SW_FO_SN<-
  sort(1-logit_to_mean(mod_Norway_1$Sol[, 2]))[c(250,5000,9750)]

prop3SW_FO_SN<-
  sort(logit_to_mean(mod_Norway_1$Sol[, 2]))[c(250,5000,9750)]


prop3SW_ds_NN<- 
  sort(totals_NN_T[,3] / (totals_NN_T[,3] + totals_NN_T[,2]))[c(250,5000,9750)]

prop2SW_ds_NN<- 
  sort(totals_NN_T[,2] / (totals_NN_T[,3] + totals_NN_T[,2]))[c(250,5000,9750)]

prop3SW_ds_SN<- 
  sort(totals_SN_T[,3] / (totals_SN_T[,3] + totals_SN_T[,2]))[c(250,5000,9750)]

prop2SW_ds_SN<- 
  sort(totals_SN_T[,2] / (totals_SN_T[,3] + totals_SN_T[,2]))[c(250,5000,9750)]


median_estimates<- rbind.data.frame(prop3SW_ds_SN [2],  prop3SW_FO_SN [2], 
                                    prop2SW_ds_SN [2],  prop2SW_FO_SN [2], 
                                    prop3SW_ds_NN [2],  prop3SW_FO_NN [2], 
                                    prop2SW_ds_NN [2],  prop2SW_FO_NN [2])

credible_intervals<- rbind(prop3SW_ds_SN [c(1,3)],  prop3SW_FO_SN [c(1,3)], 
                           prop2SW_ds_SN [c(1,3)],  prop2SW_FO_SN [c(1,3)], 
                           prop3SW_ds_NN [c(1,3)],  prop3SW_FO_NN [c(1,3)],  
                           prop2SW_ds_NN [c(1,3)],  prop2SW_FO_NN [c(1,3)])

split_analysis_plotting<- cbind.data.frame(median_estimates, credible_intervals)

colnames(split_analysis_plotting)<- c("median", "LCI", "UCI")

split_analysis_plotting$region<- 
  as.factor(c( "SN", "SN" , "SN", "SN", "NN", "NN" , "NN", "NN"))


split_analysis_plotting$habitat<- as.factor(c("coastal", "faroes" , 
                                              "coastal", "faroes",
                                              "coastal", "faroes" , 
                                              "coastal", "faroes"))

split_analysis_plotting$sea_age<- as.factor(c("3SW", "3SW" , "2SW", "2SW",
                                              "3SW", "3SW" , "2SW", "2SW"))

split_analysis_plotting$sea_age_region<- 
  factor(paste(split_analysis_plotting$sea_age, split_analysis_plotting$region,
               sep = "_"))

summary(split_analysis_plotting)

# 6.8 Plot----
ggplot(split_analysis_plotting, 
       aes(x=sea_age_region, y=median, color = region, fill = habitat))+
  
  geom_point(position=position_dodge(width=0.3), size=4, aes(shape = habitat))+
  
  geom_errorbar(aes(ymin=LCI, ymax=UCI),
                position = position_dodge(0.3), width = 0)+
  
  scale_shape_manual(values = c(1, 19))+
  xlab("")+
  ylab("")+ 
  theme_bw()+
  theme(plot.margin = unit(c(1,1,1,1), units = , "cm"),
        legend.position = "none",
        text = element_text(size = 20),
        axis.text = element_text( colour = "black"))+
  scale_x_discrete(labels=c("2SW\nNorthern Norway",
                            "2SW\nSouthern Norway",
                            "3SW\nNorthern Norway",
                            "3SW\nSouthern Norway"))+
  scale_color_manual(values = c("black", "black", "black", 
                                "black", "black", "black"))





                    ############################################################
                    #                                                          #
                    #         7.0 Estimate 1SW and MSW occurence rates         #
                    #            around the Faroe Islands relative             #
                    #   to the 1SW and MSW pfa for each reporting group----    #
                    #                                                          #
                    ############################################################

# 7.1 EBWS----
prop_EBWS_MSW_2<- prop_EBWS_MSW[-sample(1), ]

# 1SW----
ER_1SW_EBWS<-  (110*((1-prop_EBWS_MSW_2)/totals_EBWS_1SW_totals))*1000
ER_1SW_EBWS_PE<- sort(ER_1SW_EBWS[, 1])[c(250,5000,9750)]

# MSW----
ER_MSW_EBWS<-  (110*((prop_EBWS_MSW_2)/totals_EBWS_MSW_totals))*1000
ER_MSW_EBWS_PE<- sort(ER_MSW_EBWS[, 1])[c(250,5000,9750)]


# 7.2 EFNK----
prop_EFNK_MSW_2<- prop_EFNK_MSW[-sample(1), ]

# 1SW----
ER_1SW_EFNK<-  (10*((1-prop_EFNK_MSW_2)/totals_EFNK_1SW_totals))*1000
ER_1SW_EFNK_PE<- sort(ER_1SW_EFNK[, 1])[c(250,5000,9750)]

# MSW----
ER_MSW_EFNK<-  (10*((prop_EFNK_MSW_2)/totals_EFNK_MSW_totals))*1000
ER_MSW_EFNK_PE<- sort(ER_MSW_EFNK[, 1])[c(250,5000,9750)]


# 7.3 IUK----
prop_IUK_MSW_2<- prop_IUK_MSW[-sample(1), ]

# 1SW----
ER_1SW_IUK<-  (454*((1-prop_IUK_MSW_2)/totals_IUK_1SW_totals))*1000
ER_1SW_IUK_PE<- sort(ER_1SW_IUK[, 1])[c(250,5000,9750)]

# MSW----
ER_MSW_IUK<-  (454*((prop_IUK_MSW_2)/totals_IUK_MSW_totals))*1000
ER_MSW_IUK_PE<- sort(ER_MSW_IUK[, 1])[c(250,5000,9750)]


# 7.4 NN----
prop_NN_MSW_2<- prop_NN_MSW[-sample(1), ]

# 1SW----
ER_1SW_NN<-  (115*((1-prop_NN_MSW_2)/totals_NN_1SW_totals))*1000
ER_1SW_NN_PE<- sort(ER_1SW_NN[, 1])[c(250,5000,9750)]

# MSW----
ER_MSW_NN<-  (115*((prop_NN_MSW_2)/totals_NN_MSW_totals))*1000
ER_MSW_NN_PE<- sort(ER_MSW_NN[, 1])[c(250,5000,9750)]


# 7.5 SN----
prop_SN_MSW_2<- prop_SN_MSW[-sample(1), ]

# 1SW----
ER_1SW_SN<-  (817*((1-prop_SN_MSW_2)/totals_SN_1SW_totals))*1000
ER_1SW_SN_PE<- sort(ER_1SW_SN[, 1])[c(250,5000,9750)]

# MSW----
ER_MSW_SN<-  (817*((prop_SN_MSW_2)/totals_SN_MSW_totals))*1000
ER_MSW_SN_PE<- sort(ER_MSW_SN[, 1])[c(250,5000,9750)]


# 7.6 Teno----
prop_T_MSW_2<- prop_T_MSW[-sample(1), ]

# 1SW----
ER_1SW_T<-  (24*((1-prop_T_MSW_2)/totals_T_1SW_totals))*1000
ER_1SW_T_PE<- sort(ER_1SW_T[, 1])[c(250,5000,9750)]

# MSW----
ER_MSW_T<-  (24*((prop_T_MSW_2)/totals_T_MSW_totals))*1000
ER_MSW_T_PE<- sort(ER_MSW_T[, 1])[c(250,5000,9750)]

# 7.7 Build dataframe for plotting----
exploitation_rates<- rbind.data.frame(ER_1SW_EBWS_PE, ER_MSW_EBWS_PE, 
                                      ER_1SW_EFNK_PE, ER_MSW_EFNK_PE,
                                      ER_1SW_IUK_PE,  ER_MSW_IUK_PE,
                                      ER_1SW_NN_PE,   ER_MSW_NN_PE,
                                      ER_1SW_SN_PE,   ER_MSW_SN_PE,
                                      ER_1SW_T_PE,    ER_MSW_T_PE)

colnames(exploitation_rates)<- c("LCI", "median", "UCI")

exploitation_rates$reporting_group <- 
  as.factor(c(rep("EBWS", 2), rep("EFNK", 2), rep("IUK", 2), 
              rep("NN", 2), rep("SN", 2), rep("Teno", 2)))

exploitation_rates$sea_age<- as.factor(rep(c("1SW", "MSW"), 6))

exploitation_rates$reporting_group<- 
  factor(exploitation_rates$reporting_group,
         levels = c("IUK", "SN", "NN", "Teno", "EFNK", "EBWS"))

# 7.8 Plot----
fig_2b<- ggplot(exploitation_rates, 
       aes(x=as.factor(reporting_group), y=median, group = reporting_group))+
  geom_point(position= position_dodge(width = .9), size = 4, 
             aes(shape = sea_age, colour = reporting_group)) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI, colour = reporting_group), 
                position="dodge", size = 1, width = 0) +
  facet_zoom(ylim = c(0, .035), zoom.data = ifelse(a <= .035, NA, FALSE))+
  xlab("")+
  ylab("")+ 
  theme_bw()+
  scale_x_discrete(labels=c("Ireland & UK",
                            "Southern Norway",
                            "Northern Norway", 
                            "Teno", 
                            "Eastern Finnmark\nNorth Kola", 
                            "Eastern Barents\nWhite Sea")) +
                      scale_shape_manual(values = c(1, 19)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01))+
  theme(legend.position = "none",
        text = element_text(size = 20),
        axis.text = element_text( colour = "black")) +
  scale_color_manual(values = c("black", "black", "black", 
                                "black", "black", "black"))




                    ############################################################
                    #                                                          #
                    #                       8.0 Maps                           #
                    #                                                          #
                    ############################################################

setwd("C:\\Users\\osuronan\\Dropbox\\Research\\Faroese salmon\\Data\\Data for analysis")

# 8.1 Load required dataframes---- 
data_1<- read_csv("faroese_fisheries_data_for_analysis.csv")

pop_assignments_Final<- read_delim("pop_assignments_Final.csv", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)

river_coordinates<- read.delim("C:/Users/osuronan/Dropbox/Research/Faroese salmon/Data/Data for analysis/faroese diet-LH adaptation data files/Reporting_groups_Bourret_at_al.txt")

coordinates_ozerov<- read.csv("C:/Users/osuronan/Dropbox/Research/Faroese salmon/Data/Data for analysis/coordinates_ozerov.txt", sep="")

baseline_coordinates<- read_csv("C:/Users/osuronan/Dropbox/Research/Faroese salmon/Faroese population assignment work/Figures and tables/Tables/baseline_coordinates.txt")


# 8.2 Separate longline location data by fishing season (autumn/winter)
# and season (92/93 or 93/94)----
levels(as.factor(data_1$best_estimate))

table(data_1$fishing_season, data_1$season, useNA = "always")

autumn_trawls_92_93<- 
  subset(data_1, season == "92/93" & fishing_season == "autumn")

winter_trawls_92_93<- 
  subset(data_1, season == "92/93" & fishing_season == "winter")

autumn_trawls_93_94<- 
  subset(data_1, season == "93/94" & fishing_season == "autumn")

winter_trawls_93_94<- 
  subset(data_1, season == "93/94" & fishing_season == "winter")


# 8.3 Load map data for longline locations----
world<- ne_countries(scale = "medium", returnclass = "sf")

# 8.4 Map longline locations, separated by fishing season and time of year----
ggplot(data = world)+
  geom_sf()+
  coord_sf(xlim = c(-15, 10), ylim = c(58, 67), expand = T) +
  xlab("\nLongitude")+
  ylab("Latitude\n")+
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(face="bold", size=14, colour = "black"),
        axis.text.y = element_text(face="bold", size=14, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"))+
  
  geom_point(data = autumn_trawls_92_93,
             aes(longitude, latitude, color = "#F8766D"), size = 4, shape = 16) +
  
  geom_point(data = winter_trawls_92_93,
             aes(longitude, latitude, color = "#F8766D"), size = 4, shape = 17) +
  
  geom_point(data = autumn_trawls_93_94,
             aes(longitude, latitude, color = "#00BFC4"), size = 4, shape = 16) +
  
  geom_point(data = winter_trawls_93_94,
             aes(longitude, latitude, color = "#00BFC4"), size = 4, shape = 17)


# 8.5 Map with all populations---
names(pop_assignments_Final)[1]<- "uniqueID3"

match_index<- which(pop_assignments_Final$uniqueID3 %in% data_1$uniqueID3) 
all(data_1$uniqueID3 %in% pop_assignments_Final$uniqueID3)

populations<- pop_assignments_Final[match_index, ]
populations_2<- populations[c(1, 2)]
names(populations_2)[2]<- "population_best_estimate"

table(populations_2$population_best_estimate)

length(unique(populations_2$uniqueID3))
length(unique(data_1$uniqueID3))

table(populations_2$uniqueID3 %in% data_1$uniqueID3)

data_2<- merge(data_1, populations_2, by = "uniqueID3")

# 8.5.1 Required river coordinates----
table(data_2$population_best_estimate)
table(river_coordinates$Population_ID)
table(coordinates_ozerov$population_ID)

river_coordinates_2<- river_coordinates[c(1, 6, 7)]

river_coordinates_3<-
  subset(river_coordinates_2, 
           Population_ID == "Blackwater" |
           Population_ID == "Dart" |
           Population_ID == "Dionard" |
           Population_ID == "Emtsa" |
           Population_ID == "Gaula" |
           Population_ID == "Lardaiselva" |
           Population_ID == "Lebyazhya" |
           Population_ID == "Moy" |
           Population_ID == "NorthEsk" |
           Population_ID == "Numedalsiagen" |
           Population_ID == "Pongoma" |
           Population_ID == "Tana" |
           Population_ID == "Tuloma" |
           Population_ID == "Tweed" |
           Population_ID == "Yapoma")

names(river_coordinates_3)[1]<- "population_ID"

# 8.5.2 Add reporting groups----
river_coordinates_3$RG<- 
  ifelse(river_coordinates_3$population_ID == "Blackwater" |  
         river_coordinates_3$population_ID == "Moy" |
         river_coordinates_3$population_ID == "Dionard" |
         river_coordinates_3$population_ID == "NorthEsk" |
         river_coordinates_3$population_ID == "Tweed" |
         river_coordinates_3$population_ID == "Dart",
         
         "England_Ireland", 
         
 ifelse(river_coordinates_3$population_ID == "Emtsa" |
        river_coordinates_3$population_ID == "Pongoma" |
        river_coordinates_3$population_ID == "Lebyazhya" |
        river_coordinates_3$population_ID == "Yapoma",
        
        "Eastern_Barents_White_Sea",
                
 ifelse(river_coordinates_3$population_ID == "Gaula" |    
        river_coordinates_3$population_ID == "Numedalsiagen" | 
        river_coordinates_3$population_ID == "Lardaiselva",
        
        "Southern_Norway",
                       
 ifelse(river_coordinates_3$population_ID == "Tana",
        
        "Teno_system",
                              
ifelse(river_coordinates_3$population_ID == "Tuloma",
       
        "EasternFinnmark_Kola_P", NA)))))

# 8.5.3 Remove unnecessary character from coordinates----  
coordinates_ozerov[, 3]<- gsub("Â", "" , coordinates_ozerov[, 3])
coordinates_ozerov[, 4]<- gsub("Â", "" , coordinates_ozerov[, 4])

# 8.5.4 Convert coordinates to decimals----  
coordinates_ozerov[, 3]<- gsub("°", " " , coordinates_ozerov[, 3])
coordinates_ozerov[, 4]<- gsub("°", " " , coordinates_ozerov[, 4])

coordinates_ozerov_2<- coordinates_ozerov[-1, ]

# Following function from https://stackoverflow.com/questions/30879429/how-can-i-convert-degree-minute-sec-to-decimal-in-r 
angle2dec <- function(angle) {
  angle <- as.character(angle)
  x <- do.call(rbind, strsplit(angle, split=' '))
  x <- apply(x, 1L, function(y) {
    y <- as.numeric(y)
    y[1] + y[2]/60
  })
  return(x)
}

coordinates_ozerov_2$Latitude_2<- angle2dec(coordinates_ozerov_2$Latitude)
coordinates_ozerov_2$Longitude_2<- angle2dec(coordinates_ozerov_2$Longitude)

coordinates_ozerov_3<- coordinates_ozerov_2[c(1, 2, 7, 8)]
names(coordinates_ozerov_3)[3]<- "Latitude"
names(coordinates_ozerov_3)[4]<- "Longitude"

coordinates_ozerov_Suma<- coordinates_ozerov[1, c(1:4)]

# 8.5.5 Create a dataframe of all coordinates----
coordinates_all_rivers<- rbind.data.frame(coordinates_ozerov_Suma, 
                                          coordinates_ozerov_3, 
                                          river_coordinates_3)

rownames(coordinates_all_rivers)<- 1:nrow(coordinates_all_rivers)
coordinates_all_rivers$Latitude<- as.numeric(coordinates_all_rivers$Latitude)
coordinates_all_rivers$Longitude<- as.numeric(coordinates_all_rivers$Longitude)
coordinates_all_rivers$RG<- as.factor(as.character(coordinates_all_rivers$RG))

length(unique(data_2$population_best_estimate))
length(unique(coordinates_all_rivers$population_ID))
length(unique(paste(coordinates_all_rivers$Latitude, 
                    coordinates_all_rivers$Latitude)))

coordinates_all_rivers$coordinates<- paste(coordinates_all_rivers$Latitude,
                                           coordinates_all_rivers$Longitude)

table(coordinates_all_rivers$coordinates)
table(coordinates_all_rivers$population_ID)

coordinates_no_Suma<- subset(coordinates_all_rivers, population_ID != "Suma")
#write.csv(coordinates_all_rivers, file = "data_for_table1.csv")

summary(baseline_coordinates)
baseline_coordinates$Population<- 
  as.factor(as.character(baseline_coordinates$Population))

names(baseline_coordinates)[1]<- "RG"

baseline_coordinates$RG<-
  as.factor(as.character(baseline_coordinates$RG))

summary(baseline_coordinates)

# 8.5.6 Map all populations----
basemap(world, bathymetry = T, limits = c(-70, 40, 40, 80)) +
  xlab("\nLongitude")+
  ylab("Latitude\n")+
  theme_bw() +
  scale_color_manual(values = c("#d89000ff", "#0000ffff", "#f8766dff",
                                "#b79f00ff", "#ff6600ff", "#00ba38ff",
                                "#00b0f6ff", "#00bfc4ff", "#619cffff", 
                                "#f564e3ff"))+
  
  theme(legend.position = "none", 
        axis.text.x = element_text(face="bold", size=14, colour = "black"),
        axis.text.y = element_text(face="bold", size=14, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"))+
  
  geom_spatial_point(data = baseline_coordinates, 
                     aes(longitude, latitude, color = RG), size =5, shape = 16)




# 9.0 Total deep sea pre-fisheries abundance for reporting groups----
total_returns_markov_chains<-totals_EBWS_1SW_totals + totals_EBWS_MSW_totals + 
                             totals_EFNK_1SW_totals + totals_EFNK_MSW_totals +
                             totals_IUK_1SW_totals +  totals_IUK_MSW_totals +
                             totals_NN_1SW_totals +   totals_NN_MSW_totals + 
                             totals_SN_1SW_totals +   totals_SN_MSW_totals +
                             totals_T_1SW_totals +    totals_T_MSW_totals


total_returns_PE<- sort(total_returns_markov_chains[, 1])[c(250, 5000 ,9750)]

# Small salmon returns to North America. From Table 4.3.2.1 of 
# WGNAS Report 2021----
small_salmon<- 511.2+307.7

# Large salmon returns to North America. From Table 4.3.2.2 of 
# WGNAS Report 2021----
large_salmon<- 153.5+131.1

all_NA<- (small_salmon+large_salmon)*1000

total_returns_PE[2]/all_NA

# For every 1 North American fish at the Faroes, there are 8.3 European fish.
# Ceteris paribus for every 10 North American fish, there are 83 European fish.
# However, accounting fopr the actual 5.3% assignment rate to North America:

1/((5.3/94.7)/(10/83))
# Therefore, North American fish are 2.2 times less likley to found around 
# the Faroes than European fish----




                    ############################################################
                    #                                                          #
                    #                   10.0 Figure S2----                     #
                    #                                                          #
                    ############################################################

mean(fishery_data$fork_length[fishery_data$sea_age =="1"])

sd(fishery_data$fork_length[fishery_data$sea_age =="1"])

fishery_data$Gross_SA_for_plotting<- 
  as.factor(ifelse(fishery_data$sea_age == "1", "1SW", "MSW"))

table(fishery_data$Gross_SA_for_plotting, useNA = "always")
table(fishery_data$sex, useNA = "always")
fishery_data_all_sexes_known<- subset(fishery_data, is.na(sex) == FALSE)

table(fishery_data_all_sexes_known$sex, useNA = "always")

table(fishery_data_all_sexes_known$Gross_SA_for_plotting, 
      fishery_data_all_sexes_known$sex)

fishery_data_92_93_season<- 
  subset(fishery_data_all_sexes_known, season == "92/93")

fishery_data_93_94_season<- 
  subset(fishery_data_all_sexes_known, season == "93/94")

table(fishery_data_92_93_season$Gross_SA_for_plotting, 
      fishery_data_92_93_season$sex)

table(fishery_data_93_94_season$Gross_SA_for_plotting, 
      fishery_data_93_94_season$sex)




fork_lengths_1<-  ggplot(fishery_data_all_sexes_known,
                         aes(x=Gross_SA_for_plotting, y=fork_length, fill=sex))+
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  
  theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size = 30),
        axis.text = element_text( colour = "black", size = 20)) +
  xlab("") + ylab("")

fork_lengths_2<- ggplot(fishery_data_92_93_season, 
                        aes(x=Gross_SA_for_plotting, y=fork_length, fill=sex))+
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  
  theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size = 30),
        axis.text = element_text( colour = "black", size = 20)) +
  xlab("") + ylab("")

fork_lengths_3<-  ggplot(fishery_data_93_94_season, 
                         aes(x=Gross_SA_for_plotting, y=fork_length, fill=sex))+
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  
  theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size = 30),
        axis.text = element_text( colour = "black", size = 20)) +
  xlab("") + ylab("")

plot_layout<- rbind(c(1), c(2, 3))
grid.arrange(fork_lengths_1, fork_lengths_2, fork_lengths_3,
             layout_matrix = plot_layout)




                    ############################################################
                    #                                                          #
                    #             11.0 Supplementary analysis----              #
                    #                                                          #
                    ############################################################

# 11.1 Choose reporting groups with substantial numbers of 1SW fish----
fishery_data_1sw<- 
  fishery_data[fishery_data$sea_age == 1 & (fishery_data$best_estimate == "England_Ireland" | fishery_data$best_estimate == "Southern_Norway" ),]

# 11.2 Control for seasonal effects since fish can grow 
# across the fishing period----
table(fishery_data_1sw$best_estimate, fishery_data_1sw$fishing_season)

# 11.3 Test for differences in fork length within season and 
# between reporting groups----
fork_length_model<- lm(log(fishery_data_1sw$fork_length) ~ fishery_data_1sw$fishing_season/fishery_data_1sw$best_estimate)
summary(fork_length_model)

tbl_regression(fork_length_model)
