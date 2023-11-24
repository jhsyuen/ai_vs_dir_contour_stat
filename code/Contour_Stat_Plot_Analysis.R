
#-------------------------
# R Code for statistical analysis of multiple methods of contour generation (AI/DIR)
#
# Input:
#   
#  Excel data with MSC and MDA
#  Excel data with dose difference
#  Excel data with volume difference
#
# R_data that are not output:
#  Various Plots including data/stat_models/Q-Q
#  Text data of stat data  
#
# Output files:
#
#  Text file with stat data "Contour_Analysis.txt"
#  Plot of data
#
# R code submitted as part of manuscript submitted
# "AI compared to deformable image registration propagation techniques for contouring:" 
# "feasibility and accuracy for head and neck replanning" 
#
# 
#-------------------------


library(readxl)
library(tidyverse)
library(emmeans)
library(glmmTMB)
library(ggpubr)


setwd("../data")

#No statistical analysis of user training grade data

#UserTrain <- read_excel("Rdata_Grade3.xlsx")
#summary(UserTrain)
#str(UserTrain)
#DLCManGrade <- read_excel("Rdata_GradeContour6.xlsx")
#summary(DLCManGrade)
#str(DLCManGrade)



DSCMDA <- read_excel("Rdata_DSCMDAvF.xlsx")
DSCMDA %>% 
  group_by(Method) %>% 
  summarise(mean_DSC=mean(Dice), sd_DSC=sd(Dice),
            mean_MDA=mean(MDA), sd_MDA=sd(MDA))
summary(DSCMDA)


#---------------------------------------------
#Get contour data
contour_data <- read_xlsx("Rdata_DSCMDAvF.xlsx")


#-------------------------
#parotid - MDA

organ <- "Parotid"
metric = "MDA"
#-------------------------

#Contour_tab: Count of how much data each method has for each patient
contour_tab <- contour_data %>% 
  filter(Contour == organ) %>% 
  select(-metric) %>% 
  group_by( HN1, Method) %>%
  count() %>% 
  pivot_wider(names_from = Method, values_from = n)
contour_tab

#Contour_df: isolate required data based on method, patient, metric)
contour_df <- contour_data %>% 
  filter(Contour == organ ) %>% 
  rename(Metric = all_of(metric)) %>% 
  transmute(Method, Patient = factor(HN1), Metric) 

#run metric model based on glmmTMB
metric_mod <- glmmTMB(Metric ~ Method + (1|Patient), data = contour_df)
summary(metric_mod)

#check assumption with plots
plot(predict(metric_mod, re.form = NA), residuals(metric_mod))
qqnorm(residuals(metric_mod))

#model prediction with EMMIP
metric_predict <- emmip(metric_mod, ~ Method , CIs = TRUE, plotit = FALSE) 

ggplot(data = metric_predict, aes(Method, yvar) ) + 
  geom_point() + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0, linewidth = 3, col = "lightgrey")+
  geom_point(aes(Method, Metric, color = Patient), data = contour_df) +
  ylab(metric)+
  theme_classic()

#use EMMANS and then CONTRAST
ems <- emmeans(metric_mod, ~ Method)
ems

cont_par_mda <- contrast(ems, "trt.vs.ctrlk", ref = 2) #ref = 2 because VEL_DMP is second in the list above


plot(cont_par_mda, color="red") + 
  theme_classic() + 
  geom_vline(xintercept = 0)

cont=cont_par_mda

summary(cont)
confint(cont)

write("PAROTID MDA",file="Contour_Analysis.txt",append=FALSE)
write.table(summary(cont),file="Contour_Analysis.txt",append=TRUE)
write.table(confint(cont),file="Contour_Analysis.txt",append=TRUE)


#-------------------------
#parotid - DICE

organ <- "Parotid"
metric = "Dice"
#-------------------------

contour_df <- contour_data %>% 
  filter(Contour == organ ) %>% 
  rename(Metric = all_of(metric)) %>% 
  transmute(Method, Patient = factor(HN1), Metric) 

metric_mod <- glmmTMB(Metric ~ Method + (1|Patient), data = contour_df)
summary(metric_mod)

plot(predict(metric_mod, re.form = NA), residuals(metric_mod))
qqnorm(residuals(metric_mod))

metric_predict <- emmip(metric_mod, ~ Method , CIs = TRUE, plotit = FALSE) 

ggplot(data = metric_predict, aes(Method, yvar) ) + 
  geom_point() + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0, size = 3, col = "lightgrey")+
  geom_point(aes(Method, Metric, color = Patient), data = contour_df) +
  ylab(metric)+
  theme_classic()

ems <- emmeans(metric_mod, ~ Method)
ems

cont_par_dsc <- contrast(ems, "trt.vs.ctrlk", ref = 2) #ref = 2 because VEL_DMP is second in the list above

plot(cont_par_dsc) + 
  theme_classic() + 
  geom_vline(xintercept = 0)

cont=cont_par_dsc

summary(cont)
confint(cont)

write("PAROTID DICE",file="Contour_Analysis.txt",append=TRUE)
write.table(summary(cont),file="Contour_Analysis.txt",append=TRUE)
write.table(confint(cont),file="Contour_Analysis.txt",append=TRUE)

#-------------------------
#parotid - HD

organ <- "Parotid"
metric = "HD"
#-------------------------

contour_df <- contour_data %>% 
  filter(Contour == organ ) %>% 
  rename(Metric = all_of(metric)) %>% 
  transmute(Method, Patient = factor(HN1), Metric) 

metric_mod <- glmmTMB(Metric ~ Method + (1|Patient), data = contour_df)
summary(metric_mod)

plot(predict(metric_mod, re.form = NA), residuals(metric_mod))
qqnorm(residuals(metric_mod))

metric_predict <- emmip(metric_mod, ~ Method , CIs = TRUE, plotit = FALSE) 

ggplot(data = metric_predict, aes(Method, yvar) ) + 
  geom_point() + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0, size = 3, col = "lightgrey")+
  geom_point(aes(Method, Metric, color = Patient), data = contour_df) +
  ylab(metric)+
  theme_classic()

ems <- emmeans(metric_mod, ~ Method)
ems

cont_par_hd <- contrast(ems, "trt.vs.ctrlk", ref = 2) #ref = 2 because VEL_DMP is second in the list above

plot(cont_par_hd) + 
  theme_classic() + 
  geom_vline(xintercept = 0)

cont=cont_par_hd

summary(cont)
confint(cont)

write("PAROTID HD",file="Contour_Analysis.txt",append=TRUE)
write.table(summary(cont),file="Contour_Analysis.txt",append=TRUE)
write.table(confint(cont),file="Contour_Analysis.txt",append=TRUE)

#-------------------------
#Brainstem MDA

organ <- "Brainstem"
metric = "MDA"
#-------------------------

contour_tab <- contour_data %>% 
  filter(Contour == organ) %>% 
  select(-metric) %>% 
  group_by( HN1, Method) %>%
  count() %>% 
  pivot_wider(names_from = Method, values_from = n)
contour_tab

contour_df <- contour_data %>% 
  filter(Contour == organ ) %>% 
  rename(Metric = all_of(metric)) %>% 
  transmute(Method, Patient = factor(HN1), Metric) 

metric_mod <- glmmTMB(Metric ~ Method + (1|Patient), data = contour_df)
summary(metric_mod)

plot(predict(metric_mod, re.form = NA), residuals(metric_mod))
qqnorm(residuals(metric_mod))

metric_predict <- emmip(metric_mod, ~ Method , CIs = TRUE, plotit = FALSE) 

ggplot(data = metric_predict, aes(Method, yvar) ) + 
  geom_point() + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0, size = 3, col = "lightgrey")+
  geom_point(aes(Method, Metric, color = Patient), data = contour_df) +
  ylab(metric)+
  theme_classic()

ems <- emmeans(metric_mod, ~ Method)
ems

cont_bs_mda <- contrast(ems, "trt.vs.ctrlk", ref = 2) #ref = 2 because VEL_DMP is second in the list above

plot(cont_bs_mda) + 
  theme_classic() + 
  geom_vline(xintercept = 0)

cont=cont_bs_mda

summary(cont)
confint(cont)

write("BRAINSTEM MDA",file="Contour_Analysis.txt",append=TRUE)
write.table(summary(cont),file="Contour_Analysis.txt",append=TRUE)
write.table(confint(cont),file="Contour_Analysis.txt",append=TRUE)

#-------------------------
#Brainstem DICE

organ <- "Brainstem"
metric = "Dice"

#-------------------------



contour_tab <- contour_data %>% 
  filter(Contour == organ) %>% 
  select(-metric) %>% 
  group_by( HN1, Method) %>%
  count() %>% 
  pivot_wider(names_from = Method, values_from = n)
contour_tab

contour_df <- contour_data %>% 
  filter(Contour == organ ) %>% 
  rename(Metric = all_of(metric)) %>% 
  transmute(Method, Patient = factor(HN1), Metric) 

metric_mod <- glmmTMB(Metric ~ Method + (1|Patient), data = contour_df)
summary(metric_mod)

plot(predict(metric_mod, re.form = NA), residuals(metric_mod))
qqnorm(residuals(metric_mod))

metric_predict <- emmip(metric_mod, ~ Method , CIs = TRUE, plotit = FALSE) 

ggplot(data = metric_predict, aes(Method, yvar) ) + 
  geom_point() + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0, size = 3, col = "lightgrey")+
  geom_point(aes(Method, Metric, color = Patient), data = contour_df) +
  ylab(metric)+
  theme_classic()

ems <- emmeans(metric_mod, ~ Method)
ems

cont_bs_dsc <- contrast(ems, "trt.vs.ctrlk", ref = 2) #ref = 2 because VEL_DMP is second in the list above

plot(cont_bs_dsc) + 
  theme_classic() + 
  geom_vline(xintercept = 0)

cont=cont_bs_dsc

summary(cont)
confint(cont)

write("BRAINSTEM DICE",file="Contour_Analysis.txt",append=TRUE)
write.table(summary(cont),file="Contour_Analysis.txt",append=TRUE)
write.table(confint(cont),file="Contour_Analysis.txt",append=TRUE)

#-------------------------
#Brainstem HD

organ <- "Brainstem"
metric = "HD"

#-------------------------

contour_tab <- contour_data %>% 
  filter(Contour == organ) %>% 
  select(-metric) %>% 
  group_by( HN1, Method) %>%
  count() %>% 
  pivot_wider(names_from = Method, values_from = n)
contour_tab

contour_df <- contour_data %>% 
  filter(Contour == organ ) %>% 
  rename(Metric = all_of(metric)) %>% 
  transmute(Method, Patient = factor(HN1), Metric) 

metric_mod <- glmmTMB(Metric ~ Method + (1|Patient), data = contour_df)
summary(metric_mod)

plot(predict(metric_mod, re.form = NA), residuals(metric_mod))
qqnorm(residuals(metric_mod))

metric_predict <- emmip(metric_mod, ~ Method , CIs = TRUE, plotit = FALSE) 

ggplot(data = metric_predict, aes(Method, yvar) ) + 
  geom_point() + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0, size = 3, col = "lightgrey")+
  geom_point(aes(Method, Metric, color = Patient), data = contour_df) +
  ylab(metric)+
  theme_classic()

ems <- emmeans(metric_mod, ~ Method)
ems

cont_bs_hd <- contrast(ems, "trt.vs.ctrlk", ref = 2) #ref = 2 because VEL_DMP is second in the list above

plot(cont_cont_bs_hd) + 
  theme_classic() + 
  geom_vline(xintercept = 0)

cont=cont_bs_hd

summary(cont)
confint(cont)

write("BRAINSTEM HD",file="Contour_Analysis.txt",append=TRUE)
write.table(summary(cont),file="Contour_Analysis.txt",append=TRUE)
write.table(confint(cont),file="Contour_Analysis.txt",append=TRUE)


#-------------------------
#Parotid VOLDIFF 

organ <- "Parotid"
metric = "Volume_Diff"

#-------------------------

Voldata <- read_excel("Rdata_VolumeF.xlsx")

summary(Voldata)

Voldata %>% 
  group_by(Method) %>% 
  summarise(mean_VD=mean(Volume_Diff), sd_VD=sd(Volume_Diff))

str(Voldata)

#Vol_Tab: Count of how much data each method has for each patient
vol_tab <- Voldata %>% 
  filter(Contour == organ) %>% 
  select(-metric) %>% 
  group_by( Patient, Method) %>%
  count() %>% 
  pivot_wider(names_from = Method, values_from = n)
vol_tab


#vol_df: isolate required data based on method, patient, metric)
vol_df <- Voldata %>% 
  filter(Contour == organ ) %>% 
  rename(Metric = all_of(metric)) %>% 
  transmute(Method, Patient = factor(Patient), Metric) 


#run metric model based on glmmTMB
metric_mod <- glmmTMB(Metric ~ Method + (1|Patient), data = vol_df)
summary(metric_mod)

#check assumption with plots
plot(predict(metric_mod, re.form = NA), residuals(metric_mod))
qqnorm(residuals(metric_mod))

#model prediction with EMMIP
metric_predict <- emmip(metric_mod, ~ Method , CIs = TRUE, plotit = FALSE) 

ggplot(data = metric_predict, aes(Method, yvar) ) + 
  geom_point() + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0, size = 3, col = "lightgrey")+
  geom_point(aes(Method, Metric, color = Patient), data = vol_df) +
  ylab(metric)+
  theme_classic()

#use EMMANS and then CONTRAST
ems <- emmeans(metric_mod, ~ Method)
ems

cont_par_vd <- contrast(ems, "trt.vs.ctrlk", ref = 2) #ref = 2 because VEL_DMP is second in the list above

plot(cont_par_vd) + 
  theme_classic() + 
  geom_vline(xintercept = 0)

cont=cont_par_vd

summary(cont)
confint(cont)

write("PAROTID VOLDIFF",file="Contour_Analysis.txt",append=TRUE)
write.table(summary(cont),file="Contour_Analysis.txt",append=TRUE)
write.table(confint(cont),file="Contour_Analysis.txt",append=TRUE)


#-------------------------
#Brainstem VOLDIFF 

organ <- "Brainstem"
metric = "Volume_Diff"

#-------------------------


#Vol_Tab: Count of how much data each method has for each patient
vol_tab <- Voldata %>% 
  filter(Contour == organ) %>% 
  select(-metric) %>% 
  group_by( Patient, Method) %>%
  count() %>% 
  pivot_wider(names_from = Method, values_from = n)
vol_tab


#vol_df: isolate required data based on method, patient, metric)
vol_df <- Voldata %>% 
  filter(Contour == organ ) %>% 
  rename(Metric = all_of(metric)) %>% 
  transmute(Method, Patient = factor(Patient), Metric) 


#run metric model based on glmmTMB
metric_mod <- glmmTMB(Metric ~ Method + (1|Patient), data = vol_df)
summary(metric_mod)

#check assumption with plots
plot(predict(metric_mod, re.form = NA), residuals(metric_mod))
qqnorm(residuals(metric_mod))

#model prediction with EMMIP
metric_predict <- emmip(metric_mod, ~ Method , CIs = TRUE, plotit = FALSE) 

ggplot(data = metric_predict, aes(Method, yvar) ) + 
  geom_point() + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0, size = 3, col = "lightgrey")+
  geom_point(aes(Method, Metric, color = Patient), data = vol_df) +
  ylab(metric)+
  theme_classic()

#use EMMANS and then CONTRAST
ems <- emmeans(metric_mod, ~ Method)
ems

cont_bs_vd <- contrast(ems, "trt.vs.ctrlk", ref = 2) #ref = 2 because VEL_DMP is second in the list above

plot(cont) + 
  theme_classic() + 
  geom_vline(xintercept = 0)

cont=cont_bs_vd

summary(cont)
confint(cont)

write("BRAINSTEM VOLDIFF",file="Contour_Analysis.txt",append=TRUE)
write.table(summary(cont),file="Contour_Analysis.txt",append=TRUE)
write.table(confint(cont),file="Contour_Analysis.txt",append=TRUE)



#-------------------------------
#Dose data
#------------------------------

Dosedata <- read_excel("Rdata_DoseF.xlsx")

summary(Dosedata)

Dosedata %>% 
  group_by(Method) %>% 
  summarise(mean_DD=mean(Dose_Diff), sd_DD=sd(Dose_Diff))

str(Dosedata)


#-------------------------
#Parotid DOSEDIFF 

organ <- "Parotid"
metric = "Dose_Diff"

#-------------------------


#Vol_Tab: Count of how much data each method has for each patient
dose_tab <- Dosedata %>% 
  filter(Contour == organ) %>% 
  select(-metric) %>% 
  group_by( Patient, Method) %>%
  count() %>% 
  pivot_wider(names_from = Method, values_from = n)
dose_tab


#vol_df: isolate required data based on method, patient, metric)
dose_df <- Dosedata %>% 
  filter(Contour == organ ) %>% 
  rename(Metric = all_of(metric)) %>% 
  transmute(Method, Patient = factor(Patient), Metric) 


#run metric model based on glmmTMB
metric_mod <- glmmTMB(Metric ~ Method + (1|Patient), data = dose_df)
summary(metric_mod)

#check assumption with plots
plot(predict(metric_mod, re.form = NA), residuals(metric_mod))
qqnorm(residuals(metric_mod))

#model prediction with EMMIP
metric_predict <- emmip(metric_mod, ~ Method , CIs = TRUE, plotit = FALSE) 

ggplot(data = metric_predict, aes(Method, yvar) ) + 
  geom_point() + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0, size = 3, col = "lightgrey")+
  geom_point(aes(Method, Metric, color = Patient), data = dose_df) +
  ylab(metric)+
  theme_classic()

#use EMMANS and then CONTRAST
ems <- emmeans(metric_mod, ~ Method)
ems

cont_par_dd <- contrast(ems, "trt.vs.ctrlk", ref = 2) #ref = 2 because VEL_DMP is second in the list above

plot(cont_par_dd) + 
  theme_classic() + 
  geom_vline(xintercept = 0)

cont=cont_par_dd

summary(cont)
confint(cont)

write("Parotid DOSEDIFF",file="Contour_Analysis.txt",append=TRUE)
write.table(summary(cont),file="Contour_Analysis.txt",append=TRUE)
write.table(confint(cont),file="Contour_Analysis.txt",append=TRUE)



#-------------------------
#Brainstem DOSEDIFF 

organ <- "Brainstem"
metric = "Dose_Diff"

#-------------------------


#Vol_Tab: Count of how much data each method has for each patient
dose_tab <- Dosedata %>% 
  filter(Contour == organ) %>% 
  select(-metric) %>% 
  group_by( Patient, Method) %>%
  count() %>% 
  pivot_wider(names_from = Method, values_from = n)
dose_tab


#vol_df: isolate required data based on method, patient, metric)
dose_df <- Dosedata %>% 
  filter(Contour == organ ) %>% 
  rename(Metric = all_of(metric)) %>% 
  transmute(Method, Patient = factor(Patient), Metric) 


#run metric model based on glmmTMB
metric_mod <- glmmTMB(Metric ~ Method + (1|Patient), data = dose_df)
summary(metric_mod)

#check assumption with plots
plot(predict(metric_mod, re.form = NA), residuals(metric_mod))
qqnorm(residuals(metric_mod))

#model prediction with EMMIP
metric_predict <- emmip(metric_mod, ~ Method , CIs = TRUE, plotit = FALSE) 

ggplot(data = metric_predict, aes(Method, yvar) ) + 
  geom_point() + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0, size = 3, col = "lightgrey")+
  geom_point(aes(Method, Metric, color = Patient), data = dose_df) +
  ylab(metric)+
  theme_classic()

#use EMMANS and then CONTRAST
ems <- emmeans(metric_mod, ~ Method)
ems

cont_bs_dd <- contrast(ems, "trt.vs.ctrlk", ref = 2) #ref = 2 because VEL_DMP is second in the list above

plot(cont_bs_dd) + 
  theme_classic() + 
  geom_vline(xintercept = 0)

cont=cont_bs_dd

summary(cont)
confint(cont)

write("Brainstem DOSEDIFF",file="Contour_Analysis.txt",append=TRUE)
write.table(summary(cont),file="Contour_Analysis.txt",append=TRUE)
write.table(confint(cont),file="Contour_Analysis.txt",append=TRUE)

#-------------------------
#plot all contrast
#-------------------------

library(patchwork)

plot2 <- plot(cont_par_mda,  xlab="(ii) Parotid MDA (mm)", ylab="") + 
  theme_classic() + 
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = -1.21,linetype='dashed',col='black')+
  theme(text = element_text(size = 14))   

plot1 <- plot(cont_par_dsc,  xlab="(i) Parotid DSC", ylab="") + 
  theme_classic() + 
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = 0.17,linetype='dashed',col='black')+
  theme(text = element_text(size = 14))   

plot3 <- plot(cont_par_hd,  xlab="(iii) Parotid HD (mm)", ylab="") + 
  theme_classic() + 
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = -7.18,linetype='dashed',col='black')+
  theme(text = element_text(size = 14))   

plot4 <- plot(cont_par_vd,  xlab="(iv) Parotid Vol Diff (cc)", ylab="") + 
  theme_classic() + 
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = -0.32,linetype='dashed',col='black')+
  theme(text = element_text(size = 14))   

plot5 <- plot(cont_par_dd,  xlab="(v) Parotid Dose Diff (Gy)", ylab="") + 
  theme_classic() + 
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = 0.35,linetype='dashed',col='black')+
  theme(text = element_text(size = 14))   

plot7 <- plot(cont_bs_mda, xaxt = "n",  xlab="(vii) Brainstem MDA (mm)", ylab="") + 
  theme_classic() + 
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = -0.17,linetype='dashed',col='black')+
  theme(text = element_text(size = 14))   

plot6 <- plot(cont_bs_dsc,  xlab="(vi) Brainstem DSC", ylab="") + 
  theme_classic() + 
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = 0.15,linetype='dashed',col='black')+
  theme(text = element_text(size = 14))   

plot8 <- plot(cont_bs_hd,  xlab="(viii) Brainstem HD (mm)", ylab="") + 
  theme_classic() + 
  geom_vline(xintercept = 0)+
  theme_classic() + 
  geom_vline(xintercept = -7.5,linetype='dashed', col="black")+
  theme(text = element_text(size = 14))   

plot9 <- plot(cont_bs_vd,  xlab="(ix) Brainstem Vol Diff (cc)", ylab="") + 
  theme_classic() + 
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = -1.43,linetype='dashed',col='black')+
  theme(text = element_text(size = 14))   

plot10 <- plot(cont_bs_dd,  xlab="(x) Brainstem Dose Diff (Gy)", ylab="") + 
  theme_classic() + 
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = -1.37,linetype='dashed',col='black')+
  theme(text = element_text(size = 14))   


  
#See ?ggarrange

ggarrange(plot1, plot6, plot2, plot7, plot3,plot8, plot4, plot9, plot5, plot10,
                           nrow = 5,
                           ncol = 2) #nrow & ncol depend on how you want to 
#organize your plots

