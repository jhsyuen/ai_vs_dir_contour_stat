

#-------------------------
# R Code for plot of boxdata from  multiple methods of contour generation (AI/DIR)
#
# Input:
#  Excel data with self-reported grading and clinical relevance 
#  Excel data with MSC and MDA
#  Excel data with dose difference
#  Excel data with volume difference
#
#
# Output files:
#
#  Plots of data: MDA, DSC, Dose difference, Volume difference, Training, Clinical significance
#
# R code submitted as part of manuscript submitted
# "AI compared to deformable image registration propagation techniques for contouring:" 
# "feasibility and accuracy for head and neck replanning" 
#
# 
#-------------------------
library(tidyverse)
library(emmeans)
library(hrbrthemes)
library(viridis)
library(likert)
library(forcats)
library(readxl)

setwd("/cloud/project/R /Rdata")

ContourSigData <- read_excel("Rdata_ContourSigLikert.xlsx")

standardize_levelsContourSigData <- function(item_data) {
  levels_orderContourSigData <- c("Never Relevant", "Rarely Relevant", "Sometimes Relevant", "Very Often Relevant", "Always Relevant")
  factor(item_data, levels = levels_orderContourSigData)
}

# Apply the standardize_levels function to each column in data2
ContourSigData_standardized <- lapply(ContourSigData, standardize_levelsContourSigData)

class(ContourSigData_standardized)
ContourSigData_df <- as.data.frame(ContourSigData_standardized)

# Create the likert plot
ContourSigData_plot <- likert(ContourSigData_df)
plot(ContourSigData_plot)


  

ContourGradeData <- read_excel("Rdata_ContourGrade_Likert.xlsx")

standardize_levelsContourGradeData <- function(item_data) {
  levels_order <- c("Grade 5", "Grade 4", "Grade 3", "Grade 2", "Grade 1")
  factor(item_data, levels = levels_order)
}

# Apply the standardize_levels function to each column in data2
ContourGradeData_standardized <- lapply(ContourGradeData, standardize_levelsContourGradeData)

class(ContourGradeData_standardized)
ContourGradeData_df <- as.data.frame(ContourGradeData_standardized)

# Create the likert plot
ContourGradeData_plot <- likert(ContourGradeData_df)
plot(ContourGradeData_plot)  #plot ordered according to grading data

plot(ContourGradeData_plot, group.order=c("Orbit_Man","Orbit_DLC","Mandible_Man","Mandible_DLC","Cochlea_Man","Cochlea_DLC","Brainstem_Man","Brainstem_DLC","Parotid_Man","Parotid_DLC","BrachP_Man","BrachP_DLC","Oralcavity_Man","Oralcavity_DLC")) #plot ordered according to OAR manual then DLC for each OAR

      


DSCMDA <- read_excel("Rdata_DSCMDAvF.xlsx")

summary(DSCMDA)

DSCMA_BrachP <- filter(DSCMDA,Contour=="BrachP")
DSCMA_BrachP %>% 
  group_by(Method) %>% 
  summarise(mean(Dice[!Dice %in% res$out]),
            mean(MDA[!MDA %in% res$out])) #calculate mean without outlier


DSCMA_Brainstem <- filter(DSCMDA,Contour=="Brainstem")
DSCMA_Brainstem %>% 
  group_by(Method) %>% 
  summarise(mean(Dice[!Dice %in% res$out]),
            mean(HD[!HD %in% res$out]),
            mean(MDA[!MDA %in% res$out])) #calculate mean without outlier


DSCMA_Par <- filter(DSCMDA,Contour=="Parotid")
DSCMA_Par %>% 
  group_by(Method) %>% 
  summarise(mean(Dice[!Dice %in% res$out]),
            mean(HD[!HD %in% res$out]),
            mean(MDA[!MDA %in% res$out])) #calculate mean without outlier

DSCMA_coch <- filter(DSCMDA,Contour=="Cochlea")
DSCMA_coch %>% 
  group_by(Method) %>% 
  summarise(mean(Dice[!Dice %in% res$out]),
            mean(MDA[!MDA %in% res$out])) #calculate mean without outlier

DSCMA_oralc <- filter(DSCMDA,Contour=="Oralcavity")
DSCMA_oralc %>% 
  group_by(Method) %>% 
  summarise(mean(Dice[!Dice %in% res$out]),
           mean(MDA[!MDA %in% res$out])) #calculate mean without outlier


str(DSCMDA)

ggplot(DSCMDA, aes(x=Contour,y=MDA, fill=Method))+
  geom_boxplot()+
  geom_hline(yintercept = 0,linetype='dashed',col='blue')+
  geom_hline(yintercept = 2,linetype='dashed',col='green')+
  #geom_hline(yintercept = 4,linetype='dashed',col='red')+
  coord_cartesian(ylim = c(0, 16)) +  # Draw plot
  labs(title="Mean distance to agreement (MDA)", x="Contour",y="MDA (mm)") + 
  theme_classic() +
  theme(text = element_text(size = 20))    +
  geom_jitter(width=0.25, alpha=0.15)

  ggplot(DSCMDA, aes(x=Contour,y=HD, fill=Method))+
    geom_boxplot()+
    geom_hline(yintercept = 0,linetype='dashed',col='blue')+
    geom_hline(yintercept = 25,linetype='dashed',col='green')+
    #geom_hline(yintercept = 4,linetype='dashed',col='red')+
    coord_cartesian(ylim = c(0, 100)) +  # Draw plot
    labs(title="Hausdorff distance (HD)", x="Contour",y="HD (mm)") + 
    theme_classic() +
    theme(text = element_text(size = 20))    +
    geom_jitter(width=0.25, alpha=0.15)
  
ggplot(DSCMDA, aes(x=Contour,y=Dice, fill=Method))+
  geom_boxplot()+
  geom_hline(yintercept = 1.0,linetype='dashed',col='blue')+
  geom_hline(yintercept = 0.7,linetype='dashed',col='green')+
  #geom_hline(yintercept = 0.0,linetype='dashed',col='red')+
  labs(title="Dice similarity coefficient (DSC)", x="Contour",y="DSC") + 
  theme_classic()+
  theme(text = element_text(size = 20))   +
  geom_jitter(width=0.25, alpha=0.15)


Voldata <- read_excel("Rdata_VolumeF_ratio.xlsx")

summary(Voldata)

Voldata %>% 
  group_by(Method) %>% 
  summarise(mean_VD=mean(Volume_Diff), sd_VD=sd(Volume_Diff))

Voldata_oralc <- filter(Voldata,Contour=="OralCavity")
Voldata_oralc %>% 
  group_by(Method) %>% 
  summarise(mean(Volume_Diff[!Volume_Diff %in% res$out])) #calculate mean without outlier

Voldata_mand <- filter(Voldata,Contour=="Mandible")
Voldata_mand %>% 
  group_by(Method) %>% 
  summarise(mean(Volume_Diff[!Volume_Diff %in% res$out])) #calculate mean without outlier


Voldata_bs <- filter(Voldata,Contour=="Brainstem")
Voldata_bs %>% 
  group_by(Method) %>% 
  summarise(mean(Volume_Diff[!Volume_Diff %in% res$out])) #calculate mean without outlier


Voldata_par <- filter(Voldata,Contour=="Parotid")
Voldata_par %>% 
  group_by(Method) %>% 
  summarise(mean(Volume_Diff[!Volume_Diff %in% res$out])) #calculate mean without outlier



Voldata_Veldmp <- filter(Voldata,Method=="2, VEL_DMP")
Voldata_Veldmp %>% 
  group_by(Contour) %>% 
  summarise(mean(Volume_Diff))

Voldata_AI <- filter(Voldata,Method=="6, DLC")
Voldata_AI %>% 
  group_by(Contour) %>% 
  summarise(mean(Volume_Diff))


str(Voldata)

ggplot(Voldata, aes(x=Contour,y=Volume_Diff, fill=Method))+
  geom_boxplot()+
  geom_hline(yintercept = 1.0,linetype='dashed',col='blue')+
#  geom_hline(yintercept = 0.9,linetype='dashed',col='blue')+
#  geom_hline(yintercept = 1.1,linetype='dashed',col='blue')+
  geom_hline(yintercept = 0.8,linetype='dashed',col='green')+
  geom_hline(yintercept = 1.2,linetype='dashed',col='green')+
  coord_cartesian(ylim = c(0, 3)) +  # Draw plot
  labs(title="Volume ratio from manual contours", x="Contour",y="Volume ratio") + 
  theme_classic()+
  theme(text = element_text(size = 20))   +
  geom_jitter(width=0.25, alpha=0.15)


Dosedata <- read_excel("Rdata_DoseF_doseratio.xlsx")

summary(Dosedata)

Dosedata %>% 
  group_by(Method) %>% 
  summarise(mean_DD=mean(Dose_Diff), sd_DD=sd(Dose_Diff))

str(Dosedata)

Dosdata_Veldmp <- filter(Dosedata,Method=="2, VEL_DMP")
Dosdata_Veldmp %>% 
  group_by(Contour) %>% 
  summarise(mean(Dose_Diff[!Dose_Diff %in% res$out])) #calculate mean without outlier

Dosdata_AI <- filter(Dosedata,Method=="6, DLC")
Dosdata_AI %>% 
  group_by(Contour) %>% 
  summarise(mean(Dose_Diff[!Dose_Diff %in% res$out])) #calculate mean without outlier

Dosdata_par <- filter(Dosedata,Contour=="Parotid")
Dosdata_par %>% 
  group_by(Method) %>% 
  summarise(mean(Dose_Diff[!Dose_Diff %in% res$out])) #calculate mean without outlier

Dosdata_oc <- filter(Dosedata,Contour=="OralCavity")
Dosdata_oc %>% 
  group_by(Method) %>% 
  summarise(mean(Dose_Diff[!Dose_Diff %in% res$out])) #calculate mean without outlier

Dosdata_bs <- filter(Dosedata,Contour=="Brainstem")
Dosdata_bs %>% 
  group_by(Method) %>% 
  summarise(mean(Dose_Diff[!Dose_Diff %in% res$out])) #calculate mean without outlier


ggplot(Dosedata, aes(x=Contour,y=Dose_Diff, fill=Method))+
  geom_boxplot()+
  geom_hline(yintercept = 1,linetype='dashed',col='blue')+
  geom_hline(yintercept = 0.9,linetype='dashed',col='green')+
  geom_hline(yintercept = 1.1,linetype='dashed',col='green')+
  coord_cartesian(ylim = c(0, 2)) +  # Draw plot
  labs(title="Dose ratio from manual contours", x="Contour",y="Dose ratio") + 
  theme_classic()+
  theme(text = element_text(size = 20))   +
  geom_jitter(width=0.25, alpha=0.15)



