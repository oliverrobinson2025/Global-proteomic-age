library(tidyverse)
library( haven )  # To load SAS sas7bdat file with variable labels too  
library( lubridate ) 
library( Epi )
library( labelled ) 
library(survival)
library(htmlwidgets)
library(DT)
library(forcats)
library("glmnet")
library(table1)
library(forestplot)

##prepare covariate db
death_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP <- readRDS("/data/Epic/subprojects/Somalogic/work/Vivian/Phase2/Data/ForAnalyses/death_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP.rds")#import dataset prepared by data manager for mortality analysis
db<-death_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP[,c(1:220,7817:7827)]#take only covariate database
names(db)
remove(death_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP)

db$l_school[db$l_school=="Not specified"]<-NA
db$l_school<-as_factor(db$l_school)
db$smoke_stat[db$smoke_stat==4]<-NA
db$smoke_stat<-as.factor(db$smoke_stat)
db$hli_dietscore_c<-as.factor(db$hli_dietscore_c)


setwd("/data/Epic/subprojects/Somalogic/work/Oliver")
load("output/imputedcovariates.Rdata") ##(prepared in script 5)
imputedcovarites$idepic<-imputedcovarites$Idepic
imputedcovarites$Idepic<-NULL
db<-merge(db, imputedcovarites, all.x=T, by= "idepic")



###import clock data##
setwd("/data/Epic/subprojects/Somalogic/work/Oliver")
data_clocks=read.csv("output/data_clocks_new.csv")

clock_names=c("Oh","Organismal","Brain","Adipose", "Artery", "Immune","Heart","Intestine","Kidney","Liver","Lung","Muscle","Pancreas","Lehallier","Tanaka", "Wang", "Sathyan", "Global")


######


## Cox death regressions



######Apply prentice weighting for main analysis###########################################

db1<-db[db$cvd_t2d_coh==0 & db$death_status==1,] ### Select cases that are not in subcohort
db1$age_recr <- db1$age_exit_death - 1e-4 # apply prentice weighting to follow up time of cases that are not in subcohort
db2<-db[db$cvd_t2d_coh==1 & db$death_status==1,]### select cases in subcohort
db3<-db[db$cvd_t2d_coh==1 & db$death_status==0,]###select non-cases in subcohort
db.spec=rbind(db1,db2,db3) ###combined for prentice weighted new db for analysis
remove(db1,db2,db3)
############

###Run below instead for sensitivity analysis with exclusion of first 2 or 5 years of events  ###
#db$death_fut<-db$age_exit_death-db$age_recr
#db$pre_5_year<-ifelse(db$death_fut <= 5,1,0)
#db$pre_2_year<-ifelse(db$death_fut <= 2,1,0)

#dba<-db[db$death_status==1,]
##dba<-dba[dba$pre_5_year==0,] ###choose to exclude 2 or 5 years of events
#dba<-dba[dba$pre_2_year==0,]
#dbb<-db[db$death_status==0,]
#db.c=rbind(dba,dbb)
#remove(dba,dbb)

#db1<-db.c[db.c$cvd_t2d_coh==0 & db.c$death_status==1,]
#db1$age_recr <- db1$age_exit_death - 1e-4#prentice weighting applied to db for sensitivity analysis
#db2<-db.c[db.c$cvd_t2d_coh==1 & db.c$death_status==1,]
#db3<-db.c[db.c$cvd_t2d_coh==1 & db.c$death_status==0,]
#db.spec=rbind(db1,db2,db3)
#remove(db1,db2,db3, db.c)

################Run models####

Res_deaths <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                     SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
Res_deaths_60 <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                        SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())



for (clock_name in clock_names){
  
  ### merge dfs
  data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
  merged=merge(x=db.spec,y=data_plot,by="idepic")
  
  #####exclude ever smokers for sensitivity analysis
  #merged<-merged[merged$smoke_stat=="Never",]
  #############
  
  ### adjusted model combinations for age, sex and center/country
  
    mod_70plus <- coxph(Surv(age_recr, age_exit_death, death_status) ~  AgeGap_zscored + strata(age_cat_5, center, sex.x) + cluster(idepic),
                      data= merged)
  mod_60 <- coxph(Surv(age_recr, age_exit_death, death_status) ~  AgeGap_zscored + strata(age_cat_5, center, sex.x) + cluster(idepic),
                  data= merged[which(merged$age_recr >60),]) ### model only looking at over 60 years old
  
  
  temptib_70plus <- tibble(Clock= clock_name, Asso = summary(mod_70plus)$coef[1, 'exp(coef)'],
                           beta = summary(mod_70plus)$coef[1, 'coef'],
                           pval = summary(mod_70plus)$coef[1, 'Pr(>|z|)'],
                           SE= summary(mod_70plus)$coef[1, 'se(coef)'] ,
                           LCL =exp(confint(mod_70plus)[1]),
                           UCL =exp(confint(mod_70plus)[2]),
                           n =mod_70plus$n ,nevent=mod_70plus$nevent)
  temptib_60 <- tibble(Clock= clock_name, Asso = summary(mod_60)$coef[1, 'exp(coef)'],
                       beta = summary(mod_60)$coef[1, 'coef'],
                       pval = summary(mod_60)$coef[1, 'Pr(>|z|)'],
                       SE= summary(mod_60)$coef[1, 'se(coef)'] ,
                       LCL =exp(confint(mod_60)[1]),
                       UCL =exp(confint(mod_60)[2]),
                       n =mod_60$n ,nevent=mod_60$nevent)
  
  Res_deaths <- bind_rows(Res_deaths, temptib_70plus)
  Res_deaths_60 <- bind_rows(Res_deaths_60, temptib_60)
  
}

Res_deaths <- Res_deaths %>% mutate(FDR = p.adjust(pval, method = "BH"))
Res_deaths_60 <- Res_deaths_60 %>% mutate(FDR = p.adjust(pval, method = "BH"))


print(Res_deaths[order(Res_deaths$pval),])



###risk factor adjustments##############
covs<-c( "l_school", "bmi_c" ,  "alc_re",  "pa_mets" , "hli_dietscore_c" ,   "smoke_stat")
covs.imp<-c( "L_School.imp", "bmi_c" ,  "Alc_Re.imp",  "Pa_Mets.imp" ,  "Smoke_Stat.imp","Hli_Dietscore_C.imp" )

Res_deaths_adj <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                     SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
Res_deaths_imp <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                         SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())

for (clock_name in clock_names){
  
  ### merge dfs
  data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
  merged=merge(x=db.spec,y=data_plot,by="idepic")
 
  #####exclude ever smokers
  #merged<-merged[merged$smoke_stat=="Never",]
  #############
  
  ### adjusted model combinations for age, sex and center/country
  
  form<-as.formula(paste0("Surv(age_recr, age_exit_death, death_status) ~ AgeGap_zscored +", paste(covs, collapse = "+"),"+ strata(age_cat_5, center, sex.x) + cluster(idepic)"))
  mod_70plus <- coxph(form, data=merged)
  
    temptib_70plus <- tibble(Clock= clock_name, Asso = summary(mod_70plus)$coef[1, 'exp(coef)'],
                           beta = summary(mod_70plus)$coef[1, 'coef'],
                           pval = summary(mod_70plus)$coef[1, 'Pr(>|z|)'],
                           SE= summary(mod_70plus)$coef[1, 'se(coef)'] ,
                           LCL =exp(confint(mod_70plus)[1]),
                           UCL =exp(confint(mod_70plus)[1,2]),
                           n =mod_70plus$n ,nevent=mod_70plus$nevent)
 
  
  Res_deaths_adj <- bind_rows(Res_deaths_adj, temptib_70plus)
  
  form<-as.formula(paste0("Surv(age_recr, age_exit_death, death_status) ~ AgeGap_zscored +", paste(covs.imp, collapse = "+"),"+ strata(age_cat_5, center, sex.x) + cluster(idepic)"))
  mod_70plus <- coxph(form, data=merged)
  
  
  temptib_70plus <- tibble(Clock= clock_name, Asso = summary(mod_70plus)$coef[1, 'exp(coef)'],
                           beta = summary(mod_70plus)$coef[1, 'coef'],
                           pval = summary(mod_70plus)$coef[1, 'Pr(>|z|)'],
                           SE= summary(mod_70plus)$coef[1, 'se(coef)'] ,
                           LCL =exp(confint(mod_70plus)[1]),
                           UCL =exp(confint(mod_70plus)[1,2]),
                           n =mod_70plus$n ,nevent=mod_70plus$nevent)
  
  
  Res_deaths_imp <- bind_rows(Res_deaths_imp, temptib_70plus)
 
  
}

Res_deaths_adj <- Res_deaths_adj %>% mutate(FDR = p.adjust(pval, method = "BH"))
Res_deaths_imp <- Res_deaths_imp %>% mutate(FDR = p.adjust(pval, method = "BH"))


#####Save results####
save(Res_deaths,Res_deaths_60,Res_deaths_adj,Res_deaths_imp, file="output/MortalityResults_QC_v2_con2.Rdata" )
#save(Res_deaths,Res_deaths_60,Res_deaths_adj,Res_deaths_imp, file="output/MortalityResults_QC_v2_con2_pre5excl.Rdata" )
#save(Res_deaths,Res_deaths_60,Res_deaths_adj,Res_deaths_imp, file="output/MortalityResults_QC_v2_con2_neversmokers.Rdata" )
#save(Res_deaths,Res_deaths_60,Res_deaths_adj,Res_deaths_imp, file="output/MortalityResults_QC_v2_con2_pre2excl.Rdata" )


#####Compare prediction of different mortality models################################################################

db.spec$death_fut<-db.spec$age_exit_death-db.spec$age_recr ###use followup time as time scale

data_plot=data_clocks[which(data_clocks$Organ=="Global"),]
merged=merge(x=db.spec,y=data_plot,by="idepic")

mod_base_mortality <- coxph(Surv(death_fut,death_status) ~ age +sex.x +strata(center) ,data= merged)
summary(mod_base_mortality)

mod_global_mortality <- coxph(Surv(death_fut,death_status) ~ AgeGap_zscored+ age +sex.x +strata(center) ,data= merged)
summary(mod_global_mortality)

mod_RF_mortality <- coxph(Surv(death_fut,death_status) ~ L_School.imp+bmi_c+Alc_Re.imp+Pa_Mets.imp+Smoke_Stat.imp+Hli_Dietscore_C.imp + age +sex.x +strata(center) ,data= merged)
summary(mod_RF_mortality)

mod_global_RF_mortality <- coxph(Surv(death_fut,death_status) ~ AgeGap_zscored+L_School.imp+bmi_c+Alc_Re.imp+Pa_Mets.imp+Smoke_Stat.imp+Hli_Dietscore_C.imp + age +sex.x +strata(center) ,data= merged)
summary(mod_global_RF_mortality)

anova(mod_RF_mortality, mod_global_RF_mortality)##log-liklihood ratio test comparing nested models

CIresults<-rbind(c(mod_consensus_mortality[["concordance"]][6],
                   mod_consensus_mortality[["concordance"]][6]-1.96*mod_consensus_mortality[["concordance"]][7], 
                   mod_consensus_mortality[["concordance"]][6]+1.96*mod_consensus_mortality[["concordance"]][7]),
                 c(mod_RF_mortality[["concordance"]][6],
                   mod_RF_mortality[["concordance"]][6]-1.96*mod_RF_mortality[["concordance"]][7], 
                   mod_RF_mortality[["concordance"]][6]+1.96*mod_RF_mortality[["concordance"]][7]),
                 c(mod_consensusRF_mortality[["concordance"]][6],
                   mod_consensusRF_mortality[["concordance"]][6]-1.96*mod_consensusRF_mortality[["concordance"]][7], 
                   mod_consensusRF_mortality[["concordance"]][6]+1.96*mod_consensusRF_mortality[["concordance"]][7])) ##table of C-index results
############select organs into model

###prepare wide clock data###
organclocks<- c("Conventional","Organismal","Brain","Adipose", "Artery", "Immune","Heart","Intestine","Kidney","Liver","Lung","Muscle","Pancreas")
wideclockdb<-data_clocks[data_clocks$Organ==organclocks[1],]
wideclockdb<-wideclockdb[c("idepic" , "AgeGap_zscored")]

for(i in 2:length(organclocks)){
  clockdb<-data_clocks[data_clocks$Organ==organclocks[i],"AgeGap_zscored"]
  wideclockdb<-cbind(wideclockdb,clockdb)
}
names(wideclockdb)<-c("idepic", paste0(organclocks, "_AgeGap_zscored"))
remove(clockdb)

cororgan<-cor(wideclockdb[,-1])
pheatmap(cororgan, display_numbers=T) ###visualise correlation with organ age gaps

#### run lasso to select organ clocks##

merged=merge(x=db.spec,y=wideclockdb,by="idepic")

merged$age_cat_5<-droplevels(merged$age_cat_5)
merged$age_cat_5_int<-as.integer(merged$age_cat_5)
y=Surv( merged$age_recr, merged$age_exit_death, merged$death_status)
y2 <- stratifySurv(y, merged$age_cat_5_int) ##statify by age group

merged$female= ifelse(merged$sex=="Female",1,0) ##recode sex and centre variables
str(merged$center)
merged$center<-droplevels(merged$center)
centers<-as.vector(unique(merged$center))
centerdb<-matrix(NA, nrow(merged), length(centers))
for(i in 1:length(centers)){ centerdb[,i] <-ifelse(merged$center==centers[i],1,0)}
centerdb<-as.data.frame(centerdb)
names(centerdb)<-centers

x=as.matrix(cbind(merged[c(names(wideclockdb)[-1:-2], "female")],centerdb)) ###data in glmnet format

cv.fit <- cv.glmnet(x, y2, family = "cox",alpha=1,nfolds = 5)##tune lasso model
selcexted<-as.matrix(coef.glmnet(cv.fit, s= cv.fit$lambda.1se))## take optimum lasso model
selcexted<-selcexted[selcexted[,1]!=0,]###index of organs selected into model

###fit organ models and merge with other C-index results

merged=merge(x=db.spec,y=wideclockdb,by="idepic")
covs.imp<-c( "L_School.imp", "bmi_c" ,  "Alc_Re.imp",  "Pa_Mets.imp" ,  "Smoke_Stat.imp","Hli_Dietscore_C.imp" )

agegaps<-names(wideclockdb)[2]
form<-as.formula(paste0("Surv(death_fut,death_status) ~",paste(agegaps, collapse = "+") , "+ age +sex +", paste(covs.imp, collapse = "+"),"+ strata( center)"))
mod_all_Oh <- coxph(form, data=merged)
summary(mod_all_Oh)
res<-c(mod_all_Oh[["concordance"]][6],
       mod_all_Oh[["concordance"]][6]-1.96*mod_all_Oh[["concordance"]][7], 
       mod_all_Oh[["concordance"]][6]+1.96*mod_all_Oh[["concordance"]][7])
CIresults<-rbind(CIresults,res)

agegaps<-names(wideclockdb)[names(wideclockdb) %in% names(selcexted)]
form<-as.formula(paste0("Surv(death_fut,death_status) ~",paste(agegaps, collapse = "+") , "+ age +sex +", paste(covs.imp, collapse = "+"),"+ strata( center) "))
mod_sel_organ <- coxph(form, data=merged)
summary(mod_sel_organ)
res<-c(mod_sel_organ[["concordance"]][6],
       mod_sel_organ[["concordance"]][6]-1.96*mod_sel_organ[["concordance"]][7], 
       mod_sel_organ[["concordance"]][6]+1.96*mod_sel_organ[["concordance"]][7])
CIresults<-rbind(CIresults,res)

agegaps<-names(wideclockdb)[-1:-2]
form<-as.formula(paste0("Surv(death_fut,death_status) ~",paste(agegaps, collapse = "+") , "+age +sex +", paste(covs.imp, collapse = "+"),"+ strata( center) "))
mod_all_organ <- coxph(form, data=merged)
summary(mod_all_organ)
res<-c(mod_all_organ[["concordance"]][6],
       mod_all_organ[["concordance"]][6]-1.96*mod_all_organ[["concordance"]][7], 
       mod_all_organ[["concordance"]][6]+1.96*mod_all_organ[["concordance"]][7])
CIresults<-rbind(CIresults,res)

agegaps<-names(wideclockdb)[-1:-2]
form<-as.formula(paste0("Surv(age_recr, age_exit_death, death_status) ~",paste(agegaps, collapse = "+") , "+", paste(covs.imp, collapse = "+"),"+ strata(age_cat_5, center, sex) + cluster(idepic)"))
mod_all_organ <- coxph(form, data=merged)
res<-c(mod_all_organ[["concordance"]][6],
       mod_all_organ[["concordance"]][6]-1.96*mod_all_organ[["concordance"]][7], 
       mod_all_organ[["concordance"]][6]+1.96*mod_all_organ[["concordance"]][7])
CIresults<-rbind(CIresults,res)


###Cindex forest plot
Cindextext= paste0(round(CIresults[,1],2), " (",round(CIresults[,2],2), ", ", round(CIresults[,3],2), ")" )
labels<-c("Base model (age, sex, stratified by center)","Global age gap", "Risk factors only", "Global age gap and risk factors","Oh age gap and risk factors",   "Lasso selected organ clocks and risk factors", "All organ clocks and risk factors")

pdf(file=paste0("/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/figures/Cindex plot.pdf"), height=6, width =12)
forestplot(cbind(labels, c(0.5,Cindextext)), 
           mean = c(0.5,CIresults[,1]),
           lower = c(0.5,CIresults[,2]),
           upper = c(0.5,CIresults[,3]),
           xlab="Concordance Index",  zero=0.5,boxsize = 0.25, fn.ci_norm = fpDrawCircleCI,
           txt_gp =fpTxtGp(ticks = gpar(cex=1.5), xlab  = gpar(cex=1.5), cex=1.5))
dev.off()


####coefficients plot for lasso selcted organ model
agegaps<-names(wideclockdb)[names(wideclockdb) %in% names(selcexted)]
merged$bmi_sc<-scale(merged$bmi_c)
merged$alc_sc<-scale(merged$Alc_Re.imp)
merged$pa_sc<-scale(merged$Pa_Mets.imp)
form<-as.formula(paste0("Surv(age_recr, age_exit_death, death_status) ~",paste(agegaps, collapse = "+") , "+ L_School.imp+bmi_sc+alc_sc+pa_sc+Smoke_Stat.imp+Hli_Dietscore_C.imp+ strata(age_cat_5, center, sex) + cluster(idepic)"))
mod_selected_organ <- coxph(form, data=merged)

##model summary plots#
modeltab<-cbind(summary(mod_selected_organ)$coef[,'exp(coef)'],exp(confint(mod_selected_organ)), summary(mod_selected_organ)$coef[,'Pr(>|z|)'])
modeltab<-modeltab[-12,]
nullvector=c(1,1,1,1)
modeltab<-as.matrix(rbind(modeltab[1:7,],nullvector,modeltab[8:14,],nullvector,modeltab[15:16,],nullvector, modeltab[17:20,]))
plotnames<-c(rownames(modeltab)[1:7],"No Schooling completed (ref)", "Primary school completed","Technical/professional school", "Secondary school", "Longer education (incl. University deg.)", 
                                "BMI (z_scored)", "Alcohol consumption (z_scored)", "Physical activity (z_scored)","Never smoker (ref)", "Former smoker", "Current smoker",
             "Healthy diet score Q1 (lowest, ref)", "Healthy diet score Q2", "Healthy diet score Q3","Healthy diet score Q4","Healthy diet score Q5 (highest)")

pdf(file=paste0("/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/figures/Lasso_coeffiecent_plot.pdf"), height=12, width =6)
forestplot(plotnames, mean = modeltab[,1],lower = modeltab[,2],upper = modeltab[,3],title="Lasso selected organ clock model of all-cause mortality" ,
           xlab="Hazard ratio",  zero=1, boxsize = 0.25,
           txt_gp =fpTxtGp(ticks = gpar(cex=1), xlab  = gpar(cex=1)))
dev.off()

###add comparison of C indexes scorss all the conventional clcoks####

data_plot=data_clocks[which(data_clocks$Organ=="Oh"),]
merged=merge(x=db.spec,y=data_plot,by="idepic")
mod_adj_oh <- coxph(Surv(death_fut,death_status) ~ AgeGap_zscored+ L_School.imp+bmi_c+Alc_Re.imp+Pa_Mets.imp+Smoke_Stat.imp+Hli_Dietscore_C.imp + age +sex.x +strata(center) ,data= merged)

data_plot=data_clocks[which(data_clocks$Organ=="Tanaka"),]
merged=merge(x=db.spec,y=data_plot,by="idepic")
mod_adj_tanaka <- coxph(Surv(death_fut,death_status) ~ AgeGap_zscored+ L_School.imp+bmi_c+Alc_Re.imp+Pa_Mets.imp+Smoke_Stat.imp+Hli_Dietscore_C.imp + age +sex.x +strata(center) ,data= merged)

data_plot=data_clocks[which(data_clocks$Organ=="Lehallier"),]
merged=merge(x=db.spec,y=data_plot,by="idepic")
mod_adj_lehallier <- coxph(Surv(death_fut,death_status) ~ AgeGap_zscored+ L_School.imp+bmi_c+Alc_Re.imp+Pa_Mets.imp+Smoke_Stat.imp+Hli_Dietscore_C.imp + age +sex.x +strata(center) ,data= merged)

data_plot=data_clocks[which(data_clocks$Organ=="Sathyan"),]
merged=merge(x=db.spec,y=data_plot,by="idepic")
mod_adj_sathyan <- coxph(Surv(death_fut,death_status) ~ AgeGap_zscored+ L_School.imp+bmi_c+Alc_Re.imp+Pa_Mets.imp+Smoke_Stat.imp+Hli_Dietscore_C.imp + age +sex.x +strata(center) ,data= merged)

data_plot=data_clocks[which(data_clocks$Organ=="Wang"),]
merged=merge(x=db.spec,y=data_plot,by="idepic")
mod_adj_wang <- coxph(Surv(death_fut,death_status) ~ AgeGap_zscored+ L_School.imp+bmi_c+Alc_Re.imp+Pa_Mets.imp+Smoke_Stat.imp+Hli_Dietscore_C.imp + age +sex.x +strata(center) ,data= merged)

CIresults_coventional<-rbind(c(mod_RF_mortality[["concordance"]][6],
                   mod_RF_mortality[["concordance"]][6]-1.96*mod_RF_mortality[["concordance"]][7], 
                   mod_RF_mortality[["concordance"]][6]+1.96*mod_RF_mortality[["concordance"]][7]),
                 c(mod_adj_tanaka[["concordance"]][6],
                   mod_adj_tanaka[["concordance"]][6]-1.96*mod_adj_tanaka[["concordance"]][7], 
                   mod_adj_tanaka[["concordance"]][6]+1.96*mod_adj_tanaka[["concordance"]][7]),
                 c(mod_adj_lehallier[["concordance"]][6],
                   mod_adj_lehallier[["concordance"]][6]-1.96*mod_adj_lehallier[["concordance"]][7], 
                   mod_adj_lehallier[["concordance"]][6]+1.96*mod_adj_lehallier[["concordance"]][7]),
                 c(mod_adj_sathyan[["concordance"]][6],
                   mod_adj_sathyan[["concordance"]][6]-1.96*mod_adj_sathyan[["concordance"]][7], 
                   mod_adj_sathyan[["concordance"]][6]+1.96*mod_adj_sathyan[["concordance"]][7]),
                 c(mod_adj_oh[["concordance"]][6],
                   mod_adj_oh[["concordance"]][6]-1.96*mod_adj_oh[["concordance"]][7], 
                   mod_adj_oh[["concordance"]][6]+1.96*mod_adj_oh[["concordance"]][7]),
                 c(mod_adj_wang[["concordance"]][6],
                   mod_adj_wang[["concordance"]][6]-1.96*mod_adj_wang[["concordance"]][7], 
                   mod_adj_wang[["concordance"]][6]+1.96*mod_adj_wang[["concordance"]][7]),
                 c(mod_global_RF_mortality[["concordance"]][6],
                   mod_global_RF_mortality[["concordance"]][6]-1.96*mod_global_RF_mortality[["concordance"]][7], 
                   mod_global_RF_mortality[["concordance"]][6]+1.96*mod_global_RF_mortality[["concordance"]][7]))

##Cindex forest plot
Cindextext_conventional= paste0(round(CIresults_coventional[,1],3), " (",round(CIresults_coventional[,2],3), ", ", round(CIresults_coventional[,3],3), ")" )
labels_con<-c( "Risk factors only", "Tanaka","Lehallier",  "Sathyan", "Oh","Wang","Global")

pdf(file=paste0("/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/figures/Cindex plot coventional.pdf"), height=6, width =12)
forestplot(cbind(labels_con, c(Cindextext_conventional)), 
           mean = c(CIresults_coventional[,1]),
           lower = c(CIresults_coventional[,2]),
           upper = c(CIresults_coventional[,3]),
           xlab="Concordance Index",  zero=0.7, xticks = c(0.7,0.71,0.72,0.73,0.74,0.75,0.76),
           txt_gp =fpTxtGp(ticks = gpar(cex=1.5), xlab  = gpar(cex=1.5), cex=1.5))
dev.off()

######Output descriptive table of mortality sample#####
names(db)
table(db$cntr_f)
db$cntr_f<-droplevels(db$cntr_f)
table(db$country)
db$country<-droplevels(db$country)
db$death_fut<-db$age_exit_death-db$age_blood
summary(db$death_fut)
db$l_school[db$l_school=="Not specified"]<-NA
db$l_school<-droplevels(db$l_school)
db$smoke_stat<-droplevels(db$smoke_stat)

str(db$death_status)

table(db$t2d_status)

####
tabledb<-db[db$cvd_t2d_coh==1 ,]
tabledb<-tabledb[which(tabledb$idepic %in% data_clocks$idepic),]
tabledb<-to_factor(tabledb)
tabledb$death_status <- factor(tabledb$death_status, levels=c(0,1),labels=c("No", "Yes"))
label(tabledb$death_fut) <- "Follow up time deaths (yrs)"


write.csv(
  x = table1(~age_blood+
               sex+
               cntr_f+
               l_school+
               bmi_c+
               alc_re+
               smoke_stat+
               pa_mets+
               pa_index+
               hli_dietscore+
               hli_score+
               menopause+
               t2d_prev+
               cvd_prev+
               death_status+
               death_fut+
               cvd_1_status+
               cvd_fut+
               t2d_status+
               t2d_fut+
               cncr_mal_anyc+
               cancer_fut+
               neuro_status+
               fasting_c
             |country, data=tabledb, render.continuous="Mean (SD)"),
  file = "/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/descriptive_tables/mortality_sample.csv",
  row.names = FALSE  )


save(modeltab,plotnames,Cindextext, labels, file="output/for_export/summary_results/summary_results_for_mortality_predictivemodels.Rdata" )

#save.image(file="ClockMortality.Rdata")


