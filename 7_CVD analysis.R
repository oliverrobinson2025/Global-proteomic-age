library(tidyverse)
# devtools::install_github("SomaLogic/SomaDataIO@*release")
#library(SomaDataIO)
library( haven )  # To load SAS sas7bdat file with variable labels too  
library( lubridate ) 
library( Epi )
library( labelled ) 
library(survival)
library(htmlwidgets)
library(DT)
library(forcats)
library(haven)
library(gtsummary) 
library(labelled)

cvd_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP <- readRDS("/data/Epic/subprojects/Somalogic/work/Vivian/Phase2/Data/ForAnalyses/cvd_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP.rds")#import dataset prepared for CVD analysis
db<-cvd_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP[,c(1:220,7817:7827)]#take only covariates
remove(cvd_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP)

db$hli_smoke.f<-as.factor(db$hli_smoke)
levels(db$hli_smoke.f)<-c("Current, 16+ cig/day", "Current, 1-15 cig/day","Former, quit <= 10 years","Former, quit 10+ years", "Never")
db$l_school[db$l_school==5]<-NA
db$l_school<-as_factor(db$l_school)
db$smoke_stat[db$smoke_stat==4]<-NA
db$smoke_stat<-as.factor(db$smoke_stat)
db$hli_dietscore_c<-as.factor(db$hli_dietscore_c)


setwd("/data/Epic/subprojects/Somalogic/work/Oliver")
load("output/imputedcovariates.Rdata")
imputedcovarites$idepic<-imputedcovarites$Idepic
imputedcovarites$Idepic<-NULL
db<-merge(db, imputedcovarites, all.x=T, by= "idepic")


setwd("/data/Epic/subprojects/Somalogic/work/Oliver")
data_clocks=read.csv("output/data_clocks_new.csv")

clock_names=c("Oh","Organismal","Brain","Adipose", "Artery", "Immune","Heart","Intestine","Kidney","Liver","Lung","Muscle","Pancreas","Lehallier","Tanaka", "Wang", "Sathyan", "Global")


#################

## Cox cvd

###Prentice weighting for main analysis#######
db1<-db[db$cvd_t2d_coh==0 & db$cvd_1_status==1,]
db1$age_recr <- db1$age_exit_cvd_1 - 1e-4 #prentice weights
db2<-db[db$cvd_t2d_coh==1 & db$cvd_1_status==1,]
db3<-db[db$cvd_t2d_coh==1 & db$cvd_1_status==0,]
db.spec=rbind(db1,db2,db3)
remove(db1,db2,db3)

###run below instead for sensitivity analysis excluding first 5/2 years of events#
####look at followup time for exclusion
#db$cvd_fut<-db$age_exit_cvd_1-db$age_blood
#db$pre_5_year<-ifelse(db$cvd_fut <= 5,1,0)
#db$pre_2_year<-ifelse(db$cvd_fut <= 2,1,0)

#dba<-db[db$cvd_1_status==1,]
##dba<-dba[dba$pre_5_year==0,]
#dba<-dba[dba$pre_2_year==0,]
#dbb<-db[db$cvd_1_status==0,]
#db.c=rbind(dba,dbb)
#remove(dba,dbb)

#db1<-db.c[db.c$cvd_t2d_coh==0 & db.c$cvd_1_status==1,]
#db1$age_recr <- db1$age_exit_cvd_1 - 1e-4#prentice weighting
#db2<-db.c[db.c$cvd_t2d_coh==1 & db.c$cvd_1_status==1,]
#db3<-db.c[db.c$cvd_t2d_coh==1 & db.c$cvd_1_status==0,]
#db.spec=rbind(db1,db2,db3)
#remove(db1,db2,db3, db.c)


## Cox cvd


#################################

Res_cvd <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                     SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
Res_cvd_60 <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                        SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())



for (clock_name in clock_names){
  
  ### merge dfs
  data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
  merged=merge(x=db.spec,y=data_plot,by="idepic")
 
  #####exclude ever smokers
 # merged<-merged[merged$smoke_stat=="Never",]
  #############
  ### adjusted model combinations for age, sex and center/country
  
  mod_70plus <- coxph(Surv(age_recr, age_exit_cvd_1, cvd_1_status) ~  AgeGap_zscored + strata(age_cat_5, center, sex.x) + cluster(idepic),
                      data= merged)
  mod_60 <- coxph(Surv(age_recr, age_exit_cvd_1, cvd_1_status) ~  AgeGap_zscored + strata(age_cat_5, center, sex.x) + cluster(idepic),
                  data= merged[which(merged$age_recr >60),])
  
  
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
  
  Res_cvd <- bind_rows(Res_cvd, temptib_70plus)
  Res_cvd_60 <- bind_rows(Res_cvd_60, temptib_60)
  
}

Res_cvd <- Res_cvd %>% mutate(FDR = p.adjust(pval, method = "BH"))
Res_cvd_60 <- Res_cvd_60 %>% mutate(FDR = p.adjust(pval, method = "BH"))


#########risk factor adjusted CVD#############

covs<-c( "l_school", "bmi_c" ,  "alc_re",  "pa_mets" , "hli_dietscore_c" ,   "smoke_stat")
covs.imp<-c( "L_School.imp", "bmi_c" ,  "Alc_Re.imp",  "Pa_Mets.imp" ,  "Smoke_Stat.imp","Hli_Dietscore_C.imp" )

Res_cvd_adj <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                         SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
Res_cvd_imp <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())

for (clock_name in clock_names){
  
  ### merge dfs
  data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
  merged=merge(x=db.spec,y=data_plot,by="idepic")
  
  #####exclude ever smokers
 # merged<-merged[merged$smoke_stat=="Never",]
  #############
   
  ### adjusted model combinations for age, sex and center/country
  
  form<-as.formula(paste0("Surv(age_recr, age_exit_cvd_1, cvd_1_status) ~ AgeGap_zscored +", paste(covs, collapse = "+"),"+ strata(age_cat_5, center, sex.x) + cluster(idepic)"))
  mod_70plus <- coxph(form, data=merged)

  
  temptib_70plus <- tibble(Clock= clock_name, Asso = summary(mod_70plus)$coef[1, 'exp(coef)'],
                           beta = summary(mod_70plus)$coef[1, 'coef'],
                           pval = summary(mod_70plus)$coef[1, 'Pr(>|z|)'],
                           SE= summary(mod_70plus)$coef[1, 'se(coef)'] ,
                           LCL =exp(confint(mod_70plus)[1]),
                           UCL =exp(confint(mod_70plus)[1,2]),
                           n =mod_70plus$n ,nevent=mod_70plus$nevent)
  
  
  Res_cvd_adj <- bind_rows(Res_cvd_adj, temptib_70plus)
  
  form<-as.formula(paste0("Surv(age_recr, age_exit_cvd_1, cvd_1_status) ~ AgeGap_zscored +", paste(covs.imp, collapse = "+"),"+ strata(age_cat_5, center, sex.x) + cluster(idepic)"))
  mod_70plus <- coxph(form, data=merged)
 
  
  temptib_70plus <- tibble(Clock= clock_name, Asso = summary(mod_70plus)$coef[1, 'exp(coef)'],
                           beta = summary(mod_70plus)$coef[1, 'coef'],
                           pval = summary(mod_70plus)$coef[1, 'Pr(>|z|)'],
                           SE= summary(mod_70plus)$coef[1, 'se(coef)'] ,
                           LCL =exp(confint(mod_70plus)[1]),
                           UCL =exp(confint(mod_70plus)[1,2]),
                           n =mod_70plus$n ,nevent=mod_70plus$nevent)
  
  
  Res_cvd_imp <- bind_rows(Res_cvd_imp, temptib_70plus)
}

Res_cvd_adj <- Res_cvd_adj %>% mutate(FDR = p.adjust(pval, method = "BH"))
Res_cvd_imp <- Res_cvd_imp %>% mutate(FDR = p.adjust(pval, method = "BH"))


## Cox chd#################

########## prentice weighting for main analysis

db1<-db[db$cvd_t2d_coh==0 & db$chd_1_status==1,]
db1$age_recr <- db1$age_exit_chd_1 - 1e-4
db2<-db[db$cvd_t2d_coh==1 & db$chd_1_status==1,]
db3<-db[db$cvd_t2d_coh==1 & db$chd_1_status==0,]
db.spec=rbind(db1,db2,db3)
remove(db1,db2,db3)

#############exclude pre 5 or 2 yr as sensitivity analysis####
#db$chd_fut<-db$age_exit_chd_1-db$age_blood
#db$pre_5_year<-ifelse(db$chd_fut <= 5,1,0)
#db$pre_2_year<-ifelse(db$chd_fut <= 2,1,0)
#dba<-db[db$chd_1_status==1,]
##dba<-dba[dba$pre_5_year==0,]
#dba<-dba[dba$pre_2_year==0,]
#dbb<-db[db$chd_1_status==0,]
#db.c=rbind(dba,dbb)
#remove(dba,dbb)

#db1<-db.c[db.c$cvd_t2d_coh==0 & db.c$chd_1_status==1,]
#db1$age_recr <- db1$age_exit_chd_1 - 1e-4#prentice weighting
#db2<-db.c[db.c$cvd_t2d_coh==1 & db.c$chd_1_status==1,]
#db3<-db.c[db.c$cvd_t2d_coh==1 & db.c$chd_1_status==0,]
#db.spec=rbind(db1,db2,db3)
#remove(db1,db2,db3, db.c)
##################################

Res_chd <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                  SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
Res_chd_60 <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                     SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())



for (clock_name in clock_names){
  
  ### merge dfs
  data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
  merged=merge(x=db.spec,y=data_plot,by="idepic")
  
  #####exclude ever smokers
  #merged<-merged[merged$smoke_stat=="Never",]
  #############
  
  ### adjusted model combinations for age, sex and center/country
  
  mod_70plus <- coxph(Surv(age_recr, age_exit_chd_1, chd_1_status) ~  AgeGap_zscored + strata(age_cat_5, center, sex.x) + cluster(idepic),
                      data= merged)
  mod_60 <- coxph(Surv(age_recr, age_exit_chd_1, chd_1_status) ~  AgeGap_zscored + strata(age_cat_5, center, sex.x) + cluster(idepic),
                  data= merged[which(merged$age_recr >60),])
  
  
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
  
  Res_chd <- bind_rows(Res_chd, temptib_70plus)
  Res_chd_60 <- bind_rows(Res_chd_60, temptib_60)
  
}

Res_chd <- Res_chd %>% mutate(FDR = p.adjust(pval, method = "BH"))
Res_chd_60 <- Res_chd_60 %>% mutate(FDR = p.adjust(pval, method = "BH"))


#########risk factor adjusted chd#############

covs<-c( "l_school", "bmi_c" ,  "alc_re",  "pa_mets" , "hli_dietscore_c" ,   "smoke_stat")
covs.imp<-c( "L_School.imp", "bmi_c" ,  "Alc_Re.imp",  "Pa_Mets.imp" ,  "Smoke_Stat.imp","Hli_Dietscore_C.imp" )

Res_chd_adj <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
Res_chd_imp <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())

for (clock_name in clock_names){
  
  ### merge dfs
  data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
  merged=merge(x=db.spec,y=data_plot,by="idepic")
  
  #####exclude ever smokers
  #merged<-merged[merged$smoke_stat=="Never",]
  #############
  
  ### adjusted model combinations for age, sex and center/country
  
  form<-as.formula(paste0("Surv(age_recr, age_exit_chd_1, chd_1_status) ~ AgeGap_zscored +", paste(covs, collapse = "+"),"+ strata(age_cat_5, center, sex.x) + cluster(idepic)"))
  mod_70plus <- coxph(form, data=merged)
  
  temptib_70plus <- tibble(Clock= clock_name, Asso = summary(mod_70plus)$coef[1, 'exp(coef)'],
                           beta = summary(mod_70plus)$coef[1, 'coef'],
                           pval = summary(mod_70plus)$coef[1, 'Pr(>|z|)'],
                           SE= summary(mod_70plus)$coef[1, 'se(coef)'] ,
                           LCL =exp(confint(mod_70plus)[1]),
                           UCL =exp(confint(mod_70plus)[1,2]),
                           n =mod_70plus$n ,nevent=mod_70plus$nevent)
  
  
  Res_chd_adj <- bind_rows(Res_chd_adj, temptib_70plus)
  
  form<-as.formula(paste0("Surv(age_recr, age_exit_chd_1, chd_1_status) ~ AgeGap_zscored +", paste(covs.imp, collapse = "+"),"+ strata(age_cat_5, center, sex.x) + cluster(idepic)"))
  mod_70plus <- coxph(form, data=merged)
  
  temptib_70plus <- tibble(Clock= clock_name, Asso = summary(mod_70plus)$coef[1, 'exp(coef)'],
                           beta = summary(mod_70plus)$coef[1, 'coef'],
                           pval = summary(mod_70plus)$coef[1, 'Pr(>|z|)'],
                           SE= summary(mod_70plus)$coef[1, 'se(coef)'] ,
                           LCL =exp(confint(mod_70plus)[1]),
                           UCL =exp(confint(mod_70plus)[1,2]),
                           n =mod_70plus$n ,nevent=mod_70plus$nevent)
  
  
  Res_chd_imp <- bind_rows(Res_chd_imp, temptib_70plus)
}

Res_chd_adj <- Res_chd_adj %>% mutate(FDR = p.adjust(pval, method = "BH"))
Res_chd_imp <- Res_chd_imp %>% mutate(FDR = p.adjust(pval, method = "BH"))





## Cox stroke

##########Prentice weighting for main analysis

db1<-db[db$cvd_t2d_coh==0 & db$stroke_1_status==1,]
db1$age_recr <- db1$age_exit_stroke_1 - 1e-4
db2<-db[db$cvd_t2d_coh==1 & db$stroke_1_status==1,]
db3<-db[db$cvd_t2d_coh==1 & db$stroke_1_status==0,]
db.spec=rbind(db1,db2,db3)
remove(db1,db2,db3)

###exlu pre 5 or 2 year events for sensitivity analysis####
#db$stroke_fut<-db$age_exit_stroke_1-db$age_blood
#db$pre_5_year<-ifelse(db$stroke_fut <= 5,1,0)
#db$pre_2_year<-ifelse(db$stroke_fut <= 2,1,0)
#dba<-db[db$stroke_1_status==1,]
##dba<-dba[dba$pre_5_year==0,]
#dba<-dba[dba$pre_2_year==0,]
#dbb<-db[db$stroke_1_status==0,]
#db.c=rbind(dba,dbb)
#remove(dba,dbb)

#db1<-db.c[db.c$cvd_t2d_coh==0 & db.c$stroke_1_status==1,]
#db1$age_recr <- db1$age_exit_stroke_1 - 1e-4#prentice weighting
#db2<-db.c[db.c$cvd_t2d_coh==1 & db.c$stroke_1_status==1,]
#db3<-db.c[db.c$cvd_t2d_coh==1 & db.c$stroke_1_status==0,]
#db.spec=rbind(db1,db2,db3)
#remove(db1,db2,db3, db.c)
################################

Res_stroke <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                  SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
Res_stroke_60 <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                     SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())



for (clock_name in clock_names){
  
  ### merge dfs
  data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
  merged=merge(x=db.spec,y=data_plot,by="idepic")
  
  #####exclude ever smokers
 # merged<-merged[merged$smoke_stat=="Never",]
  #############
  
  ### adjusted model combinations for age, sex and center/country
  
  mod_70plus <- coxph(Surv(age_recr, age_exit_stroke_1, stroke_1_status) ~  AgeGap_zscored + strata(age_cat_5, center, sex.x) + cluster(idepic),
                      data= merged)
  mod_60 <- coxph(Surv(age_recr, age_exit_stroke_1, stroke_1_status) ~  AgeGap_zscored + strata(age_cat_5, center, sex.x) + cluster(idepic),
                  data= merged[which(merged$age_recr >60),])
  
  
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
  
  Res_stroke <- bind_rows(Res_stroke, temptib_70plus)
  Res_stroke_60 <- bind_rows(Res_stroke_60, temptib_60)
  
}

Res_stroke <- Res_stroke %>% mutate(FDR = p.adjust(pval, method = "BH"))
Res_stroke_60 <- Res_stroke_60 %>% mutate(FDR = p.adjust(pval, method = "BH"))

FUT=merged$age_exit_stroke_1-merged$age_recr
mod_age<- coxph(Surv(FUT, stroke_1_status) ~ scale(age_recr) + strata( center, sex.x) + cluster(idepic),
                data= merged)

agetib <- tibble(Clock= "Chronological", Asso = summary(mod_age)$coef[1, 'exp(coef)'],
                 beta = summary(mod_age)$coef[1, 'coef'],
                 pval = summary(mod_age)$coef[1, 'Pr(>|z|)'],
                 SE= summary(mod_age)$coef[1, 'se(coef)'] ,
                 LCL =exp(confint(mod_age)[1]),
                 UCL =exp(confint(mod_age)[2]),
                 n =mod_age$n ,nevent=mod_age$nevent,
                 FDR=NA)

FUT=FUT[which(merged$age_recr >60)]
mod_age_60<- coxph(Surv(FUT, stroke_1_status) ~ scale(age_recr) + strata( center, sex.x) + cluster(idepic),
                   data= merged[which(merged$age_recr >60),])

agetib_60 <- tibble(Clock= "Chronological", Asso = summary(mod_age_60)$coef[1, 'exp(coef)'],
                    beta = summary(mod_age_60)$coef[1, 'coef'],
                    pval = summary(mod_age_60)$coef[1, 'Pr(>|z|)'],
                    SE= summary(mod_age_60)$coef[1, 'se(coef)'] ,
                    LCL =exp(confint(mod_age_60)[1]),
                    UCL =exp(confint(mod_age_60)[2]),
                    n =mod_age_60$n ,nevent=mod_age_60$nevent,
                    FDR=NA)


Res_stroke <- bind_rows(Res_stroke, agetib)
Res_stroke_60 <- bind_rows(Res_stroke_60, agetib_60)

#########risk factor adjusted stroke#############

covs<-c( "l_school", "bmi_c" ,  "alc_re",  "pa_mets" , "hli_dietscore_c" ,   "smoke_stat")
covs.imp<-c( "L_School.imp", "bmi_c" ,  "Alc_Re.imp",  "Pa_Mets.imp" ,  "Smoke_Stat.imp","Hli_Dietscore_C.imp" )

Res_stroke_adj <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
Res_stroke_imp <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                         SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())

for (clock_name in clock_names){
  
  ### merge dfs
  data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
  merged=merge(x=db.spec,y=data_plot,by="idepic")
  
  #####exclude ever smokers
  #merged<-merged[merged$smoke_stat=="Never",]
  #############
  
  ### adjusted model combinations for age, sex and center/country
  
  form<-as.formula(paste0("Surv(age_recr, age_exit_stroke_1, stroke_1_status) ~ AgeGap_zscored +", paste(covs, collapse = "+"),"+ strata(age_cat_5, center, sex.x) + cluster(idepic)"))
  mod_70plus <- coxph(form, data=merged)
  
  temptib_70plus <- tibble(Clock= clock_name, Asso = summary(mod_70plus)$coef[1, 'exp(coef)'],
                           beta = summary(mod_70plus)$coef[1, 'coef'],
                           pval = summary(mod_70plus)$coef[1, 'Pr(>|z|)'],
                           SE= summary(mod_70plus)$coef[1, 'se(coef)'] ,
                           LCL =exp(confint(mod_70plus)[1]),
                           UCL =exp(confint(mod_70plus)[1,2]),
                           n =mod_70plus$n ,nevent=mod_70plus$nevent)
  
  
  Res_stroke_adj <- bind_rows(Res_stroke_adj, temptib_70plus)
  
  form<-as.formula(paste0("Surv(age_recr, age_exit_stroke_1, stroke_1_status) ~ AgeGap_zscored +", paste(covs.imp, collapse = "+"),"+ strata(age_cat_5, center, sex.x) + cluster(idepic)"))
  mod_70plus <- coxph(form,data=merged)
  
  temptib_70plus <- tibble(Clock= clock_name, Asso = summary(mod_70plus)$coef[1, 'exp(coef)'],
                           beta = summary(mod_70plus)$coef[1, 'coef'],
                           pval = summary(mod_70plus)$coef[1, 'Pr(>|z|)'],
                           SE= summary(mod_70plus)$coef[1, 'se(coef)'] ,
                           LCL =exp(confint(mod_70plus)[1]),
                           UCL =exp(confint(mod_70plus)[1,2]),
                           n =mod_70plus$n ,nevent=mod_70plus$nevent)
  
  
  Res_stroke_imp <- bind_rows(Res_stroke_imp, temptib_70plus)
}

Res_stroke_adj <- Res_stroke_adj %>% mutate(FDR = p.adjust(pval, method = "BH"))
Res_stroke_imp <- Res_stroke_imp %>% mutate(FDR = p.adjust(pval, method = "BH"))

#####
#save(Res_stroke,Res_stroke_60,Res_stroke_adj,Res_stroke_imp,
#  Res_chd,Res_chd_60,Res_chd_adj,Res_chd_imp,
#  Res_cvd,Res_cvd_60, Res_cvd_adj,Res_cvd_imp,file="output/CVDResults.Rdata" )
save(Res_stroke,Res_stroke_60,Res_stroke_adj,Res_stroke_imp,
     Res_chd,Res_chd_60,Res_chd_adj,Res_chd_imp,
     Res_cvd,Res_cvd_60, Res_cvd_adj,Res_cvd_imp,file="output/CVDResults_QC_v2_con2_pre2excl.Rdata")
#file="output/CVDResults_QC_v2_con2_neversmokers.Rdata")
#file="output/CVDResults_QC_v2_con2_pre5excl.Rdata"
#file="output/CVDResults_QC_v2_con2.Rdata"
#file="output/CVDResults_QC_v2_binder.Rdata" )
#file="output/CVDResults_QC_v2.Rdata" )
#file="output/CVDResults_QC.Rdata" )

###descriptive table 1##############
names(db)
table(db$cntr_f)
db$cntr_f<-droplevels(db$cntr_f)
table(db$country)
db$country<-droplevels(db$country)
db$cvd_fut<-db$age_exit_cvd_1-db$age_blood
summary(db$cvd_fut)
label(db$cvd_fut) <- "Follow up time CVD (yrs)"
db$chd_fut<-db$age_exit_chd_1-db$age_blood
label(db$chd_fut) <- "Follow up time CHD (yrs)"
db$stroke_fut<-db$age_exit_stroke_1-db$age_blood
label(db$stroke_fut) <- "Follow up time stroke (yrs)"
label(db$cntr_f) <- "Centre"

load(file="output/outlierIDs.Rdata")
tabledb<-db[db$cvd_t2d_coh==0, ]
#tabledb<-tabledb[-which(tabledb$subjectid %in% checkids),]
tabledb<-to_factor(tabledb)

write.csv(
  x =table1(~age_blood+
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
         cvd_1_status+
         cvd_fut+
         chd_1_status+
         chd_fut+
         stroke_1_status+
         stroke_fut+
         fasting_c
       |country, data=tabledb, render.continuous="Mean (SD)"),
file = "/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/descriptive_tables/cvd_sample.csv",
row.names = FALSE  )





save.image(file="ClockCVD.Rdata")
