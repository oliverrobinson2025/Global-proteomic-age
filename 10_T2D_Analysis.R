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



t2d_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP <- readRDS("/data/Epic/subprojects/Somalogic/work/Vivian/Phase2/Data/ForAnalyses/t2d_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP.rds")#import dataset prepared for T2D analysis
db<-t2d_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP[,c(1:220,7817:7827)]#take only covariates
names(db)
remove(t2d_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP)

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



data_clocks=read.csv("output/data_clocks_new.csv")



clock_names=c("Conventional","Organismal","Brain","Adipose", "Artery", "Immune","Heart","Intestine","Kidney","Liver","Lung","Muscle","Pancreas","Lehallier","Tanaka", "Wang", "Sathyan", "Consensus")

## Cox t2d
####look at followup time for exclusion
hist(db$t2d_fut[db$t2d_status==1])
db$pre_5_year<-ifelse(db$t2d_fut <= 5,1,0)
boxplot(db$t2d_fut~db$pre_5_year )
table(db$t2d_status,db$pre_5_year )

db$pre_2_year<-ifelse(db$t2d_fut <= 2,1,0)
boxplot(db$t2d_fut~db$pre_2_year )
table(db$t2d_status,db$pre_2_year )
##########

db1<-db[db$cvd_t2d_coh==0 & db$t2d_status==1,]
db1$age_recr <- db1$age_exit_t2d - 1e-4 #only for Prentice weighting
db2<-db[db$cvd_t2d_coh==1 & db$t2d_status==1,]
db3<-db[db$cvd_t2d_coh==1 & db$t2d_status==0,]
db.spec=rbind(db1,db2,db3)
remove(db1,db2,db3)

###for exclude pre5 /2 years
#dba<-db[db$t2d_status==1,]
##dba<-dba[dba$pre_5_year==0,]
#dba<-dba[dba$pre_2_year==0,]
#dbb<-db[db$t2d_status==0,]
#db.c=rbind(dba,dbb)
#remove(dba,dbb)

#db1<-db.c[db.c$cvd_t2d_coh==0 & db.c$t2d_status==1,]
#db1$age_recr <- db1$age_exit_t2d - 1e-4#prentice weighting
#db2<-db.c[db.c$cvd_t2d_coh==1 & db.c$t2d_status==1,]
#db3<-db.c[db.c$cvd_t2d_coh==1 & db.c$t2d_status==0,]
#db.spec=rbind(db1,db2,db3)
#remove(db1,db2,db3, db.c)

################################

db.spec=db.spec[db.spec$t2d_prev==0,]### exclude prevalent T2d

Res_t2d <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                  SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
Res_t2d_60 <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                     SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())



for (clock_name in clock_names){
  
  ### merge dfs
  data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
  merged=merge(x=db.spec,y=data_plot,by="idepic")
  
  #####exclude ever smokers
  #merged<-merged[merged$smoke_stat=="Never",]
  #############
  
  ### adjusted model combinations for age, sex and center/country
  
  mod_70plus <- coxph(Surv(age_recr, age_exit_t2d, t2d_status) ~  AgeGap_zscored + strata(age_cat_5, center, sex.x) + cluster(idepic),
                      data= merged)
  mod_60 <- coxph(Surv(age_recr, age_exit_t2d, t2d_status) ~  AgeGap_zscored + strata(age_cat_5, center, sex.x) + cluster(idepic),
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
  
  Res_t2d <- bind_rows(Res_t2d, temptib_70plus)
  Res_t2d_60 <- bind_rows(Res_t2d_60, temptib_60)
  
}

Res_t2d <- Res_t2d %>% mutate(FDR = p.adjust(pval, method = "BH"))
Res_t2d_60 <- Res_t2d_60 %>% mutate(FDR = p.adjust(pval, method = "BH"))



###########
#########risk factor adjusted t2d#############

covs<-c( "l_school", "bmi_c" ,  "alc_re",  "pa_mets" , "hli_dietscore_c" ,   "smoke_stat")
covs.imp<-c( "L_School.imp", "bmi_c" ,  "Alc_Re.imp",  "Pa_Mets.imp" ,  "Smoke_Stat.imp","Hli_Dietscore_C.imp" )

Res_t2d_adj <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                         SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())

Res_t2d_imp <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())

for (clock_name in clock_names){
  
  ### merge dfs
  data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
  merged=merge(x=db.spec,y=data_plot,by="idepic")
  
  #####exclude ever smokers
  #merged<-merged[merged$smoke_stat=="Never",]
  #############
  
  ### adjusted model combinations for age, sex and center/country
  
  form<-as.formula(paste0("Surv(age_recr, age_exit_t2d, t2d_status) ~ AgeGap_zscored +", paste(covs, collapse = "+"),"+ strata(age_cat_5, center, sex.x) + cluster(idepic)"))
  mod_70plus <- coxph(form, data=merged)
  
  
  temptib_70plus <- tibble(Clock= clock_name, Asso = summary(mod_70plus)$coef[1, 'exp(coef)'],
                           beta = summary(mod_70plus)$coef[1, 'coef'],
                           pval = summary(mod_70plus)$coef[1, 'Pr(>|z|)'],
                           SE= summary(mod_70plus)$coef[1, 'se(coef)'] ,
                           LCL =exp(confint(mod_70plus)[1]),
                           UCL =exp(confint(mod_70plus)[1,2]),
                           n =mod_70plus$n ,nevent=mod_70plus$nevent)
  
  
  Res_t2d_adj <- bind_rows(Res_t2d_adj, temptib_70plus)
  
  form<-as.formula(paste0("Surv(age_recr, age_exit_t2d, t2d_status) ~ AgeGap_zscored +", paste(covs.imp, collapse = "+"),"+ strata(age_cat_5, center, sex.x) + cluster(idepic)"))
  mod_70plus <- coxph(form, data=merged)
  
  
  temptib_70plus <- tibble(Clock= clock_name, Asso = summary(mod_70plus)$coef[1, 'exp(coef)'],
                           beta = summary(mod_70plus)$coef[1, 'coef'],
                           pval = summary(mod_70plus)$coef[1, 'Pr(>|z|)'],
                           SE= summary(mod_70plus)$coef[1, 'se(coef)'] ,
                           LCL =exp(confint(mod_70plus)[1]),
                           UCL =exp(confint(mod_70plus)[1,2]),
                           n =mod_70plus$n ,nevent=mod_70plus$nevent)
  
  
  Res_t2d_imp <- bind_rows(Res_t2d_imp, temptib_70plus)
}

Res_t2d_adj <- Res_t2d_adj %>% mutate(FDR = p.adjust(pval, method = "BH"))
Res_t2d_imp <- Res_t2d_imp %>% mutate(FDR = p.adjust(pval, method = "BH"))

###
#save(Res_t2d,Res_t2d_60,Res_t2d_adj, file="output/T2DResults.Rdata" )
#save(Res_t2d,Res_t2d_60,Res_t2d_adj,Res_t2d_imp, file="output/T2DResults_no_prev_t2d.Rdata" )
save(Res_t2d,Res_t2d_60,Res_t2d_adj,Res_t2d_imp, file="output/T2DResults_no_prev_t2d_QC_v2_con2_pre2excl.Rdata" )
#file="output/T2DResults_no_prev_t2d_QC_v2_con2_neversmokers.Rdata" )
#file="output/T2DResults_no_prev_t2d_QC_v2_con2_pre5excl.Rdata" 
#file="output/T2DResults_no_prev_t2d_QC_v2_con2.Rdata"
#file="output/T2DResults_no_prev_t2d_QC_v2.Rdata" )
#file="output/T2DResults_no_prev_t2d_QC.Rdata" )

###### descriptive table#####
names(db)
table(db$cntr_f)
db$cntr_f<-droplevels(db$cntr_f)
table(db$country)
db$country<-droplevels(db$country)
db$t2d_fut<-db$age_exit_t2d-db$age_blood
summary(db$t2d_fut)
label(db$t2d_fut) <- "Follow up time t2d (yrs)"
db$l_school[db$l_school=="Not specified"]<-NA
db$l_school<-droplevels(db$l_school)
db$smoke_stat<-droplevels(db$smoke_stat)


####

tabledb<-db[db$cvd_t2d_coh==0 & db$t2d_status==1,]
tabledb<-tabledb[which(tabledb$idepic %in% data_clocks$idepic),]
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
         t2d_status+
         t2d_fut+
         fasting_c
       |country, data=tabledb, render.continuous="Mean (SD)"),
file = "/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/descriptive_tables/t2d_sample.csv",
row.names = FALSE  )






save.image(file="ClockT2D.Rdata")
