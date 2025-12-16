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

cancer_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP <- readRDS("/data/Epic/subprojects/Somalogic/work/Vivian/Phase2/Data/ForAnalyses/cancer_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP.rds")#import dataset prepared for cancer analysis
db<-cancer_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP[,c(1:220,7817:7849)]#take only covariates
remove(cancer_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP)

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

data_clocks=data_clocks %>% mutate(Sex_F=as.factor(Sex_F))
table(data_clocks$Organ)

clock_names=c("Oh","Organismal","Brain","Adipose", "Artery", "Immune","Heart","Intestine","Kidney","Liver","Lung","Muscle","Pancreas","Lehallier","Tanaka", "Wang", "Sathyan", "Global")

##check subtypes
table(db$typ_tumo_lymp)
db$cncr_mal_h_lymp<-ifelse(db$typ_tumo_lymp== "Classical"|db$typ_tumo_lymp== "Malignant lymphoma"|db$typ_tumo_lymp== "NLPHL",1,0)
table(db$cncr_mal_h_lymp)
db$cncr_mal_nh_lymp<-ifelse(db$typ_tumo_lymp== "Mature B-NHL"|db$typ_tumo_lymp== "Mature T- or NK-NHL"|db$typ_tumo_lymp== "NHL"|db$typ_tumo_lymp== "Precursor NHL",1,0)
table(db$cncr_mal_nh_lymp)

cncr_mal_blad
table(db$typ_tumo_uadt[db$cncr_mal_uadt=="1"])

db$cncr_mal_esophag<-ifelse(db$typ_tumo_uadt== "04-Esophagus squamous cell carcinoma"|db$typ_tumo_uadt== "05-Esophagus other",1,0)
table(db$cncr_mal_esophag)
db$cncr_mal_orophag<-ifelse(db$typ_tumo_uadt== "01-UADT squamous cell carcinoma"|db$typ_tumo_uadt== "02-UADT other",1,0)
table(db$cncr_mal_orophag)

####look at followup time for exclusion of first 2 or 5 years of events in sensitivity analysis only
db$cancer_fut<-db$age_exit_cancer_1st-db$age_recr
db$pre_5_year<-ifelse(db$cancer_fut <= 5,1,0)
db$pre_2_year<-ifelse(db$cancer_fut <= 2,1,0)


####cox cancer
cancertypes=c("cncr_mal_anyc","cncr_mal_blad", "cncr_mal_brea","cncr_mal_brea_peri",
              "cncr_mal_brea_post","cncr_mal_brea_pre", "cncr_mal_clrt", "cncr_mal_clrt_colon",
              "cncr_mal_clrt_rectum", "cncr_mal_coru", "cncr_mal_glio", "cncr_mal_kidn", 
              "cncr_mal_live", "cncr_mal_lung", "cncr_mal_lymp", "cncr_mal_mela", 
              "cncr_mal_ovar", "cncr_mal_panc", "cncr_mal_pros", "cncr_mal_stom",
              "cncr_mal_thyr", "cncr_mal_uadt",
              "cncr_mal_orophag","cncr_mal_esophag", "cncr_mal_nh_lymp", "cncr_mal_h_lymp")
             
cancer_results = list()
cancer_results_60 = list()

for(cancertype in cancertypes){


db1<-db[db$cvd_t2d_coh==0 & db[,cancertype]==1,]
db1$age_recr <- db1$age_exit_cancer_1st - 1e-4 ## Prentice weighting
db2<-db[db$cvd_t2d_coh==1 & db[,cancertype]==1,]
db3<-db[db$cvd_t2d_coh==1 & db[,cancertype]==0,]
db.spec=rbind(db1,db2,db3)
remove(db1,db2,db3)
  
  ###exclude pre 5 or 2 year in sensitivity analysis
  
  #dba<-db[db[,cancertype]==1,]
  ##dba<-dba[dba$pre_5_year==0,]
  #dba<-dba[dba$pre_2_year==0,]
 # dbb<-db[db[,cancertype]==0,]
  #db.c=rbind(dba,dbb)
 # remove(dba,dbb)
  
 # db1<-db.c[db.c$cvd_t2d_coh==0 & db.c[,cancertype]==1,]
 # db1$age_recr <- db1$age_exit_cancer_1st - 1e-4#prentice weighting
 # db2<-db.c[db.c$cvd_t2d_coh==1 & db.c[,cancertype]==1,]
 #db3<-db.c[db.c$cvd_t2d_coh==1 & db.c[,cancertype]==0,]
#  db.spec=rbind(db1,db2,db3)
# remove(db1,db2,db3, db.c)
  
  #######


Res_cancer <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                          SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
Res_cancer_60 <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                             SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())


for (clock_name in clock_names){
  
  ### merge dfs
  data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
  merged=merge(x=db.spec,y=data_plot,by="idepic")
  
  #####exclude ever smokers in sensitivity analysis
 # merged<-merged[merged$smoke_stat=="Never",]
  #############
  
  merged$cncr_mal=as.vector(merged[,cancertype])
  if(cancertype == "cncr_mal_pros"){merged=merged[merged$Sex_F==0,]}
  if(cancertype == "cncr_mal_brea" 
     |cancertype == "cncr_mal_brea_peri" 
     |cancertype == "cncr_mal_brea_post" 
     |cancertype == "cncr_mal_brea_pre" 
     |cancertype == "cncr_mal_ovar" 
     |cancertype =="cncr_mal_coru"){merged=merged[merged$Sex_F==1,]}
  
  ### adjusted model combinations for age, sex and center/country
  
   mod_70plus <- coxph(Surv(age_recr, age_exit_cancer_1st, cncr_mal) ~  AgeGap_zscored + strata(age_cat_5, center, sex.x) + cluster(idepic),
                      data= merged)
  mod_60 <- coxph(Surv(age_recr, age_exit_cancer_1st, cncr_mal) ~  AgeGap_zscored + strata(age_cat_5, center, sex.x) + cluster(idepic),
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
  
  Res_cancer <- bind_rows(Res_cancer,  temptib_70plus)
  Res_cancer_60 <- bind_rows(Res_cancer_60, temptib_60,)
  
}

Res_cancer <- Res_cancer %>% mutate(FDR = p.adjust(pval, method = "BH"))
Res_cancer_60 <- Res_cancer_60 %>% mutate(FDR = p.adjust(pval, method = "BH"))


cancer_results[[cancertype]]<-Res_cancer
cancer_results_60[[cancertype]]<-Res_cancer_60
}


######adjusted cancer results
covs<-c( "l_school", "bmi_c" ,  "alc_re",  "pa_mets" , "hli_dietscore_c" ,   "smoke_stat")
covs.imp<-c( "L_School.imp", "bmi_c" ,  "Alc_Re.imp",  "Pa_Mets.imp" ,  "Smoke_Stat.imp","Hli_Dietscore_C.imp" )

cancer_results_adj = list()
cancer_results_imp = list()

for(cancertype in cancertypes){
  
  
  db1<-db[db$cvd_t2d_coh==0 & db[,cancertype]==1,]
  db1$age_recr <- db1$age_exit_cancer_1st - 1e-4
  db2<-db[db$cvd_t2d_coh==1 & db[,cancertype]==1,]
  db3<-db[db$cvd_t2d_coh==1 & db[,cancertype]==0,]
  db.spec=rbind(db1,db2,db3)
  remove(db1,db2,db3)
  
  ###exclude pre 5 year
  
  #dba<-db[db[,cancertype]==1,]
  ##dba<-dba[dba$pre_5_year==0,]
  #dba<-dba[dba$pre_2_year==0,]
  #dbb<-db[db[,cancertype]==0,]
  #db.c=rbind(dba,dbb)
  #remove(dba,dbb)
  
  #db1<-db.c[db.c$cvd_t2d_coh==0 & db.c[,cancertype]==1,]
  #db1$age_recr <- db1$age_exit_cancer_1st - 1e-4#prentice weighting
  #db2<-db.c[db.c$cvd_t2d_coh==1 & db.c[,cancertype]==1,]
  #db3<-db.c[db.c$cvd_t2d_coh==1 & db.c[,cancertype]==0,]
  #db.spec=rbind(db1,db2,db3)
  #remove(db1,db2,db3, db.c)
  
  #######
  
  
  
  Res_cancer_adj <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                       SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
  
  Res_cancer_imp <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                           SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
  
  for (clock_name in clock_names){
    
    ### merge dfs
    data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
    merged=merge(x=db.spec,y=data_plot,by="idepic")
    
    #####exclude ever smokers
    #merged<-merged[merged$smoke_stat=="Never",]
    #############
    
    merged$cncr_mal=as.vector(merged[,cancertype])
    if(cancertype == "cncr_mal_pros"){merged=merged[merged$Sex_F==0,]}
    if(cancertype == "cncr_mal_brea" 
       |cancertype == "cncr_mal_brea_peri" 
       |cancertype == "cncr_mal_brea_post" 
       |cancertype == "cncr_mal_brea_pre" 
       |cancertype == "cncr_mal_ovar" 
       |cancertype =="cncr_mal_coru"){merged=merged[merged$Sex_F==1,]}
    
    ### adjusted model combinations for age, sex and center/country
    
    form<-as.formula(paste0("Surv(age_recr, age_exit_cancer_1st, cncr_mal) ~ AgeGap_zscored +", paste(covs, collapse = "+"),"+ strata(age_cat_5, center, sex.x) + cluster(idepic)"))
    mod_70plus <- coxph(form, data=merged)
    
    temptib_70plus <- tibble(Clock= clock_name, Asso = summary(mod_70plus)$coef[1, 'exp(coef)'],
                             beta = summary(mod_70plus)$coef[1, 'coef'],
                             pval = summary(mod_70plus)$coef[1, 'Pr(>|z|)'],
                             SE= summary(mod_70plus)$coef[1, 'se(coef)'] ,
                             LCL =exp(confint(mod_70plus)[1]),
                             UCL =exp(confint(mod_70plus)[1,2]),
                             n =mod_70plus$n ,nevent=mod_70plus$nevent)
    
    Res_cancer_adj <- bind_rows(Res_cancer_adj,  temptib_70plus)
    
    form<-as.formula(paste0("Surv(age_recr, age_exit_cancer_1st, cncr_mal) ~ AgeGap_zscored +", paste(covs.imp, collapse = "+"),"+ strata(age_cat_5, center, sex.x) + cluster(idepic)"))
    mod_70plus <- coxph(form, data=merged)
    
    
    temptib_70plus <- tibble(Clock= clock_name, Asso = summary(mod_70plus)$coef[1, 'exp(coef)'],
                             beta = summary(mod_70plus)$coef[1, 'coef'],
                             pval = summary(mod_70plus)$coef[1, 'Pr(>|z|)'],
                             SE= summary(mod_70plus)$coef[1, 'se(coef)'] ,
                             LCL =exp(confint(mod_70plus)[1]),
                             UCL =exp(confint(mod_70plus)[1,2]),
                             n =mod_70plus$n ,nevent=mod_70plus$nevent)
    
    Res_cancer_imp <- bind_rows(Res_cancer_imp,  temptib_70plus)
  }
  
  Res_cancer_adj <- Res_cancer_adj %>% mutate(FDR = p.adjust(pval, method = "BH"))
  Res_cancer_imp <- Res_cancer_imp %>% mutate(FDR = p.adjust(pval, method = "BH"))
  
  cancer_results_adj[[cancertype]]<-Res_cancer_adj
  cancer_results_imp[[cancertype]]<-Res_cancer_imp
 
}

#######

#save(cancer_results,cancer_results_60, cancer_results_adj, cancer_results_imp,file="output/ClockCancerResults.Rdata" )
#save(cancer_results,cancer_results_60, cancer_results_adj, cancer_results_imp,file="output/ClockCancerResults_QC.Rdata" )
#save(cancer_results,cancer_results_60, cancer_results_adj, cancer_results_imp,file="output/ClockCancerResults_QC_v2.Rdata" )
#save(cancer_results,cancer_results_60, cancer_results_adj, cancer_results_imp,file="output/ClockCancerResults_QC_v2_binder.Rdata" )
#save(cancer_results,cancer_results_60, cancer_results_adj, cancer_results_imp,file="output/ClockCancerResults_QC_v2_con2.Rdata" )
#save(cancer_results,cancer_results_60, cancer_results_adj, cancer_results_imp,file="output/ClockCancerResults_QC_v2_con2_pre5excl.Rdata" )
#save(cancer_results,cancer_results_60, cancer_results_adj, cancer_results_imp,file="output/ClockCancerResults_QC_v2_con2_neversmokers.Rdata" )
save(cancer_results,cancer_results_60, cancer_results_adj, cancer_results_imp,file="output/ClockCancerResults_QC_v2_con2_pre2excl.Rdata" )


####table1
names(db)
table(db$cntr_f)
db$cntr_f<-droplevels(db$cntr_f)
table(db$country)
db$country<-droplevels(db$country)
db$cancer_fut<-db$age_exit_cancer_1st-db$age_recru
summary(db$cancer_fut)
label(db$cancer_fut) <- "Follow up time cancer (yrs)"
db$l_school[db$l_school=="Not specified"]<-NA
db$l_school<-droplevels(db$l_school)
db$smoke_stat<-droplevels(db$smoke_stat)


tabledb<-db[db$cvd_t2d_coh==0 & db$cncr_mal_anyc==1,]
tabledb<-tabledb[which(tabledb$idepic %in% data_clocks$idepic),]
tabledb<-to_factor(tabledb)

tablecancers=c("cncr_mal_anyc","cncr_mal_blad", "cncr_mal_brea","cncr_mal_brea_peri",
              "cncr_mal_brea_post","cncr_mal_brea_pre", "cncr_mal_clrt", "cncr_mal_clrt_colon",
              "cncr_mal_clrt_rectum", "cncr_mal_coru", "cncr_mal_glio", "cncr_mal_kidn", 
              "cncr_mal_live", "cncr_mal_lung", "cncr_mal_lymp", "cncr_mal_mela", 
              "cncr_mal_ovar", "cncr_mal_panc", "cncr_mal_pros", "cncr_mal_stom",
              "cncr_mal_thyr", "cncr_mal_uadt")
can_var<-paste(tablecancers, collapse = "+")

form<-as.formula(paste0("~age_blood+ sex+ cntr_f+ l_school+bmi_c+ alc_re+smoke_stat+  pa_mets+ pa_index+ hli_dietscore+ hli_score+ menopause+ t2d_prev+ cvd_prev+", 
          can_var ,
         "+cancer_fut+ fasting_c |country"))


write.csv(
  x =table1(form, data=tabledb, render.continuous="Mean (SD)"),
file = "/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/descriptive_tables/cancer_sample.csv",
row.names = FALSE  )





save.image(file="ClockCancer.Rdata")
