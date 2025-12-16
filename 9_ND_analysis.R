library(tidyverse)
library(survival)
library(ggplot2)
library(cowplot)
library( haven )

setwd("/data/Epic/subprojects/Somalogic/work/Oliver")
data_clocks=read.csv("output/data_clocks_new.csv")
data_clocks=data_clocks %>% mutate(Sex_F=as.factor(Sex_F))



clock_names=c("Oh","Organismal","Brain","Adipose", "Artery", "Immune","Heart","Intestine","Kidney","Liver","Lung","Muscle","Pancreas","Lehallier","Tanaka", "Wang", "Sathyan", "Global")


##readdata
DF_ND=readRDS("/data/Epic/subprojects/Somalogic/work/Jan/Masterfile/20241002_EPIC4ND_Masterfile_merged.rds")#import dataset prepared for EPIC4ND analyses
names(DF_ND)
table(DF_ND$analysis2)#AD
table(DF_ND$analysis3)#all_dem
table(DF_ND$analysis4)#PD
table(DF_ND$analysis8)#ALS

table(DF_ND$analysis3, DF_ND$analysis2)

####add additional variables
setwd("/data/Epic/subprojects/Somalogic/")
mm0  <- read_sas(data_file= "sources/Epi_Data/somalogic_2023.sas7bdat",
                 catalog_file="sources/Epi_Data/formats.sas7bcat")

overlap<-intersect(DF_ND[DF_ND$Has_protein_data==1,]$Idepic, imputed$Idepic)
table(DF_ND$Has_protein_data)

overlap<-intersect(names(DF_ND), names(mm0))
overlap<-overlap[-1]
include<-which( names(mm0) %in% overlap)
includvar<-names(mm0)[-include]

mm0$Idepic[1:10]
DF_ND<-merge(DF_ND, mm0[includvar], all.x=T, by= "Idepic")
remove(mm0)

DF_ND$Smoke_Stat <- 
  factor(DF_ND$Smoke_Stat, levels=c(1,2,3),
         labels=c("Never", "Former", "Smoker"))

summary(DF_ND$Date_Start)
summary(DF_ND$analysis1_EFU)
summary(DF_ND$Age_Start)

ageevent<-(as.numeric(DF_ND$analysis1_EFU-DF_ND$Date_Start))/365.25
summary(ageevent)

table(DF_ND$analysis2, DF_ND$Cvd_T2d_Coh )
table(is.na(DF_ND[DF_ND$Has_protein_data==1,]$Hli_Smoke), DF_ND[DF_ND$Has_protein_data==1,]$Smoke_Stat)#none missing - smoke stat has missing categories
table(is.na(DF_ND[DF_ND$Has_protein_data==1,]$Pa_Mets), DF_ND[DF_ND$Has_protein_data==1,]$Pa_Index)#none missing - index has missing categories
table(DF_ND[DF_ND$Has_protein_data==1,]$Smoke_Stat)

load(file="/data/Epic/subprojects/Somalogic/work/Oliver/output/imputedcovariates.Rdata" )
DF_ND<-merge(DF_ND, imputedcovarites, all.x=T, by= "Idepic")
table(DF_ND[DF_ND$Has_protein_data==1,]$Smoke_Stat.imp)
names(imputedcovarites)

DF_ND$L_School[DF_ND$L_School==5]<-NA
DF_ND$L_School<-as_factor(DF_ND$L_School)
DF_ND$Smoke_Stat[DF_ND$Smoke_Stat==4]<-NA
DF_ND$Smoke_Stat<-as.factor(DF_ND$Smoke_Stat)
DF_ND$Hli_Dietscore_C<-as.factor(DF_ND$Hli_Dietscore_C)
DF_ND$Diabet<-droplevels(DF_ND$Diabet)



#prepare AD dataset

covs<-c( "L_School", "Bmi_C" ,  "Alc_Re",  "Pa_Mets" ,  "Smoke_Stat", "Hli_Dietscore_C")
covs.imp<-c( "L_School.imp", "Bmi_C" ,  "Alc_Re.imp",  "Pa_Mets.imp" ,  "Smoke_Stat.imp","Hli_Dietscore_C.imp" )
table1vars<-c("Pa_Index","Hli_Dietscore","Hli_Score","Menopause","Cvd_Prev","T2D_Prev","Death_Status","Cvd_1_Status","T2D_Status", "Cncr_Mal_Anyc","Neuro_Status","Fasting_C", "Country")
covars<-unique(c(covs,covs.imp,table1vars,"Weights"))

DF_AD_cleaned<-DF_ND[!is.na(DF_ND$analysis2), c("Idepic", "Sex","Age_Start", "Age_Start_cat_5","Date_Start", "Center", "Cvd_T2d_Coh", "analysis2_EFU", "analysis2",covars)]
DF_AD_cleaned$ageevent<-(as.numeric(DF_AD_cleaned$analysis2_EFU-DF_AD_cleaned$Date_Start))/365.25
DF_AD_cleaned$ageevent<-DF_AD_cleaned$ageevent+DF_AD_cleaned$Age_Start
names(DF_AD_cleaned)<-c("idepic", "sex","age", "age_cat_5","Date_Start", "center", "cvd_t2d_coh", "analysis2_EFU", "indevent",covars, "ageevent")



###########Prentice weighting for main analysis###

DF_AD_cleaned1<-DF_AD_cleaned[DF_AD_cleaned$cvd_t2d_coh==0 & DF_AD_cleaned$indevent==1,]
DF_AD_cleaned1$age <- DF_AD_cleaned1$ageevent - 1e-4 ####only for Prentice weights
DF_AD_cleaned2<-DF_AD_cleaned[DF_AD_cleaned$cvd_t2d_coh==1 & DF_AD_cleaned$indevent==1,]
DF_AD_cleaned3<-DF_AD_cleaned[DF_AD_cleaned$cvd_t2d_coh==1 & DF_AD_cleaned$indevent==0,]
remove(DF_AD_cleaned)
DF_AD_cleaned=rbind(DF_AD_cleaned1,DF_AD_cleaned2,DF_AD_cleaned3)
remove(DF_AD_cleaned1,DF_AD_cleaned2,DF_AD_cleaned3)

########exclude pre 5 /2 year in sensitivity analysis only####
#DF_AD_cleaned$fut<-DF_AD_cleaned$ageevent-DF_AD_cleaned$age
#DF_AD_cleaned$pre_5_year<-ifelse(DF_AD_cleaned$fut <= 5,1,0)
#DF_AD_cleaned$pre_2_year<-ifelse(DF_AD_cleaned$fut <= 2,1,0)
#DF_AD_cleaneda<-DF_AD_cleaned[DF_AD_cleaned$indevent==1,]
##DF_AD_cleaneda<-DF_AD_cleaneda[DF_AD_cleaneda$pre_5_year==0,]
#DF_AD_cleaneda<-DF_AD_cleaneda[DF_AD_cleaneda$pre_2_year==0,]
#DF_AD_cleanedb<-DF_AD_cleaned[DF_AD_cleaned$indevent==0,]
#DF_AD_cleaned.c=rbind(DF_AD_cleaneda,DF_AD_cleanedb)
#remove(DF_AD_cleaneda,DF_AD_cleanedb)

#DF_AD_cleaned1<-DF_AD_cleaned.c[DF_AD_cleaned.c$cvd_t2d_coh==0 & DF_AD_cleaned.c$indevent==1,]
#DF_AD_cleaned1$age <- DF_AD_cleaned1$ageevent - 1e-4#prentice weighting
#DF_AD_cleaned2<-DF_AD_cleaned.c[DF_AD_cleaned.c$cvd_t2d_coh==1 & DF_AD_cleaned.c$indevent==1,]
#DF_AD_cleaned3<-DF_AD_cleaned.c[DF_AD_cleaned.c$cvd_t2d_coh==1 & DF_AD_cleaned.c$indevent==0,]
#remove(DF_AD_cleaned)
#DF_AD_cleaned=rbind(DF_AD_cleaned1,DF_AD_cleaned2,DF_AD_cleaned3)
#remove(DF_AD_cleaned1,DF_AD_cleaned2,DF_AD_cleaned3, DF_AD_cleaned.c)

######

cox_clocks=function(data_clocks, DF_AD_cleaned){
  
  clock_names=c("Conventional","Organismal","Brain","Adipose", "Artery", "Immune","Heart","Intestine","Kidney","Liver","Lung","Muscle","Pancreas","Lehallier","Tanaka", "Wang", "Sathyan", "Consensus")
  
  Res_ND <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                   SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
  Res_ND_60 <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
  Res_ND_adj <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                   SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
  Res_ND_imp <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                       SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric())
  
  # cox model
  for (clock_name in clock_names){
    
    ### merge dfs
    data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
    merged=merge(x=DF_AD_cleaned,y=data_plot,by="idepic")
    
    #####exclude ever smokers
    #merged<-merged[merged$Smoke_Stat==1,]
    #############
    
    ### adjusted model combinations for age, sex and center/country
    
    mod_70plus <- coxph(Surv(age, ageevent, indevent) ~  AgeGap_zscored + strata(age_cat_5, center, sex.x) + cluster(idepic), 
                       data= merged) 
    mod_60 <- coxph(Surv(age, ageevent, indevent) ~  AgeGap_zscored + strata(age_cat_5, center, sex.x) + cluster(idepic),
                    data= merged[which(merged$age >60),])
    
     
    form<-as.formula(paste0("Surv(age, ageevent, indevent) ~ AgeGap_zscored +", paste(covs, collapse = "+"),"+ strata(age_cat_5, center, sex.x) + cluster(idepic)"))
    mod_70adj <- coxph(form, data=merged)
  
    
    form<-as.formula(paste0("Surv(age, ageevent, indevent) ~ AgeGap_zscored +", paste(covs.imp, collapse = "+"),"+ strata(age_cat_5, center, sex.x) + cluster(idepic)"))
    mod_70imp <- coxph(form, data=merged)
   
    
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
    
    temptib_adj <- tibble(Clock= clock_name, Asso = summary(mod_70adj)$coef[1, 'exp(coef)'],
                         beta = summary(mod_70adj)$coef[1, 'coef'],
                         pval = summary(mod_70adj)$coef[1, 'Pr(>|z|)'],
                         SE= summary(mod_70adj)$coef[1, 'se(coef)'] ,
                         LCL =exp(confint(mod_70adj)[1]),
                         UCL =exp(confint(mod_70adj)[1,2]),
                         n =mod_70adj$n ,nevent=mod_70adj$nevent)
    
    temptib_imp <- tibble(Clock= clock_name, Asso = summary(mod_70imp)$coef[1, 'exp(coef)'],
                          beta = summary(mod_70imp)$coef[1, 'coef'],
                          pval = summary(mod_70imp)$coef[1, 'Pr(>|z|)'],
                          SE= summary(mod_70imp)$coef[1, 'se(coef)'] ,
                          LCL =exp(confint(mod_70imp)[1]),
                          UCL =exp(confint(mod_70imp)[1,2]),
                          n =mod_70imp$n ,nevent=mod_70imp$nevent)
    
    Res_ND <- bind_rows(Res_ND, temptib_70plus)
    Res_ND_60 <- bind_rows(Res_ND_60, temptib_60)
    Res_ND_adj <- bind_rows(Res_ND_adj, temptib_adj)
    Res_ND_imp <- bind_rows(Res_ND_imp, temptib_imp)
    
  }
  
  Res_ND <- Res_ND %>% mutate(FDR = p.adjust(pval, method = "BH"))
  Res_ND_60 <- Res_ND_60 %>% mutate(FDR = p.adjust(pval, method = "BH"))
  Res_ND_adj <- Res_ND_adj %>% mutate(FDR = p.adjust(pval, method = "BH"))
  Res_ND_imp <- Res_ND_imp %>% mutate(FDR = p.adjust(pval, method = "BH"))
  

  return(list(Res_ND,Res_ND_60,Res_ND_adj, Res_ND_imp))
}


### AD


Res_AD=cox_clocks(data_clocks, DF_AD_cleaned)
Res_AD_full=Res_AD[[1]]
Res_AD_60=Res_AD[[2]]
Res_AD_adj=Res_AD[[3]]
Res_AD_imp=Res_AD[[4]]

print(Res_AD_full[order(Res_AD_full$pval),])
print(Res_AD_60[order(Res_AD_60$pval),])
print(Res_AD_adj[order(Res_AD_adj$pval),])
print(Res_AD_imp[order(Res_AD_imp$pval),])

####prepare Dementia dataset####

DF_DEM_cleaned<-DF_ND[!is.na(DF_ND$analysis3), c("Idepic", "Sex","Age_Start", "Age_Start_cat_5","Date_Start", "Center", "Cvd_T2d_Coh", "analysis3_EFU", "analysis3",covars)]
DF_DEM_cleaned$ageevent<-(as.numeric(DF_DEM_cleaned$analysis3_EFU-DF_DEM_cleaned$Date_Start))/365.25
DF_DEM_cleaned$ageevent<-DF_DEM_cleaned$ageevent+DF_DEM_cleaned$Age_Start
names(DF_DEM_cleaned)<-c("idepic", "sex","age", "age_cat_5","Date_Start", "center", "cvd_t2d_coh", "analysis3_EFU", "indevent",covars, "ageevent")


##############Prentice weighting for main analysis

DF_DEM_cleaned1<-DF_DEM_cleaned[DF_DEM_cleaned$cvd_t2d_coh==0 & DF_DEM_cleaned$indevent==1,]
DF_DEM_cleaned1$age <- DF_DEM_cleaned1$ageevent - 1e-4
DF_DEM_cleaned2<-DF_DEM_cleaned[DF_DEM_cleaned$cvd_t2d_coh==1 & DF_DEM_cleaned$indevent==1,]
DF_DEM_cleaned3<-DF_DEM_cleaned[DF_DEM_cleaned$cvd_t2d_coh==1 & DF_DEM_cleaned$indevent==0,]
remove(DF_DEM_cleaned)
DF_DEM_cleaned=rbind(DF_DEM_cleaned1,DF_DEM_cleaned2,DF_DEM_cleaned3)
remove(DF_DEM_cleaned1,DF_DEM_cleaned2,DF_DEM_cleaned3)

########exclude pre 5 / 2year in sensitivity analysis only####
#DF_DEM_cleaned$fut<-DF_DEM_cleaned$ageevent-DF_DEM_cleaned$age
#DF_DEM_cleaned$pre_5_year<-ifelse(DF_DEM_cleaned$fut <= 5,1,0)
#DF_DEM_cleaned$pre_2_year<-ifelse(DF_DEM_cleaned$fut <= 2,1,0)

#DF_DEM_cleaneda<-DF_DEM_cleaned[DF_DEM_cleaned$indevent==1,]
##DF_DEM_cleaneda<-DF_DEM_cleaneda[DF_DEM_cleaneda$pre_5_year==0,]
#DF_DEM_cleaneda<-DF_DEM_cleaneda[DF_DEM_cleaneda$pre_2_year==0,]
#DF_DEM_cleanedb<-DF_DEM_cleaned[DF_DEM_cleaned$indevent==0,]
#DF_DEM_cleaned.c=rbind(DF_DEM_cleaneda,DF_DEM_cleanedb)
#remove(DF_DEM_cleaneda,DF_DEM_cleanedb)

#DF_DEM_cleaned1<-DF_DEM_cleaned.c[DF_DEM_cleaned.c$cvd_t2d_coh==0 & DF_DEM_cleaned.c$indevent==1,]
#DF_DEM_cleaned1$age <- DF_DEM_cleaned1$ageevent - 1e-4#prentice weighting
#DF_DEM_cleaned2<-DF_DEM_cleaned.c[DF_DEM_cleaned.c$cvd_t2d_coh==1 & DF_DEM_cleaned.c$indevent==1,]
#DF_DEM_cleaned3<-DF_DEM_cleaned.c[DF_DEM_cleaned.c$cvd_t2d_coh==1 & DF_DEM_cleaned.c$indevent==0,]
#remove(DF_DEM_cleaned)
#DF_DEM_cleaned=rbind(DF_DEM_cleaned1,DF_DEM_cleaned2,DF_DEM_cleaned3)
#remove(DF_DEM_cleaned1,DF_DEM_cleaned2,DF_DEM_cleaned3, DF_DEM_cleaned.c)

#analyse DEM

Res_DEM=cox_clocks(data_clocks, DF_DEM_cleaned)
Res_DEM_full=Res_DEM[[1]]
Res_DEM_60=Res_DEM[[2]]
Res_DEM_adj=Res_DEM[[3]]
Res_DEM_adj=Res_DEM[[3]]
Res_DEM_imp=Res_DEM[[4]]

print(Res_DEM_full[order(Res_DEM_full$pval),])
print(Res_DEM_60[order(Res_DEM_60$pval),])
print(Res_DEM_adj[order(Res_DEM_adj$pval),])
print(Res_DEM_imp[order(Res_DEM_imp$pval),])


#prepare PD dataset

DF_PD_cleaned<-DF_ND[!is.na(DF_ND$analysis4), c("Idepic", "Sex","Age_Start", "Age_Start_cat_5","Date_Start", "Center", "Cvd_T2d_Coh", "analysis4_EFU", "analysis4",covars)]
DF_PD_cleaned$ageevent<-(as.numeric(DF_PD_cleaned$analysis4_EFU-DF_PD_cleaned$Date_Start))/365.25
DF_PD_cleaned$ageevent<-DF_PD_cleaned$ageevent+DF_PD_cleaned$Age_Start
names(DF_PD_cleaned)<-c("idepic", "sex","age", "age_cat_5","Date_Start", "center", "cvd_t2d_coh", "analysis4_EFU", "indevent",covars, "ageevent")




##############Prentie weighting for main analysis#####

DF_PD_cleaned1<-DF_PD_cleaned[DF_PD_cleaned$cvd_t2d_coh==0 & DF_PD_cleaned$indevent==1,]
DF_PD_cleaned1$age <- DF_PD_cleaned1$ageevent - 1e-4
DF_PD_cleaned2<-DF_PD_cleaned[DF_PD_cleaned$cvd_t2d_coh==1 & DF_PD_cleaned$indevent==1,]
DF_PD_cleaned3<-DF_PD_cleaned[DF_PD_cleaned$cvd_t2d_coh==1 & DF_PD_cleaned$indevent==0,]
remove(DF_PD_cleaned)
DF_PD_cleaned=rbind(DF_PD_cleaned1,DF_PD_cleaned2,DF_PD_cleaned3)
remove(DF_PD_cleaned1,DF_PD_cleaned2,DF_PD_cleaned3)

########exclu pre 5 or 2 year for sensitivity only ####
#DF_PD_cleaned$fut<-DF_PD_cleaned$ageevent-DF_PD_cleaned$age
#DF_PD_cleaned$pre_5_year<-ifelse(DF_PD_cleaned$fut <= 5,1,0)
#DF_PD_cleaned$pre_2_year<-ifelse(DF_PD_cleaned$fut <= 2,1,0)

#DF_PD_cleaneda<-DF_PD_cleaned[DF_PD_cleaned$indevent==1,]
##DF_PD_cleaneda<-DF_PD_cleaneda[DF_PD_cleaneda$pre_5_year==0,]
#DF_PD_cleaneda<-DF_PD_cleaneda[DF_PD_cleaneda$pre_2_year==0,]
#DF_PD_cleanedb<-DF_PD_cleaned[DF_PD_cleaned$indevent==0,]
#DF_PD_cleaned.c=rbind(DF_PD_cleaneda,DF_PD_cleanedb)
#remove(DF_PD_cleaneda,DF_PD_cleanedb)

#DF_PD_cleaned1<-DF_PD_cleaned.c[DF_PD_cleaned.c$cvd_t2d_coh==0 & DF_PD_cleaned.c$indevent==1,]
#DF_PD_cleaned1$age <- DF_PD_cleaned1$ageevent - 1e-4#prentice weighting
#DF_PD_cleaned2<-DF_PD_cleaned.c[DF_PD_cleaned.c$cvd_t2d_coh==1 & DF_PD_cleaned.c$indevent==1,]
#DF_PD_cleaned3<-DF_PD_cleaned.c[DF_PD_cleaned.c$cvd_t2d_coh==1 & DF_PD_cleaned.c$indevent==0,]
#remove(DF_PD_cleaned)
#DF_PD_cleaned=rbind(DF_PD_cleaned1,DF_PD_cleaned2,DF_PD_cleaned3)
#remove(DF_PD_cleaned1,DF_PD_cleaned2,DF_PD_cleaned3, DF_PD_cleaned.c)

#analyse PD

Res_PD=cox_clocks(data_clocks, DF_PD_cleaned)
Res_PD_full=Res_PD[[1]]
Res_PD_60=Res_PD[[2]]
Res_PD_adj=Res_PD[[3]]
Res_PD_imp=Res_PD[[4]]

print(Res_PD_full[order(Res_PD_full$pval),])
print(Res_PD_60[order(Res_PD_60$pval),])
print(Res_PD_adj[order(Res_PD_adj$pval),])
print(Res_PD_imp[order(Res_PD_imp$pval),])

#############
#save(Res_PD_full,Res_PD_60,Res_PD_adj,Res_PD_imp,
#  Res_AD_full,Res_AD_60,Res_AD_adj,Res_AD_imp,
#  Res_DEM_full,Res_DEM_60, Res_DEM_adj,Res_DEM_imp,file="output/NDResultsv2.Rdata" )
save(Res_PD_full,Res_PD_60,Res_PD_adj,Res_PD_imp,
     Res_AD_full,Res_AD_60,Res_AD_adj,Res_AD_imp,
     Res_DEM_full,Res_DEM_60, Res_DEM_adj,Res_DEM_imp,file="output/NDResultsv2_QC_con2_pre2excl.Rdata" )
#file="output/NDResultsv2_QC_con2_neversmokers.Rdata" )
#file="output/NDResultsv2_QC_con2_pre5excl.Rdata" )
#file="output/NDResultsv2_QC_con2.Rdata"


##tables
###ND descriptive tables##############


form<-as.formula(~age+
                   sex+
                   center+
                   L_School+
                   Bmi_C+
                   Alc_Re+
                   Smoke_Stat+
                   Pa_Mets+
                   Pa_Index+
                   Hli_Dietscore+
                   Hli_Score+
                   Menopause+
                   Cvd_Prev+
                   T2D_Prev +
                   Death_Status+
                   Cvd_1_Status+
                   T2D_Status+
                   Cncr_Mal_Anyc+
                   Neuro_Status+
                   indevent+
                   FUT+
                   Fasting_C
                 |Country)

tabledb<-DF_AD_cleaned[DF_AD_cleaned$cvd_t2d_coh==1,]
tabledb<-tabledb[which(tabledb$idepic %in% data_clocks$idepic),]
tabledb<-to_factor(tabledb)
tabledb$indevent<-as.factor(tabledb$indevent)
tabledb$FUT<-tabledb$ageevent-tabledb$age
label(tabledb$FUT) <- "Follow up time (yrs)"
write.csv(
  x =table1(form, data=tabledb, render.continuous="Mean (SD)"),
  file = "/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/descriptive_tables/AD_comparison_cohort.csv",
  row.names = FALSE  )


tabledb<-DF_DEM_cleaned[DF_DEM_cleaned$cvd_t2d_coh==0,]
tabledb<-tabledb[which(tabledb$idepic %in% data_clocks$idepic),]
tabledb<-to_factor(tabledb)
tabledb$indevent<-as.factor(tabledb$indevent)
tabledb$FUT<-tabledb$ageevent-tabledb$age
label(tabledb$FUT) <- "Follow up time (yrs)"
write.csv(
  x =table1(form, data=tabledb, render.continuous="Mean (SD)"),
  file = "/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/descriptive_tables/Dem_sample.csv",
  row.names = FALSE  )


tabledb<-DF_PD_cleaned[DF_PD_cleaned$cvd_t2d_coh==0,]
tabledb<-tabledb[which(tabledb$idepic %in% data_clocks$idepic),]
tabledb<-to_factor(tabledb)
tabledb$indevent<-as.factor(tabledb$indevent)
tabledb$FUT<-tabledb$ageevent-tabledb$age
label(tabledb$FUT) <- "Follow up time (yrs)"
tabledb<-tabledb[!is.na(tabledb$FUT),]
write.csv(
  x =table1(form, data=tabledb, render.continuous="Mean (SD)"),
  file = "/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/descriptive_tables/PD_sample.csv",
  row.names = FALSE  )








save.image(file="ClockNDv2.Rdata")
