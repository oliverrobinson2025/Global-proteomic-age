library(tidyverse)
library( haven )
library(gtsummary) 
library(labelled)
install.packages("mice")
library(mice)

setwd("/data/Epic/subprojects/Somalogic/work/Oliver")
data_clocks=read.csv("output/data_clocks_new.csv")
table(data_clocks$Organ)

setwd("/data/Epic/subprojects/Somalogic/")
mm0  <- read_sas(data_file= "sources/Epi_Data/somalogic_2023.sas7bdat",
                 catalog_file="sources/Epi_Data/formats.sas7bcat") ### read covariate database in R

overlap<-intersect(data_clocks$idepic, mm0$Idepic)
names(mm0)
mm0<-mm0[mm0$Idepic %in% overlap,] ###match data with subjects with proteomic age

###imputation#####
#format covariates of interest
covs<-c( "l_school", "bmi_c" ,  "alc_re",  "pa_mets" , "hli_dietscore_c" ,   "hli_smoke.f")

mm0$L_School[mm0$L_School==5]<-NA
mm0$L_School<-as_factor(mm0$L_School)

summary(mm0$Alc_Re)
str(mm0$Hli_Alcohol)
is.factor(mm0$Hli_Alcohol)

summary(mm0$Hli_Dietscore_C)
mm0$Hli_Dietscore_C<-as.factor(mm0$Hli_Dietscore_C)
summary(mm0$MRMed_Score)
mm0$MRMed_Score_C<-as.factor(mm0$MRMed_Score_C)

mm0$Smoke_Stat[mm0$Smoke_Stat==4]<-NA
mm0$Smoke_Stat<-as.factor(mm0$Smoke_Stat)

table(mm0$Sex)
mm0$Sex<-as.factor(mm0$Sex)
mm0$Cvd_T2d_Coh<-as.factor(mm0$Cvd_T2d_Coh)

####run MICE###

vars<-c("Age_Blood","Pa_Mets","Alc_Re", "Bmi_C", "Center", "Cvd_T2d_Coh", "Fasting_C", "Hli_Dietscore_C", "Smoke_Stat", "Hli_Score","L_School", "Sex", "MRMed_Score" )

toimpute<-mm0[vars]
md.pattern(toimpute, rotate.names = T)
Napattern<-toimpute
Napattern<-ifelse(is.na(toimpute),1,0)
summary(table(toimpute$Center, Napattern[,6]))

imputedEport <- mice(toimpute, m=2, maxit = 50, method = 'cart', seed = 500)
completeimportEport <- mice::complete(imputedEport,1)
md.pattern(completeimportEport, rotate.names = T)

imputedcovarites<-completeimportEport[c("Alc_Re",  "Center",  "Fasting_C", "Hli_Dietscore_C", "Smoke_Stat", "Hli_Score","L_School", "MRMed_Score","Pa_Mets")]

names(imputedcovarites)<-paste0(names(imputedcovarites),".imp")
imputedcovarites$Idepic<-mm0$Idepic

save(imputedcovarites, file="/data/Epic/subprojects/Somalogic/work/Oliver/output/imputedcovariates.Rdata" )
