library(tidyverse)
library( haven )
library(gtsummary) 
library(labelled)
install.packages("mice")
library(mice)

setwd("/data/Epic/subprojects/Somalogic/work/Oliver")
data_clocks=read.csv("output/data_clocks_new.csv")
data_clocks=data_clocks %>% mutate(Sex_F=as.factor(Sex_F))
clock_names=c("Conventional","Organismal","Brain","Adipose", "Artery", "Immune","Heart","Intestine","Kidney","Liver","Lung","Muscle","Pancreas","Lehallier","Tanaka", "Wang", "Sathyan", "Consensus")
table(data_clocks$Organ)

setwd("/data/Epic/subprojects/Somalogic/")
mm0  <- read_sas(data_file= "sources/Epi_Data/somalogic_2023.sas7bdat",
                 catalog_file="sources/Epi_Data/formats.sas7bcat")

overlap<-intersect(data_clocks$idepic, mm0$Idepic)
names(mm0)
mm0<-mm0[mm0$Idepic %in% overlap,]

###imputation
covs<-c( "l_school", "bmi_c" ,  "alc_re",  "pa_mets" , "hli_dietscore_c" ,   "hli_smoke.f")


table(mm0$L_School)#5=NA
mm0$L_School[mm0$L_School==5]<-NA
str(mm0$L_School)
mm0$L_School<-as_factor(mm0$L_School)

summary(mm0$Bmi_C)#5=NA
str(mm0$Bmi_C)

summary(mm0$Alc_Re)#26=NA
str(mm0$Hli_Alcohol)
summary(mm0$Hli_Alcohol)###434 NAs

summary(mm0$Pa_Mets)#313=NA
is.factor(mm0$Hli_Alcohol)

summary(mm0$Hli_Dietscore_C)#434
mm0$Hli_Dietscore_C<-as.factor(mm0$Hli_Dietscore_C)
summary(mm0$MRMed_Score)#26
boxplot(mm0$MRMed_Score~mm0$Hli_Dietscore_C)
mm0$MRMed_Score_C<-as.factor(mm0$MRMed_Score_C)

summary(mm0$Hli_Smoke)#571
table(mm0$Hli_Smoke)#571
table(mm0$Smoke_Stat)#103 "4=NA"
mm0$Smoke_Stat[mm0$Smoke_Stat==4]<-NA
mm0$Smoke_Stat<-as.factor(mm0$Smoke_Stat)
summary(mm0$Hli_Score)#791

table(mm0$Sex)
mm0$Sex<-as.factor(mm0$Sex)
mm0$Cvd_T2d_Coh<-as.factor(mm0$Cvd_T2d_Coh)
tapply(mm0$Hli_Score,mm0$Cvd_T2d_Coh, summary )

tapply(mm0$Hli_Score,mm0$Cvd_T2d_Coh, summary )

mm0$Fasting_C<-as.factor(mm0$Fasting_C)

####
mm0$idepic<-mm0$Idepic
data_plot=data_clocks[which(data_clocks$Organ=="Consensus"),]
merged=merge(x=mm0[which(mm0$Cvd_T2d_Coh==1),],y=data_plot,by="idepic")

mod1<-lm( AgeGap_zscored ~ Age_Blood+ Center+ Sex +MRMed_Score_C+L_School, data=merged)
summary(mod1)

####


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

imputedcovarites$Alc_Re[1:10]
mm0$Alc_Re[1:10]

names(imputedcovarites)<-paste0(names(imputedcovarites),".imp")
imputedcovarites$Idepic<-mm0$Idepic

save(imputedcovarites, file="/data/Epic/subprojects/Somalogic/work/Oliver/output/imputedcovariates.Rdata" )
