rm(list=ls())
setwd("/data/Epic/subprojects/Somalogic/work/Oliver")

library(tidyverse)
devtools::install_github("SomaLogic/SomaDataIO@*release")
library(SomaDataIO)
library( haven )  # To load SAS sas7bdat file with variable labels too  
library( lubridate ) 
library( Epi )
library( labelled ) 
library(survival)
library(htmlwidgets)
library(DT)
library(forcats)
remotes::install_github("privefl/bigutilsr")
library(bigutilsr) ### for the detection of outliers based on local outlier factor and Tukey's rule (see https://privefl.github.io/blog/detecting-outlier-samples-in-pca/)


####import co variate data ####
pathWD <- "work/Vivian/MultiEndpoints"
is_seq <- function(.x) grepl("^seq\\.[0-9]{4}", .x) # regex for analytes


mm0       <- read_sas(data_file= "/data/Epic/subprojects/Somalogic/sources/Epi_Data/somalogic_2023.sas7bdat",
                      catalog_file="/data/Epic/subprojects/Somalogic/sources/Epi_Data/formats.sas7bcat")

##### import somalogic data ####
suff0      <-".hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP"
suff1      <-".hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP"# ".hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck" #   

soma_adat0 <- read_adat(paste0("/data/Epic/subprojects/Depot_Somalogic/sources/Raw_Data/SS-2230719.SS-2215915.SS-2230718.SS-2230720_v4.1_CitratePlasma", suff0, ".adat"))
soma_adat1 <- read_adat(paste0("/data/Epic/subprojects/Depot_Somalogic/sources/Raw_Data/SS-2342103_SS-2342071_SS-2340734_v4.1_CitratePlasma", suff1, ".adat"))

soma_adat0$batch=1
soma_adat1$batch=2
soma_adat  <- bind_rows(soma_adat0, soma_adat1)
soma_red   <- soma_adat %>% select(1:33,7630)
soma_red$ID=rownames(soma_red)

rm(soma_adat0, soma_adat1)

soma_all <- inner_join(soma_red  %>% mutate(tojoin=SubjectID), mm0 %>% mutate(tojoin=Idepic_Aliq), by="tojoin")
names(soma_all) <- tolower(names(soma_all))

table(soma_all$sampletype)

######## display #####


soma_adat[1:5,1:5]


#### Exclusion of "flagged" samples (those for which the normalization scaling factors are not within the acceptable range) #####
#### and exclusion of the Somalogic inhouse QCs, Cals, and EPIC QCs (QC#A1)
FlaggedToEliminate        <- soma_adat %>% filter(RowCheck =="FLAG") %>% select(SubjectID) %>% pull 
soma_adat                 <- soma_adat %>% filter(RowCheck !="FLAG")


#### Exclusion of outliers, based on PCA and LOF (https://privefl.github.io/blog/detecting-outlier-samples-in-pca/)
peakTable                <- data.frame(soma_adat  %>% select_if(grepl("seq",names(.))))
pca                      <- prcomp(peakTable, scale. = TRUE, rank. = 10)
U                        <- pca$x
llof                     <- LOF(U)
OutliersToEliminate      <- soma_adat %>% slice(which(llof > tukey_mc_up(llof))) %>% select(SubjectID) %>% pull

soma_adat <- soma_adat %>% slice(-which(llof > tukey_mc_up(llof)))

Eliminated               <- list(Flagged  = FlaggedToEliminate, ##no SMP24 24 SMP 247 247
                                 Outliers = OutliersToEliminate) ##no SMP15 16 SMP 9 9

#####prepare protein and covariate data as specified for organ age package####

overlap<-intersect(mm0$Idepic_Aliq, soma_adat$SubjectID)

soma_adat2<-soma_adat[soma_adat$SubjectID %in% overlap,c(23,34:7629)]
colnames(soma_adat2)[c(1:5,7590:7597)]

new_names=substring(colnames(soma_adat2)[2:7597],5)
new_names2=gsub("[.]","-",new_names)

colnames(soma_adat2)=c("ID",new_names2)
summary(soma_adat2$`10000-28`)

mm1<-mm0[mm0$Idepic_Aliq %in% overlap,]

id_order<-match(mm1$Idepic_Aliq, soma_adat2$ID)
soma_adat2<-soma_adat2[id_order,]

mm1$Idepic_Aliq[1:10]
soma_adat2$ID[1:10]


write.csv(soma_adat2,"Input_prot_clocks_all.csv",row.names = F)

names(mm1)
names(soma_all)
soma_all<-soma_all[id_order,]
soma_all$subjectid[1:10]

table(soma_all$Sex)
print_labels(soma_all$sex)
l=lapply(soma_all, is.labelled)
l

soma_all %>% 
  select(idepic,subjectid,age_blood,sex) %>% 
  mutate(Sex_F=ifelse(sex==2,1,0)) -> Input_prot_clocks_age_sex

colnames(Input_prot_clocks_age_sex)=c("idepic","ID","Age","sex","Sex_F")

table(Input_prot_clocks_age_sex$Sex_F)

write.csv(Input_prot_clocks_age_sex,"Input_prot_clocks_all_age_sex.csv",row.names = F)

##prepare lifted 5K protein data for other clocks####
####lift version to 5k version##
is_intact_attr(soma_adat)
attr(soma_adat, "Header.Meta")$HEADER$StudyMatrix |> as.character()
from_space <- getSignalSpace(soma_adat)
from_space
is.null(checkSomaScanVersion(from_space))
soma_adat_5k <- lift_adat(soma_adat, bridge = "7k_to_5k")
is_lifted(soma_adat_5k) 
is.soma_adat(soma_adat_5k)
getSignalSpace(soma_adat_5k)
getSomaScanVersion(soma_adat_5k)



####
overlap<-intersect(mm0$Idepic_Aliq, soma_adat_5k$SubjectID)
soma_adat_5k
soma_adat2<-soma_adat[soma_adat_5k$SubjectID %in% overlap,c(23,34:7629)]
colnames(soma_adat2)[c(1:5,7590:7597)]

new_names=substring(colnames(soma_adat2)[2:7597],5)
new_names2=gsub("[.]","-",new_names)

colnames(soma_adat2)=c("ID",new_names2)
summary(soma_adat2$`10000-28`)

mm1<-mm0[mm0$Idepic_Aliq %in% overlap,]

id_order<-match(mm1$Idepic_Aliq, soma_adat2$ID)
soma_adat2<-soma_adat2[id_order,]

mm1$Idepic_Aliq[1:10]
soma_adat2$ID[1:10]


write.csv(soma_adat2,"Input_prot_clocks_all_5K.csv",row.names = F)


