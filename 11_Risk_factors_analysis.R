library(haven)
library(tidyverse)
library(gtsummary) 
library(labelled)
library(table1)
library(pheatmap)

CrossSectional_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP <- readRDS("/data/Epic/subprojects/Somalogic/work/Vivian/Phase2/Data/ForAnalyses/CrossSectional_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP.rds")
db<-CrossSectional_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP[,c(1:220,7817:7824)]
remove(CrossSectional_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP)

db$hli_alcohol.f<-as.factor(db$hli_alcohol)
table(db$hli_bmic.f)
levels(db$hli_alcohol.f)<-c(">60 (g/d)", ">24-60 (g/d)", ">12-24 (g/d)", ">6-12 (g/d)", "0-6 (g/d)")
db$hli_dietscore_c.f<-as.factor(db$hli_dietscore_c)
levels(db$hli_dietscore_c.f)<-c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)")
db$hli_smoke.f<-as.factor(db$hli_smoke)
levels(db$hli_smoke.f)<-c("Current, 16+ cig/day", "Current, 1-15 cig/day","Former, quit <= 10 years","Former, quit 10+ years", "Never")
db$hli_bmic.f<-as.factor(db$hli_bmic)
levels(db$hli_bmic.f)<-c(">=30","26-<30","24-<26","22-<24","<22")
db$hli_pamets.f<-as.factor(db$hli_pamets)
levels(db$hli_pamets.f)<-c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)")
db$hli_score_c.f<-as.factor(db$hli_score_c )
table(db$hli_score_c.f)


summary(db$hli_score)
vQuin = quantile(db$hli_score, c(0:5/5), na.rm = T)
db$hli_score.q = cut(db$hli_score,  vQuin,  include.lowest = T,  labels = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"))
db$hli_score.q<-as.factor(db$hli_score.q)
table(db$hli_score.q)

db$hli_score.sc<-scale(db$hli_score)
db$alc_sc<-scale(db$alc_re)
db$bmi_sc<-scale(db$bmi_c)

db$currentsmoker<-ifelse(db$hli_smoke<=1,1,0)

db$pa_mets.sc<-scale(db$pa_mets)
db$hli_dietscore.sc<-scale(db$hli_dietscore)
db$highedu<-ifelse(db$l_school=="Longer education (incl. University deg.)",1,0)


######
setwd("/data/Epic/subprojects/Somalogic/work/Oliver")
data_clocks=read.csv("output/data_clocks_new.csv")



###########

clock_names=c("Oh","Organismal","Brain","Adipose", "Artery", "Immune","Heart","Intestine","Kidney","Liver","Lung","Muscle","Pancreas","Lehallier","Tanaka", "Wang", "Sathyan", "Global")




####

riskfactors=c("alc_sc", "bmi_sc" ,  "currentsmoker",  "pa_mets.sc" ,  "hli_dietscore.sc" , "highedu" )
covs<-c("age_blood","center", "female")
covs.adj<-c("age_blood","center", "female", "alc_sc", "bmi_sc" ,  "currentsmoker",  "pa_mets.sc" ,  "hli_dietscore.sc" , "highedu")


RF_results = list()
RF_results_adj = list()

for(riskfactor in riskfactors){
  
  RF_tib <- tibble(Clock = character(), beta = numeric(),  pval = numeric(),
                       SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric())
  RF_tib_adj <- tibble(Clock = character(), beta = numeric(),  pval = numeric(),
                   SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric())
  
  
  for (clock_name in clock_names){
    
  data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
  merged=merge(x=db[which(db$cvd_t2d_coh==1),],y=data_plot,by="idepic")
 
  form<-as.formula(paste0("AgeGap_zscored ~", riskfactor, "+", paste(covs, collapse = "+")))
  lmloop<-lm(form, data=merged)
  
  temptib_RF <- tibble(Clock= clock_name, beta = summary(lmloop)$coef[2, 'Estimate'],
                           pval = summary(lmloop)$coef[2, 'Pr(>|t|)'],
                           SE= summary(lmloop)$coef[2, 'Std. Error'] ,
                           LCL =confint(lmloop)[2,1],
                           UCL =confint(lmloop)[2,2],
                           n =nobs(lmloop) )

  RF_tib <- bind_rows(RF_tib,  temptib_RF)
  
  form<-as.formula(paste0("AgeGap_zscored ~", riskfactor, "+", paste(covs.adj, collapse = "+")))
  lmloop<-lm(form, data=merged)
  
  temptib_RF <- tibble(Clock= clock_name, beta = summary(lmloop)$coef[2, 'Estimate'],
                       pval = summary(lmloop)$coef[2, 'Pr(>|t|)'],
                       SE= summary(lmloop)$coef[2, 'Std. Error'] ,
                       LCL =confint(lmloop)[2,1],
                       UCL =confint(lmloop)[2,2],
                       n =nobs(lmloop) )
  
    RF_tib_adj<- bind_rows(RF_tib_adj,  temptib_RF)
  
  RF_results[[riskfactor]]<-RF_tib
  RF_results_adj[[riskfactor]]<-RF_tib_adj
  }
  }
  

############figures
conventclocks<-c("Tanaka","Lehallier", "Sathyan","Oh","Wang", "Global")
library(forestplot)

riskfactor="hli_score.sc"
plottable=RF_results[[riskfactor]]
plottable=plottable[plottable$Clock %in% conventclocks,]
plottable=plottable[match(conventclocks,plottable$Clock ),]

pdf(file=paste0("/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/figures/HLi_con_clocks.pdf"), height=10, width =12)
forestplot(plottable$Clock, #cbind(plottable$Clock, paste("p=",signif(plottable$pval,2))  )
           mean = plottable$beta,
           lower = plottable$LCL,
           upper = plottable$UCL,title=paste0("Healthy Lifestyle Index", " (n=", plottable$n[1],")"), lwd.ci=3,
           xlab="Agegap zscore per SD difference in Healthy Lifestyle Index",  zero=0,  , line.margin = .2, xticks=seq(-0.2,0.1,0.05),
           txt_gp =fpTxtGp( gpar(cex=1.5), xlab = gpar(cex=1.5), ticks  = gpar(cex=1.5)))
dev.off()

#######riskfactor heatmap
riskfactors=c("alc_sc", "bmi_sc" ,  "currentsmoker",  "pa_mets.sc" ,  "hli_dietscore.sc", "highedu")          


rfheatmap<- RF_results_adj[[riskfactors[1]]]
rfheatmap<-rfheatmap["beta"]
names(rfheatmap)<-riskfactors[1]

for(i in 2:length(riskfactors)){
  plottable<-RF_results_adj[[riskfactors[i]]]
  plottable<-plottable["beta"]
  names(plottable)<-riskfactors[i]
  rfheatmap<-cbind(rfheatmap,plottable)
}

rownames(rfheatmap)<-clock_names
names(rfheatmap)<-c("Alcohol consumption", "BMI" ,  "Current smoker",  "Physical activity (METs)" ,  "Healthy diet score" , "High education level") 
rfheatmap<- rfheatmap[conventclocks,]
rownames(rfheatmap)<-conventclocks

rfheatmappval<- RF_results_adj[[riskfactors[1]]]
rfheatmappval<-rfheatmappval["pval"]
names(rfheatmappval)<-riskfactors[1]

for(i in 2:length(riskfactors)){
  plottable<-RF_results_adj[[riskfactors[i]]]
  plottable<-plottable["pval"]
  names(plottable)<-riskfactors[i]
  rfheatmappval<-cbind(rfheatmappval,plottable)
}

rownames(rfheatmappval)<-clock_names
rfheatmappval<- rfheatmappval[conventclocks,]

test<-c(as.matrix(rfheatmappval))
FDR<-p.adjust(test, method="BH")
FDrcutoff<-cbind(test,FDR)
FDRval<-max(FDrcutoff[,1][FDrcutoff[,2]<0.05])

star<-matrix(ifelse(rfheatmappval < 0.05, "*", ""),nrow(rfheatmappval),ncol(rfheatmappval))
star[rfheatmappval < FDRval ]<- "**"


rg<-max(abs(rfheatmap))

pdf(file=paste0("/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/figures/risk_factor_heatmap_con_clocks.pdf"), height=6, width =8)
pheatmap(t(rfheatmap),  cluster_rows=F, cluster_cols=F,display_numbers = t(star),
         fontsize_number=15,
         breaks=seq(-rg,rg, length.out=50),  
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

####categorical risk factors#####
riskfactors.f=c("hli_dietscore_c.f", "hli_smoke.f", "hli_alcohol.f"  ,"hli_bmic.f" ,  "hli_pamets.f"  , "hli_score.q")          
covs<-c("age_blood","center", "female")

RF_results.f = list()

for(riskfactor in riskfactors.f){
  
  RF_tib <- tibble(Clock = character(), beta2 = numeric(),beta3 = numeric(),beta4 = numeric(),beta5 = numeric(), 
                   pval2 = numeric(),pval3 = numeric(),pval4 = numeric(),pval5 = numeric(),
                  LCL2 = numeric(),LCL3 = numeric(),LCL4 = numeric(),LCL5 = numeric(),
                  UCL2 = numeric(),UCL3 = numeric(),UCL4 = numeric(),UCL5 = numeric(),
                  n = numeric())
  
  for (clock_name in clock_names){
    
    data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
    merged=merge(x=db[which(db$cvd_t2d_coh==1),],y=data_plot,by="idepic")
    
    form<-as.formula(paste0("AgeGap_zscored ~", riskfactor, "+", paste(covs, collapse = "+")))
    lmloop<-lm(form, data=merged)
    
    temptib_RF <- tibble(Clock= clock_name, 
                         beta2 = summary(lmloop)$coef[2, 'Estimate'],
                         beta3 = summary(lmloop)$coef[3, 'Estimate'],
                         beta4 = summary(lmloop)$coef[4, 'Estimate'],
                         beta5 = summary(lmloop)$coef[5, 'Estimate'],
                         pval2 = summary(lmloop)$coef[2, 'Pr(>|t|)'],
                         pval3 = summary(lmloop)$coef[3, 'Pr(>|t|)'],
                         pval4 = summary(lmloop)$coef[4, 'Pr(>|t|)'],
                         pval5 = summary(lmloop)$coef[5, 'Pr(>|t|)'],
                         LCL2 =confint(lmloop)[2,1],
                         LCL3 =confint(lmloop)[3,1],
                         LCL4 =confint(lmloop)[4,1],
                         LCL5 =confint(lmloop)[5,1],
                         UCL2 =confint(lmloop)[2,2],
                         UCL3 =confint(lmloop)[3,2],
                         UCL4 =confint(lmloop)[4,2],
                         UCL5 =confint(lmloop)[5,2],
                         n =nobs(lmloop) )
    
    RF_tib <- bind_rows(RF_tib,  temptib_RF)
    
    RF_results.f[[riskfactor]]<-RF_tib
  }
}

#cat forest plot
par(mfrow=c(2,3))
riskfactor=riskfactors.f[4]
plottable=RF_results.f[[riskfactor]]
plottable=plottable[plottable$Clock %in% conventclocks,]


forestplot(c("Oh",conventclocks[2:5], "Consensus"), 
           mean = cbind(0,plottable$beta2,plottable$beta3,plottable$beta4,plottable$beta5),
           lower = cbind(0,plottable$LCL2,plottable$LCL3,plottable$LCL4,plottable$LCL5),
           upper = cbind(0,plottable$UCL2,plottable$UCL3,plottable$UCL4,plottable$UCL5),
           title=paste0(riskfactor, " (n=", plottable$n[1],")"),
           xlab="agegap_zscore per RF",  zero=0, boxsize = 0.1,line.margin = .2,
           txt_gp =fpTxtGp(ticks = gpar(cex=1), xlab  = gpar(cex=1)))

#by clock
clockname="Global"
riskfactornames= c("Healthy Diet Score", "Smoking Status" , "Alcohol Consumption","Body Mass Index", "Physical activity","Healthy Lifestyle Index", "Level of Schooling","Smoking Status" , "Alcohol Consumption","Body Mass Index")
i=6 
riskfactor=riskfactors.f[i]
plottable=RF_results.f[[riskfactor]]
plottable=plottable[plottable$Clock == clockname,]

riskfactor
catnames<-levels(db$hli_score.q)#change manually for each plot

pdf(file=paste0("/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/figures/",riskfactor, "_quintiles_Global.pdf"), height=10, width =12)
forestplot(cbind(catnames, paste("n=",as.numeric(table(merged[riskfactors.f[i]])))),
           mean = rbind(0,plottable$beta2,plottable$beta3,plottable$beta4,plottable$beta5),
           lower = rbind(0,plottable$LCL2,plottable$LCL3,plottable$LCL4,plottable$LCL5),
           upper = rbind(0,plottable$UCL2,plottable$UCL3,plottable$UCL4,plottable$UCL5),
           title=riskfactornames[i],lwd.ci=3,
           xlab=paste("Difference in",clockname, "age gap zscore"),  zero=0, ,xticks=seq(-0.5,0.2,0.1),
           txt_gp =fpTxtGp( gpar(cex=1.5), xlab = gpar(cex=1.5), ticks  = gpar(cex=1.5)))
dev.off()

###subcohort table 1##############
names(db)
table(db$cntr_f)
db$cntr_f<-droplevels(db$cntr_f)
table(db$country)
db$country<-droplevels(db$country)
db$cvd_fut<-db$age_exit_cvd_1-db$age_recr
summary(db$cvd_fut)
label(db$cvd_fut) <- "Follow up time CVD (yrs)"
db$death_fut<-db$age_exit_death-db$age_recr
summary(db$death_fut)
label(db$death_fut) <- "Follow up time deaths (yrs)"
db$cancer_fut<-db$age_exit_cancer_1st-db$age_recr
summary(db$cancer_fut)
label(db$cancer_fut) <- "Follow up time cancer (yrs)"
db$t2d_fut<-db$age_exit_t2d-db$age_recr
summary(db$t2d_fut)
label(db$t2d_fut) <- "Follow up time type 2 diabetes (yrs)"
db$l_school[db$l_school=="Not specified"]<-NA
db$l_school<-droplevels(db$l_school)
db$smoke_stat<-droplevels(db$smoke_stat)
db$death_status <- factor(db$death_status, levels=c(0,1),labels=c("No", "Yes"))
db$t2d_status <- factor(db$t2d_status, levels=c(0,1),labels=c("No", "Yes"))

table(db$t2d_status)

####

tabledb<-db[db$cvd_t2d_coh==0,]
tabledb<-tabledb[which(tabledb$idepic %in% data_clocks$idepic),]
tabledb<-to_factor(tabledb)
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
       chd_1_status+
       cvd_fut+
       t2d_status+
       t2d_fut+
       cncr_mal_anyc+
       cancer_fut+
       neuro_status+
       fasting_c
       |country, data=tabledb, render.continuous="Mean (SD)"),
       file = "/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/descriptive_tables/main_subcohort.csv",
       row.names = FALSE  )




save(RF_results_adj,RF_results.f,RF_results, file="output/for_export/summary_results/summary_results_for_riskfactors.Rdata" )
save.image(file="cx_riskfactors.Rdata")
