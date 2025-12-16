library(forestplot)
library(haven)
library(tidyverse)
library(pheatmap)
library(grid)
library(forestploter)
###import results#####

load("/data/Epic/subprojects/Somalogic/work/Oliver/output/MortalityResults_QC_v2_con2.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/CVDResults_QC_v2_con2.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/T2DResults_no_prev_t2d_QC_v2_con2.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/ClockCancerResults_QC_v2_con2.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/NDResultsv2_QC_con2.Rdata")

all_results<-cancer_results
all_results[["AD"]]<-Res_AD_full
all_results[["PD"]]<-Res_PD_full
all_results[["ALS"]]<-Res_ALS_full
all_results[["DEM"]]<-Res_DEM_full
all_results[["cvd"]]<-Res_cvd
all_results[["chd"]]<-Res_chd
all_results[["stroke"]]<-Res_stroke
all_results[["t2d"]]<-Res_t2d
all_results[["deaths"]]<-Res_deaths

all_results_adj<-cancer_results_adj
all_results_adj[["AD"]]<-Res_AD_adj
all_results_adj[["PD"]]<-Res_PD_adj
all_results_adj[["ALS"]]<-Res_ALS_adj
all_results_adj[["DEM"]]<-Res_DEM_adj
all_results_adj[["cvd"]]<-Res_cvd_adj
all_results_adj[["chd"]]<-Res_chd_adj
all_results_adj[["stroke"]]<-Res_stroke_adj
all_results_adj[["t2d"]]<-Res_t2d_adj
all_results_adj[["deaths"]]<-Res_deaths_adj

all_results_imp<-cancer_results_imp
all_results_imp[["AD"]]<-Res_AD_imp
all_results_imp[["PD"]]<-Res_PD_imp
all_results_imp[["ALS"]]<-Res_ALS_imp
all_results_imp[["DEM"]]<-Res_DEM_imp
all_results_imp[["cvd"]]<-Res_cvd_imp
all_results_imp[["chd"]]<-Res_chd_imp
all_results_imp[["stroke"]]<-Res_stroke_imp
all_results_imp[["t2d"]]<-Res_t2d_imp
all_results_imp[["deaths"]]<-Res_deaths_imp

###make Results by clock ####
clock_names=c("Oh","Organismal","Brain","Adipose", "Artery", "Immune","Heart","Intestine","Kidney","Liver","Lung","Muscle","Pancreas","Lehallier","Tanaka", "Wang", "Sathyan", "Global")

outcomenames<-c("All cause dementia", "Alzheimer's disease","Parkinson's disease","Amyotrophic lateral sclerosis","Cardiovascular Disease","Coronary heart disease","Stroke","Type 2 diabetes","All-cause mortality","Any cancer","Bladder cancer","Breast cancer","Peri-menopausal breast cancer","Post-menopausal breast cancer","Pre-menopausal breast cancer","Colon-rectum cancer","Colon cancer","Rectum cancer","Endometrial cancer","Glioma cancer","Kidney cancer","Liver cancer","Lung cancer","Lymphoma","Melanoma cancer","Ovary cancer","Pancreas cancer","Prostate cancer","Stomach cancer","Thyroid cancer","Upper aero-digestive tract cancer", "UADT (non-Oesophageal) ","UADT (Oesophageal) Cancer", "Non-Hodgkins Lymphoma", "Hodgkins Lymphoma")
outcomes=c("DEM", "AD","PD","ALS","cvd","chd","stroke", "t2d","deaths","cncr_mal_anyc","cncr_mal_blad", "cncr_mal_brea","cncr_mal_brea_peri","cncr_mal_brea_post","cncr_mal_brea_pre", "cncr_mal_clrt", "cncr_mal_clrt_colon", "cncr_mal_clrt_rectum", "cncr_mal_coru", "cncr_mal_glio", "cncr_mal_kidn", "cncr_mal_live", "cncr_mal_lung", "cncr_mal_lymp", "cncr_mal_mela", "cncr_mal_ovar", "cncr_mal_panc", "cncr_mal_pros", "cncr_mal_stom", "cncr_mal_thyr", "cncr_mal_uadt","cncr_mal_orophag","cncr_mal_esophag", "cncr_mal_nh_lymp", "cncr_mal_h_lymp")

all_clock_results = list()

for (clock_name in clock_names){
  
  Res_clock <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric(), outcome= character())
  
  for(outcome in outcomes){
    
    plottable=all_results[[outcome]]
    plottable$outcome<-outcome
    Res_clock <-rbind(Res_clock, plottable[plottable$Clock==clock_name,])
  }
  
  Res_clock$outcomenames<-outcomenames
  Res_clock$FDR<-p.adjust(Res_clock$pval, method = "BH")
  
  all_clock_results[[clock_name]]<-Res_clock
  
}

Conventional<-all_clock_results[["Oh"]] ### add results for Oh clock that are reffered to as "Conventional" for organ age presnetation
Conventional$Clock<-"Conventional"
all_clock_results[["Conventional"]]<-Conventional

##adjusted all clock results

all_clock_results_adj = list()

for (clock_name in clock_names){
  
  Res_clock <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric(), outcome= character())
  
  for(outcome in outcomes){
    
    plottable=all_results_adj[[outcome]]
    plottable$outcome<-outcome
    Res_clock <-rbind(Res_clock, plottable[plottable$Clock==clock_name,])
  }
  
  Res_clock$outcomenames<-outcomenames
  Res_clock$FDR<-p.adjust(Res_clock$pval, method = "BH")
  
  all_clock_results_adj[[clock_name]]<-Res_clock
  
}

Conventional<-all_clock_results_adj[["Oh"]]
Conventional$Clock<-"Conventional"
all_clock_results_adj[["Conventional"]]<-Conventional

#imputed all clcok results

all_clock_results_imp = list()

for (clock_name in clock_names){
  
  Res_clock <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric(), outcome= character())
  
  for(outcome in outcomes){
    
    plottable=all_results_imp[[outcome]]
    plottable$outcome<-outcome
    Res_clock <-rbind(Res_clock, plottable[plottable$Clock==clock_name,])
  }
  
  Res_clock$outcomenames<-outcomenames
  Res_clock$FDR<-p.adjust(Res_clock$pval, method = "BH")
  
  all_clock_results_imp[[clock_name]]<-Res_clock
  
}

Conventional<-all_clock_results_imp[["Oh"]]
Conventional$Clock<-"Conventional"
all_clock_results_imp[["Conventional"]]<-Conventional
remove(Conventional)

####
####clock performance########
setwd("/data/Epic/subprojects/Somalogic/work/Oliver")

data_clocks=read.csv("output/data_clocks_new.csv")

table(data_clocks$Organ)


clockperformance<-matrix(NA, length(clock_names), 8)
clockperformance[,1]<-clock_names

for (i in 1:length(clock_names)){
  
  data_plot=data_clocks[which(data_clocks$Organ==clock_names[i]),]
  gg_test=cor.test(data_plot$Age,  data_plot$Predicted_Age)
  clockperformance[i,2]<-round(gg_test$estimate,2)
  clockperformance[i,3]<-signif(gg_test$p.value,2)
  clockperformance[i,4]<-round(mean(abs(data_plot$Predicted_Age- data_plot$Age)),2)
  clockperformance[i,5]<-round(mean(data_plot$AgeGap),2)
  clockperformance[i,6]<-round(sd(data_plot$AgeGap),2)
  clockperformance[i,7]<-length(data_plot$AgeGap_zscored[abs(data_plot$AgeGap_zscored)>5])
  clockperformance[i,8]<-as.numeric(gg_test$estimate)
}

clockperformance<-as.data.frame(clockperformance)
names(clockperformance)<-c("Clock", "r", "pval", "MAE", "mean(Agegap)", "sd(Agegap)", "N_outliers(5SD)", "r_num")


for (clock_name in clock_names){
  
  
  data_plot=data_clocks[which(data_clocks$Organ==clock_name),]
  #gg_min=floor(min(c(data_plot$Age,data_plot$Predicted_Age)))
  #gg_max=ceiling(max(c(data_plot$Age,data_plot$Predicted_Age)))
  gg_min=15
  gg_max=95
  
  gg_test=cor.test(data_plot$Age,
                   data_plot$Predicted_Age)
  
  mae<-mean(abs(data_plot$Age-data_plot$Predicted_Age))
  
  p=ggplot(data_plot, aes(x = Age, y = Predicted_Age, color=AgeGap_zscored)) +
    geom_point() +
    lims(x=c(gg_min,gg_max), y=c(gg_min,gg_max)) +
    ggtitle(paste0(clock_name," clock (r=",round(gg_test$estimate,2),", MAE=",signif(mae,2),")")) +
    xlab("Chronological age") +
    ylab("Predicted age") +
    geom_abline(intercept = 0, colour="grey") +
    scale_colour_gradient2()
  
  pdf(file=paste0("/data/Epic/subprojects/Somalogic/work/Oliver/output/figures/",clock_name, "_correlation_plot.pdf"), height=4, width =6)
  
    print(p)
    
    dev.off()
}



########################clock cor heat map#####

AAmatrix<-as.data.frame(cbind(data_clocks$AgeGap_zscored[data_clocks$Organ=="Tanaka"],
                              data_clocks$AgeGap_zscored[data_clocks$Organ=="Oh"],
                              data_clocks$AgeGap_zscored[data_clocks$Organ=="Sathyan"],
                              data_clocks$AgeGap_zscored[data_clocks$Organ=="Lehallier"],
                              data_clocks$AgeGap_zscored[data_clocks$Organ=="Wang"]))
names(AAmatrix)<-c("Tanaka", "Oh","Sathyan", "Lehallier","Wang" )
corAA<-cor(AAmatrix)


pdf(file=paste0("/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/figures/Con_clock_corrheatmap.pdf"), height=4, width =6)
pheatmap(corAA, display_numbers=T,na_col="white", cluster_rows=F, cluster_cols=F,border_color=NA,fontsize_number=15, fontsize = 15  )
dev.off()

###venn diagram clocks###
load(file="all_model_aptermers.Rdata")
install.packages("VennDiagram")
library(VennDiagram)

x=list(  Tanaka, Wang, Sathyan, Lehallier, Oh)



venn.plot <- venn.diagram(
  x ,
  filename = NULL,
  col = "black",
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  alpha = 0.50,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  margin = 0.05
)
grid.draw(venn.plot)



###make heatmap of conclock-disease asociations####
outcomenames.sh<-c("All cause dementia", "Alzheimer's disease","Parkinson's disease","Coronary heart disease","Stroke","Type 2 diabetes","All-cause mortality","Any cancer","Bladder cancer","Breast cancer","Colon cancer","Rectum cancer","Endometrial cancer","Glioma cancer","Kidney cancer","Liver cancer","Lung cancer","Lymphoma","Melanoma cancer","Ovary cancer","Pancreas cancer","Prostate cancer","Stomach cancer","Thyroid cancer","Upper aero-digestive tract cancer")
con_clock_names=c("Tanaka","Lehallier", "Sathyan","Oh","Wang", "Global")

clockheatmap<- as.data.frame(all_clock_results_imp[[con_clock_names[1]]])
rownames(clockheatmap)<-clockheatmap$outcomenames
clockheatmap<-clockheatmap["beta"]
names(clockheatmap)<-con_clock_names[1]

for(i in 2:length(con_clock_names)){
  plottable<-all_clock_results_imp[[con_clock_names[i]]]
  plottable<-plottable["beta"]
  names(plottable)<-con_clock_names[i]
  clockheatmap<-cbind(clockheatmap,plottable)
}



clockheatmap<-clockheatmap[outcomenames.sh,]
names(clockheatmap) <-con_clock_names

clockheatpval<- as.data.frame(all_clock_results_imp[[con_clock_names[1]]])
rownames(clockheatpval)<-clockheatpval$outcomenames
clockheatpval<-clockheatpval["pval"]
names(clockheatpval)<-con_clock_names[1]

for(i in 2:length(con_clock_names)){
  plottable<-all_clock_results_imp[[con_clock_names[i]]]
  plottable<-plottable["pval"]
  names(plottable)<-con_clock_names[i]
  clockheatpval<-cbind(clockheatpval,plottable)
}

clockheatpval<-clockheatpval[outcomenames.sh,]

test<-c(as.matrix(clockheatpval))
FDR<-p.adjust(test, method="BH")
FDrcutoff<-cbind(test,FDR)
FDRval<-max(FDrcutoff[,1][FDrcutoff[,2]<0.05])

star<-matrix(ifelse(clockheatpval < 0.05, "*", ""),nrow(clockheatpval),ncol(clockheatpval))
star[clockheatpval < FDRval ]<- "**"


rg<-max(abs(clockheatmap))

pdf(file=paste0("/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/figures/Con_clock_heatmap_plot.pdf"), height=6, width =10)
pheatmap((clockheatmap),  cluster_rows=T, cluster_cols=T,display_numbers = (star),
         fontsize_number=10,
         breaks=seq(-rg,rg, length.out=50),  
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

hitmatix<-matrix(ifelse(clockheatpval < FDRval & clockheatmap >0, 1, 0),nrow(clockheatpval),ncol(clockheatpval))


##convention clock performance
Conventional_clockperformance<-clockperformance[clockperformance$Clock %in% con_clock_names,]

hits<-rep(NA,length(con_clock_names))
for (i in 1:length(con_clock_names)){hits[i]<-sum(hitmatix[,i])}
diseasehts<-cbind(c("Tanaka","Lehallier", "Sathyan","Oh","Wang", "Consensus"), hits)
par(mar= c(5.1, 8.1, 4.1, 2.1))
barplot(as.numeric(diseasehts[,2]),
        horiz=T,xlab="Number of positive sigificant disease associations" ,
        names.arg=diseasehts[,1], las=2 )
par(mar= c(5.1, 4.1, 4.1, 2.1))


###figure by outcome conventional clocks and adjusted ####
# Set-up theme
tm <- forest_theme(base_size = 10,
                   refline_gp = gpar("solid"),
                   ci_pch = c(15, 18),#style for the CI
                   ci_col = c("blue", "red"),#color for the CI
                   ci_lwd = 1.6, #size for the CI line
                   title_gp = gpar(cex=1,fontface="bold"),
                   legend_name = "",
                   legend_value = c("Base model", "Risk factor adjusted"),
                   legend_gp = gpar(fontsize = 8.5,cex = 1.1,
                                    lwd=1.8),
                   legend_position = "bottom",
                   core=list(bg_params=list(fill=c("white"))),#set the color of background
                   vertline_lty = c("dashed", "dotted"),
                   vertline_col = c("#d6604d", "#bababa"))




####

for(i in 1:length(outcomes)){

plottableA=all_results[[outcomes[i]]]
plottableA=plottableA[plottableA$Clock %in% con_clock_names,]
plottableA=plottableA[match(con_clock_names,plottableA$Clock ),]
plottableB=all_results_imp[[outcomes[i]]]
plottableB=plottableB[plottableB$Clock %in% con_clock_names,]
plottableB=plottableB[match(con_clock_names,plottableB$Clock ),]
plottable=cbind(plottableA, plottableB)
names(plottable)[11:20]<-paste0(names(plottable)[1:10], "_adj")
plottable<-plottable[,c(-10,-20)]

plottable<-plottable %>% mutate(
  Asso = as.numeric(Asso),
  LCL = as.numeric(LCL),
  UCL = as.numeric(UCL),
  Asso_adj = as.numeric(Asso_adj),
  LCL_adj = as.numeric(LCL_adj),
  UCL_adj = as.numeric(UCL_adj),
  
  #reserve 2 digits for HR
  `HR (95% CI)` = paste(sprintf("%.2f (%.2f, %.2f)", 
                                Asso, LCL, UCL),
                        sprintf("%.2f (%.2f, %.2f)", 
                                Asso_adj, LCL_adj, UCL_adj),
                        sep = "\n"),
  
  #reserve 2 digits for P value and also I use Scientific notation
  #If p-value is lower than the threshold, * is added
  pval = ifelse(pval<FDRval, ##FDRval from conclcoks
                paste0(sprintf("%.2e", pval),"*"), 
                sprintf("%.2e", pval)),
  pval_adj = ifelse(pval_adj<FDRval, 
                    paste0(sprintf("%.2e", pval_adj),"*"), 
                    sprintf("%.2e", pval_adj)),
  `P value` = paste(pval,
                    pval_adj,
                    sep = "\n"),
  Clock = paste(Clock," ",sep = "\n"))

plottable$` ` <- paste(rep(" ", 20), collapse = " ")##create blank space for forest plot
xaxis_high<-signif(max(append(plottable$UCL,plottable$UCL_adj)),3)##calculate the max limit for X-axis


p <- forest(plottable[,c(1, 19, 20, 21)],
            est = list(plottable$Asso,
                       plottable$Asso_adj),
            lower = list(plottable$LCL,
                         plottable$LCL_adj),
            upper = list(plottable$UCL,
                         plottable$UCL_adj),
            ci_column = 4,
            ref_line = 1,
            nudge_y = 0.4,
            sizes = 0.9,
            theme = tm,
            title = paste0(outcomenames[i], 
                           " (n=", plottable$n[1],
                           " nevent=",plottable$nevent[1],")"),
            xlab="Hazard ratio per age gap zscore",
            xlim = c(0.8,xaxis_high),
            ticks_at = seq(0.8,xaxis_high,0.2))

p

ggplot2::ggsave(filename = paste0("output/for_export/figures/",outcomenames[i],"_con_clock.png"), plot = p,
               dpi = 300,
                width = 7.5, height = 7.5, units = "in")

grid.newpage()
}

# figure 4 adjusted and base----------------------------------------------------------------
clock_name="Sathyan"

plottableA=all_clock_results[[clock_name]]
plottableA=plottableA[plottableA$outcomenames %in% outcomenames.sh,]
plottableB=all_clock_results_imp[[clock_name]]
plottableB=plottableB[plottableB$outcomenames %in% outcomenames.sh,]
plottable=cbind(plottableA, plottableB)
names(plottable)[13:24]<-paste0(names(plottable)[1:12], "_adj")
plottable$pval_fdr<-paste0("p= ", signif(plottable$pval,2), "*")## following lines redundant as FDRval from heatmap used instead
plottable$FDR<-p.adjust(plottable$pval, method = "BH")##
plottable$pval_fdr[plottable$FDR>0.05]<-paste0("p= ",signif(plottable$pval,2))[plottable$FDR>0.05]##
plottable$pval_fdr_adj<-paste0("p= ",signif(plottable$pval_adj,2), "*")##
plottable$FDR_adj<-p.adjust(plottable$pval_adj, method = "BH")##
plottable$pval_fdr_adj[plottable$FDR_adj>0.05]<-paste0("p= ",signif(plottable$pval_adj,2))[plottable$FDR_adj>0.05]##

plottable=plottable[order(-plottable$Asso),]


# Generate point estimation and 95% CI. 
#Paste two CIs together and separate by line break.
#add new variable of N/n_event
##load plottable first
plottable<-plottable %>% mutate(
  Asso = as.numeric(Asso),
  LCL = as.numeric(LCL),
  UCL = as.numeric(UCL),
  Asso_adj = as.numeric(Asso_adj),
  LCL_adj = as.numeric(LCL_adj),
  UCL_adj = as.numeric(UCL_adj),
  
  #reserve 2 digits for HR
  `HR (95% CI)` = paste(sprintf("%.2f (%.2f, %.2f)", 
                                Asso, LCL, UCL),
                        sprintf("%.2f (%.2f, %.2f)", 
                                Asso_adj, LCL_adj, UCL_adj),
                        sep = "\n"),
  
  #reserve 2 digits for P value and also I use Scientific notation
  #If p-value is lower than the threshold, * is added
  pval = ifelse(pval<FDRval, 
                paste0(sprintf("%.2e", pval),"*"), 
                sprintf("%.2e", pval)),
  pval_adj = ifelse(pval_adj<FDRval, 
                    paste0(sprintf("%.2e", pval_adj),"*"), 
                    sprintf("%.2e", pval_adj)),
  `P value` = paste(pval,
                    pval_adj,
                    sep = "\n"),
  Disease = paste(outcomenames," ",sep = "\n"),
  `N (N event)` = paste0(n," (",nevent,") ",sep = "\n"))
# `N/n_event` = paste0(n,"/",nevent," ",sep = "\n"))

plottable$` ` <- paste(rep(" ", nrow(plottable)), collapse = " ")##create blank space for forest plot
xaxis_high<-signif(max(append(plottable$UCL,plottable$UCL_adj)),3)##calculate the max limit for X-axis



p <- forest(plottable[,c(29, 30, 27, 28, 31)], #the columns should be outcome, n/n_event, HR, P, and blank column
            est = list(plottable$Asso,
                       plottable$Asso_adj),
            lower = list(plottable$LCL,
                         plottable$LCL_adj),
            upper = list(plottable$UCL,
                         plottable$UCL_adj),
            ci_column = 5,
            ref_line = 1,
            nudge_y = 0.4,
            sizes = 0.9,
            theme = tm,
            title = paste0(clock_name, " clock"),
            xlab="Hazard ratio per age gap zscore",
            xlim = c(0.5,xaxis_high),
            ticks_at = seq(0.5,xaxis_high,0.5))

#p

ggplot2::ggsave(filename = paste0("output/for_export/figures/",clock_name, "_all_outcomes.png"), plot = p,
               dpi = 500,
               width = 10, height = 15, units = "in")




########compare with CPRD#####
load("output/KuanCPRD.Rdata" )
Res_compare<-merge(all_clock_results_imp[["Consensus"]], KuanCPRD, by = "outcomenames")

table(Res_compare$cluster)


par(mar= c(5.1, 8.1, 4.1, 2.1))

pdf(file=paste0("/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/figures/GlobalKuanBoxplot.pdf"), height=8, width =8)
boxplot(as.numeric(Res_compare$beta)~Res_compare$cluster, col= c(1,3,4,5) ,ylab="Log hazard per Global proteomic age gap zscore in EPIC",
        xlab="Age cluster",cex.lab=1.5, cex.axis=1.5)
text(3.2,0.39,"p value = 0.002", pos=1, cex=1.6)
dev.off()

summary(aov(as.numeric(Res_compare$beta)~Res_compare$cluster))

t.test(as.numeric(Res_compare$beta[Res_compare$cluster==3|Res_compare$cluster==4])~Res_compare$cluster[Res_compare$cluster==3|Res_compare$cluster==4])

plot(density(as.numeric(log(Res_compare$Beta_G))))
str(Res_compare)
Res_compare$outcomenames

table(Res_compare$cluster)

cor.test(Res_compare$Asso, as.numeric(Res_compare$Beta_G), method= "spearman")

pos<-rep(4,22)
pos[c(4,13,22,9)]<-1
pos[c(2,10,6)]<-2

pdf(file=paste0("/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/figures/GlobalKuanScatterplot_verticalline.pdf"), height=10, width =10)
plot(Res_compare$beta, as.numeric(Res_compare$Beta_G), pch=19, #xlim=c(-0.1,0.23),ylim=c(-0.1,0.23),
     col=Res_compare$cluster,xlab="Log hazard per Global proteomic age gap zscore in EPIC",
     ylab="Rate of disease onset increase with age from UK health records",cex.lab=1.8, cex.axis=1.8)
text(Res_compare$beta, as.numeric(Res_compare$Beta_G), labels= Res_compare$outcomenames,  pos=pos)
legend("bottomright",pch=19,legend = c("Age cluster 1","Age cluster 3", "Age cluster 4", "Age cluster 5" ), col= c(1,3,4,5))
abline(lm(as.numeric(Res_compare$Beta_G)~Res_compare$beta ), col="blue", lty=2 )
text(0.2, 0.2, "Correlation = 0.49, p value = 0.02",col="blue", pos=1)
dev.off()

par(mar= c(5.1, 4.1, 4.1, 2.1))


####Organ clcok perfromance##########
organclocks<-c("Conventional","Organismal","Brain","Adipose", "Artery", "Immune","Heart","Intestine","Kidney","Liver","Lung","Muscle","Pancreas")

Organ_clockperformance<-clockperformance[clockperformance$Clock %in% organclocks,]
Organ_clockperformance$r_num<-as.numeric(Organ_clockperformance$r_num)

par(mar= c(5.1, 8.1, 4.1, 2.1))
pdf(file=paste0("/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/figures/OrganCorrBarplot.pdf"), height=8, width =6)
barplot(Organ_clockperformance$r_num[order(Organ_clockperformance$r_num)],
        horiz=T,xlab="Correlation with chronological age in EPIC" ,
        names.arg=Organ_clockperformance$Clock[order(Organ_clockperformance$r_num)], las=2 )
dev.off()
par(mar= c(5.1, 4.1, 4.1, 2.1))

####organ correlation heat map####
###prepare wide data###
wideclockdb<-data_clocks[data_clocks$Organ==organclocks[1],]
wideclockdb<-wideclockdb[c("idepic" , "AgeGap_zscored")]

for(i in 2:length(organclocks)){
  clockdb<-data_clocks[data_clocks$Organ==organclocks[i],"AgeGap_zscored"]
  wideclockdb<-cbind(wideclockdb,clockdb)
}
names(wideclockdb)<-c("idepic", organclocks)
remove(clockdb)

cororgan<-cor(wideclockdb[,-1])

pdf(file=paste0("/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/figures/OrganCorrheatmap.pdf"), height=4, width =6)
pheatmap(cororgan, display_numbers=T)
dev.off()
####

###make heatmap of organ-disease asociations####


organaheatmap<- as.data.frame(all_clock_results_imp[[organclocks[1]]])
rownames(organaheatmap)<-organaheatmap$outcomenames
organaheatmap<-organaheatmap["beta"]
names(organaheatmap)<-organclocks[1]

for(i in 2:length(organclocks)){
  plottable<-all_clock_results_imp[[organclocks[i]]]
  plottable<-plottable["beta"]
  names(plottable)<-organclocks[i]
  organaheatmap<-cbind(organaheatmap,plottable)
}



organaheatmap<-organaheatmap[outcomenames.sh,]

organaheatpval<- as.data.frame(all_clock_results_imp[[organclocks[1]]])
rownames(organaheatpval)<-organaheatpval$outcomenames
organaheatpval<-organaheatpval["pval"]
names(organaheatpval)<-organclocks[1]

for(i in 2:length(organclocks)){
  plottable<-all_clock_results_imp[[organclocks[i]]]
  plottable<-plottable["pval"]
  names(plottable)<-organclocks[i]
  organaheatpval<-cbind(organaheatpval,plottable)
}

organaheatpval<-organaheatpval[outcomenames.sh,]

test<-c(as.matrix(organaheatpval))
FDR<-p.adjust(test, method="BH")
FDrcutoff<-cbind(test,FDR)
FDRval_org<-max(FDrcutoff[,1][FDrcutoff[,2]<0.05])

star<-matrix(ifelse(organaheatpval < 0.05, "*", ""),nrow(organaheatpval),ncol(organaheatpval))
star[organaheatpval < FDRval_org ]<- "**"

rg<-max(abs(organaheatmap))

pdf(file=paste0("/data/Epic/subprojects/Somalogic/work/Oliver/output/for_export/figures/Organ_clock_heatmap_plot.pdf"), height=6, width =10)
pheatmap((organaheatmap),  cluster_rows=T, cluster_cols=T,display_numbers = (star),
         fontsize_number=10,
         breaks=seq(-rg,rg, length.out=50),  
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

####num hits###
organhitmatix<-matrix(ifelse(organaheatpval < FDRval_org & organaheatmap >0, 1, 0),nrow(organaheatmap),ncol(organaheatmap))


orgnanhits<-rep(NA,length(organclocks))
for (i in 1:length(organclocks)){orgnanhits[i]<-sum(organhitmatix[,i])}
diseasehtsorg<-cbind(names(organaheatmap), orgnanhits)
diseasehtsorg=diseasehtsorg[order(diseasehtsorg[,2]),]
par(mar= c(5.1, 8.1, 4.1, 2.1))
barplot(as.numeric(diseasehtsorg[,2]),
        horiz=T,xlab="Number of positive significant disease associations" ,
        names.arg=diseasehtsorg[,1], las=2 )
par(mar= c(5.1, 4.1, 4.1, 2.1))


# figure 6 -organ clocks by outcome
# Set-up theme
tm <- forest_theme(base_size = 8,
                   refline_gp = gpar("solid"),
                   ci_pch = 16,#style for the CI
                   ci_col = "black",
                   ci_lwd = 1.6, #size for the CI line
                   title_gp = gpar(cex=1,fontface="bold"),
                   legend_name = "",
                   legend_gp = gpar(fontsize = 8.5,cex = 1.1,
                                    lwd=1.8),
                   legend_position = "top",
                   core=list(bg_params=list(fill=c("white"))))#set the color of background

for(i in length(outcome)){

plottable=all_results_imp[[outcomes[i]]]
plottable=plottable[plottable$Clock %in% organclocks,]
plottable$pval_fdr<-paste0("p= ", signif(plottable$pval,2), "*")
plottable$pval_fdr[plottable$pval>FDRval_org]<-paste0("p= ",signif(plottable$pval,2))[plottable$pval>FDRval_org]
plottable=plottable[order(-plottable$Asso),]
plottable$HRest=paste0(round(plottable$Asso,2), " (",round(plottable$LCL,2), ", ", round(plottable$UCL,2), ")" )



plottable<-plottable %>% mutate(
  Asso = as.numeric(Asso),
  LCL = as.numeric(LCL),
  UCL = as.numeric(UCL),
  
  #reserve 2 digits for HR
  `HR (95% CI)` = paste(sprintf("%.2f (%.2f, %.2f)", 
                                Asso, LCL, UCL),
                        sep = "\n"),
  
  #reserve 2 digits for P value and also I use Scientific notation
  #If p-value is lower than the threshold, * is added
  pval = ifelse(pval<FDRval_org, 
                paste0(sprintf("%.2e", pval),"*"), 
                sprintf("%.2e", pval)),
  `P value` = paste(pval))

plottable$` ` <- paste(rep(" ", nrow(plottable)), collapse = " ")##create blank space for forest plot
xaxis_high<-signif(max(plottable$UCL),2)##calculate the max limit for X-axis



p <- forest(plottable[,c(1, 13, 14, 15)], #the columns should be clock, HR, P, and blank column
            est = list(plottable$Asso),
            lower = list(plottable$LCL),
            upper = list(plottable$UCL),
            ci_column = 4,
            ref_line = 1,
            nudge_y = 0.4,
            sizes = 0.6,
            theme = tm,
            xlab="Hazard ratio per age gap zscore",
            title = paste0(outcomenames[i],
                           " (n=", plottable$n[1],
                           " nevent=",plottable$nevent[1],")"),
            xlim = c(0.8,xaxis_high),
            ticks_at = seq(0.8,xaxis_high,0.2))
p

ggplot2::ggsave(filename = paste0("output/for_export/figures/",outcomenames[i],"_Organ_clock.png"), plot = p,
                dpi = 300,
                width = 7.5, height = 7.5, units = "in")

grid.newpage()
}
#############Sensitivity exclude first 5 year##############################################
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/MortalityResults_QC_v2_con2_pre5excl.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/CVDResults_QC_v2_con2_pre5excl.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/T2DResults_no_prev_t2d_QC_v2_con2_pre5excl.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/ClockCancerResults_QC_v2_con2_pre5excl.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/NDResultsv2_QC_con2_pre5excl.Rdata")

all_results_5eclu<-cancer_results
all_results_5eclu[["AD"]]<-Res_AD_full
all_results_5eclu[["PD"]]<-Res_PD_full
all_results_5eclu[["ALS"]]<-Res_ALS_full
all_results_5eclu[["DEM"]]<-Res_DEM_full
all_results_5eclu[["cvd"]]<-Res_cvd
all_results_5eclu[["chd"]]<-Res_chd
all_results_5eclu[["stroke"]]<-Res_stroke
all_results_5eclu[["t2d"]]<-Res_t2d
all_results_5eclu[["deaths"]]<-Res_deaths



all_results_5eclu_imp<-cancer_results_imp
all_results_5eclu_imp[["AD"]]<-Res_AD_imp
all_results_5eclu_imp[["PD"]]<-Res_PD_imp
all_results_5eclu_imp[["ALS"]]<-Res_ALS_imp
all_results_5eclu_imp[["DEM"]]<-Res_DEM_imp
all_results_5eclu_imp[["cvd"]]<-Res_cvd_imp
all_results_5eclu_imp[["chd"]]<-Res_chd_imp
all_results_5eclu_imp[["stroke"]]<-Res_stroke_imp
all_results_5eclu_imp[["t2d"]]<-Res_t2d_imp
all_results_5eclu_imp[["deaths"]]<-Res_deaths_imp

###make Results by clock ####
all_clock_results_5eclu = list()

for (clock_name in clock_names){
  
  Res_clock <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric(), outcome= character())
  
  for(outcome in outcomes){
    
    plottable=all_results_5eclu[[outcome]]
    plottable$outcome<-outcome
    Res_clock <-rbind(Res_clock, plottable[plottable$Clock==clock_name,])
  }
  
  Res_clock$outcomenames<-outcomenames
  Res_clock$FDR<-p.adjust(Res_clock$pval, method = "BH")
  
  all_clock_results_5eclu[[clock_name]]<-Res_clock
  
}

Oh<-all_clock_results_5eclu[["Conventional"]]
Oh$Clock<-"Oh"
all_clock_results_5eclu[["Oh"]]<-Oh
remove(Oh)


#imputed all clcok results

all_clock_results_imp_5eclu = list()

for (clock_name in clock_names){
  
  Res_clock <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric(), outcome= character())
  
  for(outcome in outcomes){
    
    plottable=all_results_5eclu_imp[[outcome]]
    plottable$outcome<-outcome
    Res_clock <-rbind(Res_clock, plottable[plottable$Clock==clock_name,])
  }
  
  Res_clock$outcomenames<-outcomenames
  Res_clock$FDR<-p.adjust(Res_clock$pval, method = "BH")
  
  all_clock_results_imp_5eclu[[clock_name]]<-Res_clock
  
}



####
#############Sensitvity exlude first 2 year##############################################
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/MortalityResults_QC_v2_con2_pre2excl.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/CVDResults_QC_v2_con2_pre2excl.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/T2DResults_no_prev_t2d_QC_v2_con2_pre2excl.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/ClockCancerResults_QC_v2_con2_pre2excl.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/NDResultsv2_QC_con2_pre2excl.Rdata")

all_results_2eclu<-cancer_results
all_results_2eclu[["AD"]]<-Res_AD_full
all_results_2eclu[["PD"]]<-Res_PD_full
all_results_2eclu[["ALS"]]<-Res_ALS_full
all_results_2eclu[["DEM"]]<-Res_DEM_full
all_results_2eclu[["cvd"]]<-Res_cvd
all_results_2eclu[["chd"]]<-Res_chd
all_results_2eclu[["stroke"]]<-Res_stroke
all_results_2eclu[["t2d"]]<-Res_t2d
all_results_2eclu[["deaths"]]<-Res_deaths



all_results_2eclu_imp<-cancer_results_imp
all_results_2eclu_imp[["AD"]]<-Res_AD_imp
all_results_2eclu_imp[["PD"]]<-Res_PD_imp
all_results_2eclu_imp[["ALS"]]<-Res_ALS_imp
all_results_2eclu_imp[["DEM"]]<-Res_DEM_imp
all_results_2eclu_imp[["cvd"]]<-Res_cvd_imp
all_results_2eclu_imp[["chd"]]<-Res_chd_imp
all_results_2eclu_imp[["stroke"]]<-Res_stroke_imp
all_results_2eclu_imp[["t2d"]]<-Res_t2d_imp
all_results_2eclu_imp[["deaths"]]<-Res_deaths_imp

###make Results by clock ####
all_clock_results_2eclu = list()

for (clock_name in clock_names){
  
  Res_clock <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric(), outcome= character())
  
  for(outcome in outcomes){
    
    plottable=all_results_2eclu[[outcome]]
    plottable$outcome<-outcome
    Res_clock <-rbind(Res_clock, plottable[plottable$Clock==clock_name,])
  }
  
  Res_clock$outcomenames<-outcomenames
  Res_clock$FDR<-p.adjust(Res_clock$pval, method = "BH")
  
  all_clock_results_2eclu[[clock_name]]<-Res_clock
  
}



#imputed all clcok results

all_clock_results_imp_2eclu = list()

for (clock_name in clock_names){
  
  Res_clock <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric(), outcome= character())
  
  for(outcome in outcomes){
    
    plottable=all_results_2eclu_imp[[outcome]]
    plottable$outcome<-outcome
    Res_clock <-rbind(Res_clock, plottable[plottable$Clock==clock_name,])
  }
  
  Res_clock$outcomenames<-outcomenames
  Res_clock$FDR<-p.adjust(Res_clock$pval, method = "BH")
  
  all_clock_results_imp_2eclu[[clock_name]]<-Res_clock
  
}


########
#############Sensitivity exclude ever smokers##############################################
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/MortalityResults_QC_v2_con2_neversmokers.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/CVDResults_QC_v2_con2_neversmokers.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/T2DResults_no_prev_t2d_QC_v2_con2_neversmokers.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/ClockCancerResults_QC_v2_con2_neversmokers.Rdata")
load("/data/Epic/subprojects/Somalogic/work/Oliver/output/NDResultsv2_QC_con2_neversmokers.Rdata")

all_results_neversmokers<-cancer_results
all_results_neversmokers[["AD"]]<-Res_AD_full
all_results_neversmokers[["PD"]]<-Res_PD_full
all_results_neversmokers[["ALS"]]<-Res_ALS_full
all_results_neversmokers[["DEM"]]<-Res_DEM_full
all_results_neversmokers[["cvd"]]<-Res_cvd
all_results_neversmokers[["chd"]]<-Res_chd
all_results_neversmokers[["stroke"]]<-Res_stroke
all_results_neversmokers[["t2d"]]<-Res_t2d
all_results_neversmokers[["deaths"]]<-Res_deaths



all_results_neversmokers_imp<-cancer_results_imp
all_results_neversmokers_imp[["AD"]]<-Res_AD_imp
all_results_neversmokers_imp[["PD"]]<-Res_PD_imp
all_results_neversmokers_imp[["ALS"]]<-Res_ALS_imp
all_results_neversmokers_imp[["DEM"]]<-Res_DEM_imp
all_results_neversmokers_imp[["cvd"]]<-Res_cvd_imp
all_results_neversmokers_imp[["chd"]]<-Res_chd_imp
all_results_neversmokers_imp[["stroke"]]<-Res_stroke_imp
all_results_neversmokers_imp[["t2d"]]<-Res_t2d_imp
all_results_neversmokers_imp[["deaths"]]<-Res_deaths_imp

###make Results by clock ####
all_clock_results_neversmokers = list()

for (clock_name in clock_names){
  
  Res_clock <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric(), outcome= character())
  
  for(outcome in outcomes){
    
    plottable=all_results_neversmokers[[outcome]]
    plottable$outcome<-outcome
    Res_clock <-rbind(Res_clock, plottable[plottable$Clock==clock_name,])
  }
  
  Res_clock$outcomenames<-outcomenames
  Res_clock$FDR<-p.adjust(Res_clock$pval, method = "BH")
  
  all_clock_results_neversmokers[[clock_name]]<-Res_clock
  
}



#imputed all clcok results

all_clock_results_imp_neversmokers = list()

for (clock_name in clock_names){
  
  Res_clock <- tibble(Clock = character(), beta = numeric(), Asso = numeric(), pval = numeric(),
                      SE = numeric(),LCL = numeric(),UCL = numeric(),n = numeric(), nevent = numeric(), outcome= character())
  
  for(outcome in outcomes){
    
    plottable=all_results_neversmokers_imp[[outcome]]
    plottable$outcome<-outcome
    Res_clock <-rbind(Res_clock, plottable[plottable$Clock==clock_name,])
  }
  
  Res_clock$outcomenames<-outcomenames
  Res_clock$FDR<-p.adjust(Res_clock$pval, method = "BH")
  
  all_clock_results_imp_neversmokers[[clock_name]]<-Res_clock
  
}

Oh<-all_clock_results_imp_neversmokers[["Conventional"]]
Oh$Clock<-"Oh"
all_clock_results_imp_neversmokers[["Oh"]]<-Oh
remove(Oh)

####

####set up plot looking at adjusted, all included and 5yr excluded side by side 
tm <- forest_theme(base_size = 10,
                   refline_gp = gpar("solid"),
                   ci_pch = c(15, 18),#style for the CI
                   ci_col = c("blue", "red"),#color for the CI
                   ci_lwd = 1.6, #size for the CI line
                   title_gp = gpar(cex=1,fontface="bold"),
                   legend_name = "",
                   legend_value = c("Whole population", "First 5 years of events excluded"),
                   legend_gp = gpar(fontsize = 8.5,cex = 1.1,
                                    lwd=1.8),
                   legend_position = "bottom",
                   core=list(bg_params=list(fill=c("white"))),#set the color of background
                   vertline_lty = c("dashed", "dotted"),
                   vertline_col = c("#d6604d", "#bababa"))



clock_name="Consensus"

plottableA=all_clock_results_imp[[clock_name]]
plottableA=plottableA[plottableA$outcomenames %in% outcomenames.sh,]
plottableB=all_clock_results_imp_5eclu[[clock_name]]
plottableB=plottableB[plottableB$outcomenames %in% outcomenames.sh,]
plottable=cbind(plottableA, plottableB)
names(plottable)[13:24]<-paste0(names(plottable)[1:12], "_adj")
plottable$pval_fdr<-paste0("p= ", signif(plottable$pval,2), "*")
plottable$FDR<-p.adjust(plottable$pval, method = "BH")
plottable$pval_fdr[plottable$FDR>0.05]<-paste0("p= ",signif(plottable$pval,2))[plottable$FDR>0.05]
plottable$pval_fdr_adj<-paste0("p= ",signif(plottable$pval_adj,2), "*")
plottable$FDR_adj<-p.adjust(plottable$pval_adj, method = "BH")
plottable$pval_fdr_adj[plottable$FDR_adj>0.05]<-paste0("p= ",signif(plottable$pval_adj,2))[plottable$FDR_adj>0.05]
plottable=plottable[order(-plottable$Asso),]


# Generate point estimation and 95% CI. 
#Paste two CIs together and separate by line break.
#add new variable of N/n_event
##load plottable first
plottable<-plottable %>% mutate(
  Asso = as.numeric(Asso),
  LCL = as.numeric(LCL),
  UCL = as.numeric(UCL),
  Asso_adj = as.numeric(Asso_adj),
  LCL_adj = as.numeric(LCL_adj),
  UCL_adj = as.numeric(UCL_adj),
  
  #reserve 2 digits for HR
  `HR (95% CI)` = paste(sprintf("%.2f (%.2f, %.2f)", 
                                Asso, LCL, UCL),
                        sprintf("%.2f (%.2f, %.2f)", 
                                Asso_adj, LCL_adj, UCL_adj),
                        sep = "\n"),
  
  #reserve 2 digits for P value and also I use Scientific notation
  #If p-value is lower than the threshold, * is added
  pval = ifelse(pval<FDRval, 
                paste0(sprintf("%.2e", pval),"*"), 
                sprintf("%.2e", pval)),
  pval_adj = ifelse(pval_adj<FDRval, 
                    paste0(sprintf("%.2e", pval_adj),"*"), 
                    sprintf("%.2e", pval_adj)),
  `P value` = paste(pval,
                    pval_adj,
                    sep = "\n"),
  Disease = paste(outcomenames," ",sep = "\n"),
  `N (N event)` = paste(
    (paste0(n," (",nevent,") ")),
    (paste0(n_adj," (",nevent_adj,") ")),
    sep = "\n"))



plottable$` ` <- paste(rep(" ", nrow(plottable)), collapse = " ")##create blank space for forest plot
xaxis_high<-signif(max(append(plottable$UCL,plottable$UCL_adj)),3)##calculate the max limit for X-axis



p <- forest(plottable[,c(29, 30, 27, 28, 31)], #the columns should be outcome, n/n_event, HR, P, and blank column
            est = list(plottable$Asso,
                       plottable$Asso_adj),
            lower = list(plottable$LCL,
                         plottable$LCL_adj),
            upper = list(plottable$UCL,
                         plottable$UCL_adj),
            ci_column = 5,
            ref_line = 1,
            nudge_y = 0.4,
            sizes = 0.9,
            theme = tm,
            title = paste0("Global clock"),
            xlab="Hazard ratio per age gap zscore",
            xlim = c(0.5,xaxis_high),
            ticks_at = seq(0.5,xaxis_high,0.5))

p

ggplot2::ggsave(filename = paste0("output/for_export/figures/",clock_name, "_all_outcomes_5years_excluded.png"), plot = p,
   dpi = 500,
   width = 10, height = 15, units = "in")


save(all_clock_results_imp, all_clock_results_adj, all_clock_results, all_clock_results_imp_5eclu, all_clock_results_imp_2eclu, all_clock_results_5eclu, all_clock_results_2eclu, all_clock_results_imp_neversmokers, all_clock_results_neversmokers,
     file="output/for_export/summary_results/summary_results.Rdata")

#save.image(file="FiguresQC_v2_con2.Rdata")

           
