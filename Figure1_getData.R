# Figure 1
# generate data from TCGA and TRACERx (1A-B)
library(reshape)
library(ggplot2)
library(ggpubr)

createTableMode <- FALSE
createPlotMode  <- TRUE

# Figure 1AB
{
#### make summary table with info about GD and LOH
if(createTableMode)
{
cancer.types    <- c("LUAD","LUSC")

mutTableLoc <- "/camp/lab/swantonc/working/Saioa/TCGA_data_process/Release_dir/Mut_Tables_v5_310818_GenomeDuplication/Allmut_tables_withtiming_310818/"
copyNumLoc  <- "/camp/lab/swantonc/working/Saioa/ASCAT_TCGA_liftedhg38/"

cancer.total.mat <- matrix(nrow=0, ncol=8)
colnames(cancer.total.mat) <- c("GD", "LOH_prop","n_mutations","cancertype","LOH_mutations", "median_LOH_mutations","haploidLOH_prop","muts_in_LOH")

for (cancer in 1:length(cancer.types)){
  
  print (paste0(cancer, "/", length(cancer.types),": ",cancer.types[cancer]))
  Muttable <- readRDS(paste0(mutTableLoc,cancer.types[cancer],"_310818.RDS"))
  patients <- unique(Muttable$patient)
  
  ################################################################################
  ### make table summarizing results per patient
  
  cancer.summary.mat <- matrix(nrow = length(patients), ncol=8)
  
  rownames(cancer.summary.mat) <-  patients
  colnames(cancer.summary.mat) <- c("GD", "LOH_prop","n_mutations","cancertype","LOH_mutations", "median_LOH_mutations","haploidLOH_prop","muts_in_LOH")
  
  for (p in 1:length(patients)){
    
    print(paste0(p,"/", length(patients))) 
    cancer.summary.mat[p,4] <- cancer.types[cancer]
    
    details.patient <- Muttable[Muttable$patient==patients[p],]
    if (ncol(details.patient)==1) details.patient = t(details.p)
    
    cancer.summary.mat[p,3] <- nrow(details.patient) 
    
    ##mutations in LOH  
    mutLOHn <- nrow(details.patient[details.patient$minor_raw<0.25,])
    cancer.summary.mat[p,8] <- mutLOHn
    
    ### GD
    cancer.summary.mat[p,1] <- as.character(details.patient$GD.status[1])
    
    ## LOH proportion
    copynumberfiles <- grep(patients[p], list.files(copyNumLoc))
    copynumberfile <- list.files(copyNumLoc)[copynumberfiles[1]]
    
    copynumber <- read.table(paste0(copyNumLoc,copynumberfile), header=T)
    totallength <- copynumber[,4]-copynumber[,3]
    copynumber <- cbind(copynumber, totallength)
    
    LOH.sum <- sum(copynumber[copynumber$nBraw<0.25,9])
    total.sum <- sum(copynumber$totallength)
    LOH.prop <- LOH.sum/total.sum
    cancer.summary.mat[p,2] <- LOH.prop
    
    LOHhaploid.sum <- sum(copynumber[copynumber$nBraw<0.25 & copynumber$nMajor==1,9])
    LOHhaploid.prop <- LOHhaploid.sum/total.sum
    cancer.summary.mat[p,7] <- LOHhaploid.prop
    
    #####
    cancer.summary.mat[p,5] <- as.numeric(cancer.summary.mat[p,2])*as.numeric(cancer.summary.mat[p,3] )
  }
  
  cancer.summary.mat[,6] <- median(as.numeric(cancer.summary.mat[,5]))
  
  cancer.total.mat <- rbind(cancer.total.mat, cancer.summary.mat)
  
  #mean haploid LOH
  cancer.GD <- cancer.total.mat[which(cancer.total.mat[,1]=="GD"),]
  mean(as.numeric(cancer.GD[,7]))
  cancer.nGD <- cancer.total.mat[which(cancer.total.mat[,1]=="nGD"),]
  mean(as.numeric(cancer.nGD[,7]))
  
}

write.table(cancer.total.mat, "~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/MullersRatchetAnalysis/Data/LOH_GD_TCGA_Figure1.txt", quote=F)
}
if(createPlotMode)
{
  cancer.types    <- c("LUAD","LUSC")
  
  cancer.total.mat <- read.table("~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/MullersRatchetAnalysis/Data/LOH_GD_TCGA_Figure1.txt",header=TRUE)
  #### LOH and haploid LOH proportions (Fig 1A,B)
  
  for (cancer in 1:length(cancer.types)){
    
    #paste0(pdf("Fig2AB_LOHproportions_",cancer.types[cancer],".pdf", width = 4, height = 3))
    
    cancer.summary <- data.frame(cancer.total.mat)
    cancer.summary[,c(2:3,5:7)] <- apply(cancer.summary[,c(2:3,5:7)], 2, as.numeric)
    cancer.summary[,1] <- as.character(cancer.summary[,1])
    cancer.summary[,1] <- gsub("GD_early","GD",cancer.summary[,1])
    cancer.summary[,1] <- gsub("GD_late","GD",cancer.summary[,1])
    
    cancer.summary[,2] <- as.numeric(as.character(cancer.summary[,2]))
    cancer.summary <- cancer.summary[cancer.summary$cancertype == cancer.types[cancer],]
    
    par(mar=c(5,15,5,15))
    bp <- barplot(cbind(table(cancer.summary$GD))/nrow(cancer.summary)
                  ,col=c(rgb(39, 49, 111, maxColorValue=255, alpha=255)
                         ,rgb(141, 146, 147, maxColorValue=255, alpha=255))
                  ,beside = FALSE,las=1
                  ,border=FALSE)
    text(x = bp[1],y = 0.05,labels =table(cancer.summary$GD)[1] ,col='white')
    text(x = bp[1],y = (table(cancer.summary$GD)/nrow(cancer.summary))[1]+0.05,labels =table(cancer.summary$GD)[2] ,col='white')
    
    
    
    cancer.summary.m = melt(cancer.summary, measure.vars = c("LOH_prop","haploidLOH_prop"), id.vars = c("GD"))
    
    print (ggplot(data=cancer.summary.m, aes(x=variable, y=value, fill=GD)) +
             geom_boxplot()+ 
             ggtitle(paste0("TCGA-",cancer.types[cancer])) +
             scale_fill_manual("", values = c("GD" = rgb(39, 49, 111, maxColorValue=255, alpha=255), "nGD" = rgb(141, 146, 147, maxColorValue=255, alpha=255))) +
             stat_compare_means(method = "t.test", label = "p.signif") + 
             facet_wrap(~variable, scale="free", nrow = 1) +
             ylim(0,0.8)+
             theme(
               axis.text.x = element_blank(),
               axis.title.x = element_blank(), 
               axis.title.y = element_blank()))
    
    
    #dev.off()
  }

}

#generate data from TRACERx
if(createTableMode)
{
  mutTableLoc <- "/camp/lab/swantonc/working/Saioa/Tracerx100_110918/"
  DataFolder <- "~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/MullersRatchetAnalysis/Data/"
  
  cancer.types<- c("LUAD","LUSC")
  
  cancer.total.mat <- matrix(nrow=0, ncol=8)
  colnames(cancer.total.mat) <- c("GD", "LOH_prop","n_mutations","cancertype","LOH_mutations", "median_LOH_mutations","haploidLOH_prop","muts_in_LOH")
  
  load(paste(DataFolder,"gd_status.20160818.RData",sep=""))
  GDstatus <- gd_status
  gd_status <- gsub("Subclonal GD", "GD", gd_status)
  gd_status <- gsub("Clonal GD", "GD", gd_status)
  gd_status <- gsub("Not GD", "nGD", gd_status)
  GDstatus <- gsub("Subclonal GD", "scGD", GDstatus)
  GDstatus <- gsub("Clonal GD", "clGD", GDstatus)
  GDstatus <- gsub("Not GD", "nGD", GDstatus)
  
  
  cancer.types<- c("TX100_LUAD","TX100_LUSC")
  
  for (cancer in 1:length(cancer.types)){
    print (paste0(cancer, "/", length(cancer.types),": ",cancer.types[cancer]))
    
    Muttable <- readRDS(paste0(mutTableLoc,"mutTableAll.20160818_withGD_",cancer.types[cancer],".RDS"))
    
    patients <- unique(Muttable$SampleID)
    
    ################################################################################
    ### make table summarizing results per patient
    
    cancer.summary.mat <- matrix(nrow = length(patients), ncol=8)
    
    rownames(cancer.summary.mat) <-  patients
    colnames(cancer.summary.mat) <- c("GD", "LOH_prop","n_mutations","cancertype","LOH_mutations", "median_LOH_mutations","haploidLOH_prop","muts_in_LOH")
    
    for (p in 1:length(patients)){
      
      print(paste0(p,"/", length(patients))) 
      ### cancer type
      cancer.summary.mat[p,4] <- paste0("TX100_",cancer.types[cancer])
      
      details.patient <- Muttable[Muttable$SampleID==patients[p],]
      if (ncol(details.patient)==1) details.patient = t(details.p)
      
      cancer.summary.mat[p,3] <- nrow(details.patient) 
      
      ##mutations in LOH  
      
      details.patient.sub <- details.patient[-which(is.na(details.patient$MinorCPN)),]
      MinorCPN_region <- strsplit(gsub("R[1-9]:","",details.patient.sub$MinorCPN), ";")
      MinorCPN_region <- lapply(MinorCPN_region, as.numeric)
      
      mutLOH <- which(lapply(MinorCPN_region, function(x) if (length(x[x==0]) > 0) {TRUE} else {FALSE}) == TRUE)
      mutLOHn <- nrow(details.patient.sub[mutLOH,])
      
      cancer.summary.mat[p,8] <- mutLOHn
      
      ### GD
      cancer.summary.mat[p,1] <- gd_status[which(names(gd_status)==patients[p])]
      
      ## LOH proportion
      load(paste(DataFolder,"/tracerx.ascat.seg.129samples.r002.2.20160818.RData",sep=""))
      
      copynumber <- tracerx.ascat.seg[grep(patients[p], tracerx.ascat.seg$sample),]
      
      totallength <- copynumber[,4]-copynumber[,3]
      copynumber <- cbind(copynumber, totallength)
      
      regions <- unique(copynumber$sample)
      LOH.prop.vec <-NULL
      LOHhaploid.prop.vec <- NULL
      
      for (r in 1:length(regions)){
        copynumber.region <- copynumber[which(copynumber$sample==regions[r]),]
        
        LOH.sum.region <- sum(copynumber.region[copynumber.region$nBraw<0.25,15])
        total.sum.region <- sum(copynumber.region$totallength)
        LOH.prop.region <- LOH.sum.region/total.sum.region
        LOH.prop.vec <- c(LOH.prop.vec, LOH.prop.region)
        
        LOHhaploid.sum.region <- sum(copynumber.region[copynumber.region$nBraw<0.25 & copynumber.region$nMajor==1,15])
        LOHhaploid.prop.region <- LOHhaploid.sum.region/total.sum.region
        LOHhaploid.prop.vec <- c(LOHhaploid.prop.vec, LOHhaploid.prop.region)
        
      }
      
      cancer.summary.mat[p,2] <- mean(LOH.prop.vec)
      cancer.summary.mat[p,7] <- mean(LOHhaploid.prop.vec)
      
      
      cancer.summary.mat[p,5] <- as.numeric(cancer.summary.mat[p,2])*as.numeric(cancer.summary.mat[p,3] )
    }
    
    cancer.summary.mat[,6] <- median(as.numeric(cancer.summary.mat[,5]))
    
    cancer.total.mat <- rbind(cancer.total.mat, cancer.summary.mat)
    
    #mean haploid LOH
    cancer.GD <- cancer.total.mat[which(cancer.total.mat[,1]=="GD"),]
    mean(as.numeric(cancer.GD[,7]))
    cancer.nGD <- cancer.total.mat[which(cancer.total.mat[,1]=="nGD"),]
    mean(as.numeric(cancer.nGD[,7]))
    
  }
  
  write.table(cancer.total.mat, paste(DataFolder,"/LOH_GD_TRACERx.txt",sep=""), quote=F)
  
}
if(createPlotMode)
{
  mutTableLoc <- "/camp/lab/swantonc/working/Saioa/Tracerx100_110918/"
  DataFolder <- "~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/MullersRatchetAnalysis/Data/"
  
  cancer.types<- c("TX100_LUAD","TX100_LUSC")
  
  
  load(paste(DataFolder,"gd_status.20160818.RData",sep=""))
  GDstatus <- gd_status
  gd_status <- gsub("Subclonal GD", "GD", gd_status)
  gd_status <- gsub("Clonal GD", "GD", gd_status)
  gd_status <- gsub("Not GD", "nGD", gd_status)
  GDstatus <- gsub("Subclonal GD", "scGD", GDstatus)
  GDstatus <- gsub("Clonal GD", "clGD", GDstatus)
  GDstatus <- gsub("Not GD", "nGD", GDstatus)
  cancer.total.mat <- read.table(paste(DataFolder,"/LOH_GD_TRACERx.txt",sep=""))
  
  for (cancer in 1:length(cancer.types)){
    
    #paste0(pdf("Fig2AB_LOHproportions_",cancer.types[cancer],".pdf", width = 4, height = 3))
    
    cancer.summary <- data.frame(cancer.total.mat)
    cancer.summary[,c(2:3,5:7)] <- apply(cancer.summary[,c(2:3,5:7)], 2, as.numeric)
    cancer.summary[,1] <- as.character(cancer.summary[,1])
    cancer.summary[,1] <- gsub("GD_early","GD",cancer.summary[,1])
    cancer.summary[,1] <- gsub("GD_late","GD",cancer.summary[,1])
    
    cancer.summary[,2] <- as.numeric(as.character(cancer.summary[,2]))
    cancer.summary <- cancer.summary[cancer.summary$cancertype == cancer.types[cancer],]
    
    GDstatus[rownames(cancer.summary)]
    
    bpPlot <- cbind(table(GDstatus[rownames(cancer.summary)]))/nrow(cancer.summary)
    par(mar=c(5,15,5,15))
    bp <- barplot(  bpPlot[c('clGD','scGD','nGD'),,drop=FALSE]
                    ,col=c(rgb(39, 49, 111, maxColorValue=255, alpha=255)
                           ,rgb(159, 164, 211, maxColorValue=255, alpha=255)
                           ,rgb(141, 146, 147, maxColorValue=255, alpha=255))
                    ,beside = FALSE,las=1
                    ,border=FALSE)
    text(x = bp[1],y = 0.05,labels =table(GDstatus[rownames(cancer.summary)])['clGD'] ,col='white')
    text(x = bp[1],y = table(GDstatus[rownames(cancer.summary)])['clGD']/length(GDstatus[rownames(cancer.summary)])+0.03,labels =table(GDstatus[rownames(cancer.summary)])['scGD'] ,col='white')
    text(x = bp[1],y = 0.95,labels =table(GDstatus[rownames(cancer.summary)])['nGD'] ,col='white')
    
    
    cancer.summary.m = melt(cancer.summary, measure.vars = c("LOH_prop","haploidLOH_prop"), id.vars = c("GD"))
    
    print (ggplot(data=cancer.summary.m, aes(x=variable, y=value, fill=GD)) +
             geom_boxplot()+ 
             ggtitle(paste0("TRACERx-",cancer.types[cancer])) +
             scale_fill_manual("", values = c("GD" = rgb(39, 49, 111, maxColorValue=255, alpha=255), "nGD" = rgb(141, 146, 147, maxColorValue=255, alpha=255))) +
             stat_compare_means(method = "t.test", label = "p.signif") + 
             facet_wrap(~variable, scale="free", nrow = 1) +
             ylim(0,0.8)+
             theme(
               axis.text.x = element_blank(),
               axis.title.x = element_blank(), 
               axis.title.y = element_blank()))
    
    
    #dev.off()
  }
}

}

#Figure 1C
{
# next let's plot the boxplot looking at proportion clonality
{
library(dplyr)
library(ggplot2)
library(tidyr)
library(GenomicRanges)
library(ggbeeswarm)
library(ggpubr)
library(cowplot)
library(knitr)
source("/camp/lab/swantonc/working/lime/projects/PanCan/scripts/utils.R")
source("/camp/lab/swantonc/working/lime/projects/PanCan/scripts/cn_functions.R")
#Loading CN Data
load("/camp/lab/swantonc/working/lime/projects/PanCan/analysis_data/CN_heterogeneity/CN_heterogeneity_segmentation_output.RDS_list11042018_output/CN_Het_details_by_patient_df.RData")
cn_summaries_info_patient_df = CN_Het_details_by_patient_df
patients_to_include = unique(cn_summaries_info_patient_df$patient_id)
#Loading Sample Info
load("/camp/lab/swantonc/working/lime/projects/PanCan/analysis_data/cn_summaries/sample_data_06042018.RData")
sample_info = sample_df %>%
  filter(patient_id%in%patients_to_include)
#Patient Info
load("/camp/lab/swantonc/working/lime/projects/PanCan/analysis_data/cn_summaries/patient_data_06042018.RData")
patient_info = patient_df%>%
  filter(patient_id%in%patients_to_include)
patient_info_no_cn = patient_df%>%
  filter(!patient_id%in%patients_to_include) %>%
  select(patient_id)
cn_df = cn_summaries_info_patient_df %>%
  separate_rows(regions_cn,sep=";")
cn_df_LOH = cn_df %>%
  ungroup() %>%
  filter(regions_cn%in%c("Amp_LOH","DeepLoss_LOH","Gain_LOH","Loss_LOH","Neutral_LOH")) %>%
  select(patient_id,segName,segsize_gn_percent,CN_Subclonal) %>%
  distinct()
total_amount_of_loh = cn_df_LOH %>%
  group_by(patient_id) %>%
  summarise(total_loh=sum(segsize_gn_percent))
loh_subclonal_clonal_diff = cn_df_LOH %>%
  group_by(patient_id,CN_Subclonal) %>%
  summarise(segsize_gn_percent=sum(segsize_gn_percent)) %>%
  spread(CN_Subclonal,segsize_gn_percent,fill=0) %>%
  mutate(clonal_subclonal_diff=Clonal-Subclonal) %>%
  left_join(.,patient_info,by="patient_id") %>%
  filter(data_set_name=="jamal-hanjani")
# subclonal_clonal_diff_order = loh_subclonal_clonal_diff$patient_id[order(loh_subclonal_clonal_diff$clonal_subclonal_diff)]
# loh_subclonal_clonal_diff %>%
#   ungroup() %>%
#   mutate(patient_id=factor(patient_id,levels=subclonal_clonal_diff_order),
#          clonal_greater=ifelse(clonal_subclonal_diff>0,"Clonal","Subclonal")) %>%
#   ggplot(aes(x=patient_id,y=clonal_subclonal_diff,fill=clonal_greater)) +
#   geom_bar(stat="identity") +
#   facet_wrap(~GD_predicted,scales="free_x") +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.line.x = element_blank()) +
#   scale_fill_manual(values=c("orange","grey40"))
loh_subclonal_clonal_diff %>%
  select(patient_id,Clonal,Subclonal,GD_predicted) %>%
  melt(id.vars=c("patient_id","GD_predicted")) %>%
  ggplot(aes(x=GD_predicted,y=value,fill=variable)) +
  geom_boxplot() +
  stat_compare_means(method="t.test",paired=TRUE)+
  scale_fill_manual(values=c("grey40","orange")) +
  ylab("% Genome LOH")
}
}

#Figure 1D 
{
  library(reshape)
  library(ggpubr)
  library(scatterplot3d)
  
  BRCA_subtypes_tcga <- read.table('~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/Fig1E/tcga_BRCA_subtype.txt', header = T, stringsAsFactors = F)
  BRCA_sub <- c("TNBC", "ER", "HER2")
  
  cancer.types<- c("BRCA","LUAD","LUSC","LGG","HNSC","OV","THCA","PRAD","GBM","SKCM","COAD","STAD","BLCA","LIHC","KIRC","CESC","KIRP","SARC","ESCA","PAAD","PCPG","READ","TGCT","THYM","ACC","MESO","UVM","KICH","UCS","CHOL","DLBC","UCEC")
  
  a=1
  
  summary.mat <- matrix(ncol=12, nrow=(length(BRCA_sub)+length(cancer.types)-1))
  rownames(summary.mat) <- c(BRCA_sub, cancer.types[-which(cancer.types=="BRCA")])
  colnames(summary.mat) <- c("%GD","LOH_GD","LOH_nGD","haploidLOH_GD","haploidLOH_nGD","totalhaploid_LOH","LOH_nGD/GD","LOH_GD/GD","haploidLOH_nGD/GD","haploidLOH_GD/GD","n","mutationrate")
  
  for (c in 1:length(cancer.types)){
    
    print(c)
    print('reading summary')
    cancer.summary<-read.table(paste0("~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/Fig1E/LOH_GD/",cancer.types[c],"_GD_LOH_major1_040518.txt"))
    cancer.summary2<-read.table(paste0("~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/Fig1E/LOH_GD/",cancer.types[c],"_GD_LOH.txt"))
    
    if(cancer.types[c]=="BRCA"){
      
      for (s in 1:length(BRCA_sub)){
        
        cancer.summary<-read.table(paste0("~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/Fig1E/LOH_GD/",cancer.types[c],"_GD_LOH_major1_040518.txt"))
        cancer.summary2<-read.table(paste0("~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/Fig1E/LOH_GD/",cancer.types[c],"_GD_LOH.txt"))
        
        subset_ids <- BRCA_subtypes_tcga$Id[which(BRCA_subtypes_tcga$subtype==BRCA_sub[s])]
        cancer.summary <- cancer.summary[which(rownames(cancer.summary)%in%subset_ids),]
        cancer.summary2 <- cancer.summary2[which(rownames(cancer.summary2)%in%subset_ids),]
        
        cancer.summary.m <- melt(cancer.summary, measure.vars = c("LOH_prop"), id.vars = c("GD"))
        cancer.summary.m$cancer <- rep(cancer.types[c],nrow(cancer.summary.m))
        
        cancer.summary.m2 <- melt(cancer.summary2, measure.vars = c("LOH_prop"), id.vars = c("GD"))
        
        summary.mat[s,1] <- length(which(cancer.summary.m$GD=="GD"))/nrow(cancer.summary.m)
        summary.mat[s,2] <- mean(na.omit(cancer.summary.m$value[which(cancer.summary.m$GD=="GD")]))
        summary.mat[s,3] <- mean(na.omit(cancer.summary.m$value[which(cancer.summary.m$GD=="nGD")]))
        summary.mat[s,4] <- mean(na.omit(cancer.summary.m2$value[which(cancer.summary.m2$GD=="GD")]))
        summary.mat[s,5] <- mean(na.omit(cancer.summary.m2$value[which(cancer.summary.m2$GD=="nGD")]))
        summary.mat[s,6] <- mean(na.omit(cancer.summary.m2$value))
        summary.mat[s,7] <- -(summary.mat[c,6]/summary.mat[c,1])
        
        summary.mat[s,8] <- -(summary.mat[c,2]/summary.mat[c,1])
        summary.mat[s,9] <- -(summary.mat[c,5]/summary.mat[c,1])
        summary.mat[s,10] <- -(summary.mat[c,4]/summary.mat[c,1])
        summary.mat[s,11] <- -nrow(cancer.summary.m)
        
        
        if (a==1) cancer.summary.merge <- cancer.summary.m
        if (a!=1) cancer.summary.merge <- rbind(cancer.summary.merge, cancer.summary.m)
        
        
        a=a+1
      }
      
      
      
    }
    
    cancer.summary.m <- melt(cancer.summary, measure.vars = c("LOH_prop"), id.vars = c("GD"))
    cancer.summary.m$cancer <- rep(cancer.types[c],nrow(cancer.summary.m))
    
    cancer.summary.m2 <- melt(cancer.summary2, measure.vars = c("LOH_prop"), id.vars = c("GD"))
    
    summary.mat[c+2,1] <- length(which(cancer.summary.m$GD=="GD"))/nrow(cancer.summary.m)
    summary.mat[c+2,2] <- mean(na.omit(cancer.summary.m$value[which(cancer.summary.m$GD=="GD")]))
    summary.mat[c+2,3] <- mean(na.omit(cancer.summary.m$value[which(cancer.summary.m$GD=="nGD")]))
    summary.mat[c+2,4] <- mean(na.omit(cancer.summary.m2$value[which(cancer.summary.m2$GD=="GD")]))
    summary.mat[c+2,5] <- mean(na.omit(cancer.summary.m2$value[which(cancer.summary.m2$GD=="nGD")]))
    summary.mat[c+2,6] <- mean(na.omit(cancer.summary.m2$value))
    summary.mat[c+2,7] <- -(summary.mat[c,6]/summary.mat[c,1])
    
    summary.mat[c+2,8] <- -(summary.mat[c,2]/summary.mat[c,1])
    summary.mat[c+2,9] <- -(summary.mat[c,5]/summary.mat[c,1])
    summary.mat[c+2,10] <- -(summary.mat[c,4]/summary.mat[c,1])
    summary.mat[c+2,11] <- -nrow(cancer.summary.m)
    
    
    if (a==1) cancer.summary.merge <- cancer.summary.m
    if (a!=1) cancer.summary.merge <- rbind(cancer.summary.merge, cancer.summary.m)
    
    
    a=a+1
  }
  
  summary.mat = as.data.frame(summary.mat)
  
  #pdf("/Users/saioa/Desktop/FigS6.pdf")
  scatter.smooth(x=summary.mat$`%GD`, y=summary.mat$haploidLOH_nGD,  col="red", pch=20, cex=2, ylab = "prop. haploid LOH in nWGD tumors", xlab="prop. WGD tumors", main="")
  text(x=summary.mat$`%GD`, y=summary.mat$haploidLOH_nGD, labels=rownames(summary.mat), pos=3)
  #dev.off()
}
