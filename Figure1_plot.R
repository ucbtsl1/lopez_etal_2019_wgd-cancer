#### make LOH/GD summary table for TCGA-LUSC and LUAD samples
#### plots for Figures 1A,B

library(reshape)
library(ggplot2)
library(ggpubr)



cancer.types<- c("LUAD","LUSC")
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


#

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(GenomicRanges)
library(ggbeeswarm)
library(ggpubr)
library(cowplot)
library(knitr)
source("/camp/lab/swantonc/working/lime/tct_lime/projects/tx_analysis/tracerx-exome-analysis/copy_num/segments/nicolai_functions.R")
#Loading TX Data
load("/camp/lab/swantonc/working/lime/tct_lime/projects/tx100_data/analysis_data/copy_number/gd_status.20160818.RData")
tx100_gd_df = data.frame(cbind(names(gd_status),NA,(gd_status),"TX100"))
colnames(tx100_gd_df) = c("ltx_id","X2","GD_Status","TX100")
load("/camp/lab/swantonc/working/mcgrann/projects/RFolders/ASCATrevised/tracerx.ascat.seg.129samples.r002.2.20160818.RData")
#tracerx.ascat.seg
regionsToUse <- read.table("/camp/lab/swantonc/working/mcgrann/projects/RFolders/ASCATrevised/ok_samples.cnv.txt",stringsAsFactors=FALSE)[,1]
tracerx.ascat.seg_df = tracerx.ascat.seg %>%
  filter(sample%in%gsub("LTX0","LTX",regionsToUse))
region_count = table(substr(unique(tracerx.ascat.seg_df$sample),1,8))
at_least2_regions = names(region_count[region_count>1])
tracerx.ascat.seg_df = tracerx.ascat.seg_df %>%
  filter(substr(sample,1,8)%in%at_least2_regions)
tracerx.ascat.seg_consen = lapply((unique(substr(tracerx.ascat.seg_df$sample,1,8))),function(sample_name){
  print(sample_name)
  tracerx.ascat.seg_sample = tracerx.ascat.seg_df %>%
    filter(grepl(sample_name,sample))
  tracerx.ascat.seg_consen_cnTotal = min.cons.fn(tracerx.ascat.seg_sample,colToUse = "cnTotal")
  tracerx.ascat.seg_consen_cnTotal = tracerx.ascat.seg_consen_cnTotal %>%
    data.table() %>%
    melt.data.table(.,id.vars=c("Chr","Start","End"),variable.name="sample_id",value.name="cnTotal")
  tracerx.ascat.seg_consen_nAraw = min.cons.fn(tracerx.ascat.seg_sample,colToUse = "nAraw")
  tracerx.ascat.seg_consen_nAraw = tracerx.ascat.seg_consen_nAraw %>%
    data.table() %>%
    melt.data.table(.,id.vars=c("Chr","Start","End"),variable.name="sample_id",value.name="nAraw")
  tracerx.ascat.seg_consen_nBraw = min.cons.fn(tracerx.ascat.seg_sample,colToUse = "nBraw")
  tracerx.ascat.seg_consen_nBraw = tracerx.ascat.seg_consen_nBraw %>%
    data.table() %>%
    melt.data.table(.,id.vars=c("Chr","Start","End"),variable.name="sample_id",value.name="nBraw")
  tracerx.ascat.seg_consen_df = full_join(tracerx.ascat.seg_consen_cnTotal,tracerx.ascat.seg_consen_nAraw,by=c("Chr","Start","End","sample_id")) %>%
    full_join(.,tracerx.ascat.seg_consen_nBraw,by=c("Chr","Start","End","sample_id"))
})
tracerx.ascat.seg_consen_df = bind_rows(tracerx.ascat.seg_consen) %>%
  mutate(patient_id=substr(sample_id,1,8)) %>%
  mutate(seg_name=paste0(Chr,":",Start,"-",End))
tracerx.ascat.seg_by_patient_df = tracerx.ascat.seg_consen_df %>%
  mutate(within_samp_cl_LOH_present=ifelse(nBraw<0.5&nAraw>=1,"hasLOH","noLOH")) %>%
  mutate(sample_LOH=paste(sample_id,within_samp_cl_LOH_present,sep="-")) %>%
  group_by(patient_id,seg_name) %>%
  summarise(LOH_status=paste(sort(unique(within_samp_cl_LOH_present)),collapse=";"),
            within_samp_cl_LOH_present_list = paste(within_samp_cl_LOH_present,collapse=";"),
            nAraw_list = paste(signif(nAraw,2),collapse=";"),
            nBraw_list = paste(signif(nBraw,2),collapse=";"),
            num_regions=length(unique(sample_id)))
tracerx.ascat.seg_anno_df = tracerx.ascat.seg_consen_df %>%
  distinct(seg_name,Chr,Start,End,patient_id) %>%
  mutate(seg_size=End-Start) %>%
  left_join(.,tracerx.ascat.seg_by_patient_df,by=c("patient_id","seg_name")) %>%
  mutate(LOH_status=ifelse(LOH_status=="hasLOH;noLOH","hasSubclonalLOH",LOH_status)) %>%
  group_by(patient_id,LOH_status) %>%
  summarise(total_seg_size=sum(seg_size)) %>%
  spread(LOH_status,total_seg_size,fill=0) %>%
  mutate(hasLOH_pc_gn=hasLOH/(hasLOH+hasSubclonalLOH+noLOH)) %>%
  mutate(hasSubclonalLOH_pc_gn=hasSubclonalLOH/(hasLOH+hasSubclonalLOH+noLOH)) %>%
  mutate(noLOH_pc_gn=noLOH/(hasLOH+hasSubclonalLOH+noLOH))
prop_genome_loh_plot = tracerx.ascat.seg_anno_df %>%
  melt(id.vars=c("patient_id"),measure.vars=c("hasLOH_pc_gn","hasSubclonalLOH_pc_gn")) %>%
  mutate(ltx_id=substr(patient_id,3,8)) %>%
  left_join(.,tx100_gd_df,by="ltx_id") %>%
  mutate(variable=gsub("hasLOH_pc_gn","Clonal LOH",variable)) %>%
  mutate(variable=gsub("hasSubclonalLOH_pc_gn","Subclonal LOH",variable)) %>%
  ggplot(aes(x=variable,y=value)) +
  geom_boxplot() +
  facet_wrap(~GD_Status) +
  theme_cowplot() +
  ylab("Prop. Genome")
#ggsave(plot=prop_genome_loh_plot,file="/camp/lab/swantonc/working/lime/tct_lime/projects/Selection/plots/prop_genome_loh_plot.pdf")
clonalGD_prop_genome_loh_plot = tracerx.ascat.seg_anno_df %>%
  melt(id.vars=c("patient_id"),measure.vars=c("hasLOH_pc_gn","hasSubclonalLOH_pc_gn")) %>%
  mutate(ltx_id=substr(patient_id,3,8)) %>%
  left_join(.,tx100_gd_df,by="ltx_id") %>%
  filter(GD_Status=="Clonal GD") %>%
  mutate(variable=gsub("hasLOH_pc_gn","Clonal LOH",variable)) %>%
  mutate(variable=gsub("hasSubclonalLOH_pc_gn","Subclonal LOH",variable)) %>%
  ggplot(aes(x=variable,y=value)) +
  geom_boxplot() +
  facet_wrap(~GD_Status) +
  theme_cowplot() +
  ylab("Prop. Genome")
#ggsave(plot=prop_genome_loh_plot,file="/camp/lab/swantonc/working/lime/tct_lime/projects/Selection/plots/clonalGD_prop_genome_loh_plot.pdf")
#P-val and Effect Size
GRPA = tracerx.ascat.seg_anno_df %>%
  melt(id.vars=c("patient_id"),measure.vars=c("hasLOH_pc_gn","hasSubclonalLOH_pc_gn")) %>%
  mutate(ltx_id=substr(patient_id,3,8)) %>%
  left_join(.,tx100_gd_df,by="ltx_id") %>%
  mutate(variable=gsub("hasLOH_pc_gn","Clonal LOH",variable)) %>%
  mutate(variable=gsub("hasSubclonalLOH_pc_gn","Subclonal LOH",variable)) %>%
  filter(GD_Status=="Subclonal GD") %>%
  filter(variable=="Clonal LOH") %>%
  pull(value)
GRPB = tracerx.ascat.seg_anno_df %>%
  melt(id.vars=c("patient_id"),measure.vars=c("hasLOH_pc_gn","hasSubclonalLOH_pc_gn")) %>%
  mutate(ltx_id=substr(patient_id,3,8)) %>%
  left_join(.,tx100_gd_df,by="ltx_id") %>%
  mutate(variable=gsub("hasLOH_pc_gn","Clonal LOH",variable)) %>%
  mutate(variable=gsub("hasSubclonalLOH_pc_gn","Subclonal LOH",variable)) %>%
  filter(GD_Status=="Subclonal GD") %>%
  filter(variable=="Subclonal LOH") %>%
  pull(value)
median(GRPA)
median(GRPB)
t.test(GRPA,GRPB,exact=FALSE,paired=T)$p.value
abs(mean(GRPA)-mean(GRPB))/sd(c(GRPA,GRPB))

# final figure

library(reshape)
library(ggpubr)
library(scatterplot3d)
dataFolder <- "~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/MullersRatchetAnalysis/Data/"


BRCA_subtypes_tcga <- read.table(paste(dataFolder,'/tcga_BRCA_subtype.txt',sep=""), header = T, stringsAsFactors = F)
BRCA_sub <- c("TNBC", "ER", "HER2")

cancer.types<- c("BRCA","LUAD","LUSC","LGG","HNSC","OV","THCA","PRAD","GBM","SKCM","COAD","STAD","BLCA","LIHC","KIRC","CESC","KIRP","SARC","ESCA","PAAD","PCPG","READ","TGCT","THYM","ACC","MESO","UVM","KICH","UCS","CHOL","DLBC","UCEC")

a=1

summary.mat <- matrix(ncol=12, nrow=(length(BRCA_sub)+length(cancer.types)-1))
rownames(summary.mat) <- c(BRCA_sub, cancer.types[-which(cancer.types=="BRCA")])
colnames(summary.mat) <- c("%GD","LOH_GD","LOH_nGD","haploidLOH_GD","haploidLOH_nGD","totalhaploid_LOH","LOH_nGD/GD","LOH_GD/GD","haploidLOH_nGD/GD","haploidLOH_GD/GD","n","mutationrate")

for (c in 1:length(cancer.types)){
  
  print(c)
  print('reading summary')
  cancer.summary<-read.table(paste0(dataFolder,cancer.types[c],"_GD_LOH_major1_040518.txt"))
  cancer.summary2<-read.table(paste0(dataFolder,cancer.types[c],"_GD_LOH.txt"))
  
  if(cancer.types[c]=="BRCA"){
    
    for (s in 1:length(BRCA_sub)){
      
      cancer.summary<-read.table(paste0(dataFolder,cancer.types[c],"_GD_LOH_major1_040518.txt"))
      cancer.summary2<-read.table(paste0(dataFolder,cancer.types[c],"_GD_LOH.txt"))
      
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
cor.test(summary.mat$`%GD`,summary.mat$haploidLOH_nGD)
plot(summary.mat$`%GD`,summary.mat$haploidLOH_nGD)

0.6468971*0.6468971
#dev.off()






