#### make LOH/GD summary table for TRACERx100 samples
#### plots for Figures 2A,B


library(reshape)
library(ggplot2)
library(ggpubr)

cancer.types<- c("LUAD","LUSC")

cancer.total.mat <- matrix(nrow=0, ncol=8)
colnames(cancer.total.mat) <- c("GD", "LOH_prop","n_mutations","cancertype","LOH_mutations", "median_LOH_mutations","haploidLOH_prop","muts_in_LOH")

load("~/TX100_gd_status.RData")
gd_status <- gsub("Subclonal GD", "GD", gd_status)
gd_status <- gsub("Clonal GD", "GD", gd_status)
gd_status <- gsub("Not GD", "nGD", gd_status)


for (cancer in 1:length(cancer.types)){
  print (paste0(cancer, "/", length(cancer.types),": ",cancer.types[cancer]))
  
  Muttable <- readRDS(paste0("~/TX100_Muttables/",cancer.types[cancer],"_muttable.RDS"))
  
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
    load("~/ASCAT_TX100/tracerx.ascat.seg.129samples.r002.2.20160818.RData")
    
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
  
write.table(cancer.total.mat, "~/LOH_GD_TRACERx.txt", quote=F)



#### LOH and haploid LOH proportions (Fig 2A,B)

for (cancer in 1:length(cancer.types)){
  
  paste0(pdf("Fig2AB_LOHproportions_",cancer.types[cancer],"_tracerx.pdf", width = 4, height = 3))
  
  cancer.summary <- data.frame(cancer.total.mat)
  cancer.summary[,c(2:3,5:7)] <- apply(cancer.summary[,c(2:3,5:7)], 2, as.numeric)
  
  cancer.summary[,2] <- as.numeric(as.character(cancer.summary[,2]))
  cancer.summary[,4] <- as.character(cancer.summary[,4])
  
  cancer.summary <- cancer.summary[cancer.summary$cancertype == paste0("TX100_",cancer.types[cancer])]
  cancer.summary.m = melt(cancer.summary, measure.vars = c("LOH_prop","haploidLOH_prop"), id.vars = c("GD"))
  
  print (ggplot(data=cancer.summary.m, aes(x=variable, y=value, fill=GD)) +
           geom_boxplot()+ 
           ggtitle(paste0("TX100-",cancer.types[cancer])) +
           scale_fill_manual("", values = c("GD" = "tan3", "nGD" = "gray27")) +
           stat_compare_means(method = "t.test", label = "p.signif") + 
           facet_wrap(~variable, scale="free", nrow = 1) +
           ylim(0,0.8)+
           theme(axis.text.x = element_blank(),
                 axis.title.x = element_blank(), 
                 axis.title.y = element_blank()))
  

dev.off()
}

