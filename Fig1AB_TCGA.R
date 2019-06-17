#### make LOH/GD summary table for TCGA-LUSC and LUAD samples
#### plots for Figures 2A,B


library(reshape)
library(ggplot2)
library(ggpubr)

cancer.types<- c("LUAD","LUSC")

#### make summary table with info about GD and LOH

cancer.total.mat <- matrix(nrow=0, ncol=8)
colnames(cancer.total.mat) <- c("GD", "LOH_prop","n_mutations","cancertype","LOH_mutations", "median_LOH_mutations","haploidLOH_prop","muts_in_LOH")

for (cancer in 1:length(cancer.types)){
  
  print (paste0(cancer, "/", length(cancer.types),": ",cancer.types[cancer]))
  Muttable <- readRDS(paste0("~/TCGA_Muttables/",cancer.types[cancer],"_muttable.RDS"))
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
    copynumberfiles <- grep(patients[p], list.files("~/ASCAT_TCGA/"))
    copynumberfile <- list.files("~/ASCAT_TCGA/")[copynumberfiles[1]]
    
    copynumber <- read.table(paste0("~/ASCAT_TCGA/",copynumberfile), header=T)
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

write.table(cancer.total.mat, "~/LOH_GD_TCGA.txt", quote=F)





#### LOH and haploid LOH proportions (Fig 2A,B)

for (cancer in 1:length(cancer.types)){
  
  paste0(pdf("Fig2AB_LOHproportions_",cancer.types[cancer],".pdf", width = 4, height = 3))

  cancer.summary <- data.frame(cancer.total.mat)
  cancer.summary[,c(2:3,5:7)] <- apply(cancer.summary[,c(2:3,5:7)], 2, as.numeric)
  
  cancer.summary[,2] <- as.numeric(as.character(cancer.summary[,2]))
  cancer.summary <- cancer.summary[cancer.summary$cancertype == cancer.types[cancer],]
  cancer.summary.m = melt(cancer.summary, measure.vars = c("LOH_prop","haploidLOH_prop"), id.vars = c("GD"))
  
  print (ggplot(data=cancer.summary.m, aes(x=variable, y=value, fill=GD)) +
           geom_boxplot()+ 
           ggtitle(paste0("TCGA-",cancer.types[cancer])) +
           scale_fill_manual("", values = c("GD" = "tan3", "nGD" = "gray27")) +
           stat_compare_means(method = "t.test", label = "p.signif") + 
           facet_wrap(~variable, scale="free", nrow = 1) +
           ylim(0,0.8)+
           theme(
                 axis.text.x = element_blank(),
                 axis.title.x = element_blank(), 
                 axis.title.y = element_blank()))
  

dev.off()
}




