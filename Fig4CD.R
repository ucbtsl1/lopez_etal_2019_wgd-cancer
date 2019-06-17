
library(dplyr)
library(GenomicRanges)
library(Homo.sapiens)
library(BSgenome)
library(rtracklayer)
library(gtools)
library(grid)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(DescTools)
library(gdata)

source('GDsims_functions.R')
source('essentialgenes_perarm_functions.R')


##### use essential Chr Arm Score from Davoli et al. 2013

davoli <- read.xls("s6_davoli.xlsx")
davoli$Arm<-as.character(davoli$Arm)


#### haploid LOH

dir.ascat.tcga <-"/ASCAT_TCGA/"

Ploidy <- read.table(paste0(dir.ascat.tcga, "combined.acf.ploidy.txt"), header=T, stringsAsFactors = F)
all.patients <- list.files(dir.ascat.tcga)
all.patients <- substr(all.patients,1,12)
all.patients <- all.patients[-c(1,2)]

patients <- Ploidy$ID[Ploidy$cancer_type%in%c("LUSC","LUAD")]

preGD.vec <- NULL
postGD.vec <- NULL
haploidLOH.vec <- NULL
haploidLOH.vec.GD <- NULL

n=1

patients <- patients[1:length(patients)]
for (patient in patients){
  
  print(paste0(n,"/",length(patients)))
  ASCATfile <- paste0(dir.ascat.tcga,patient,".segments.raw.txt")
  
  if(file.exists(ASCATfile)){
    seg <- read.table(ASCATfile, header=T)
    
    sub.mat.copy <- seg
    sub.mat.copy <- filter(sub.mat.copy, chr%in%1:22)
    
    sub.mat.copy.arm <- get.seg.mat.arm(sub.mat.copy)
    
    sub.mat.minor   <- cbind(as.character(sub.mat.copy.arm$sample)
                             ,as.character(sub.mat.copy.arm$chr)
                             ,sub.mat.copy.arm$startpos
                             ,sub.mat.copy.arm$endpos
                             ,apply(cbind(round(sub.mat.copy.arm$nAraw),round(sub.mat.copy.arm$nBraw)),1,min))
    
    colnames(sub.mat.minor) <- c("Sample","Chrom","Start","End","val") 
    
    sub.mat.cn      <- cbind(as.character(sub.mat.copy.arm$sample)
                             ,as.character(sub.mat.copy.arm$chr)
                             ,sub.mat.copy.arm$startpos
                             ,sub.mat.copy.arm$endpos
                             ,apply(cbind(sub.mat.copy.arm$nMajor, sub.mat.copy.arm$nMinor),1,sum))
    
    colnames(sub.mat.cn) <- c("Sample","Chrom","Start","End","val") 
    
    TCGA.barcode <- unique(sub.mat.copy$sample)
    
    GD.pval              <- genome.doub.sig.tcga(sample=TCGA.barcode,seg.mat.minor=sub.mat.minor,seg.mat.copy=sub.mat.cn,number.of.sim=10000)
    
    ploidy <- Ploidy$ploidy[Ploidy$ID==patient]
    
    GD.status            <- fun.GD.status(GD.pval=GD.pval,ploidy.val=round(ploidy))
    print(GD.status)
    
    obshaploid <- find_haploid_LOH(sample=TCGA.barcode,seg.mat.minor=sub.mat.minor,seg.mat.copy=sub.mat.cn,number.of.sim=10000)
    
    haploidLOH <- obshaploid[[1]]
    haploidLOH.vec <- c(haploidLOH.vec, haploidLOH)
    
    if (GD.status == "GD") {
      post <- obshaploid[[3]]
      haploidLOH.GD <- obshaploid[[1]]
      
      postGD.vec <- c(postGD.vec, post )
      haploidLOH.vec.GD <- c(haploidLOH.vec.GD, haploidLOH.GD)
    }
    
   
    n=n+1
  }
}


postGD <- table(postGD.vec)
haploidLOHvec <- table(haploidLOH.vec)



for (n in 1:length(postGD)){
  
  if (length(grep(".5", names(postGD)[n]))==0) names(postGD)[n]<-paste0(names(postGD)[n],"p")
  if (length(grep(".5", names(postGD)[n]))>0) names(postGD)[n]<-gsub("\\.5","q",names(postGD)[n])
  
}
postGD <- postGD[mixedorder(names(postGD))]
postGD <- as.data.frame(postGD)
colnames(postGD) <- c("Arm", "postGD")

for (n in 1:length(haploidLOHvec)){
  
  if (length(grep(".5", names(haploidLOHvec)[n]))==0) names(haploidLOHvec)[n]<-paste0(names(haploidLOHvec)[n],"p")
  if (length(grep(".5", names(haploidLOHvec)[n]))>0) names(haploidLOHvec)[n]<-gsub("\\.5","q",names(haploidLOHvec)[n])
  
}
haploidLOHvec <- haploidLOHvec[mixedorder(names(haploidLOHvec))]
haploidLOHvec <- as.data.frame(haploidLOHvec)
colnames(haploidLOHvec) <- c("Arm", "n_haploidLOH")



davoli_LOH <- merge(davoli, merge(postGD, haploidLOHvec, all=TRUE), all=TRUE)


score.type <- 'CharmEssential_score'


gg1 <- ggscatter(davoli_LOH, x = "n_haploidLOH", y = score.type, 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "n_haploidLOH", ylab = score.type,
          label = davoli_LOH$Arm)


gg2 <- ggscatter(davoli_LOH, x = "postGD", y = score.type, 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "postGD", ylab = score.type,
          label = davoli_LOH$Arm)

