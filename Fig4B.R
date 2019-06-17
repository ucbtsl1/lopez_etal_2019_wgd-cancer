library(dplyr)
library(lattice)
library(akima)
library(reshape2)
library(karyoploteR)
library(plyr)
library(ggpubr)

source('GDsims_functions.R')

dir.ascat.tcga <-"/ASCAT_TCGA/"

Ploidy <- read.table(paste0(dir.ascat.tcga, "combined.acf.ploidy.txt"), header=T, stringsAsFactors = F)

all.patients <- list.files(dir.ascat.tcga)
all.patients <- substr(all.patients,1,12)
all.patients <- all.patients[-c(1,2)]

lung.patients <- Ploidy$ID[Ploidy$cancer_type%in%c("LUSC","LUAD")]


p.val.vec <- NULL
initialLOH.vec <- NULL
gains.vec <- NULL
losses.vec <- NULL
obsLOH <- NULL
simLOH.vec <- NULL
gains.losses.summary <- matrix(ncol=9, nrow=0)


colnames(gains.losses.summary) <- c('TCGA_barcode', 'gains', 'losses', 'losses.preGD', 'losses.postGD', 'obshaploid',
                                    'mean_nhap_simulated', 'pvalue', 'ploidy')

n<-1

for (patient in lung.patients){
  
  tryCatch({
    
    print(paste0(n,"/",length(lung.patients)))
    
    ASCATfile <- paste0(dir.ascat.tcga,patient,".segments.raw.txt")
    
    if(file.exists(ASCATfile)){
      
      sub.mat.copy <- read.table(ASCATfile, header=T, stringsAsFactors = F)
      sub.mat.copy <- filter(sub.mat.copy, chr%in%1:22)
      
      ###A: for chr arm level
      sub.mat.copy.arm <- get.seg.mat.arm(sub.mat.copy)
      
      sub.mat.minor   <- cbind(sub.mat.copy.arm$sample
                               ,as.character(sub.mat.copy.arm$chr)
                               ,sub.mat.copy.arm$startpos
                               ,sub.mat.copy.arm$endpos
                               ,apply(cbind(round(sub.mat.copy.arm$nAraw),round(sub.mat.copy.arm$nBraw)),1,min))
      
      colnames(sub.mat.minor) <- c("Sample","Chrom","Start","End","val") 
      
      sub.mat.cn      <- cbind(sub.mat.copy.arm$sample
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
      
      if (GD.status == "GD"){
        
        gains_losses <- CN_parameters.tcga(TCGA.barcode, sub.mat.minor, sub.mat.cn, number.of.sim = 10000)
        
        ##only do if we observe initial LOH
        
        if (gains_losses[[1]]>0) {
          
          sims <- simHaploidLoss(gains_losses[[1]],gains_losses[[2]],gains_losses[[3]],gains_losses[[4]],sims=1000)
          
          simLOH <- sims[[1]]
          pvalsims<- sims[[2]]
          simpreGD<- sims[[3]]
          p_simpreGD<- sims[[4]]
          simpostGD<- sims[[5]]
          p_simpostGD<- sims[[6]]
          
          p.val.vec <- c(p.val.vec,pvalsims)
          initialLOH.vec <- c(initialLOH.vec, gains_losses[[1]])
          gains.vec <- c(gains.vec, gains_losses[[2]])
          losses.vec <- c(losses.vec, gains_losses[[3]])
          obsLOH <- c(obsLOH, gains_losses[[4]])
          simLOH.vec <- c(simLOH.vec, simLOH)
          
          
          gains_losses <- c(TCGA.barcode, gains_losses, sims, ploidy)
          gains.losses.summary <- rbind(gains.losses.summary, gains_losses)
          
        }
      }
      n<- n+1
    }
    
  } , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


totalresults <- cbind(p.val.vec, initialLOH.vec, gains.vec, losses.vec, obsLOH, simLOH.vec)

#######

totalresults$difference <- totalresults$simLOH.vec-totalresults$obsLOH

significant <- totalresults$p.val.vec<0.05
significant[which(significant==TRUE)] = "blue"
significant[which(significant==FALSE)] = "darkgoldenrod4"
totalresults$significant <- significant

significant2 <- totalresults$p.val.vec<0.05
totalresults$border <- !significant2

totalresults <- totalresults[order(as.numeric(totalresults$difference), decreasing=T),]



barplot(totalresults$difference, col = totalresults$significant,  ylim=c(-4,8),ylab = "Haploid LOH Sims-observed", border=totalresults$significant)
abline(h=0)
legend('topright', c('p<0.05', 'ns'), fill= c('blue', 'darkgoldenrod4'))
