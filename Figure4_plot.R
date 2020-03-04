# first locate the TCGA chromosome data
DataFolder <- "~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/MullersRatchetAnalysis/Data/"

load("~/Downloads/RFolders/txscratch/Data/tcga.luad.trxpipe.ascat.seg.395samples.20160818.RData")
load("~/Downloads/RFolders/txscratch/Data/tcga.lusc.trxpipe.ascat.seg.346samples.20160818.RData")
load("~/Downloads/RFolders/ASCATrevised/tracerx.ascat.seg.129samples.r002.2.20160818.RData")
#tracerx.ascat.seg
regionsToUse <- read.table("~/Downloads/RFolders/ASCATrevised/ok_samples.cnv.txt",stringsAsFactors=FALSE)[,1]
tracerx.ascat.seg_df = tracerx.ascat.seg %>%
filter(sample%in%gsub("LTX0","LTX",regionsToUse))



luscMutTable <- load("~/Downloads/RFolders/txscratch/Data/mutTableTCGALUSC.20160511.RData")
luscMutTable <- get(luscMutTable)

samplesToUseLUSC <- intersect(luscMutTable$SampleID,substr(tcga.lusc.trxpipe.ascat.seg$sample,3,8))

luadMutTable <- load("~/Downloads/RFolders/txscratch/Data/mutTableTCGALUAD.20160909.RData")
luadMutTable <- get(luadMutTable)

samplesToUseLUAD <- intersect(luadMutTable$SampleID,substr(tcga.luad.trxpipe.ascat.seg$sample,3,8))


get.seg.mat.arm <- function (seg.mat.copy)
{
  load("~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/MullersRatchetAnalysis/Data/centromere.RData") ### centromere positions
  seg.mat.copy$PQ <- NA
  
  seg.mat.copy.arm  <- c()
  for (chr in c(1:22))
  {
    chr.mat.copy <- seg.mat.copy[seg.mat.copy$chr==chr,]
    
    #label those that finish before and after the centromere with obvious labelling
    chrPs <- chr.mat.copy$endpos<centromere[centromere[,1]==chr,3]
    if(TRUE%in%chrPs)
    {
      chr.mat.copy[chr.mat.copy$endpos<centromere[centromere[,1]==chr,3],]$PQ <- paste(chr,"p",sep="")
      
    }
    
    chrQs <- chr.mat.copy$startpos>centromere[centromere[,1]==chr,3]
    if(TRUE%in%chrQs)
    {
      chr.mat.copy[chr.mat.copy$startpos>centromere[centromere[,1]==chr,3],]$PQ <- paste(chr,"q",sep="")
    }
    
    chr.mat.copy.done <- chr.mat.copy[!is.na(chr.mat.copy$PQ),,drop=FALSE]
    # work out how many cases do we need to do more
    chr.mat.missing <- chr.mat.copy[is.na(chr.mat.copy$PQ),,drop=FALSE]
    
    if(nrow(chr.mat.missing)==0)
    {
      seg.mat.copy.arm <- rbind(seg.mat.copy.arm,chr.mat.copy.done)
      
    }
    
    if(nrow(chr.mat.missing)>=1)
    {
      chr.mat.missing.p <- chr.mat.missing.q <-chr.mat.missing 
      
      chr.mat.missing.p$endpos <- as.numeric(centromere[centromere[,1]==chr,3])-1
      chr.mat.missing.q$startpos <- as.numeric(centromere[centromere[,1]==chr,3])+1
      
      chr.mat.missing.p[chr.mat.missing.p$endpos<centromere[centromere[,1]==chr,3],]$PQ <- paste(chr,"p",sep="")
      chr.mat.missing.q[chr.mat.missing.q$startpos>centromere[centromere[,1]==chr,3],]$PQ <- paste(chr,"q",sep="")
      
      chr.mat.copy.compl <- rbind(chr.mat.copy.done,chr.mat.missing.p,chr.mat.missing.q)
      
      chr.mat.copy.compl <- chr.mat.copy.compl[order(chr.mat.copy.compl$sample,chr.mat.copy.compl$startpos),,drop=FALSE]
      
      seg.mat.copy.arm <- rbind(seg.mat.copy.arm,chr.mat.copy.compl)
      
    }
    
  }
  
  return(seg.mat.copy.arm)
}

#restrict to the correct samples where we have everything
tcga.luad.trxpipe.ascat.seg <- tcga.luad.trxpipe.ascat.seg[substr(tcga.luad.trxpipe.ascat.seg$sample,3,8)%in%samplesToUseLUAD,,drop=FALSE]
tcga.lusc.trxpipe.ascat.seg <- tcga.lusc.trxpipe.ascat.seg[substr(tcga.lusc.trxpipe.ascat.seg$sample,3,8)%in%samplesToUseLUSC,,drop=FALSE]

# count haploid LOH



tcga.seg <- rbind(tcga.luad.trxpipe.ascat.seg,tcga.lusc.trxpipe.ascat.seg,tracerx.ascat.seg_df)
#tcga.seg <- tcga.luad.trxpipe.ascat.seg
tcga.seg <- tcga.seg[tcga.seg$chr%in%1:22,,drop=FALSE]
tmp.seg <- get.seg.mat.arm(seg.mat.copy = tcga.seg)
tmp.seg[tmp.seg$PQ=='21p',]$PQ <- '21q'

davoli <- read.xls(paste(DataFolder,"1-s2.0-S0092867413012877-mmc6.xlsx",sep=""),skip=1,stringsAsFactors=FALSE)
davoli <- davoli[1:39,]
davoli$Arm<-as.character(davoli$Arm)


unique(tmp.seg$PQ)
unique(davoli$Arm)

# now let's count haploid LOH for each arm
chrArmHaploidMedian <- c()
chrArmHaploidMean <- c()
chrArmHaploidSum <- c()
chrArmHaploidSegs  <- c()
#chrArmPostGDLoss   <- c()
chrArmpostWGDMedian <- c()
chrArmpostWGDMean <- c()
chrArmpostWGDSum <- c()
chrArmpostWGDSegs <- c()


chrSize <- c()
for (chrArm in unique(tmp.seg$PQ))
{
  chrArmSeg  <- tmp.seg[tmp.seg$PQ==chrArm,]
  chrArmSize <- max(as.numeric(chrArmSeg$endpos))-min(as.numeric(chrArmSeg$startpos))
  chrSize<- c(chrSize,chrArmSize)
  haploidLOH_samples <- rep(0,length(unique(chrArmSeg$sample)))
  postGD_samples <- rep(0,length(unique(chrArmSeg$sample)))
  
  names(haploidLOH_samples) <- unique(chrArmSeg$sample)
  names(postGD_samples)     <- names(haploidLOH_samples)
  chrArmSegHap <- chrArmSeg[chrArmSeg$nAraw<=1&chrArmSeg$nBraw<=0.75,,drop=FALSE]
  chrArmSegWGDloss <- chrArmSeg[chrArmSeg$nBraw<1.5&chrArmSeg$nBraw>=0.75&chrArmSeg$nAraw>=1.5&(chrArmSeg$nAraw+chrArmSeg$nBraw)<=chrArmSeg$Ploidy,,drop=FALSE]
  
  for (sample in unique(chrArmSegHap$sample))
  {
    specSampleHap <- chrArmSegHap[chrArmSegHap$sample%in%sample,,drop=FALSE]
    sum(specSampleHap$endpos-specSampleHap$startpos)/chrArmSize
    
    haploidLOH_samples[sample] <- sum(as.numeric(specSampleHap$endpos)-as.numeric(specSampleHap$startpos))/chrArmSize
    
    
    
  }
  
  for (sample in unique(chrArmSegWGDloss$sample))
  {
    specSampleWGDLoss <- chrArmSegWGDloss[chrArmSegWGDloss$sample%in%sample,,drop=FALSE]
    sum(specSampleWGDLoss$endpos-specSampleWGDLoss$startpos)/chrArmSize
    
    postGD_samples[sample] <- sum(as.numeric(specSampleWGDLoss$endpos)-as.numeric(specSampleWGDLoss$startpos))/chrArmSize
    
  }
  
  
  chrArmHaploidMedian <- c(chrArmHaploidMedian,median(haploidLOH_samples))
  chrArmHaploidMean <- c(chrArmHaploidMean,mean(haploidLOH_samples))
  chrArmHaploidSum <- c(chrArmHaploidSum,sum(haploidLOH_samples))
  chrArmHaploidSegs <- c(chrArmHaploidSegs,length(haploidLOH_samples[which(haploidLOH_samples>0)]))
  
  chrArmpostWGDMedian <- c(chrArmpostWGDMedian,median(postGD_samples))
  chrArmpostWGDMean <- c(chrArmpostWGDMean,mean(postGD_samples))
  chrArmpostWGDSum <- c(chrArmpostWGDSum,sum(postGD_samples))
  chrArmpostWGDSegs <- c(chrArmpostWGDSegs,length(postGD_samples[which(postGD_samples>0)]))
  
}

names(chrArmHaploidMean) <- unique(tmp.seg$PQ)
names(chrArmHaploidSum) <- unique(tmp.seg$PQ)
names(chrSize) <-  unique(tmp.seg$PQ)

DavoliScore <- davoli$CharmEssential_score
GenomicMeasure <- chrArmpostWGDMean
funPlotDavoliScore <- function(DavoliScore
                               ,GenomicMeasure
                               ,cexSize=1.5
                               ,y_lab="Essential Score"
                               ,x_lab="Haploid Mean"
                               ,chrNames= c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p",
                               "9q","10p","10q","11p","11q","12p","12q","13q","14q","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q")
                               )
{
  
  ScoreToUse <- GenomicMeasure
  #ScoreToUse <- chrArmHaploidSum
  cexSize   <- 1.5
  crtest <- cor.test(ScoreToUse,as.numeric(as.character(DavoliScore)))
  
  #par(mar=c(3,2,3,2))
  plot(ScoreToUse,as.numeric(as.character(DavoliScore))
       ,las=1#
       #,cex=chrSize/max(chrSize)+1
       ,cex=cexSize
       ,ylim=c(0,round(max(as.numeric(as.character(DavoliScore)))))
       ,xlab=x_lab
       ,ylab=y_lab
       ,pch=16,col='gray'
  )
  abline(lm(as.numeric(as.character(DavoliScore))~ScoreToUse), col="black",lwd=3)
  
  x <- ScoreToUse
  y <- as.numeric(as.character(DavoliScore))
  lm.out <- lm(as.numeric(as.character(DavoliScore))~(ScoreToUse))
  lm.out <- lm(y ~ x)
  newx = seq(-0.1,0.12,by = 0.01)
  conf_interval <- predict(lm.out, newdata=data.frame(x=newx), interval="confidence",
                           level = 0.95)
  # lines(newx, conf_interval[,2], col="blue", lty=2)
  # lines(newx, conf_interval[,3], col="blue", lty=2)
  
  
  polygon(x = c(newx, rev(newx)),
          y = c(conf_interval[,2], 
                rev(conf_interval[,3])),
          col =  adjustcolor("gray", alpha.f = 0.5), border = NA)
  
  points(ScoreToUse,as.numeric(as.character(DavoliScore)),col='#636363',pch=16,cex=cexSize)
  
  # matlines(newx, conf_interval[,2:3], col = "blue", lty=2)
  
  
  text(x,y
       ,labels = chrNames,pos = 1,cex=0.75,offset = 0.5)
  
  text(x = max(ScoreToUse),y = as.numeric(max(DavoliScore)),labels = c(paste('r=',signif(crtest$estimate,3),'\np=',signif(crtest$p.value,2))),adj = 1)
  
}

pdf("~/Desktop/R_Figures/Figure4_plots.pdf",width=5,height=5,useDingbats = FALSE)
par(mar=c(5,5,5,5))
funPlotDavoliScore(DavoliScore = davoli$CharmEssential_score
                   ,GenomicMeasure = chrArmHaploidMean
                   ,y_lab="Essential Score"
                   ,x_lab = "Mean proportion haploid")

funPlotDavoliScore(DavoliScore = davoli$CharmEssential_score
                   ,GenomicMeasure = chrArmpostWGDMean
                   ,x_lab = "Mean proportion post duplication loss")

#print(gg1)
dev.off()

# next let's also do the simulations

library(dplyr)
library(lattice)
library(akima)
library(reshape2)
library(karyoploteR)
library(plyr)
library(ggpubr)

source('~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/Fig4CD/GDsims_functions.R')


p.val.vec <- NULL
initialLOH.vec <- NULL
gains.vec <- NULL
losses.vec <- NULL
obsLOH <- NULL
simLOH.vec <- NULL
ploidy.vec <- NULL
names.vec <- NULL
gains.losses.summary <- matrix(ncol=9, nrow=0)


colnames(gains.losses.summary) <- c('TCGA_barcode', 'gains', 'losses', 'losses.preGD', 'losses.postGD', 'obshaploid',
                                    'mean_nhap_simulated', 'pvalue', 'ploidy')

n<-1


for (patient in unique(tmp.seg$sample))
  {
  
    print(paste0(n,"/",length(unique(tmp.seg$sample))))
    

      sub.mat.copy <- tmp.seg[tmp.seg$sample%in%patient,,drop=FALSE]
      
      ###A: for chr arm level
      sub.mat.copy.arm <- get.seg.mat.arm(sub.mat.copy)
      
      sub.mat.minor   <- cbind(sub.mat.copy.arm$sample
                               ,ifelse(gsub("p","",(gsub("q",".5",(as.character(sub.mat.copy.arm$PQ)))))%in%"21","21.5",gsub("p","",(gsub("q",".5",(as.character(sub.mat.copy.arm$PQ))))))
                               ,sub.mat.copy.arm$startpos
                               ,sub.mat.copy.arm$endpos
                               ,apply(cbind(round(sub.mat.copy.arm$nAraw),round(sub.mat.copy.arm$nBraw)),1,min))
      
      colnames(sub.mat.minor) <- c("Sample","Chrom","Start","End","val") 
      
      sub.mat.cn      <- cbind(sub.mat.copy.arm$sample
                               ,ifelse(gsub("p","",(gsub("q",".5",(as.character(sub.mat.copy.arm$PQ)))))%in%"21","21.5",gsub("p","",(gsub("q",".5",(as.character(sub.mat.copy.arm$PQ))))))
                               ,sub.mat.copy.arm$startpos
                               ,sub.mat.copy.arm$endpos
                               ,apply(cbind(sub.mat.copy.arm$nMajor, sub.mat.copy.arm$nMinor),1,sum))
      
      colnames(sub.mat.cn) <- c("Sample","Chrom","Start","End","val") 
      
      
      TCGA.barcode <- unique(sub.mat.copy$sample)
      #GD.pval              <- genome.doub.sig.tcga(sample=TCGA.barcode,seg.mat.minor=sub.mat.minor,seg.mat.copy=sub.mat.cn,number.of.sim=10000)
      
      ploidy <- sub.mat.copy$Ploidy[1]
      #GD.status            <- fun.GD.status(GD.pval=GD.pval,ploidy.val=round(ploidy))
      #print(GD.status)
      
#
        
        gains_losses <- CN_parameters.tcga(TCGA.barcode, sub.mat.minor, sub.mat.cn, number.of.sim = 10000)
        
        ##only do if we observe initial LOH
        
        if (gains_losses[[1]]>0&as.numeric(ploidy>2.5)) {
          
          sims <- simHaploidLoss(gains_losses[[1]],gains_losses[[2]],gains_losses[[3]],gains_losses[[4]],sims=1000)
          
          simLOH <- sims[[1]]
          pvalsims<- sims[[2]]
          # simpreGD<- sims[[3]]
          # p_simpreGD<- sims[[4]]
          # simpostGD<- sims[[5]]
          # p_simpostGD<- sims[[6]]
          
          p.val.vec <- c(p.val.vec,pvalsims)
          initialLOH.vec <- c(initialLOH.vec, gains_losses[[1]])
          gains.vec <- c(gains.vec, gains_losses[[2]])
          losses.vec <- c(losses.vec, gains_losses[[3]])
          obsLOH <- c(obsLOH, gains_losses[[4]])
          simLOH.vec <- c(simLOH.vec, simLOH)
          ploidy.vec <- c(ploidy.vec,ploidy)
          names.vec <- c(names.vec,TCGA.barcode)
          
          gains_losses <- c(TCGA.barcode, gains_losses, sims, ploidy)
          gains.losses.summary <- rbind(gains.losses.summary, gains_losses)
          
        }
        n<- n+1
        
      }


# ploidy.table <- tmp.seg[!duplicated(tmp.seg$sample),,drop=FALSE]
# ploidy.vec   <- ploidy.table$Ploidy
# names(ploidy.vec)  <- ploidy.table$sample
# 
# unlist(gains.losses.summary[,1])



totalresults <- cbind(p.val.vec, initialLOH.vec, gains.vec, losses.vec, obsLOH, simLOH.vec,ploidy.vec,names.vec)
totalresults <- data.frame(totalresults,stringsAsFactors = F)
#######

totalresults$difference <- as.numeric(totalresults$simLOH.vec)-as.numeric(totalresults$obsLOH)

significant <- as.numeric(totalresults$p.val.vec)<0.05
significant[which(significant==TRUE)] = rgb(red = 246,green = 182,blue = 186,maxColorValue = 255)
significant[which(significant==FALSE)] = "#4D4C4C" #darkgoldenrod4"
totalresults$significant <- significant

significant2 <- totalresults$p.val.vec<0.05
totalresults$border <- !significant2

totalresults <- totalresults[order(as.numeric(totalresults$difference), decreasing=T),]

write.table(totalresults,file="~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/MullersRatchetAnalysis/Data/Simulations.txt",sep="\t"
            ,quote=FALSE
            ,col.names=NA)
pdf("~/Desktop/R_Figures/simPloidy.pdf",width=8,height=5,useDingbats = FALSE)
par(mfrow=c(2,1))
par(mar=c(0.2,3,0.5,2))
barplot(totalresults$difference, col = totalresults$significant,  ylim=c(-4,8),ylab = "Haploid LOH Sims-observed", border=totalresults$significant,las=1)
abline(h=0)
par(mar=c(6,3,0.2,2))
barplot(-as.numeric(totalresults$ploidy.vec),ylim=c(-5,-2), col = '#D4D4D4',las=1,border='#D4D4D4',axes=FALSE,ylab='Ploidy')#,  #ylim=c(-2,-6),ylab = "Ploidy", border='gray')
axis(side = 2,at = -5:-2,las=1,labels = rev(c(2:5)))
dev.off()
#abline(h=0)



