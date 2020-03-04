library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

DataFolder <- "~/Google_Drive/Manuscripts/Mullers_Ratchet_paper/FinalSubmission/Rscripts/MullersRatchetAnalysis/Data/"

regions <- c("D_HCT001_BS",
             "D_HCT001_SU_T1-R1",
             "D_HCT001_SU_T10-R1",            
             "D_HCT001_SU_T11-R1","D_HCT001_SU_T12-R1",                 
             "D_HCT001_SU_T13-R1","D_HCT001_SU_T14-R1",                 
             "D_HCT001_SU_T2-R1","D_HCT001_SU_T3-R1",                  
             "D_HCT001_SU_T4-R1","D_HCT001_SU_T5-R1",                  
             "D_HCT001_SU_T6-R1","D_HCT001_SU_T7-R1",                  
             "D_HCT001_SU_T8-R1","D_HCT001_SU_T9-R1")

mut.table <- read.table(paste(DataFolder,"HCT001.Exome.SNV.xls",sep="")
                        ,sep="\t"
                        ,stringsAsFactors=FALSE
                        ,header=TRUE)

muttable.sub <- mut.table[,c(1,grep("Varscan.D_HCT001", colnames(mut.table)))]


myplots <- list()
a <- 1

###### common mutations

commonmutations <- muttable.sub[which(rowSums(muttable.sub[,2:15])==14),1]

timing.mat <- matrix(nrow=length(commonmutations), ncol=length(regions))
rownames(timing.mat) <- commonmutations
colnames(timing.mat) <- regions

for(c in 1:length(regions)){
  
  TCGA.earlyLate <- read.table(paste0(DataFolder,regions[c],"_earlyORlate.txt"))
  timing.mat[,c] <- as.character(TCGA.earlyLate[match(rownames(timing.mat), TCGA.earlyLate$key),ncol(TCGA.earlyLate)])
  
}

timing.mat <- as.data.frame(timing.mat)
colnames(timing.mat) <- c("T1","T4_DC14","T50_TC13","T50_TC16","T50_TC17","T50_TC17b","T50_TC35","T4_DC25","T4_TC13","T4_TC16","T4_TC17","T4_TC17b","T4_TC35","T50_DC14","T50_DC25") 
timing.mat <- timing.mat[,-c(1,2,8,14,15)]


timing.mat2 <- timing.mat %>% 
  lapply(table) %>% 
  lapply(as.data.frame) %>% 
  Map(cbind,var = names(timing.mat),.) %>% 
  rbind_all() %>% 
  group_by(var) %>% 
  mutate(prop = Freq / sum(Freq))

timing.mat3 <- select(timing.mat2, var, Var1, prop)
mean(unlist(timing.mat3[which(timing.mat3[,2]=="early"),3]))
sd(unlist(timing.mat3[which(timing.mat3[,2]=="early"),3]))

myplots[[a]] <- ggplot(timing.mat3, aes(factor(var), prop, fill = Var1)) + 
  ggtitle("common mutations") +
  geom_bar(stat="identity", position = "fill") + 
  scale_fill_brewer(palette = "Set1", name= "timing") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", y="Proportion")


###### private mutations

muts_not_indiploids <- muttable.sub[which(rowSums(muttable.sub[,c(2,3,9,10)])==0),1]


timing.mat <- matrix(nrow=length(muts_not_indiploids), ncol=length(regions))
rownames(timing.mat) <- muts_not_indiploids
colnames(timing.mat) <- regions

for(c in 1:length(regions)){
  
  TCGA.earlyLate <- read.table(paste0(DataFolder,regions[c],"_earlyORlate.txt"))
  timing.mat[,c] <- as.character(TCGA.earlyLate[match(rownames(timing.mat), TCGA.earlyLate$key),ncol(TCGA.earlyLate)])
  
}

timing.mat <- as.data.frame(timing.mat)
colnames(timing.mat) <- c("T1","T4_DC14","T50_TC13","T50_TC16","T50_TC17","T50_TC17b","T50_TC35","T4_DC25","T4_TC13","T4_TC16","T4_TC17","T4_TC17b","T4_TC35","T50_DC14","T50_DC25") 
timing.mat <- timing.mat[,-c(1,2,8,14,15)]

timing.mat2 <- timing.mat %>% 
  lapply(table) %>% 
  lapply(as.data.frame) %>% 
  Map(cbind,var = names(timing.mat),.) %>% 
  rbind_all() %>% 
  group_by(var) %>% 
  mutate(prop = Freq / sum(Freq))

timing.mat3 <- select(timing.mat2, var, Var1, prop)
mean(unlist(timing.mat3[which(timing.mat3[,2]=="late"),3]))
sd(unlist(timing.mat3[which(timing.mat3[,2]=="late"),3]))

myplots[[a+1]] <- ggplot(timing.mat3, aes(factor(var), prop, fill = Var1)) + 
  ggtitle("tetraploid private mutations") +
  geom_bar(stat="identity", position = "fill") + 
  scale_fill_brewer(palette = "Set1", name= "timing") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", y="Proportion")

pdf ("Fig3B.pdf")
do.call(grid.arrange, c(myplots, nrow=2))
dev.off()


### snake plot "mutations_in_LOH"

library(dplyr)

options(scipen = 999)

### read tables with LOH/GD info for TCGA and TRACERx - split by early/late

cancermat <- read.table(paste(DataFolder,"LOH_GD_TCGA.txt",sep=""))
cancermat.tx100 <- read.table(paste(DataFolder,"LOH_GD_TRACERX100.txt",sep=""))

### only LUSC and LUAD

cancermat <- filter(cancermat, cancertype%in%c("LUSC", "LUAD"))
cancermat$cancertypeGD <- paste0("TCGA-",cancermat$cancertype, "_", cancermat$GD)

cancermat.tx100$cancertype <- gsub("_","-",cancermat.tx100$cancertype )
cancermat.tx100$cancertypeGD <- paste0(cancermat.tx100$cancertype, "_", cancermat.tx100$GD)

cancermat <- rbind(cancermat, cancermat.tx100)
cancermat <- cancermat[order(cancermat$cancertypeGD),]
cancers <- unique(cancermat$cancertypeGD)


median.vec <- NULL

for (c in 1:length(cancers)){
  
  median <- median(as.numeric(cancermat[cancermat$cancertypeGD == cancers[c],8]))
  median <- rep(median, nrow(cancermat[cancermat$cancertypeGD == cancers[c],]))
  median.vec <- c(median.vec, median)
}

cancermat$median.mutsLOH <- median.vec


###

### manually order
manual.order <- c("TX100-LUSC_nGD","TX100-LUSC_GD_early","TX100-LUSC_GD_late","TCGA-LUSC_nGD","TCGA-LUSC_GD_early","TCGA-LUSC_GD_late","TX100-LUAD_nGD","TX100-LUAD_GD_early","TX100-LUAD_GD_late","TCGA-LUAD_nGD","TCGA-LUAD_GD_early","TCGA-LUAD_GD_late")
cancermat$order <- match(cancermat$cancertypeGD, manual.order)

cancermat <- cancermat[order(cancermat$order, cancermat$muts_in_LOH),]

sampleTable <- cancermat

cancer.order <- unique(sampleTable$cancertypeGD)
numberOfCancerTypes <- length(unique(sampleTable$cancertypeGD))

x <- c(0.05,(numberOfCancerTypes)+0.2)

ylab <- paste('muts_inLOH',sep="")
ylim <- c(-0.15,3.5)


pdf("~/Fig3c.pdf", height = 7, width = 10)
par(mar=c(12,4.1,4.1,2.1))

plot(1
     ,type="n"
     ,axes=F
     ,xlim=c(min(x),max(x))
     ,ylim= ylim
     ,xlab = ''
     ,ylab = ylab
     ,xaxt='n'
     ,xaxs="i"
     ,yaxs="i"
     ,pch=16
     ,main="")

axis(side = 2,at = 0:4,c(1,10,100,1000,10000) ,las=2 )

#plot the white/grey rect
for (i in seq(1,length(cancer.order),by=6))
{
  xleft   <- i-0.90
  xright  <- i+2.10
  ybottom <- -4.5
  ytop    <- 4.05
  
  rect(xleft,ybottom,xright,ytop,col='grey90',border=FALSE)
  
}

box()
for (x in c(-4,-3,-2,-1,0,1,2,3,4))
{
  abline(h=x,lty='dashed',col='gray',)
  
}
axis(side=1,at=seq(from=0.625,by=1,length.out = length(cancer.order)),labels = cancer.order,las=2,tick = FALSE, font = 2)


# start with the first cancerType
k <- 1
for (cancerType in cancer.order)
{
  cancerTypeTable <-sampleTable[sampleTable$cancertypeGD%in%cancerType,,drop=FALSE]
  
  x.start <- k-0.55
  x.end   <- k-0.20
  y.start  <- log(median(as.numeric(cancerTypeTable$muts_in_LOH)),10)
  
  y.end    <- y.start
  segments(x.start,y.start,x.end,y.end,lwd=4,col='orange')
  
  
  start <- k -0.75
  end   <- k  
  length.out <- nrow(cancerTypeTable)
  xseq <- seq(start,end,length.out = length.out)
  yseq <- log(as.numeric(cancerTypeTable$muts_in_LOH),10)
  
  
  if (cancerType%in%c("TX100-LUAD_nGD", "TX100-LUSC_nGD","TCGA-LUAD_nGD", "TCGA-LUSC_nGD")) points(xseq,yseq,pch=16,col="black")  
  if (cancerType%in%c("TX100-LUAD_GD_early", "TX100-LUSC_GD_early","TCGA-LUAD_GD_early", "TCGA-LUSC_GD_early")) points(xseq,yseq,pch=16,col="brown2")
  if (cancerType%in%c("TX100-LUAD_GD_late", "TX100-LUSC_GD_late","TCGA-LUAD_GD_late", "TCGA-LUSC_GD_late")) points(xseq,yseq,pch=16,col="deepskyblue4")
  
  
  k <- k+1
  
  
}

dev.off()