# Figure 2D#####
# part 1 ####
outPutLoc <- "/camp/lab/swantonc/working/mcgrann/projects/pfSim_Rver10/output_20190829_essential_mutRate///"
tableLength <- 6
removeNA <- TRUE
NameToUse <- 'GDmutRate1.'
GDsummary <- read.table(paste(outPutLoc,'/',NameToUse,'combinedSummaryAll.txt',sep=""),stringsAsFactors=FALSE,sep="\t")
GDsummary <- do.call(rbind,sapply(GDsummary,strsplit,split=" "))
LocalSPs  <- unique(GDsummary[GDsummary[,1]=='localsp',2])
GDcosts   <- unique(GDsummary[GDsummary[,1]=='GDcost',2])
localUPs   <- unique(GDsummary[GDsummary[,1]=='local.U_p',2])
localSEs   <- unique(GDsummary[GDsummary[,1]=='local.se',2])
localUEs   <- unique(GDsummary[GDsummary[,1]=='local.U_e',2])
columnOfInterest1 <- 'GDcost'
columnOfInterest2 <- 'local.U_p'


if(columnOfInterest1=='GDcost')
{
  y_column <- GDcosts
}
if(columnOfInterest1=='localsp')
{
  y_column <- LocalSPs
}
if(columnOfInterest1=='local.U_p')
{
  y_column <- localUPs
}
if(columnOfInterest1=='local.se')
{
  y_column <- localSEs
}
if(columnOfInterest1=='local.U_e')
{
  y_column <- localUEs
}



if(columnOfInterest2=='GDcost')
{
  x_column <- GDcosts
}
if(columnOfInterest2=='localsp')
{
  x_column <- LocalSPs
}
if(columnOfInterest2=='local.U_p')
{
  x_column <- localUPs
}
if(columnOfInterest2=='local.se')
{
  x_column <- localSEs
}
if(columnOfInterest2=='local.U_e')
{
  x_column <- localUEs
}


GDsummary            <- cbind(GDsummary,rep(1:(nrow(GDsummary)/tableLength), each = tableLength)  )
GDmatrix            <- matrix(NA,ncol=length(x_column),nrow=length(y_column))
GDmatrix_nruns      <- GDmatrix


#check how many runs we have of the different x and y columns
colnames(GDmatrix) <- as.character(sort(as.numeric(x_column)))
rownames(GDmatrix) <- as.character(sort(as.numeric(y_column)))
colnames(GDmatrix_nruns) <- as.character(sort(as.numeric(x_column)))
rownames(GDmatrix_nruns) <- as.character(sort(as.numeric(y_column)))


for (GDcost in y_column)
{
  IDs <- GDsummary[GDsummary[,1]==columnOfInterest1,3][(as.numeric(GDsummary[GDsummary[,1]==columnOfInterest1,2])==as.numeric(GDcost))]
  tmpSummary <- GDsummary[GDsummary[,3]%in%IDs,]
  
  for (localSP in x_column)
  {
    
    SP_IDs <- tmpSummary[tmpSummary[,1]==columnOfInterest2,3][as.numeric(tmpSummary[tmpSummary[,1]==columnOfInterest2,2])==as.numeric(localSP)]
    SP_summary <- tmpSummary[tmpSummary[,3]%in%SP_IDs,]
    GDmatrix_nruns[as.character(GDcost),as.character(localSP)] <- length(unique(SP_summary[,3]))
    GDmatrix[as.character(GDcost),as.character(localSP)] <-     mean(as.numeric(SP_summary[SP_summary[,1]=='propWGD',2]))
    
    
    
  }
}

if(removeNA)
{
  GDmatrix <- GDmatrix[,!is.na(colSums(GDmatrix)),drop=FALSE]
  GDmatrix <- GDmatrix[!is.na(rowSums(GDmatrix)),,drop=FALSE]
}

if(removeNA)
{
  GDmatrix_nruns <- GDmatrix_nruns[,colnames(GDmatrix),drop=FALSE]
  GDmatrix_nruns <- GDmatrix_nruns[rownames(GDmatrix),,drop=FALSE]
}



pdf(paste("~/Desktop/R_Figures/mutRateGD_scatter_sp_",LocalSPs,".pdf",sep=""),width=5,height=5,useDingbats = FALSE)
par(mar=c(5,5,5,5))
some.colors2<-colorRampPalette(c("aliceblue","blue4"))(1000)

plot(GDmatrix,as.numeric(colnames(GDmatrix))
     ,xlab=c('GD proportion')
     ,ylab=c('#deleterious alterations \nper cell per generation')
     ,cex=3.5
     ,main=c('GD proportion vs. mutation rate')
     ,pch=16
     ,col=some.colors2[c(round(GDmatrix*1000))]
     ,las=1)
text(x = 0,y = 1,labels = paste("sp=",LocalSPs),adj=0)
text(x = 0,y = 0.92,labels = paste("swgd=",GDcosts),adj=0)
lines(lowess(c(GDmatrix),as.numeric(colnames(GDmatrix)),f=1.5), col="red",lwd=2.5)

dev.off()
GDmatrix_part1 <- GDmatrix
LocalSPs_part1 <- LocalSPs
GDcosts_part1  <- GDcosts

# part 2 ####

outPutLoc <- "/camp/lab/swantonc/working/mcgrann/projects/pfSim_Rver10/output_201909_essential_mutRate/"
tableLength <- 6
removeNA <- TRUE
NameToUse <- 'GDmutRate.'
GDsummary <- read.table(paste(outPutLoc,'/',NameToUse,'combinedSummaryAll.txt',sep=""),stringsAsFactors=FALSE,sep="\t")
GDsummary <- do.call(rbind,sapply(GDsummary,strsplit,split=" "))
LocalSPs  <- unique(GDsummary[GDsummary[,1]=='localsp',2])
GDcosts   <- unique(GDsummary[GDsummary[,1]=='GDcost',2])
localUPs   <- unique(GDsummary[GDsummary[,1]=='local.U_p',2])
localSEs   <- unique(GDsummary[GDsummary[,1]=='local.se',2])
localUEs   <- unique(GDsummary[GDsummary[,1]=='local.U_e',2])
columnOfInterest1 <- 'GDcost'
columnOfInterest2 <- 'local.U_p'


if(columnOfInterest1=='GDcost')
{
  y_column <- GDcosts
}
if(columnOfInterest1=='localsp')
{
  y_column <- LocalSPs
}
if(columnOfInterest1=='local.U_p')
{
  y_column <- localUPs
}
if(columnOfInterest1=='local.se')
{
  y_column <- localSEs
}
if(columnOfInterest1=='local.U_e')
{
  y_column <- localUEs
}



if(columnOfInterest2=='GDcost')
{
  x_column <- GDcosts
}
if(columnOfInterest2=='localsp')
{
  x_column <- LocalSPs
}
if(columnOfInterest2=='local.U_p')
{
  x_column <- localUPs
}
if(columnOfInterest2=='local.se')
{
  x_column <- localSEs
}
if(columnOfInterest2=='local.U_e')
{
  x_column <- localUEs
}



GDsummary            <- cbind(GDsummary,rep(1:(nrow(GDsummary)/tableLength), each = tableLength)  )
GDmatrix            <- matrix(NA,ncol=length(x_column),nrow=length(y_column))
GDmatrix_nruns      <- GDmatrix


#check how many runs we have of the different x and y columns
colnames(GDmatrix) <- as.character(sort(as.numeric(x_column)))
rownames(GDmatrix) <- as.character(sort(as.numeric(y_column)))
colnames(GDmatrix_nruns) <- as.character(sort(as.numeric(x_column)))
rownames(GDmatrix_nruns) <- as.character(sort(as.numeric(y_column)))


for (GDcost in y_column)
{
  IDs <- GDsummary[GDsummary[,1]==columnOfInterest1,3][(as.numeric(GDsummary[GDsummary[,1]==columnOfInterest1,2])==as.numeric(GDcost))]
  tmpSummary <- GDsummary[GDsummary[,3]%in%IDs,]
  
  for (localSP in x_column)
  {
    
    SP_IDs <- tmpSummary[tmpSummary[,1]==columnOfInterest2,3][as.numeric(tmpSummary[tmpSummary[,1]==columnOfInterest2,2])==as.numeric(localSP)]
    SP_summary <- tmpSummary[tmpSummary[,3]%in%SP_IDs,]
    GDmatrix_nruns[as.character(GDcost),as.character(localSP)] <- length(unique(SP_summary[,3]))
    GDmatrix[as.character(GDcost),as.character(localSP)] <-     mean(as.numeric(SP_summary[SP_summary[,1]=='propWGD',2]))
    
    
    
  }
}

if(removeNA)
{
  GDmatrix <- GDmatrix[,!is.na(colSums(GDmatrix)),drop=FALSE]
  GDmatrix <- GDmatrix[!is.na(rowSums(GDmatrix)),,drop=FALSE]
}

if(removeNA)
{
  GDmatrix_nruns <- GDmatrix_nruns[,colnames(GDmatrix),drop=FALSE]
  GDmatrix_nruns <- GDmatrix_nruns[rownames(GDmatrix),,drop=FALSE]
}


pdf(paste("~/Desktop/R_Figures/mutRateGD_scatter_sp_",LocalSPs,".pdf",sep=""),width=5,height=5,useDingbats = FALSE)
par(mar=c(5,5,5,5))
some.colors2<-colorRampPalette(c("aliceblue","blue4"))(1000)

plot(GDmatrix,as.numeric(colnames(GDmatrix))
     ,xlab=c('GD proportion')
     ,ylab=c('#deleterious alterations \nper cell per generation')
     ,cex=3.5
     ,main=c('GD proportion vs. mutation rate')
     ,pch=16
     ,col=some.colors2[c(round(GDmatrix*1000))]
     ,las=1)
text(x = 0,y = 1,labels = paste("sp=",LocalSPs),adj=0)
text(x = 0,y = 0.92,labels = paste("swgd=",GDcosts),adj=0)
#text(x = 0,y = 0.84,labels = paste("#generations=1000"),adj=0)

lines(lowess(c(GDmatrix),as.numeric(colnames(GDmatrix)),f=1.5), col="red",lwd=2.5)

dev.off()

#  combined plot ####
pdf(paste("~/Desktop/R_Figures/mutRateGD_scatter_sp_",LocalSPs,".pdf",sep=""),width=5,height=5,useDingbats = FALSE)
par(mfrow=c(1,2))
par(mar=c(5,5,5,0.5))
some.colors2<-colorRampPalette(c("aliceblue","blue4"))(1000)

plot(GDmatrix_part1,as.numeric(colnames(GDmatrix_part1))
     ,xlab=c('GD proportion')
     ,ylab=c('#deleterious alterations \nper cell per generation')
     ,cex=3.5
     #,main=c('GD proportion vs. mutation rate')
     ,pch=16
     ,col=some.colors2[c(round(GDmatrix_part1*1000))]
     ,las=1)
text(x = 0,y = 1,labels = paste("sp=",LocalSPs_part1),adj=0)
text(x = 0,y = 0.92,labels = paste("swgd=",GDcosts_part1),adj=0)

lines(lowess(c(GDmatrix_part1),as.numeric(colnames(GDmatrix_part1)),f=1.5), col="red",lwd=2.5)

par(mar=c(5,0.5,5,5))
plot(GDmatrix,as.numeric(colnames(GDmatrix))
     ,xlab=c('GD proportion')
     #,ylab=c('#deleterious alterations \nper cell per generation')
     ,cex=3.5
    # ,main=c('GD proportion vs. mutation rate')
     ,pch=16
     ,col=some.colors2[c(round(GDmatrix*1000))]
     ,las=1
     ,yaxt='n')
text(x = 0,y = 1,labels = paste("sp=",LocalSPs),adj=0)
text(x = 0,y = 0.92,labels = paste("swgd=",GDcosts),adj=0)
#text(x = 0,y = 0.84,labels = paste("#generations=1000"),adj=0)


lines(lowess(c(GDmatrix),as.numeric(colnames(GDmatrix)),f=1.5), col="red",lwd=2.5)
# 
# l <- loess.sd(c(GDmatrix),as.numeric(colnames(GDmatrix)),nsigma = 1.96,span=1)
# lines(l$x, l$y)
# lines(l$x, l$upper, lty=2)
# lines(l$x, l$lower, lty=2)


dev.off()



