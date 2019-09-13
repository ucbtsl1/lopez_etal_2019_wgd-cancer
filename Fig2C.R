# FIGURE2C

# choose which columns to plot
columnOfInterest1= 'GDcost'
columnOfInterest2= 'localsp'
# specify max fitness cost
maxSP  <- 0.1
GDCOST.RANGE=NULL
SP.RANGE=NULL
removeNA=TRUE

# load in data
GDsummary1 <- read.table("/SAN/mcgranahanlab/general/mcgrann/GDscripts/output_20190828_essential_combinedSummaryAll.txt",stringsAsFactors = FALSE,sep="\t")
GDsummary2 <- read.table("/camp/lab/swantonc/working/mcgrann/projects/pfSim_Rver10/output_20190826_essential_v2_combinedSummaryAll.txt",stringsAsFactors = FALSE,sep="\t")
GDsummary3 <- read.table("/camp/lab/swantonc/working/mcgrann/projects/pfSim_Rver10/output_20190826_essential_v3_combinedSummaryAll.txt",stringsAsFactors = FALSE,sep="\t")

GDsummary <- rbind(GDsummary1,GDsummary2,GDsummary3)
GDsummary <- do.call(rbind,sapply(GDsummary,strsplit,split=" "))
LocalSPs  <- unique(GDsummary[GDsummary[,1]=='localsp',2])
GDcosts   <- unique(GDsummary[GDsummary[,1]=='GDcost',2])
localUPs   <- unique(GDsummary[GDsummary[,1]=='local.U_p',2])
localSEs   <- unique(GDsummary[GDsummary[,1]=='local.se',2])
localUEs   <- unique(GDsummary[GDsummary[,1]=='local.U_e',2])

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
  x_column <- x_column[as.numeric(x_column)<maxSP]
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


# for the sake of argument fix UP, SE and UE, and consider GDcost versus LocalSPs
#GDsummary[,3]       <- nrow(GDsummary)/tableLength
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

#plot(as.numeric(colnames(GDmatrix)))

pdf("~/Desktop/R_Figures/Fig3C.pdf",width=8,height=6)
plot.new()
par(mar=c(6,6,3,10),new=T)
some.colors2<-colorRampPalette(c("aliceblue","blue4"))(500000)

image(t(GDmatrix)
      #, xlab = 'passenger deleteriousness'
      #,xlim=c(min(as.numeric(x_column)),max(as.numeric(x_column)))
      #,ylim=c(min(as.numeric(y_column)),max(as.numeric(y_column)))
      ,ylab='GD cost'
      ,zlim = c(0,1)
      , axes=F,col=some.colors2
      , xlab='')

box(lty=1)
axis(side = 2,las=1,labels = rownames(GDmatrix),at=seq(0,1,length.out = length(rownames(GDmatrix))))
axis(side = 1,las=2,labels = colnames(GDmatrix),at=seq(0,1,length.out = length(colnames(GDmatrix))))


#1-as.numeric(1/(1+as.numeric(rownames(GDmatrix))))

par(las=2,mgp=c(0,1,0),mar=c(0,0.5,6,0.95),fig=c(0.87,0.93,0.2,1), new=T)
colindex<-t(matrix(seq(0,100,length.out=100),ncol=1,nrow=100))
image(0,1:100,colindex,xaxt="n",yaxt="n",xlab="",ylab="",col=some.colors2)
box(lty=1)
labels.vec=round(signif(min(colindex)+(max(colindex)-min(colindex))*seq(0,1,length.out=6),2), digits = 2)
axis(4,at=seq(1,100,length.out=6),labels=labels.vec,las=2,tick=FALSE,cex.axis=1,line=-0.8)

dev.off()
