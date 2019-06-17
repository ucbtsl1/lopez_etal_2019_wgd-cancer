### snake plot "mutations_in_LOH"

  library(dplyr)

  options(scipen = 999)

  ### read tables with LOH/GD info for TCGA and TRACERx - split by early/late
  
  cancermat <- read.table("~/LOH_GD_TCGA_earlylate.txt")
  cancermat.tx100 <- read.table("~/LOH_GD_TRACERX100_earlylate.txt")
  
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
