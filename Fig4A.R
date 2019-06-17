library(magrittr)
library(ggplot2)
library(gridExtra)

dataset <- c("tracerx_all","tcga_lusc","tcga_luad")
par(mfrow=c(1,3))

for (d in dataset){
      if (d == "tracerx_all"){
        output.dir = "~/dNdScv_output_TRACERx/all"
        mle = paste0(output.dir,"/all_combinedmat_EDITED.txt")
      } 
  
      if (d == "tcga_lusc") {
        output.dir = "~/dNdScv_output_TCGA/LUSC"
        mle = paste0(output.dir,"/LUSC_combinedmat_EDITED.txt")
      }
        
      if (d == "tcga_luad") {
        output.dir = "~/dNdScv_output_TCGA/LUAD"
        mle = paste0(output.dir,"/LUAD_combinedmat_EDITED.txt")
      } 
        
      mle = read.table(mle, sep=";")
        

      myplots <- list()
      a <- 1
=      
      geneset <- c("essential_genes","cancer_genes", "all_genes")
      
      mle.filtered.combined <- data.frame()
      nmutations.combined <- matrix(ncol=2)
      colnames(nmutations.combined) <- c("type", "label")
      
      for (g in 1:length(geneset)){
      
            mle$genelist <- as.vector(mle$genelist)
            
            filter.genelist <- geneset[g]
            
            filter.timing <- c("early", "late","all")
            filter.LOH <- c("LOH", "noLOH")
            filter.mutation.type <- "wtru"
            
            mle.filtered <-  mle[mle$genelist==filter.genelist & 
                                   mle$mutation.type==filter.mutation.type &
                                   mle$LOH%in%filter.LOH &
                                   mle$timing%in%filter.timing
                                 ,]
            
            rowstoremove2 <- which(mle.filtered[,2]=="allGD")
            if  (length(rowstoremove2>0)) {mle.filtered <- mle.filtered[-rowstoremove2,]}
            
            rowstoremove2b <- which(mle.filtered[,2]=="nonGD")
            if  (length(rowstoremove2b>0)) {mle.filtered <- mle.filtered[-rowstoremove2b,]}
            
            rowstoremove3 <- which(mle.filtered[,2]=="GD" & mle.filtered[,3]=="all")
            if  (length(rowstoremove3>0)) {mle.filtered <- mle.filtered[-rowstoremove3,]}
            
    
            x = c(4,3,2,1)
            y = rep(-1, 4)
            
            mle.filtered$x <- x
            mle.filtered$y <- y

            
            rowstoremove <- which(mle.filtered[,2]=="nonGD")
            if  (length(rowstoremove>0)) {mle.filtered <- mle.filtered[-rowstoremove,]}
  
            if (g==1) mle.filtered.combined <- mle.filtered
            if (g>1) mle.filtered.combined <- rbind(mle.filtered.combined, mle.filtered)

            }
      
          mle.filtered.combined$mle[which(mle.filtered.combined$mle<1)] <- -1/mle.filtered.combined$mle[which(mle.filtered.combined$mle<1)]
          dotchart(mle.filtered.combined$mle)
          abline(v=1, col="red")
    
    
          mle.filtered.combined$mle.log <- log(mle.filtered.combined$mle)
          mle.filtered.combined$high.log <- log(mle.filtered.combined$high)
          mle.filtered.combined$low.log <- log(mle.filtered.combined$low)


          dotchart(mle.filtered.combined$mle.log, pt.cex=2, main=d, pch=20, xlim = c(-2, 3.5))

          x=seq(1,nrow(mle.filtered.combined))

          arrows(mle.filtered.combined$high.log,x, mle.filtered.combined$low.log, x, length=0.05, angle=90, code=3)

          abline(v=0, col="red")


}


pdf(paste0("~/Fig4.pdf"), width = 12, height = 4)

do.call(grid.arrange, c(myplots, ncol=3))

dev.off()
