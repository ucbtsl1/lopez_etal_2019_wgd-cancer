cmdArgs    <- commandArgs()
cancer <- cmdArgs[6]
cancer  <- as.character(cancer)
options(scipen=999)

filename <- paste("~/TCGA_muttables_TRACERxprocessed/",cancer,".RDS",sep="")

### 1. Run dNdS
### Driver discovery (positive selection) in cancer exomes/genomes
### load libraries

library("seqinr")
library("Biostrings")
library("MASS")
library("GenomicRanges")
library("dndscv")
library("car")
library("dplyr")
library("magrittr")
library("gridExtra")

martinc.function <- "~/dndscv_SL.R"
source(martinc.function)

### input files 

cosmic.genes <- as.matrix(read.csv("~/Census_allMon Sep  3 13_55_10 2018.csv"))[,1]
essential_genes <- as.matrix(read.table("~/essential_genes_1734_list.txt")) ### list of esential genes from Blomen et al., 2015

toremove_fromessentials <- which(essential_genes%in%cosmic.genes)
if (length(toremove_fromessentials) > 0) essential_genes <- essential_genes[-toremove_fromessentials]

combined.mat=matrix(ncol=4, nrow=0)

geneset <- c("essential_genes","cancer_genes","all_genes")
timing <-c("all","late", "early")
LOHfilter <-c("noLOH", "LOH", "allLOH")


for (g in 1:length(geneset)){
  tryCatch({
  if (geneset[g] == "cancer_genes"){
    target_genes <- cosmic.genes
  }

  if (geneset[g] == "essential_genes"){
    target_genes <- essential_genes
  }
  
 
  details.raw<-readRDS(filename)
  details.raw <- details.raw[!is.na(details.raw$start),,drop=FALSE]
  details.raw <- details.raw[!is.na(details.raw$Reference_Base),,drop=FALSE]
  details.raw <- details.raw[!is.na(details.raw$Alternate_Base),,drop=FALSE]
  details.raw <- details.raw[!is.na(details.raw$mut.multi),,drop=FALSE]
 
  for (t in 1:length(timing)){
    tryCatch({
  
      for (l in 1:length(LOHfilter)){
      tryCatch({
   
      details <- details.raw
      details$mut.multi <- as.numeric(as.character(details$mut.multi))
      details$major_raw <- as.numeric(as.character(details$major_raw))
      details$minor_raw <- as.numeric(as.character(details$minor_raw))
   
   
      if (timing[t]=="early") details <- details[details$mut.multi>=1.75&details$major_raw>=1.75,,drop=FALSE]
      if (timing[t]=="late") details <- details[details$mut.multi<=1.25&details$major_raw>=1.75,,drop=FALSE]
   
      if (LOHfilter[l]=="LOH") details <- details[details$minor_raw<0.25,,drop=FALSE]
      if (LOHfilter[l]=="noLOH") details <- details[details$minor_raw>0.25,,drop=FALSE]

      mutations <- select(details, Patient, chr, start, ref, var, GD.status)
   
      colnames(mutations) <- c("sampleID","chr","pos","ref","mut", "GD")

      mutations_GenomeDoubled <- mutations[mutations$GD%in%'GD',,drop=FALSE]
      mutations_noGDfilter <- mutations
   
   
        if (geneset[g]=="all_genes")  {
        dndsout_GenomeDoubled <- dndscv_sl(mutations_GenomeDoubled)
        dndsout_noGDfilter <- dndscv_sl(mutations_noGDfilter)
    
        } else {
        dndsout_GenomeDoubled <- dndscv_sl(mutations_GenomeDoubled, gene_list=target_genes)
        dndsout_noGDfilter <- dndscv_sl(mutations_noGDfilter, gene_list=target_genes)
    
        }
      
   #### printDoubled scores per gene
   
   output.dir=paste("~/dNdScv_output_TCGA/",cancer,sep="")
   if (!file.exists(output.dir)){
     dir.create(output.dir)
   }
   
   #### print Doubled scores per gene
   sel_cvGD <- dndsout_GenomeDoubled$sel_cv
   signif_genesGD <- sel_cvGD[, c("gene_name","n_syn","n_mis","n_non","n_spl","wmis_cv","wnon_cv","wspl_cv","pmis_cv","ptrunc_cv","pallsubs_cv","qmis_cv","qtrunc_cv","qallsubs_cv")]
   colnames(signif_genesGD) <- c("gene_name","n_syn","n_mis","n_non","n_spl","wmis_cv","wnon_cv","wspl_cv","pmis_cv","ptrunc_cv","pallsubs_cv","qmis_cv","qtrunc_cv","qallsubs_cv")
   rownames(signif_genesGD) = NULL
   write.table(signif_genesGD, paste(output.dir,"/",cancer,"_GD_MLEpergene_",geneset[g],"_",timing[t],"_",LOHfilter[l],".txt", sep=""),  quote=F, row.names = F )

   #### print all scores per gene
   sel_allGD <- dndsout_noGDfilter$sel_cv
   signif_genesallGD <- sel_allGD[, c("gene_name","n_syn","n_mis","n_non","n_spl","wmis_cv","wnon_cv","wspl_cv","pmis_cv","ptrunc_cv","pallsubs_cv","qmis_cv","qtrunc_cv","qallsubs_cv")]
   colnames(signif_genesallGD) <- c("gene_name","n_syn","n_mis","n_non","n_spl","wmis_cv","wnon_cv","wspl_cv","pmis_cv","ptrunc_cv","pallsubs_cv","qmis_cv","qtrunc_cv","qallsubs_cv")
   rownames(signif_genesallGD) = NULL
   write.table(signif_genesallGD, paste(output.dir,"/",cancer,"_allGD_MLEpergene_",geneset[g],"_",timing[t],"_",LOHfilter[l],".txt", sep=""),  quote=F, row.names = F )

   ## make global output
   GDoutput.mat=as.matrix(dndsout_GenomeDoubled$globaldnds)
   rownames(GDoutput.mat)=c(paste("wmis_GD_", geneset[g],"_", timing[t],"_", LOHfilter[l], sep=""),
                            paste("wnom_GD", geneset[g],"_", timing[t],"_", LOHfilter[l], sep=""),
                            paste("wspl_GD", geneset[g],"_", timing[t],"_", LOHfilter[l], sep=""),
                            paste("wtru_GD", geneset[g],"_", timing[t],"_", LOHfilter[l], sep=""),
                            paste("wall_GD", geneset[g],"_", timing[t],"_", LOHfilter[l], sep=""))
   
   
   allGDoutput.mat=as.matrix(dndsout_noGDfilter$globaldnds)
   rownames(allGDoutput.mat)=c(paste("wmis_allGD", geneset[g],"_", timing[t],"_", LOHfilter[l], sep=""),
                              paste("wnom_allGD", geneset[g],"_", timing[t],"_", LOHfilter[l], sep=""),
                              paste("wspl_allGD", geneset[g],"_", timing[t],"_", LOHfilter[l], sep=""),
                              paste("wtru_allGD", geneset[g],"_", timing[t],"_", LOHfilter[l], sep=""),
                              paste("wall_allGD", geneset[g],"_", timing[t],"_", LOHfilter[l], sep=""))
   
   
   combined.mat=rbind(combined.mat, GDoutput.mat)
   combined.mat=rbind(combined.mat, allGDoutput.mat)
   
    } , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}
    } , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}
    } , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    } 

  colnames(combined.mat) <- c("type", "mle", "low", "high")
  write.table(combined.mat, paste(output.dir,"/",cancer,"_combinedmat.txt", sep=""), quote=F, append=TRUE)

  
  
  ##### reformat combined tables

  output.extended <- matrix(ncol=8, nrow=nrow(combined.mat))
  colnames(output.extended)=c("genelist", "doubled", "timing", "LOH", "mutation.type", "mle", "high", "low")
  output.extended[,5:8]=combined.mat[,1:4]
  output.extended[grep("essential_genes", rownames(combined.mat)),1]="essential_genes"
  output.extended[grep("all_genes", rownames(combined.mat)),1]="all_genes"
  output.extended[grep("cancer_genes", rownames(combined.mat)),1]="cancer_genes"
  output.extended[grep("_GD", rownames(combined.mat)),2]="GD"
  output.extended[grep("_allGD", rownames(combined.mat)),2]="allGD"
  output.extended[grep("early", rownames(combined.mat)),3]="early"
  output.extended[grep("late", rownames(combined.mat)),3]="late"
  output.extended[grep("genes_all_", rownames(combined.mat)),3]="all"
  output.extended[grep("_LOH", rownames(combined.mat)),4]="LOH"
  output.extended[grep("_noLOH", rownames(combined.mat)),4]="noLOH"
  output.extended[grep("_allLOH", rownames(combined.mat)),4]="all"
  rownames(output.extended)=rownames(combined.mat)
  
  write.table(output.extended, paste(output.dir,"/",cancer,"_combinedmat_EDITED.txt", sep=""), sep=";", quote=F)

 
  
