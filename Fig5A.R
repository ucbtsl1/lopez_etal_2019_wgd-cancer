library(ggplot2)
library(dplyr)
library(ggrepel)

cancer = "LUSC"

significance.threshold <- 0.05

cancergenes <- as.matrix(read.csv("~/Census_allMon Sep  3 13_55_10 2018.csv"))[,1]

### read dNdScv outputs (TRACERx + TCGA combined)

## early mutations in LOH
earlyLOH <- read.table(paste0("~/",cancer,"_allGD_MLEpergene_all_genes_early_LOH.txt"), header=T)
earlyLOH.sub <- select(earlyLOH, gene_name, wnon_cv, qtrunc_cv)
names(earlyLOH.sub) <- c("gene_name", "wnon_cv_earlyLOH", "qtrunc_cv_earlyLOH")

## all mutations in noLOH
noLOH <- read.table(paste0("~/",cancer,"_allGD_MLEpergene_all_genes_all_noLOH.txt"), header=T)
noLOH.sub <- select(noLOH, gene_name, wnon_cv, qtrunc_cv)
names(noLOH.sub) <- c("gene_name", "wnon_cv_noLOH", "qtrunc_cv_noLOH")

## all mutations (i.e. not split by LOH/noLOH)
allLOH <- read.table(paste0("~/",cancer,"_allGD_MLEpergene_all_genes_all_allLOH.txt"), header=T)
allLOH.sub <- select(allLOH, gene_name, wnon_cv, qtrunc_cv)
names(allLOH.sub) <- c("gene_name", "wnon_cv_all", "qtrunc_cv_all")

truncmutations.sub <- merge(earlyLOH.sub, allLOH.sub, all=TRUE)
truncmutations <- merge(truncmutations.sub, noLOH.sub, all=TRUE)

truncmutations.significant <- filter(truncmutations, qtrunc_cv_all<significance.threshold | qtrunc_cv_earlyLOH<significance.threshold)

truncmutations.significant$cosmic <- truncmutations.significant$gene_name%in%cancergenes

where.significant <- NULL

for (i in 1:nrow(truncmutations.significant)){
  if (truncmutations.significant$qtrunc_cv_earlyLOH[i] <= significance.threshold & truncmutations.significant$qtrunc_cv_all[i] > significance.threshold) signif.i <- "earlyLOH"
  else if (truncmutations.significant$qtrunc_cv_earlyLOH[i] > significance.threshold & truncmutations.significant$qtrunc_cv_all[i] <= significance.threshold) signif.i <- "allLOH"
  else if (truncmutations.significant$qtrunc_cv_earlyLOH[i] <= significance.threshold & truncmutations.significant$qtrunc_cv_all[i] <= significance.threshold) signif.i <- "both"
  where.significant <- c(where.significant, signif.i)
}

truncmutations.significant$wheresignificant <- where.significant
limit <- max(c(truncmutations.significant$wnon_cv_earlyLOH, truncmutations.significant$wnon_cv_noLOH))

truncmutations.significant$color.wheresignificant <- truncmutations.significant$wheresignificant
truncmutations.significant$color.wheresignificant <- gsub("allLOH","orange",truncmutations.significant$color.wheresignificant)
truncmutations.significant$color.wheresignificant <- gsub("both","grey",truncmutations.significant$color.wheresignificant)
truncmutations.significant$color.wheresignificant <- gsub("earlyLOH","cyan",truncmutations.significant$color.wheresignificant)

color.vec <- sort(unique(truncmutations.significant$color.wheresignificant), decreasing = T)

pdf(paste0("~/",cancer,"_Fig5A.pdf"), width = 6, height = 5)

  ggplot(truncmutations.significant, aes(x=wnon_cv_earlyLOH, y=wnon_cv_noLOH, color=cosmic)) + 
    labs(title = cancer) +
    theme(plot.title = element_text(size=20, hjust = 0.5)) +
    geom_point(shape=21, color="black", fill="darkgrey") + 
    geom_abline(intercept = 0, linetype = 2) +
    scale_colour_manual(values=c("red","black")) +
    geom_label_repel(aes(label=truncmutations.significant$gene_name, fontface=2, fill=wheresignificant), nudge_x = 0.25, nudge_y = 0.2) +
    scale_fill_manual(values=color.vec) +
    coord_cartesian(xlim = c(0,limit+10), ylim = c(0,limit+10))


dev.off()