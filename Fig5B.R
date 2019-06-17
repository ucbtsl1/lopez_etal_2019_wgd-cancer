library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(gplots)


cancer.types<- c("SKCM","LUSC","LUAD","COAD","HNSC","BRCA", "BLCA","CESC","LIHC","KIRC","KIRP","ESCA","SARC","PRAD","GBM","LGG","PAAD","THCA","READ","ACC","OV","KICH","CHOL","PCPG","THYM","UCS","UVM","MESO","TGCT","DLBC","UCEC","LAML","STAD")
cancer.types<- sort(cancer.types)

cancergenes.mat <- read.csv("~/Census_allMon Sep  3 13_55_10 2018.csv")
cancergenes <- as.character(unique(cancergenes.mat$Gene.Symbol))
  
mat.list <- list()

significant.count <- matrix(nrow=length(cancer.types), ncol=3)
rownames(significant.count) <- cancer.types
colnames(significant.count) <- c("all", "LOH","whichgenes")


a=1
for (c in 1:length(cancer.types)){
  tryCatch({
print(c)
significance.threshold <- 0.05

####
#read dNdScv outputs for TCGA-MC3 calls
####
earlyLOH <- read.table(paste0("~/",cancer.types[c],"_allGD_MLEpergene_all_genes_early_LOH.txt"), header=T)
earlyLOH.sub <- select(earlyLOH, gene_name, wnon_cv, qtrunc_cv)
names(earlyLOH.sub) <- c("gene_name", "wnon_cv_earlyLOH", "q_earlyLOH")

noLOH <- read.table(paste0("~/",cancer.types[c],"_allGD_MLEpergene_all_genes_all_noLOH.txt"), header=T)
noLOH.sub <- select(noLOH, gene_name, wnon_cv, qtrunc_cv)
names(noLOH.sub) <- c("gene_name", "wnon_cv_noLOH", "q_noLOH")

allLOH <- read.table(paste0("~/",cancer.types[c],"_allGD_MLEpergene_all_genes_all_allLOH.txt"), header=T)
allLOH.sub <- select(allLOH, gene_name, wnon_cv, qtrunc_cv)
names(noLOH.sub) <- c("gene_name", "wnon_cv_all", "q_all")

truncmutations <- merge(earlyLOH.sub, noLOH.sub, all=TRUE)
truncmutations <- merge(truncmutations, noLOH.sub, all=TRUE)

significant.count[c,1] <- length(which(truncmutations$q_all<significance.threshold))
significant.count[c,2] <- length(which(truncmutations$q_earlyLOH<significance.threshold & truncmutations$q_all>significance.threshold))
significant.count[c,3] <- paste(as.character(truncmutations$gene_name[which(truncmutations$q_earlyLOH<significance.threshold & truncmutations$q_all>significance.threshold)]), collapse=";")

truncmutations.significant <- filter(truncmutations, q_earlyLOH<=significance.threshold, q_all>significance.threshold, q_noLOH>significance.threshold)

truncmutations.significant <- truncmutations.significant[order(truncmutations.significant$q_earlyLOH, truncmutations.significant$q_noLOH, truncmutations.significant$q_all),]


names(truncmutations.significant) <- c("gene", paste0(cancer.types[c],"_cv_earlyLOH"), paste0(cancer.types[c],"_q_earlyLOH"), paste0(cancer.types[c],"_cv_noLOH"), paste0(cancer.types[c],"_q_noLOH"), paste0(cancer.types[c],"_cv_all"), paste0(cancer.types[c],"_q_all"))

if (a==1) truncmutations.merged<-truncmutations.significant
else (truncmutations.merged <- merge(truncmutations.merged, truncmutations.significant, all=TRUE))

a=a+1

} , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


significant.count <- significant.count[-which(significant.count[,2]==0),]
significant.count <- significant.count[order(significant.count[,2], decreasing = T),]

cancers.tokeep.withnew <- rownames(significant.count)

genestokeep <- as.character(truncmutations.merged$gene)


###

a=1
for (c in 1:length(cancers.tokeep.withnew)){
  tryCatch({
    print(c)
     
    ####
    earlyLOH <- read.table(paste0("mount/camp/lab/swantonc/working/Saioa/TCGA_data_MC3_process/dNdS/",cancers.tokeep.withnew[c],"/",cancers.tokeep.withnew[c],"_allGD_MLEpergene_all_genes_early_LOH.txt"), header=T)
    earlyLOH.sub <- select(earlyLOH, gene_name, wnon_cv, qtrunc_cv)
    names(earlyLOH.sub) <- c("gene_name", "wnon_cv_earlyLOH", "q_earlyLOH")
    
    noLOH <- read.table(paste0("mount/camp/lab/swantonc/working/Saioa/TCGA_data_MC3_process/dNdS/",cancers.tokeep.withnew[c],"/",cancers.tokeep.withnew[c],"_allGD_MLEpergene_all_genes_all_noLOH.txt"), header=T)
    noLOH.sub <- select(noLOH, gene_name, wnon_cv, qtrunc_cv)
    names(noLOH.sub) <- c("gene_name", "wnon_cv_noLOH", "q_noLOH")
    
    allLOH <- read.table(paste0("mount/camp/lab/swantonc/working/Saioa/TCGA_data_MC3_process/dNdS/",cancers.tokeep.withnew[c],"/",cancers.tokeep.withnew[c],"_allGD_MLEpergene_all_genes_all_allLOH.txt"), header=T)
    allLOH.sub <- select(allLOH, gene_name, wnon_cv, qtrunc_cv)
    names(noLOH.sub) <- c("gene_name", "wnon_cv_all", "q_all")
    
    truncmutations <- merge(earlyLOH.sub, noLOH.sub, all=TRUE)
    truncmutations <- merge(truncmutations, noLOH.sub, all=TRUE)
    truncmutations.significant <- filter(truncmutations, gene_name%in%genestokeep)
    

    truncmutations.significant <- truncmutations.significant[order(truncmutations.significant$q_earlyLOH, truncmutations.significant$q_noLOH, truncmutations.significant$q_all),]

    names(truncmutations.significant) <- c("gene", paste0(cancers.tokeep.withnew[c],"_cv_earlyLOH"), paste0(cancers.tokeep.withnew[c],"_q_earlyLOH"), paste0(cancers.tokeep.withnew[c],"_cv_noLOH"), paste0(cancers.tokeep.withnew[c],"_q_noLOH"), paste0(cancers.tokeep.withnew[c],"_cv_all"), paste0(cancers.tokeep.withnew[c],"_q_all"))
    
    if (a==1) truncmutations.merged<-truncmutations.significant
    else (truncmutations.merged <- merge(truncmutations.merged, truncmutations.significant, all=TRUE))
    
    a=a+1
    
  } , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


cols = 2:ncol(truncmutations.merged)
cols.q = seq(3,ncol(truncmutations.merged),2)
cols.cv=seq(2,ncol(truncmutations.merged),2)

truncmutations.merged[,cols] = apply(truncmutations.merged[,cols], 2, function(x) as.numeric(as.character(x)));

### order by count: how many times is significant

truncmutations.merged$count = apply(truncmutations.merged[,cols.q], 1, function(x) length(which(x<0.05)))

### order by sum

truncmutations.significant.mat <- as.matrix(truncmutations.merged[,c(cols.cv, ncol(truncmutations.merged))])
rownames(truncmutations.significant.mat) <- truncmutations.merged[,1]
truncmutations.significant.mat<-truncmutations.significant.mat[order(truncmutations.significant.mat[,ncol(truncmutations.significant.mat)]),]

truncmutations.significant.mat<-truncmutations.significant.mat[,-ncol(truncmutations.significant.mat)]

##### remove the "all" column
cols_to_remove <- seq(3,ncol(truncmutations.significant.mat),3)
truncmutations.significant.mat <- truncmutations.significant.mat[,-cols_to_remove]
grep.vec=NULL

reverse.cancers.tokeep.withnew <- rev(cancers.tokeep.withnew)
for (ci in 1:length(reverse.cancers.tokeep.withnew)){
grep.ci <- grep(reverse.cancers.tokeep.withnew[ci], colnames(truncmutations.significant.mat))
grep.vec <- c(grep.ci, grep.vec)
  
}
###

truncmutations.significant.mat <- truncmutations.significant.mat[,grep.vec]

color.totalgenes=rownames(truncmutations.significant.mat)
color.totalgenes[!color.totalgenes%in%cancergenes]="red"
color.totalgenes[color.totalgenes%in%cancergenes]="black"

truncmutations.significant.mat.raw <- truncmutations.significant.mat
truncmutations.significant.mat <- truncmutations.significant.mat.raw

####### NOW color by gene specific in cancer type

truncmutations.significant.mat[truncmutations.significant.mat == 0] <- 0.0000000000001




pdf("Fig5B.pdf", width = 8, height = 5)


heatmap.mat = t(log(abs(truncmutations.significant.mat)))

quan.vals = c(0,1)


layout(matrix(c(1,1,1,1,1,1,1,2,
                3,3,3,3,3,3,3,4,
                3,3,3,3,3,3,3,4,
                3,3,3,4,3,3,3,4), nrow = 4, ncol = 8, byrow = TRUE))


#### plot 1 gene counts

par(mar=c(0.5,4,1,0))
x<-barplot(t(significant.count), col=c("grey","green"), ylim=c(0,110), xlab = "", las=1)
y<-as.numeric(t(significant.count[,1]))
text(x,y+20,labels=significant.count[,2])

#### empty plot 
plot(0,type='n',axes=FALSE,ann=FALSE)

cap <- 3.92
heatmap.mat[heatmap.mat>cap] <- cap
colscale=c(0, cap)

some.colors2<-colorRampPalette(c("white","blue4"))(500000)
par(mar=c(6,6,3,2))
image(1:dim(heatmap.mat)[1],1:dim(heatmap.mat)[2],heatmap.mat,xaxt="n",yaxt="n",xlab="",ylab="",col=some.colors2,zlim=colscale)


mtext(rownames(truncmutations.significant.mat), side=2, line=1, adj=1, at=seq(1,length(color.totalgenes),1),las=1,cex=0.5, col=color.totalgenes, font=2)
axis(1,seq(1,ncol(truncmutations.significant.mat),1),labels=rep(c("LOH","noLOH"), length(cancers.tokeep.withnew)),cex.axis=0.7,las=2, font=2)


### add stars to cancer-specific known cosmic genes

for (j in 1:ncol(truncmutations.significant.mat)){
  typecancer <- strsplit(colnames(truncmutations.significant.mat)[j], "_")[[1]][1]
  
  cancergenes_specific <-unlist(strsplit(significant.count[match(typecancer, rownames(significant.count)),3],";"))
  
  for (r in 1:nrow(truncmutations.significant.mat)){
    if (!rownames(truncmutations.significant.mat)[r]%in%cancergenes_specific) truncmutations.significant.mat[r,j] = truncmutations.significant.mat[r,j]*-1
  }
}


 truncmutations.significant.mat[is.na(truncmutations.significant.mat)] <- 0

 for (c in seq(1,ncol(truncmutations.significant.mat),2)){
 for (r in 1:nrow(truncmutations.significant.mat)){
   if (truncmutations.significant.mat[r,c]>0)  text(c,r, labels="*", col="orange", cex=2, font=2)
 }
 }

 #######
 
abline(v=seq(0.5,ncol(truncmutations.significant.mat),2))

## for legend
par(las=2,mgp=c(0,1,0),mar=c(0,1,10,0.5),fig=c(0.89,0.94,0.05,0.95), new=T)
colindex<-t(matrix(seq(colscale[1],colscale[2],length.out=100),ncol=1,nrow=100))
image(0,1:100,colindex,xaxt="n",yaxt="n",xlab="",ylab="",col=some.colors2,zlim=colscale)
box(lty=1)
labels.vec=round(signif(min(colindex)+(max(colindex)-min(colindex))*seq(0,1,length.out=6),2), digits = 2)

labels.vec <- c("1","2","5","11","22",">50")
axis(4,at=seq(1,100,length.out=6),labels=labels.vec,las=2,tick=FALSE,cex.axis=2,line=-0.8)

dev.off()

