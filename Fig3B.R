library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

regions <- c("D_HCT001_BS",
             "D_HCT001_SU_T1-R1",
             "D_HCT001_SU_T10-R1",            
             "D_HCT001_SU_T11-R1","D_HCT001_SU_T12-R1",                 
             "D_HCT001_SU_T13-R1","D_HCT001_SU_T14-R1",                 
             "D_HCT001_SU_T2-R1","D_HCT001_SU_T3-R1",                  
             "D_HCT001_SU_T4-R1","D_HCT001_SU_T5-R1",                  
             "D_HCT001_SU_T6-R1","D_HCT001_SU_T7-R1",                  
             "D_HCT001_SU_T8-R1","D_HCT001_SU_T9-R1")

mut.table <- read.table("~/HCT001.Exome.SNV.xls"
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
  
TCGA.earlyLate <- read.table(paste0("~/PyClone_singleRegion_phylo/",regions[c],"_earlyORlate.txt"))
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
  
  TCGA.earlyLate <- read.table(paste0("~/PyClone_singleRegion_phylo/",regions[c],"_earlyORlate.txt"))
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
