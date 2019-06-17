library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(GenomicRanges)
library(ggbeeswarm)
library(ggpubr)
library(cowplot)
library(knitr)


#Loading CN Data
cn_summaries_info_patient_df <- readRDS("~/LOHsegments.RDS")
patients_to_include = unique(cn_summaries_info_patient_df$patient_id)

sample_info <- readRDS("~/sample_info.RDS")
patient_info <- readRDS("~/patient_info.RDS")

cn_df = cn_summaries_info_patient_df %>%
  separate_rows(regions_cn,sep=";")

cn_df_LOH = cn_df %>%
  ungroup() %>%
  filter(regions_cn%in%c("Amp_LOH","DeepLoss_LOH","Gain_LOH","Loss_LOH","Neutral_LOH")) %>%
  select(patient_id,segName,segsize_gn_percent,CN_Subclonal) %>%
  distinct()

total_amount_of_loh = cn_df_LOH %>%
  group_by(patient_id) %>%
  summarise(total_loh=sum(segsize_gn_percent))

loh_subclonal_clonal_diff = cn_df_LOH %>%
  group_by(patient_id,CN_Subclonal) %>%
  summarise(segsize_gn_percent=sum(segsize_gn_percent)) %>%
  spread(CN_Subclonal,segsize_gn_percent,fill=0) %>%
  mutate(clonal_subclonal_diff=Clonal-Subclonal) %>%
  left_join(.,patient_info,by="patient_id")

subclonal_clonal_diff_order = loh_subclonal_clonal_diff$patient_id[order(loh_subclonal_clonal_diff$clonal_subclonal_diff)]


pdf("~/Figure2D.pdf", height = 6, width = 12)

loh_subclonal_clonal_diff %>%
  ungroup() %>%
  mutate(patient_id=factor(patient_id,levels=subclonal_clonal_diff_order),
         clonal_greater=ifelse(clonal_subclonal_diff>0,"Clonal","Subclonal")) %>%
  ggplot(aes(x=patient_id,y=clonal_subclonal_diff,fill=clonal_greater)) +
  geom_bar(stat="identity") +
  facet_wrap(~GD_predicted,scales="free_x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  scale_fill_manual(values=c("darkolivegreen4","deepskyblue1"))


loh_subclonal_clonal_diff %>%
  select(patient_id,Clonal,Subclonal,GD_predicted) %>%
  melt(id.vars=c("patient_id","GD_predicted")) %>%
  ggplot(aes(x=GD_predicted,y=value,fill=variable)) +
  geom_boxplot() +
  stat_compare_means(method="t.test",paired=TRUE, label="p.signif")+
  scale_fill_manual(values=c("darkolivegreen4","deepskyblue1")) +
  ylab("% Genome LOH")

###calculate mean
loh_subclonal_clonal_diff_GD <- loh_subclonal_clonal_diff[which(loh_subclonal_clonal_diff$GD_predicted=="GD"),]
meanclonal <- mean(loh_subclonal_clonal_diff_GD$Clonal)
meansubclonal <- mean(loh_subclonal_clonal_diff_GD$Subclonal)

percent_loh_clonal_subclonal_plot = loh_subclonal_clonal_diff %>%
  filter(GD_predicted!="SubclonalGD") %>%
  select(patient_id,Clonal,Subclonal,GD_predicted) %>%
  melt(id.vars=c("patient_id","GD_predicted")) %>%
  ggplot(aes(x=GD_predicted,y=value,fill=GD_predicted)) +
  geom_boxplot() +
  stat_compare_means(method="t.test",paired=FALSE,comparisons = list(c("GD","nGD")))+
  ylab("% Genome LOH") +
  facet_wrap(~variable)

percent_loh_added_up_plot = loh_subclonal_clonal_diff %>%
  mutate(percent_LOH_for_comparison=ifelse(GD_predicted=="GD",Clonal,Clonal+Subclonal)) %>%
  mutate(GD_Label=ifelse(GD_predicted=="GD","GD Clonal","nGD Clonal+Subclonal")) %>%
  filter(GD_predicted%in%c("GD","nGD")) %>%
  ggplot(aes(x=GD_Label,y=percent_LOH_for_comparison,fill=GD_predicted)) +
  geom_boxplot() +
  stat_compare_means(method="t.test",paired=FALSE) +
  ylab("% Genome LOH") +
  ylim(0,1)


plot_grid(percent_loh_clonal_subclonal_plot,percent_loh_added_up_plot,rel_widths = c(2,1))

dev.off()