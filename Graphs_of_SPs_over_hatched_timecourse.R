# Checking genes of interest across post-hatch dataset

# VST heatmaps
deseq2VST <- vst(ddsObj, blind = T)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
#res.dat <- as.data.frame(resAsh)
#res.dat <- column_to_rownames(res.dat, 'Target_ID')
#sigGenes <- rownames(res.dat[res.dat$padj <= .05 & abs(res.dat$log2FoldChange) > 1,])
#deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

# Convert the VST counts to long format for ggplot2
library(reshape2)

# First compare wide vs long version
deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

library(tidyverse)

VST_genes <- full_join(sampleinfo, deseq2VST_long, by = join_by("SampleName" == "variable")) #Add timepoint groups
VST_Ave <- aggregate(VST_genes$value, by = list(VST_genes$timepoint, VST_genes$Gene), FUN=mean)
colnames(VST_Ave)[3] <- 'VST_value'
VST_Ave_SP <- VST_Ave %>% mutate(SP = case_when(Group.2 %in% c("TMUE_3000011680.1", 
                                                         "TMUE_2000006651.1", 
                                                         "TMUE_2000006858.1", 
                                                         "TMUE_3000013335.1", 
                                                         "TMUE_3000013624.1") ~ '8w_upreg'))
colnames(VST_Ave_SP)
VST_Ave_SP$Group.1 <- factor(VST_Ave_SP$Group.1, levels = c("L1_Unhatched", "L1_Hatched", "L1_3hr", "L1_24hr", "L2"))
ggplot(VST_Ave_SP, aes(x=Group.1, y=VST_value, group=Group.2)) +
  geom_line(aes(color = SP))
  geom_point(aes(shape=SP, color = SP))

  
  
# Same but with TPMs
All_TSVs_raw <- rownames_to_column(All_TSVs_raw, "Gene")
TPM_long <- melt(All_TSVs_raw, id.vars=c("Gene"))
TPM_genes <- full_join(sampleinfo, TPM_long, by = join_by("time_sample" == "variable")) #Add timepoint groups
TPM_Ave <- aggregate(TPM_genes$value, by = list(TPM_genes$timepoint, TPM_genes$Gene), FUN=mean)
colnames(TPM_Ave)[3] <- 'VST_value'
TPM_Ave_SP <- TPM_Ave %>% mutate(SP = case_when(Group.2 %in% c("TMUE_3000011680.1", 
                                                               "TMUE_2000006651.1", 
                                                               "TMUE_2000006858.1", 
                                                               "TMUE_3000013335.1", 
                                                               "TMUE_3000013624.1") ~ '8w_upreg'))
colnames(TPM_Ave_SP)
TPM_Ave_SP$Group.1 <- factor(TPM_Ave_SP$Group.1, levels = c("L1_Unhatched", "L1_Hatched", "L1_3hr", "L1_24hr", "L2"))
ggplot(TPM_Ave_SP, aes(x=Group.1, y=VST_value, group=Group.2)) +
  geom_line(aes(color = SP))
geom_point(aes(shape=SP, color = SP))

subset(TPM_Ave_SP)
TPM_Ave_SP_Only <- subset(TPM_Ave_SP, SP == "8w_upreg")
ggplot(TPM_Ave_SP_Only, aes(x=Group.1, y=VST_value, group=Group.2)) +
  geom_line(aes(color = Group.2))
