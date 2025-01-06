# DEseq2 from Kallisto outputs
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(biomaRt)
library(rhdf5)
library(tximport)
library(DESeq2)
library(tidyverse)
library(apeglm)
library(EnhancedVolcano)

#Read in locally (when drive is down)
setwd("~/Desktop/Early")
sampleinfo <- read.csv("KallistoOutputs/SampleInfo3.csv") #Has removed samples with low VST ave

# read in some sample information
setwd("/Volumes/duque_correa/csci_duque_correa/Lab_Members/Paul/Hatching_RNAseq")
sampleinfo <- read.csv("KallistoOutputs/SampleInfo3.csv") #Has removed samples with low VST ave
sampleinfo$timepoint <- as.factor(sampleinfo$timepoint)
sampleinfo$time_sample

# Name the files to read in
files <- file.path("KallistoOutputs/Early", sampleinfo$time_sample, "abundance.h5")
files <- set_names(files, sampleinfo$SampleName)

# To make the tx2gene here I borrowed from Sleuth by making a list of the gene_IDs and transcript IDs listed by Kallisto from biomart
mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)
tx2gene <- getBM(mart = mart, filters = "species_id_1010", 
                 value = "trmuriprjeb126", 
                 attributes = c("external_gene_id", "wbps_transcript_id"))
tx2gene <- dplyr::rename(tx2gene, TXNAME = wbps_transcript_id, GENEID = external_gene_id)

# import kallisto files
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, txOut = TRUE)
saveRDS(txi, "Hatching_txi.rds")

# prepare raw counts
rawCounts <- round(txi$counts, 0)
dim(rawCounts)

keep <- rowSums(rawCounts) > 5
table(keep, useNA = "always")

filtCounts <- rawCounts[keep,]
dim(filtCounts)

# set up the model 
simple.model <- as.formula(~ timepoint)

sampleinfo <- mutate(sampleinfo, timepoint = fct_relevel(timepoint, "L1_Hatched")) # This is updated depending on the test
model.matrix(simple.model, data = sampleinfo)
cbind(sampleinfo, model.matrix(simple.model, data = sampleinfo))

# Import the things into DESEq2 ----
ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData  = sampleinfo,
                                       design = ~ timepoint)
# subset the low count genes ----
keep <- rowSums(counts(ddsObj.raw)) > 5 
ddsObj.filt <- ddsObj.raw[keep ,  ]
# DESeq2 Manual ----
ddsObj <- DESeq(ddsObj.filt)
results.simple <- results(ddsObj, alpha=0.05)
# Sig DEGs 
sum(results.simple$padj < 0.05, na.rm = TRUE)

# Count DEGs for a specific output - in a simple model you can just sum the results
resultsNames(ddsObj)
results.Unhatch_v_Hatch <- results(ddsObj, name = "timepoint_L1_Hatched_vs_L1_Unhatched", alpha = 0.05)
sum(results.Unhatch_v_Hatch$padj < 0.05, na.rm = TRUE)

# step 1 - estimate size factors ----
ddsObj <- estimateSizeFactors(ddsObj)
res <- results(ddsObj, name="timepoint_L1_Hatched_vs_L1_Unhatched")
res <- results(ddsObj, contrast=c("timepoint", "L1_Hatched","L1_3hr"))

# because we are interested in Unhatched vs Hatched, we set 'coef=2'
resApe <- lfcShrink(ddsObj, coef=2, type="apeglm")
resNorm <- lfcShrink(ddsObj, coef=2, type="normal")
resAsh <- lfcShrink(ddsObj, coef=2, type="ashr")
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
DESeq2::plotMA(resApe, xlim=xlim, ylim=ylim, main="apeglm")
DESeq2::plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
DESeq2::plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

#idx <- identify(res$baseMean, res$log2FoldChange)
#rownames(res)[idx]

write.csv(as.data.frame(resAsh), 
          file="timepoint_hatching_results.csv")

###################################################
########## Venn Diagram of Hatched ################
###################################################

# Here we do pairwise comparisons for gene lists, then make a venn of the gene lists.

# Count DEGs for a specific output - in a simple model you can just sum the results
resultsNames(ddsObj)

results.Unhatch_v_Hatch <- results(ddsObj, name = "timepoint_L1_Hatched_vs_L1_Unhatched", alpha = 0.05)
sum(results.Unhatch_v_Hatch$padj < 0.05, na.rm = TRUE)
write.csv(as.data.frame(results.Unhatch_v_Hatch), file="results.Unhatch_v_Hatch.csv")

results.Hatch_v_3hr <- results(ddsObj, name = "timepoint_L1_3hr_vs_L1_Hatched", alpha = 0.05)
sum(results.Hatch_v_3hr$padj < 0.05, na.rm = TRUE)
write.csv(as.data.frame(results.Hatch_v_3hr), file="results.3hr_v_Hatch.csv")
# Looked at these in excel and filtered
Hatched_Upreg_3hr <- read_excel("Hatched_Upreg_v3hr.xlsx")
Hatched_Upreg_3hr <- left_join(Hatched_Upreg_3hr, genes, by = join_by(Target_ID == target_id))
write.csv(Hatched_Upreg_3hr, file="Hatched_Upred_v3hr_GO.csv")

X3hr_upreg_vHatched <- left_join(X3hr_Upreg_vHatched, genes, by = join_by(Target_ID == target_id))
write.csv(X3hr_upreg_vHatched, file="3hr_Upreg_vHatched_GO.csv")