# DEseq2 from Kallisto outputs
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("rhdf5")
BiocManager::install("tximportData")
BiocManager::install("devtools")    # only if devtools not yet installed
BiocManager::install("biomaRt")
BiocManager::install("limma")
BiocManager::install("AnnotationHub")
BiocManager::install('apeglm')
BiocManager::install('EnhancedVolcano')

library(biomaRt)
library(rhdf5)
library(tximport)
library(DESeq2)
library(tidyverse)
library(apeglm)
library(EnhancedVolcano)

# read in some sample information
setwd("/Volumes/duque_correa/csci_duque_correa/Lab_Members/Paul/Hatching_RNAseq")
sampleinfo <- read.csv("KallistoOutputs/SampleInfo.csv")
sampleinfo$timepoint <- as.factor(sampleinfo$timepoint)
sampleinfo$time_sample <- as.factor(sampleinfo$time_sample)

# Adding bonus groups based on PCA results for ggvenn
group_time <- c("6W_A", "8W_A", "8W_A", "6W_A", "6W_A","6W_A", "6W_B","8W_A", "8W_B","8W_A","8W_B")
group_time <- c("A", "C", "C", "A", "A", "A", "B","C", "D","C","D")
sampleinfo <- cbind.data.frame(sampleinfo, group_time)

sampleinfo$time_sample
files <- file.path("KallistoOutputs", sampleinfo$time_sample, "abundance.h5")
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

# because of alphabetical assignment our beta 0 / 1 are the wrong way around
sampleinfo <- mutate(sampleinfo, timepoint = fct_relevel(timepoint, "Egg_6w"))
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
results.6W_v_8W <- results(ddsObj, name = "timepoint_Egg_8w_vs_Egg_6w", alpha = 0.05)
sum(results.6W_v_8W$padj < 0.05, na.rm = TRUE)

# step 1 - estimate size factors ----
ddsObj <- estimateSizeFactors(ddsObj)
res <- results(ddsObj, name="timepoint_Egg_8w_vs_Egg_6w")
res <- results(ddsObj, contrast=c("timepoint", "Egg_6w","Egg_8w"))

# because we are interested in 6W vs 8W, we set 'coef=2'
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

write.csv(as.data.frame(resApe), 
          file="timepoint_hatching_results.csv")


##########################################################
############Global plots #################################
##########################################################

res.dat <- as.data.frame(resAsh)
res.dat <- rownames_to_column(res.dat, "Target_ID")

EnhancedVolcano(res.dat,
                lab = res.dat$Target_ID,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano - global',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)

res.Global.list <- filter(res.dat, log2FoldChange > 0.5 | log2FoldChange < 0.5, padj <  0.01)
write.csv(res.Global.list, "res_global.csv")

##########################################################
############GO terms curated #############################
##########################################################

mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)
genes <- getBM(mart = mart, filters = "species_id_1010", 
               value = "trmuriprjeb126", 
               attributes = c("wbps_gene_id", "wbps_transcript_id", "external_gene_id", "go_name_1006"))
genes <- dplyr::rename(genes, target_id = wbps_transcript_id,
                       ens_gene = wbps_gene_id, ext_gene = external_gene_id, go_name = go_name_1006)
res.dat <- as.data.frame(resAsh)
res.dat <- rownames_to_column(res.dat, "Target_ID")
res.GO <- full_join(genes, res.dat, by = join_by(target_id == Target_ID))
res.GO.list <- filter(res.GO, log2FoldChange > 1 | log2FoldChange < 1, padj <  0.01)
write.csv(res.GO.list, "res_GO.csv")
#res.GO.list <- filter(res.GO, log2FoldChange > 0.5 | log2FoldChange < 0.5, padj <  10e-32)
res.GO.list <- distinct(res.GO.list, target_id, log2FoldChange, padj)



res.GO.extracellular <- filter(res.GO, go_name == "extracellular region")
EnhancedVolcano(res.GO.extracellular,
                lab = res.GO.extracellular$target_id,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano - GO - Extracellular',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)


res.GO.serineinhibitor <- filter(res.GO, go_name == "serine-type endopeptidase inhibitor activity")
EnhancedVolcano(res.GO.serineinhibitor,
                lab = res.GO.serineinhibitor$target_id,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano - global',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)

res.GO.serineinhibitor <- filter(res.GO, go_name == "serine-type endopeptidase inhibitor activity")
EnhancedVolcano(res.GO.serineinhibitor,
                lab = res.GO.serineinhibitor$target_id,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano - global',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)


##########################################################
############Eichenburger secretome proteins ##############
##########################################################

res.dat <- as.data.frame(resAsh)
res.dat <- rownames_to_column(res.dat, "Target_ID")

genes_and_transcripts <- arrange(genes, target_id, ext_gene)
genes_and_transcripts <- distinct(genes_and_transcripts, target_id, ext_gene)

setwd("/Volumes/duque_correa/csci_duque_correa/Paul/Hatching_RNAseq")
library(readxl)
Other_names_Tmuris <- read_excel("MEROPS/Other_names_Tmuris.xlsx")
Other_names_Tmuris <- distinct(Other_names_Tmuris, New, Old)

Eichenberger_proteins <- read_csv("Eichenberger_proteins.csv")
Eichenberger_proteins <- Eichenberger_proteins %>%
  inner_join(Other_names_Tmuris, by = join_by(Secreted == Old), keep = FALSE)

Secreted <- left_join(Eichenberger_proteins, genes_and_transcripts, by = join_by(New == ext_gene))
Secreted <- Secreted[,2:4]

res.Secreted <-  full_join(Secreted, res.dat, by = join_by(target_id == Target_ID))
res.Secreted.minimal <-  left_join(Secreted, res.dat, by = join_by(target_id == Target_ID))
res.Secreted$Presence <- !is.na(res.Secreted$New) & res.Secreted$New != ""

EnhancedVolcano(res.Secreted.minimal,
                lab = res.Secreted.minimal$target_id,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Secreted - curated types',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  res.Secreted$Presence == TRUE,'red', 'grey')

names(keyvals)[keyvals == 'grey'] <- 'Absent'
names(keyvals)[keyvals == 'red'] <- 'Present'

EnhancedVolcano(res.Secreted,
                lab = res.Secreted$New,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = res.Secreted$New[which(names(keyvals) %in% c('Present'))],
                xlab = bquote(~Log[2]~ 'fold change'),
                #title = 'Custom colour over-ride',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 4,
                labSize = 1,
                shape = c(16, 16, 16, 16),
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'none',
                legendLabSize = 15,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                arrowheads = FALSE,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black')

res.Secreted.list <- filter(res.Secreted.minimal, log2FoldChange > 0.5 | log2FoldChange < 0.5, padj <  10e-32)
write.csv(res.Secreted.list, "Secreted_list.csv")


##########################################################
############MEROPS ANNOTED FROM 50 HELMINTHS #############
##########################################################

# MANUALLY ADDED MEROPS CLANS PROPERLY! NAMED MEROPS_updated
MEROPS <- read_csv("MEROPS_annotated.csv", col_types = cols(...1 = col_skip()))
MEROPS <- MEROPS %>% 
  left_join(tx2gene, by = join_by(Gene_ID == GENEID))

res.dat <- as.data.frame(resAsh)
res.dat <- rownames_to_column(res.dat, "Target_ID")
res.MEROPS.minimal <- left_join(MEROPS, res.dat, by = join_by(TXNAME == Target_ID))
res.MEROPS <- full_join(MEROPS, res.dat, by = join_by(TXNAME == Target_ID))
res.MEROPS$Presence <- !is.na(res.MEROPS$Type) & res.MEROPS$Type != ""

EnhancedVolcano(res.MEROPS.minimal,
                lab = res.MEROPS.minimal$Family,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Secreted - curated types',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  res.MEROPS$Presence == TRUE,'red', 'grey')

names(keyvals)[keyvals == 'grey'] <- 'Absent'
names(keyvals)[keyvals == 'red'] <- 'Present'

EnhancedVolcano(res.MEROPS,
                lab = res.MEROPS$TXNAME,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = res.MEROPS$TXNAME[which(names(keyvals) %in% c('Present'))],
                xlab = bquote(~Log[2]~ 'fold change'),
                #title = 'Custom colour over-ride',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 4,
                labSize = 1,
                shape = c(16, 16, 16, 16),
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'none',
                legendLabSize = 15,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                arrowheads = FALSE,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black')

res.MEROPS.list <- filter(res.MEROPS.minimal, log2FoldChange > 0.5 | log2FoldChange < 0.5, padj <  0.01)
write.csv(res.MEROPS.list, "MEROPS_list.csv")

