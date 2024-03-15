#########################################
### HATCHING DATA - GOULDING 2024 #######
#########################################
# Post Kallisto Output - importing data for DESeq2 analysis

#################
### PACKAGES ####
#################
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("rhdf5")
#BiocManager::install("tximportData")
#BiocManager::install("devtools")    # only if devtools not yet installed
#BiocManager::install("biomaRt")
#BiocManager::install("limma")
#BiocManager::install('apeglm')
#BiocManager::install("AnnotationHub")
#BiocManager::install('EnhancedVolcano')

library(biomaRt)
library(rhdf5)
library(tximport)
library(DESeq2)
library(tidyverse)
library(apeglm)
library(EnhancedVolcano)

####################
### DATA IMPORT ####
####################

setwd("/Users/pma37/Desktop/Hatching_RNAseq")
library(readr)
sampleinfo <- read_csv("s2c.csv")
sampleinfo$timepoint <- as.factor(sampleinfo$timepoint)
sampleinfo$time_sample <- as.factor(sampleinfo$time_sample)
sampleinfo$time_sample

files <- file.path("KallistoOutputsCombinedReps", sampleinfo$sample_number, "abundance.h5")
files <- set_names(files, sampleinfo$sample_number)

# To make the tx2gene here I borrowed from Sleuth by making a list of the gene_IDs and transcript IDs listed by Kallisto from biomart
mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)
tx2gene <- getBM(mart = mart, filters = "species_id_1010", 
                 value = "trmuriprjeb126", 
                 attributes = c("external_gene_id", "wbps_transcript_id"))
tx2gene <- dplyr::rename(tx2gene, TXNAME = wbps_transcript_id, GENEID = external_gene_id)

# import kallisto files
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, txOut = TRUE)


#######################
### FILTER SAMPLES ####
#######################

# prepare filtered counts >5 rowsum
rawCounts <- round(txi$counts, 0)
dim(rawCounts)
keep <- rowSums(rawCounts) > 5
table(keep, useNA = "always")
filtCounts <- rawCounts[keep,]
dim(filtCounts)

# Start exploration
summary(filtCounts)
boxplot(filtCounts, main= "Raw Counts", las = 2)
plot(rowMeans(filtCounts), rowSds(filtCounts),
     main="Raw Counts sd vs mean",
     xlim=c(0,10000),
     ylim=c(0,5000))

# VST version of normalised counts - all look good!
vst_counts <- vst(filtCounts)
boxplot(vst_counts,
        xlab="",
        ylab="VSTcounts",
        las = 2)
abline(h=median(vst_counts), col="blue")
plot(rowMeans(vst_counts), rowSds(vst_counts),
     main="VST Counts sd vs mean")

###################
### RUN DESEQ2 ####
###################

# set up the model 
simple.model <- as.formula(~ timepoint)
# because of alphabetical assignment our beta 0 / 1 are the wrong way around
sampleinfo <- mutate(sampleinfo, timepoint = fct_relevel(timepoint, "6W"))
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
results.6W_v_8W <- results(ddsObj, name = "timepoint_8W_vs_6W", alpha = 0.05)
sum(results.6W_v_8W$padj < 0.05, na.rm = TRUE)

# step 1 - estimate size factors ----
ddsObj <- estimateSizeFactors(ddsObj)
res <- results(ddsObj, name="timepoint_8W_vs_6W")
res <- results(ddsObj, contrast=c("timepoint", "6W","8W"))

# because we are interested in 6W vs 8W, we set 'coef=2'
resApe <- lfcShrink(ddsObj, coef=2, type="apeglm")
resNorm <- lfcShrink(ddsObj, coef=2, type="normal")
resAsh <- lfcShrink(ddsObj, coef=2, type="ashr")
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
DESeq2::plotMA(resApe, xlim=xlim, ylim=ylim, main="apeglm")
DESeq2::plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
DESeq2::plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

write.csv(as.data.frame(resAsh), 
          file="timepoint_hatching_results.csv")


#########################
### MAKE THE PCA PLOT ###
#########################

# PCA from VST heatmaps
deseq2VST <- vst(ddsObj, blind = T)
plotPCA(deseq2VST, intgroup = c('timepoint'))

###################################
### MAKE THE DENDRIGRAM CLUSTER ###
###################################

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
res.dat <- as.data.frame(resAsh)
sigGenes <- rownames(res.dat[res.dat$padj <= .05 & abs(res.dat$log2FoldChange) > 1,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

# Convert the VST counts to long format for ggplot2
library(reshape2)

# First compare wide vs long version
deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

head(deseq2VST_wide)
head(deseq2VST_long)

library(tidyverse)

# Now overwrite our original data frame with the long format
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

deseq2VST <- deseq2VST %>%
  mutate(
    timepoint = case_when(
      variable == 1 ~ '6W',
      variable == 2 ~ '8W',
      variable == 3 ~ '8W',
      variable == 4 ~ '6W',
      variable == 5 ~ '6W',
      variable == 6 ~ '6W',
      variable == 7 ~ '6W',
      variable == 8 ~ '8W',
      variable == 9 ~ '8W',
      variable == 10 ~ '8W',
      variable == 11 ~ '8W',
    )
  )

# Convert the significant genes back to a matrix for clustering
deseq2VSTMatrix <- dcast(deseq2VST, Gene ~ variable)
rownames(deseq2VSTMatrix) <- deseq2VSTMatrix$Gene
deseq2VSTMatrix$Gene <- NULL

# Compute a distance calculation on both dimensions of the matrix
distanceGene <- dist(deseq2VSTMatrix)
distanceSample <- dist(t(deseq2VSTMatrix))

# Cluster based on the distance calculations
clusterGene <- hclust(distanceGene, method="average")
clusterSample <- hclust(distanceSample, method="average")

# Construct a dendogram for samples
#install.packages("ggdendro")
library(ggdendro)
library(viridis)
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

# Re-factor samples for ggplot2
deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])

# Construct the heatmap. note that at this point we have only clustered the samples NOT the genes
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

# Combine the dendrogram and the heatmap
#install.packages("gridExtra")
library(gridExtra)
grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))

# Load in libraries necessary for modifying plots
#install.packages("gtable")
library(gtable)
library(grid)

# Modify the ggplot objects
library(ggplot2)
sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0))
heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))

# Convert both grid based objects to grobs
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)

# Check the widths of each grob
sampleDendrogramGrob$widths
heatmapGrob$widths

# Add in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

# Make sure every width between the two grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)

# Arrange the grobs into a plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))

# Draw the plot
grid.draw(finalGrob)

# Re-order the sample data to match the clustering we did
sampleinfo$timepoint <- factor(sampleinfo$timepoint, levels = c("6W", "8W"))
samples <- ggplot(sampleinfo, aes(x=timepoint, y=1, fill=timepoint)) + geom_tile() + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)) + scale_fill_manual(values=c("#E5954D", "#4194A6")) + theme_void()

# Convert the clinical plot to a grob
sampleGrob <- ggplotGrob(samples)

# Make sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleGrob$widths <- as.list(maxWidth)

finalGrob <- arrangeGrob(sampleDendrogramGrob, sampleGrob, heatmapGrob, ncol=1, heights=c(2,1,5))
grid.draw(finalGrob)

################################################################################
################# Step 1: create dendrogram for genes ##########################
library(ggplot2)
library(grid)
library(gtable)
library(viridis)
library(gridExtra)
library(ggdendro)

# we already did the clustering for genes in the tutorial, get the data to make a dendrogram with ggplot
geneModel <- as.dendrogram(clusterGene)
geneDendrogramData <- segment(dendro_data(geneModel, type = "rectangle"))

# construct the dendrogram in ggplot
geneDendrogram <- ggplot(geneDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() + scale_y_reverse(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) + theme_dendro()

################################################################################
################# Step 2: Re-arrange the heatmap cells #########################

# re-factor genes for ggplot2
deseq2VST$Gene <- factor(deseq2VST$Gene, levels=clusterGene$labels[clusterGene$order])

# recreate the heatmap with this new factoring
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())

################################################################################
################# Step 3: convert to everything to grobs #######################

# note! before this step as mentioned you might need to alter the expand parameters in the plot scales for all the plots we do that here

# convert the heatmap to a grob
heatmapGrob <- ggplotGrob(heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)))

# convert the dendrogram to a grob
# note! we flipped the axis above so the x-axis is now displayed as what we would think of as the y-axis
geneDendrogramGrob <- ggplotGrob(geneDendrogram + scale_x_discrete(expand=c(0, 0)))

# we already have a sample Dendrogram, but here it is again
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0)))

################################################################################
######### Step 4: align the gene dendrograms to match the heatmap ##############

# check that the both the heatmap and gene dendrogram have the same number of vertical elements
length(heatmapGrob$heights) == length(geneDendrogramGrob$heights)

# make sure every height between the two grobs is the same
maxHeight <- unit.pmax(geneDendrogramGrob$heights, heatmapGrob$heights)
geneDendrogramGrob$heights <- as.list(maxHeight)
heatmapGrob$heights <- as.list(maxHeight)

################################################################################
# Step 4b: we have a new heatmap so we need to re-align the horizontal elements #

# repeat the steps in the tutorial

# check the widths of each grob
sampleDendrogramGrob$widths
heatmapGrob$widths
sampleGrob$widths

# add in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

# make sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleGrob$widths <- as.list(maxWidth)

################################################################################
############### Step 5: create a blank panel ###################################

# we can use grid graphics for this
blankPanel <- grid.rect(gp=gpar(col="white"))

################################################################################
############### Step 6: Arrange the final result ###############################

# arrange all the plots together
finalGrob_v2 <- arrangeGrob(blankPanel, sampleDendrogramGrob, geneDendrogramGrob, heatmapGrob, ncol=2, nrow=2, widths=c(1,5), heights=c(2,6))

# draw the final result
grid.draw(finalGrob_v2)


##########################################################
############Global plots #################################
##########################################################

res.dat <- as.data.frame(resAsh)
res.dat <- rownames_to_column(res.dat, "Target_ID")

EnhancedVolcano(res.dat,
                lab = res.dat$Target_ID,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Volcano - global',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)

res.Global.list <- filter(res.dat, log2FoldChange > 0.1 | log2FoldChange < 0.1, padj <  0.05)
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
res.GO.list <- distinct(res.GO.list, target_id, log2FoldChange, padj)

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
                pCutoff = 0.05,
                FCcutoff = 1,
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
                y = 'padj',
                selectLab = res.Secreted$New[which(names(keyvals) %in% c('Present'))],
                xlab = bquote(~Log[2]~ 'fold change'),
                #title = 'Custom colour over-ride',
                pCutoff = 0.05,
                FCcutoff = 1,
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
