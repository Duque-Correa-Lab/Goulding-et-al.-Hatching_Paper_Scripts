
# VST heatmaps
deseq2VST <- vst(ddsObj, blind = T)
plotPCA(deseq2VST, intgroup = c('timepoint'))
# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
res.dat <- column_to_rownames(res.dat, 'Target_ID')
sigGenes <- rownames(res.dat[res.dat$padj <= .05 & abs(res.dat$log2FoldChange) > 3,])
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
install.packages("ggdendro")
library(ggdendro)
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

# Re-factor samples for ggplot2
deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])

# Construct the heatmap. note that at this point we have only clustered the samples NOT the genes
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

# Combine the dendrogram and the heatmap
install.packages("gridExtra")
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

