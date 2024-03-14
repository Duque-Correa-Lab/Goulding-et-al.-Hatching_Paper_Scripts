# Name the files to read in - .tsv
setwd("~/Desktop/Early")
TSV_names <- paste0(sampleinfo$SampleID, ".tsv")
TSVs <- file.path("Clust_Analysis", sampleinfo$SampleID, TSV_names)
TSVs <- set_names(TSVs, sampleinfo$SampleID)
txi_tsv <- tximport(TSVs, type = "kallisto", tx2gene = tx2gene, txOut = TRUE)
txi_tsv_table <- as.data.frame(txi_tsv)
?tximport
# prepare raw counts
rawCounts <- round(txi_tsv$counts, 0)
dim(rawCounts)
keep <- rowSums(rawCounts) > 5
table(keep, useNA = "always")
filtCounts <- rawCounts[keep,]
dim(filtCounts)

All_TSVs_raw <- as.data.frame(txi_tsv$counts)
All_TSVs <- as.data.frame(filtCounts)
write.csv(All_TSVs_raw, "TPMs_All.csv")
sampleinfo$SampleID


# Same as above but for est_counts
setwd("~/Desktop/Early")
TSV_names <- paste0(sampleinfo$time_sample, ".tsv")
TSVs <- file.path("Clust_Analysis","Raw_data", sampleinfo$time_sample, TSV_names)
TSVs <- set_names(TSVs, sampleinfo$time_sample)
txi_tsv <- tximport(TSVs, type = "kallisto", tx2gene = tx2gene, txOut = TRUE)

# prepare raw counts
rawCounts <- round(txi_tsv$abundance, 0)
dim(rawCounts)
keep <- rowSums(rawCounts) > 5
table(keep, useNA = "always")
filtCounts <- rawCounts[keep,]
dim(filtCounts)

All_TSVs_raw <- as.data.frame(txi_tsv$counts)
All_TSVs <- as.data.frame(filtCounts)
write.csv(All_TSVs_raw, "TPMs_All_Abundance.csv")