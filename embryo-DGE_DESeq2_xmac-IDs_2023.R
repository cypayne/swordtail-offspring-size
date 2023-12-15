## Project: Offspring size evolution
## Payne et al
##
## DESeq2 DGE of X. malinche, X. birchmanni, and X. cortezi late-stage embryos,
## pseudoaligned against Xbirchmanni pseudotranscriptome generated with X. maculatus genome
## RNAseq collected in 2023 from nonpregnant and pregnant mothers --> see metadata in ovary_embryo_rnaseq_samples_FebAug23_combined.txt
##
## add annotations to outfile with:
## Scripts/refs/add-gtf-annots-col.py Scripts/refs/Xiphophorus_maculatus_LG.Xipmac4.4.2.81_unique-tx_w-mito.gtf Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/ovary-xmacID-combined2023_dge_lfc-shr_all.csv
##
## cyp XII-2023

## load libraries
library(tximportData)
library(tximport)
library(GenomicFeatures)
library(readr)
library(DESeq2)
library(PCAtools)
library(ggplot2)

## read in kallisto sample files
dir <- "Data/embryo_ovary_dge_combined_2023/ovary_embryo_kallisto_output_FebAug23_combined/kallisto_posttrim_xbir-allhaps-xmac-txtome"
outdir <- "Data/embryo_ovary_dge_combined_2023/embryo_xmac-gtf_dge/"
of_header   <- "embryo-xmacID-combined2023_dge"
samples <- read.table("Data/embryo_ovary_dge_combined_2023/ovary_embryo_kallisto_output_FebAug23_combined/ovary_embryo_rnaseq_samples_FebAug23_combined.txt", header = TRUE)

## subset samples
tissue <- "embryo"
samples <- samples[samples$tissue == tissue,]
samples

# grab file paths to those samples
files <- file.path(dir, paste(samples$file_basename,"kallisto",sep="_"), "abundance.h5")
names(files) <- paste0(samples$sample_ID)
files

## load transcript annotations
txdb <- makeTxDbFromGFF(file="Scripts/refs/Xiphophorus_maculatus_LG.Xipmac4.4.2.81_unique-tx_w-mito.gtf", format="gtf")

# create tx2gene table, match transcript name to geneid
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# get gene-level count information from kallisto counts
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)
head(txi$counts)


#################### DGE analysis ####################

## instantiate the DESeqDataSet with the design formula: ~ species + origin + collection
dds     <- DESeqDataSetFromTximport(txi,
                                    colData = samples,
                                    design = ~ species + origin + collection)

## set Xbirchmanni as reference level
dds$species <- relevel(dds$species, ref = "Xbirchmanni")

## determine the best fit for data
## by plotting or quantitatively with median absolute residual (smaller = better)
# parametric: 0.11318
dds <- DESeq(dds, fitType="parametric")
plotDispEsts(dds, main="Dispersion plot with parametric fit")
residual <- mcols(dds)$dispGeneEst - mcols(dds)$dispFit
absres<-abs(residual)
summary(absres)
# versus local: 0.069
dds <- DESeq(dds, fitType="local")
plotDispEsts(dds, main="Dispersion plot with local fit")
residual <- mcols(dds)$dispGeneEst - mcols(dds)$dispFit
absres<-abs(residual)
summary(absres)

# chose local fit, run DESeq
dds <- DESeq(dds, fitType="local")
resultsNames(dds)

# save dds object
saveRDS(dds, paste0(outdir,of_header,"_dds.rds"))

## save vst object for heatmap
vst <- vst(dds,blind=FALSE)
saveRDS(vst, paste0(outdir,of_header,"_vst.rds"))

## load dds
dds <- readRDS("Data/embryo_ovary_dge_combined_2023/embryo_xmac-gtf_dge/embryo-xmacID-combined2023_dge_dds.rds")

# compare expression between species (all "late" stage)
res.malvbirch_late <- lfcShrink(dds, coef="species_Xmalinche_vs_Xbirchmanni", type="ashr")
res.corvbirch_late <- lfcShrink(dds, coef="species_Xcortezi_vs_Xbirchmanni", type="ashr")
res.malvcor_late <- lfcShrink(dds, contrast=list("species_Xmalinche_vs_Xbirchmanni","species_Xcortezi_vs_Xbirchmanni"), type="ashr")


#################### START PRODUCE OUT TABLES ####################

## output results
# collect all results DESeq2 object variables starting with "res."res_list <- lapply(ls(pattern = "^res\\."), get)
names(res_list)<- (ls(pattern = "^res\\."))
names(res_list)
# loop through all results objects, output to csv files
for (i in 1:length(res_list)) {

  # create output file with all values for current results object
  res <- res_list[[i]]
  # reorder results object by adj p-val (increasing)
  res_reorder <- res[order(res$padj), ]
  # merge results values with the normalized counts data for all samples
  res_out <- merge(as.data.frame(res_reorder), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  # rename first column
  names(res_out)[1] <- "Gene"
  # output this to unique file
  write.csv(res_out, file = paste0(outdir,of_header,"_lfc-shr_",names(res_list[i]), ".csv", sep = ""))

  # create output file with log-fold change and adjusted p-values from all created res. objects
  # grab only LFC and padj columns from results object
  res_slim <- data.frame(res_out$Gene,res_out$log2FoldChange,res_out$padj)
  names(res_slim) <- c( "Gene", paste("LFC_",names(res_list[i]),sep=""), paste("padj_",names(res_list[i]),sep="") )

  # if you're on the first object, initialize resdata.shr
  if (i == 1) {
    resdata.shr <- res_slim
  }
  # otherwise merge new results values with master table of LFC and padj values
  else {
    resdata.shr <- merge(resdata.shr,as.data.frame(res_slim), by="Gene",sort=FALSE)
  }
}

# merge and output results values with the normalized counts data for all samples
#row.names(resdata.shr) <- resdata.shr$Gene
resdata.shr <- merge(resdata.shr, as.data.frame(counts(dds, normalized=TRUE)), by.x="Gene", by.y="row.names",sort=FALSE)
write.csv(resdata.shr, file=paste(outdir,of_header,"_lfc-shr_all.csv",sep=""))

dim(res_out) # 19176


### Summary of expression patterns
embryo_dge <- read.csv("Data/embryo_ovary_dge_combined_2023/embryo_xmac-gtf_dge/embryo-xmacID-combined2023_dge_lfc-shr_all.csv_with-annots.csv",header=T)
embryo_dge_expressed_all <- subset(embryo_dge, !is.na(padj_res.malvbirch_late) & !is.na(padj_res.malvcor_late) & !is.na(padj_res.corvbirch_late)) # 854 not expressed in all species, 18322 expressed
embryo_dge_expressed_one <- subset(embryo_dge, !is.na(padj_res.malvbirch_late) | !is.na(padj_res.malvcor_late) | !is.na(padj_res.corvbirch_late)) # 854 not expressed in any species, 18322 expressed in at least one

embryo_padj0.05_malvbir_late <- subset(embryo_dge, padj_res.malvbirch_late < 0.05) # 315
embryo_padj0.05_malvcor_late <- subset(embryo_dge, padj_res.malvcor_late < 0.05) # 4630
embryo_padj0.05_corvbir_late <- subset(embryo_dge, padj_res.corvbirch_late < 0.05) # 5237
length(intersect(embryo_padj0.05_malvcor_late$Gene,embryo_padj0.05_corvbir_late$Gene)) # 3741 shared for the xcor comparisons
