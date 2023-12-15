## Project: Offspring size evolution
## Payne et al
##
## DESeq2 DGE of X. malinche, X. birchmanni, and X. cortezi early- and late-pregnancy ovary tissue,
## pseudoaligned against Xbirchmanni pseudotranscriptome generated with X. maculatus genome
## RNAseq collected in 2023 from nonpregnant and pregnant mothers --> see metadata in ovary_embryo_rnaseq_samples_FebAug23_combined.txt
## Decided to code "some still yolking - 0/2" stages as stage 0
##
## add annotations to outfile with:
## Scripts/refs/add-gtf-annots-col.py Scripts/refs/Xiphophorus_maculatus_LG.Xipmac4.4.2.81_unique-tx_w-mito.gtf Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/ovary-xmacID-combined2023_dge_lfc-shr_all.csv
##
## note: dropping two samples (1M2O and 1B2O) because they are malinche stage 15 and birchmanni stage 5, with
##       outlier expression patterns (resembling late stage?), but want to be safe
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

## read in kallisto abundance files
dir <- "Data/embryo_ovary_dge_combined_2023/ovary_embryo_kallisto_output_FebAug23_combined/kallisto_posttrim_xbir-allhaps-xmac-txtome"
outdir <- "Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/"
of_header   <- "ovary-xmacID-combined2023_dge"
samples <- read.table("Data/embryo_ovary_dge_combined_2023/ovary_embryo_kallisto_output_FebAug23_combined/ovary_embryo_rnaseq_samples_FebAug23_combined.txt", header = TRUE)

## subset samples
tissue <- "ovary"
samples <- samples[samples$tissue == tissue,]

# drop samples 1M2O and 1B2O (outliers)
samples <- samples[samples$sample_ID != "1M2O" & samples$sample_ID != "1B2O",]
samples

# grab file paths to those samples
files <- file.path(dir, paste(samples$file_basename,"kallisto",sep="_"), "abundance.h5")
names(files) <- paste0(samples$sample_ID)
files

## load transcript annotations
#txdb <- makeTxDbFromGFF(file="Scripts/refs/Xiphophorus_maculatus_LG.Xipmac4.4.2.81_unique-tx_w-mito.gtf", format="gtf")
txdb <- makeTxDbFromGFF(file="Data/refs/Xiphophorus_maculatus_LG.Xipmac4.4.2.81_unique-tx_w-mito.gtf", format="gtf")

# create tx2gene table, match transcript name to geneid
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# get gene-level count information from kallisto counts
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)
head(txi$counts)


#################### DGE analysis ####################

## instantiate the DESeqDataSet with the design formula: ~ stage_group + origin + collection
dds     <- DESeqDataSetFromTximport(txi,
                                    colData = samples,
                                    design = ~ stage_group + origin + collection)

## determine the best fit for data
## by plotting or quantitatively with median absolute residual (smaller = better)
# parametric: 0.0437
dds <- DESeq(dds, fitType="parametric")
plotDispEsts(dds, main="Dispersion plot with parametric fit")
residual <- mcols(dds)$dispGeneEst - mcols(dds)$dispFit
absres<-abs(residual)
summary(absres)
# versus local: 0.03661
dds <- DESeq(dds, fitType="local")
plotDispEsts(dds, main="Dispersion plot with local fit")
residual <- mcols(dds)$dispGeneEst - mcols(dds)$dispFit
absres<-abs(residual)
summary(absres)

# chose local fit, run DESeq
dds <- DESeq(dds, fitType="local")
resultsNames(dds)

## save dds object
saveRDS(dds, paste0(outdir,of_header,"_dds.rds"))

## save vst object for heatmap
vst <- vst(dds,blind=FALSE)
saveRDS(vst, paste0(outdir,of_header,"_vst.rds"))

## load dds
dds <- readRDS("Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/ovary-xmacID-combined2023_dge_dds.rds")

## compare expression between species - early stage
res.malvbirch_early <- lfcShrink(dds, coef="stage_group_mal_early_vs_bir_early", type="ashr")
res.corvbirch_early <- lfcShrink(dds, coef="stage_group_cor_early_vs_bir_early", type="ashr")
res.malvcor_early <- lfcShrink(dds, contrast=list("stage_group_mal_early_vs_bir_early","stage_group_cor_early_vs_bir_early"), type="ashr")

## compare expression between species - late stage
res.malvbirch_late <- lfcShrink(dds, contrast=list("stage_group_mal_late_vs_bir_early","stage_group_bir_late_vs_bir_early"), type="ashr")
res.corvbirch_late <- lfcShrink(dds, contrast=list("stage_group_cor_late_vs_bir_early","stage_group_bir_late_vs_bir_early"), type="ashr")
res.malvcor_late <- lfcShrink(dds, contrast=list("stage_group_mal_late_vs_bir_early","stage_group_cor_late_vs_bir_early"), type="ashr")

## compare expression between stage categories - within species
res.latevearly_mal <- lfcShrink(dds, contrast=list("stage_group_mal_late_vs_bir_early","stage_group_mal_early_vs_bir_early"), type="ashr")
res.latevearly_cor <- lfcShrink(dds, contrast=list("stage_group_cor_late_vs_bir_early","stage_group_cor_early_vs_bir_early"), type="ashr")
res.latevearly_birch <- lfcShrink(dds, coef="stage_group_bir_late_vs_bir_early", type="ashr")


#################### START PRODUCE OUT TABLES ####################

## output results
# collect all results DESeq2 object variables starting with "res."
res_list <- lapply(ls(pattern = "^res\\."), get)
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

  # if on the first object, initialize resdata.shr
  if (i == 1) {
    resdata.shr <- res_slim
  }
  # otherwise merge new results values with master table of LFC and padj values
  else {
    resdata.shr <- merge(resdata.shr,as.data.frame(res_slim), by="Gene",sort=FALSE)
  }
}

# merge and output results values with the normalized counts data for all samples
resdata.shr <- merge(resdata.shr, as.data.frame(counts(dds, normalized=TRUE)), by.x="Gene", by.y="row.names",sort=FALSE)
write.csv(resdata.shr, file=paste(outdir,of_header,"_lfc-shr_all.csv",sep=""))

dim(res_out) # 19176


### Summary of expression patterns
ovary_dge <- read.csv("Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/ovary-xmacID-combined2023_dge_lfc-shr_all.csv_with-annots.csv",header=T)
ovary_dge_expressed_all <- subset(ovary_dge, !is.na(padj_res.malvbirch_early) & !is.na(padj_res.malvbirch_late) & !is.na(padj_res.malvcor_late) & !is.na(padj_res.malvcor_late) & !is.na(padj_res.corvbirch_late) & !is.na(padj_res.corvbirch_late)) # 1416 not expressed in all species, 17760 expressed
ovary_dge_expressed_one <- subset(ovary_dge, !is.na(padj_res.malvbirch_early) | !is.na(padj_res.malvbirch_late) | !is.na(padj_res.malvcor_late) | !is.na(padj_res.malvcor_late) | !is.na(padj_res.corvbirch_late) | !is.na(padj_res.corvbirch_late)) # 687 not expressed in at least one species, 18489 expressed
ovary_dge_expressed_mal <- subset(ovary_dge, !is.na(padj_res.latevearly_mal)) # 1416 not expressed in mal, 17760 expressed
ovary_dge_expressed_bir <- subset(ovary_dge, !is.na(padj_res.latevearly_birch)) # 2504 not expressed bir, 16672 expressed
ovary_dge_expressed_cor <- subset(ovary_dge, !is.na(padj_res.latevearly_cor)) # 1779 not expressed in cor, 17397 expressed

ovary_padj0.05_malvbir_early <- subset(ovary_dge, padj_res.malvbirch_early < 0.05) # 2202
ovary_padj0.05_malvbir_late <- subset(ovary_dge, padj_res.malvbirch_late < 0.05) # 4934
dim(intersect(ovary_padj0.05_malvbir_early,ovary_padj0.05_malvbir_late)) # 1672 shared
(2202-1672) + (4934-1672) # 3792 change expression from early to late pregnancy
ovary_padj0.05_malvcor_early <- subset(ovary_dge, padj_res.malvcor_early < 0.05) # 7957
ovary_padj0.05_malvcor_late <- subset(ovary_dge, padj_res.malvcor_late < 0.05) # 5639
dim(intersect(ovary_padj0.05_malvcor_early,ovary_padj0.05_malvcor_late)) # 4066 shared
(7957-4066) + (5639-4066) # 5464 change expression from early to late pregnancy
ovary_padj0.05_corvbir_early <- subset(ovary_dge, padj_res.corvbirch_early < 0.05) # 5393
ovary_padj0.05_corvbir_late <- subset(ovary_dge, padj_res.corvbirch_late < 0.05) # 5425
dim(intersect(ovary_padj0.05_corvbir_early,ovary_padj0.05_corvbir_late)) # 3178 shared
(5393-3178) + (5425-3178) # 4462 change expression from early to late pregnancy

n_expressed <- 17760

## get the genes that are sig different within species between late and stage ovaries
ovary_padj0.05_latevearly_mal <- subset(ovary_dge, padj_res.latevearly_mal < 0.05) # 310 /17760
ovary_padj0.05_latevearly_bir <- subset(ovary_dge, padj_res.latevearly_birch < 0.05) # 63 /16672
ovary_padj0.05_latevearly_cor <- subset(ovary_dge, padj_res.latevearly_cor < 0.05) # 130 /17397
