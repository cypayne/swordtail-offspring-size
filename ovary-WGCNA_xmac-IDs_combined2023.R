## Project: Offspring size evolution
## Payne et al
##
## WGCNA analysis of ovary gene expression data between Xmal, Xbir, Xcor
## Hierarchical clustering of genes into modules based on
## expression profile across samples, and correlation of
## module expression trend with traits of interest (i.e.
## genotype and mother origin treatment)
##
##     Need:
##           kallisto output (abundance.h5 file per sample)
##           GFF/GTF with annotated transcripts
##           sample file
##
## cyp XII-2023

## Useful tutorials:
## https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/

## R version 4.1.2 (2021-11-01)
sessionInfo()

## load libraries
library(tximportData)
library(tximport)
library(GenomicFeatures)
library(readr)
library(DESeq2)
library(rhdf5)
library(WGCNA)

## set out directory
outdir <- "Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/"
out_header <- "ovary-xmac_combined-FebAug23_"

######  Process input data - DGE analysis with DESeq2 ######

# specify tissue
tissue <- "ovary"

# specify kallisto directory
dir <- "Data/embryo_ovary_dge_combined_2023/ovary_embryo_kallisto_output_FebAug23_combined/kallisto_posttrim_xbir-allhaps-xmac-txtome"

## read in sample files
samples <- read.table("Data/embryo_ovary_dge_combined_2023/ovary_embryo_kallisto_output_FebAug23_combined/ovary_embryo_rnaseq_samples_FebAug23_combined.txt", header = TRUE)

## subset data
# subset tissue
samples <- samples[samples$tissue == tissue, ]
# drop outlier samples 1M2O and 1B2O
samples <- samples[samples$sample_ID != "1M2O" & samples$sample_ID != "1B2O",]
samples


###### Load VST object from DESeq2 ######

vst <- readRDS("Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/ovary-xmacID-combined2023_dge_vst.rds")


###### Running WGCNA ######

## set parameter
options(stringsAsFactors = FALSE)

## load vst counts
vsd <- assay(vst)

## transpose, put into dataframe
datExpr0 <- as.data.frame(t(vsd))


### filter genes and samples

# check for missing entries, weights below a threshold,
# and zero-variance genes, returns list of genes/samples
# that pass these filters
gsg <- goodSamplesGenes(datExpr0, minFraction = 1/2, verbose = 3) # 261 dropped (<14/28 samples or 0 variance) --> 18816
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");

# plot sample tree to look for outliers
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# remove sample outliers by choosing height for branch cut

# plot a line to show the cutoff
# this step is procedural, no samples were dropped
abline(h = 15, col = "red")
# determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==0) # change the cluster value, in this case only 1 cluster=0
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# didn't remove anything


## make sample dendogram with Trait heat map

# one-hot encode your traits:
library("caret")
traits <-samples[, c("sample_ID","stage_group","origin","collection")]
rownames(traits)<-samples$sample
traits <- traits[, -1]
traits

dmy <- dummyVars(" ~ .", data = traits)
one_hotTraits <- data.frame(predict(dmy, newdata = traits))
allTraits<-one_hotTraits
allTraits
datTraits <- allTraits

# re-cluster samples
sampleTree2 <- hclust(dist(datExpr), method = "average")
# convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(datTraits[, c("stage_groupbir_early","stage_groupbir_late","stage_groupcor_early","stage_groupcor_late","stage_groupmal_early","stage_groupmal_late","originlab","collection2023Aug")], signed = FALSE)
# plot the sample dendrogram and the colors underneath.
pdf(file=paste0(outdir,out_header,'sample-tree.pdf'),height=7,width=9.5)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# save objects for subsequent easy loading:
save(datExpr, datTraits, file = paste0(outdir,out_header,"dataInput.RData"))

# to reload these objects, uncomment and run:
#lnames <- load(file = paste0(outdir,out_header,"dataInput.RData"))


## choose soft-thresholding power
## test a set of soft-thresholding powers
## choose lowest power for which scale-free topoplogy index (SFT.R.sq) reaches 0.9 or inflection point of plot
## https://support.bioconductor.org/p/87024/
powers <- c(c(1:10), seq(from = 12, to=30, by=2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# plot the results:
pdf(file=paste0(outdir,out_header,'pickSoftThreshold.pdf'),height=7,width=9.5)
par(mfrow = c(1,2));
cex1 = 0.9;

# plot 1: Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# plot 2: Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# opted for WGCNA guide recommendation -
# if working with <20 samples, recommended soft-threshold for unsigned network: 9
# since working with 28 samples, opted for soft-threshold of 8 as recommended for 20-30 samples


## cluster genes into modules

## we'll use a single block (i.e. all 19k genes in one go, use maxBlockSize=20000)
# maxBlockSize = 20000, blockSizePenaltyPower = Inf : if maxBlockSize > total tx
# maxBlockSize = 10000, blockSizePenaltyPower = 5 : otherwise
net <- blockwiseModules(datExpr, power = 8,
                       TOMType = "unsigned", minModuleSize = 20,
                       maxBlockSize = 20000, blockSizePenaltyPower = Inf,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = FALSE, pamRespectsDendro = FALSE,
                       loadTOMs = FALSE,
                       saveTOMFileBase = "ovary-xmac_combined-FebAug23_TOM",
                       verbose = 3)

# convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# plot the dendrogram and the module colors underneath
pdf(file=paste0(outdir,out_header,'WGCNA_blockwiseModule-dendrogram_onehot.pdf'),height=8.5,width=12)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# save modules as an Rdata object
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = paste0(outdir,out_header,"networkConstruction-auto.RData"))


## correlate modules with traits of interest

# load the expression and trait data saved in the first part
lnames = load(file = paste0(outdir,out_header,"dataInput.RData"))
lnames
# load network data saved in the second part
lnames = load(file = paste0(outdir,out_header,"networkConstruction-auto.RData"))
lnames

# get module eigenvalues and calculate correlation with each trait
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# calculate correlations
moduleTraitCor = cor(MEs, cbind(datTraits$stage_groupbir_early,datTraits$stage_groupbir_late,datTraits$stage_groupcor_early,datTraits$stage_groupcor_late,datTraits$stage_groupmal_early,datTraits$stage_groupmal_late,datTraits$originlab,datTraits$collection2023Aug), use = "p")

# calculate p-values
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


## plot correlations
textMatrix =  paste(signif(moduleTraitCor, 2), sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf(file=paste0(outdir,out_header,"WGCNA_module-trait-heatmap.pdf"),height=10,width=8)
par(mar = c(6, 8.5, 3, 3))
# Plot heatmap of ME-trait correlations
hmcolors <- colorRampPalette(c("deepskyblue3", "white", "gold"))(n = 100)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c("Xbir_early","Xbir_late","Xcor_early","Xcor_late","Xmal_early","Xmal_late","origin","collection"),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = hmcolors,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


## subset MEs with significant trait correlations

# subset species-stage related modules:
subset(moduleTraitPvalue,moduleTraitPvalue[,1]<0.05) # xbir early
subset(moduleTraitPvalue,moduleTraitPvalue[,2]<0.05) # xbir late
subset(moduleTraitPvalue,moduleTraitPvalue[,3]<0.05) # xcor early
subset(moduleTraitPvalue,moduleTraitPvalue[,4]<0.05) # xcor late
subset(moduleTraitPvalue,moduleTraitPvalue[,5]<0.05) # xmal early
subset(moduleTraitPvalue,moduleTraitPvalue[,6]<0.05) # xmal late
# subset origin related modules:
subset(moduleTraitPvalue,moduleTraitPvalue[,7]<0.05) # origin (lab)
# subset collection related modules:
subset(moduleTraitPvalue,moduleTraitPvalue[,8]<0.05) # collection (Aug23)


## output significant ME pvalues by trait
colnames(moduleTraitPvalue) <- c("Xbir_early","Xbir_late","Xcor_early","Xcor_late","Xmal_early","Xmal_late","origin","collection")
write.csv(moduleTraitPvalue,paste0(outdir,out_header,"MEtraitpvals.csv"))
write.csv(subset(moduleTraitPvalue,moduleTraitPvalue[,1]<0.05),paste0(outdir,out_header,"WGCNA_MEtraitpvals_sig-xbir-early.csv"),quote=F)
write.csv(subset(moduleTraitPvalue,moduleTraitPvalue[,2]<0.05),paste0(outdir,out_header,"WGCNA_MEtraitpvals_sig-xbir-late.csv"),quote=F)
write.csv(subset(moduleTraitPvalue,moduleTraitPvalue[,3]<0.05),paste0(outdir,out_header,"WGCNA_MEtraitpvals_sig-xcor-early.csv"),quote=F)
write.csv(subset(moduleTraitPvalue,moduleTraitPvalue[,4]<0.05),paste0(outdir,out_header,"WGCNA_MEtraitpvals_sig-xcor-late.csv"),quote=F)
write.csv(subset(moduleTraitPvalue,moduleTraitPvalue[,5]<0.05),paste0(outdir,out_header,"WGCNA_MEtraitpvals_sig-xmal-early.csv"),quote=F)
write.csv(subset(moduleTraitPvalue,moduleTraitPvalue[,6]<0.05),paste0(outdir,out_header,"WGCNA_MEtraitpvals_sig-xmal-late.csv"),quote=F)
write.csv(subset(moduleTraitPvalue,moduleTraitPvalue[,7]<0.05),paste0(outdir,out_header,"WGCNA_MEtraitpvals_sig-origin-lab.csv"),quote=F)
write.csv(subset(moduleTraitPvalue,moduleTraitPvalue[,8]<0.05),paste0(outdir,out_header,"WGCNA_MEtraitpvals_sig-collection-Aug23.csv"),quote=F)


## plot heatmap for modules with at least one significant trait relationship
sigmoduleTraitCor <- subset(moduleTraitCor,moduleTraitPvalue[,1] <0.05 | moduleTraitPvalue[,2] <0.05 | moduleTraitPvalue[,3] <0.05 | moduleTraitPvalue[,4] <0.05 | moduleTraitPvalue[,5] <0.05 | moduleTraitPvalue[,6] <0.05 | moduleTraitPvalue[,7] <0.05 | moduleTraitPvalue[,8] <0.05)
sigmoduleTraitPvalue<- subset(moduleTraitPvalue,moduleTraitPvalue[,1] <0.05 | moduleTraitPvalue[,2] <0.05 | moduleTraitPvalue[,3] <0.05 | moduleTraitPvalue[,4] <0.05 | moduleTraitPvalue[,5] <0.05 | moduleTraitPvalue[,6] <0.05 | moduleTraitPvalue[,7] <0.05 | moduleTraitPvalue[,8] <0.05)
sigMEs <-subset(names(MEs),moduleTraitPvalue[,1] <0.05 | moduleTraitPvalue[,2] <0.05 | moduleTraitPvalue[,3] <0.05 | moduleTraitPvalue[,4] <0.05 | moduleTraitPvalue[,5] <0.05 | moduleTraitPvalue[,6] <0.05 | moduleTraitPvalue[,7] <0.05 | moduleTraitPvalue[,8] <0.05)
textMatrix =  paste(signif(sigmoduleTraitCor, 2), sep = "");
dim(textMatrix) = dim(sigmoduleTraitCor)
hmcolors <- colorRampPalette(c("deepskyblue3", "white", "gold"))(n = 100)

# plot heatmap
pdf(file=paste0(outdir,out_header,'WGCNA_module-trait-heatmap_sigMEs.pdf'),height=7.5,width=5.5)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = sigmoduleTraitCor,
               xLabels = c("Xbir_early","Xbir_late","Xcor_early","Xcor_late","Xmal_early","Xmal_late","origin","collection"),
               yLabels = sigMEs,
               ySymbols = sigMEs,
               colorLabels = FALSE,
               #               colors = blueWhiteRed(50),
               colors = hmcolors,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


## plot heatmap for species-stage-associated modules
sigmoduleTraitCor <- subset(moduleTraitCor,moduleTraitPvalue[,7] > 0.05 & moduleTraitPvalue[,8] > 0.05)
sigmoduleTraitPvalue<- subset(moduleTraitPvalue,moduleTraitPvalue[,4] <0.05)
sigMEs <-subset(names(MEs),moduleTraitPvalue[,4] <0.05)
textMatrix =  paste(signif(sigmoduleTraitCor, 2), sep = "");
dim(textMatrix) = dim(sigmoduleTraitCor)
hmcolors <- colorRampPalette(c("deepskyblue3", "white", "gold"))(n = 100)

# plot heatmap
pdf(file='TT-brain-WGCNA_module-trait-heatmap_sigTempMEs.pdf',height=3,width=5)
par(mar = c(4, 9, 2, 1));
labeledHeatmap(Matrix = sigmoduleTraitCor,
               xLabels = c("Xbir","Xcor","Xmal","early_stage","late_stage","lab_origin","wild_origin"),
               yLabels = sigMEs,
               ySymbols = sigMEs,
               colorLabels = FALSE,
               #               colors = blueWhiteRed(50),
               colors = hmcolors,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


## output module membership for all modules with significant trait relationship

# make outfiles for all significant modules, in bulk
for( MEname in sigMEs ) {
  MEname <- substring(MEname, 3)
  print(MEname)
  ME<-names(datExpr)[moduleColors==MEname]
  write.csv(ME,paste0(outdir,out_header,MEname,"_genes.csv"))
}

# save sigMEs object for future use
save(sigMEs, file = paste0(outdir,out_header,"sigMEs.RData"))
