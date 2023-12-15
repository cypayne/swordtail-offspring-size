##### FIGURES
##### Recent evolution of large offspring size and post-fertilization nutrient provisioning in swordtails
##### Payne et al
##### 12/2023

## load required libraries
library(ggplot2)
library(visreg)
library(dplyr)
library(emmeans)

## set global color variables
malcol <- "#00BFC4"
bircol <- "#F8766D"
corcol <- "#FFDC74"
birxmal_f1col <- "#c03ae2"
malxbir_f1col <- "#8f31de"
malxbir_f2col <- "#BD2FFF"
hybrid_col <- "#E162FF"


## Figure 1A: 5 species Xiphophorus fry standard length
fry_size_data <- read.table("Data/newborn_fry_size_data.csv",header=T,sep=',')
xipho_sp_size <- subset(fry_size_data, species_site %in% c("XbirCOAC","XmalCHIC","XcorPTHC","XpygPTHC","XvarCOAC"))
# average fry size within a brood
xipho_sp_size_summed <- as.data.frame(
  xipho_sp_size %>%
    group_by(brood_ID) %>%
    reframe(sl_mm_mean = mean(sl_mm,na.rm=T),species=first(species)))
# plot raw data
species_order <- c("Xmal","Xbir","Xcor","Xpyg","Xvar")
color_list <- c(malcol,bircol,corcol,"#FBAC87","#B3E561")
text_col <- "black"
xipho_sp_fry_size_plot <-
  ggplot(xipho_sp_size_summed, aes(x=factor(species,levels=species_order), y=sl_mm_mean, fill=factor(species,levels=species_order), color=factor(species,levels=species_order))) +
  scale_color_manual(values=color_list) +
  scale_fill_manual(values=color_list) +
  geom_jitter(alpha=0.6,width=0.1, cex=0.7) +
  stat_summary(fun.data = mean_se, width=0.4, size=0.8, geom = "crossbar", fill=NA) +
  geom_violin(color=NA, alpha=0.2, trim=FALSE) +
  theme_minimal() +
  theme(
    legend.position="none",
    legend.title=element_blank(),
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.1),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  ylim(6,15) +
  xlab("Xiphophorus species") +
  ylab("standard length (mm)")
xipho_sp_fry_size_plot
ggsave(xipho_sp_fry_size_plot,filename='Figures/Fig1A_Xipho-5species-fry-std-length.pdf',height=8.5,width=8.5,bg = "transparent")


## Figure 1B: Hybrid fry standard length
fry_size_data <- read.table("Data/newborn_fry_size_data.csv",header=T,sep=',')
mal_bir_hyb_size <- subset(fry_size_data, species_site %in% c("XbirCOAC","XmalCHIC","Xmal_xbir_F1","Xmal_xbir_F2","Xmal_xbirCALL"))
# average fry size within a brood
mal_bir_hyb_size_summed <- as.data.frame(
  mal_bir_hyb_size %>%
    group_by(brood_ID) %>%
    reframe(sl_mm_mean = mean(sl_mm,na.rm=T),species_site=first(species_site)))
# plot raw data
group_order <- c('XmalCHIC','XbirCOAC','Xmal_xbir_F1','Xmal_xbir_F2','Xmal_xbirCALL')
color_list <- c(malcol, bircol, malxbir_f1col, malxbir_f2col, hybrid_col)
text_col <- "black"
sl <- ggplot(mal_bir_hyb_size_summed, aes(x=factor(species_site,levels=group_order), y=sl_mm_mean, fill=factor(species_site,levels=group_order), color=factor(species_site,levels=group_order))) +
  scale_color_manual(values=color_list) +
  scale_fill_manual(values=color_list) +
  geom_jitter(alpha=0.6,width=0.1, cex=0.7) +
  stat_summary(fun.data = mean_se, width=0.4, size=0.8, geom = "crossbar", fill=NA) +
  geom_violin(color=NA, alpha=0.2, trim=FALSE) +
  theme_minimal() +
  theme(
    legend.position="none",
    legend.title=element_blank(),
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.1),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  ylim(6,15) +
  scale_x_discrete(labels = c('Xmal','Xbir',expression('F'[1]),expression('F'[2]),'CALL')) +
  xlab("Xiphophorus species") +
  ylab("standard length (mm)")
sl
ggsave(sl,filename='Figures/Fig1B_Xipho-Xmal-Xbirch-f1-f2-CALL-fry_std-length.pdf',bg = "transparent",width=8.5,height=8.5)


## Figure 2A: Xmalinche (CHIC) and Xbirchmanni (COAC) embryo size over development partial residual plots
combined_embryo_data <- read.table("Data/all-pops_combined_embryo_datasets.csv",header=T,sep=',')
combined_subset <- subset(combined_embryo_data,(popcoll=="CHIC_V-2022" | popcoll=="CHIC_VIII_20" | popcoll=="COAC_VIII_20" | popcoll=="COAC_II-2023" | popcoll=="COAC_IX-2022" | popcoll=="PTHC_II-2023") & stage!="0" & stage!="2" & stage!="5" & stage!="?" & stage!="-")
# average embryo size for embryos of the same stage within a brood
combined_subset_summed <- as.data.frame(
  combined_subset %>%
    group_by(brood_ID,species,population,collection,stage) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T),brood_size=mean(brood_size,na.rm=T),mother_std_length=first(mother_std_length)))
# choose model
combined_full_model <- lm(embryo_dry_weight_g_mean~species+stage+collection+brood_size+mother_std_length,data=combined_subset_summed)
step(combined_full_model) # embryo_dry_weight_g ~ species + stage + collection + mother_std_length
combined_fit <- lm(embryo_dry_weight_g_mean ~ stage + species + collection + mother_std_length,data=combined_subset_summed)
# calculate partial residuals
combined_stage_partial_residuals <- visreg(combined_fit, "stage",by="species",plot=F)$res
# reorder levels
combined_stage_partial_residuals$species <- factor(combined_stage_partial_residuals$species, levels = c("Xmalinche","Xbirchmanni","Xcortezi"))
# plot developmental profile
color_list <- c(malcol,bircol,corcol)
text_col = "black"
combined_stage_partial_residuals <- subset(combined_stage_partial_residuals,species!="Xcortezi")
embryo_dev_profiles <- ggplot(combined_stage_partial_residuals, aes(x=stage, y=visregRes, fill=species, color=species)) +
  geom_jitter(aes(col=species),alpha=0.6,width=0.1, cex=0.7) +
  scale_color_manual(values=color_list,breaks=c('Xmalinche', 'Xbirchmanni', 'Xcortezi'),labels=c('Xmal', 'Xbir', 'Xcor')) +
  stat_summary(fun.data = mean_se, width=0.7, size=0.8, geom = "crossbar", fill=NA) +
  theme_minimal() +
  theme(
    legend.position=c(0.97,1.03),
    legend.title=element_blank(),
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.1),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  ylim(0.0018,0.0051) +
  xlab("stage") +
  ylab("partial residuals of embryo dry weight (g)")
embryo_dev_profiles
ggsave(embryo_dev_profiles,filename='Figures/Fig2A_partial-residuals_dev-profile_all-stages_Xmal-Xbir_black-text_mean-se.pdf',height=8.5,width=13,bg = "transparent")

## Supp Fig S5: Developmental profiles with raw data
combined_subset_summed_xmal_xbir <- subset(combined_subset_summed,species!="Xcortezi")
raw_embryo_dev_profiles <- ggplot(combined_subset_summed_xmal_xbir, aes(x=stage, y=embryo_dry_weight_g_mean, fill=species, color=species)) +
  geom_jitter(aes(col=species),alpha=0.6,width=0.1, cex=0.7) +
  scale_color_manual(values=color_list,breaks=c('Xmalinche', 'Xbirchmanni', 'Xcortezi'),labels=c('Xmal', 'Xbir', 'Xcor')) +
  stat_summary(fun.data = mean_se, width=0.7, size=0.8, geom = "crossbar", fill=NA) +
  theme_minimal() +
  theme(
    legend.position=c(0.97,1.03),
    legend.title=element_blank(),
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.1),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  ylim(0.0018,0.0051) +
  xlab("stage") +
  ylab("embryo dry weight (g)")
raw_embryo_dev_profiles
ggsave(raw_embryo_dev_profiles,filename='Figures/FigSX_raw_dev-profile_all-stages_Xmal-Xbir_black-text_mean-se.pdf',height=8.5,width=10.5,bg = "transparent")


## Figure 2B: Comparison of stage 0 dry weight from CHIC, COAC, and PTHC
stage0_data <- read.table("Data/CHIC_COAC_PTHC_fully-yolked-stage0_dry-weights_mother-length.csv",header=T,sep=',')
stage0_data$population[stage0_data$population=="CHIC-DOWN"] <- "CHIC"
stage0_data$population[stage0_data$population=="CHIC-UP"] <- "CHIC"
# average egg size within a brood
stage0_data_summed <- as.data.frame(
  stage0_data %>%
    group_by(brood_ID,species,population,collection) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T),brood_size=mean(brood_size,na.rm=T),mother_std_length=first(mother_std_length)))
# choose model
stage0_fit <- lm(embryo_dry_weight_g_mean~species+collection+population:collection+brood_size+mother_std_length,stage0_data_summed)
step(stage0_fit) # brood_size + mother_std_length
# Statistical comparison of means
stage0_fit_lm <- lm(embryo_dry_weight_g_mean~species+brood_size+mother_std_length,stage0_data_summed)
# compare estimated marginal means (EMMs), adjusted pairwise comparisons with Tukey post-hoc test
emmeans(stage0_fit_lm,"species")
# contrast                 estimate       SE df t.ratio p.value
# Xbirchmanni - Xcortezi   7.51e-04 0.000640 48   1.173  0.4746
# Xbirchmanni - Xmalinche -7.02e-05 0.000421 48  -0.167  0.9848
# Xcortezi - Xmalinche    -8.22e-04 0.000924 48  -0.889  0.6495

# calculate partial residuals accounting for covariates
stage0_partial_residuals <- visreg(stage0_fit_lm, "species",plot=F)$res

# plot without X. cortezi
stage0_size_violin <- ggplot(stage0_partial_residuals, aes(x=species, y=visregRes, col=species, fill=species)) + theme_bw() +
  scale_color_manual(values=color_list,breaks=c("Xmalinche","Xbirchmanni")) +
  scale_fill_manual(values=color_list,breaks=c("Xmalinche","Xbirchmanni")) +
  geom_jitter(alpha=0.6,width=0.1, cex=0.7) +
  stat_summary(fun.data = mean_se, width=0.6, size=0.8, geom = "crossbar", fill=NA) +
  geom_violin(color=NA, alpha=0.2, trim=FALSE) +
  theme_minimal() +
  theme(
    legend.position="none",
    legend.title=element_blank(),
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.2),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  ylim(0.0000,0.0059) +
  xlab("species") +
  ylab("partial residuals of yolked egg dry weight (g)") +
  scale_x_discrete(limits=c("Xmalinche","Xbirchmanni"),labels=c("Xmal","Xbir"))
stage0_size_violin
ggsave(stage0_size_violin,filename='Figures/Fig2B_COAC-CHIC-egg-size_violin.pdf',bg = "transparent",width=6,height=8.5)

## Supp Fig 17: Plot egg size with X. cortezi
text_col <- "black"
species_order <- c("Xmalinche","Xbirchmanni", "Xcortezi")
color_list <- c(malcol,bircol,corcol)
stage0_size_violin <- ggplot(stage0_partial_residuals, aes(x=species, y=visregRes, col=species, fill=species)) + theme_bw() +
  #geom_boxplot(aes(col=species),position="identity",alpha=0,lwd=1.5) +
  scale_color_manual(values=color_list,breaks=c("Xmalinche","Xbirchmanni","Xcortezi")) +
  scale_fill_manual(values=color_list,breaks=c("Xmalinche","Xbirchmanni","Xcortezi")) +
  geom_jitter(alpha=0.6,width=0.1, cex=0.7) +
  stat_summary(fun.data = mean_se, width=0.4, size=0.8, geom = "crossbar", fill=NA) +
  geom_violin(color=NA, alpha=0.2, trim=FALSE) +
  theme_minimal() +
  theme(
    legend.position="none",
    legend.title=element_blank(),
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.2),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  ylim(0.0001,0.0059) +
  xlab("species") +
  ylab("partial residuals of yolked egg dry weight (g)") +
  scale_x_discrete(limits=c("Xmalinche","Xbirchmanni","Xcortezi"),labels=c("Xmal","Xbir","Xcor"))
stage0_size_violin
ggsave(stage0_size_violin,filename='Figures/Fig2B_COAC-CHIC-PTHC-egg-size_violin.pdf',bg = "transparent",width=6,height=8.5)


## Figure 2C: Partial residual plot of embryo size for lab-reared Xmalinche, Xbirchmanni, and F1s of both cross directions from roof tanks
roof_chic_coac_f1 <- read.csv("Data/IV-2023_roof-tank_CHIC-COAC-F1_embryo-weights.csv",header=T,sep=",")
roof_chic_coac_f1$stage <- as.factor(roof_chic_coac_f1$stage)
# drop early stages (have <=2 samples)
roof_chic_coac_f1_subset <- subset(roof_chic_coac_f1,stage!="0" & stage!="5" & stage!="10" & stage!="15" )
# remove one outlier - likely erroneous data recording
roof_chic_coac_f1_subset <- subset(roof_chic_coac_f1_subset,embryo_dry_weight_g < 0.008)
# choose model
roof_chic_coac_f1_full_model <- lm(embryo_dry_weight_g~species+stage+brood_size+mother_std_length+brood_ID,data=roof_chic_coac_f1_subset)
step(roof_chic_coac_f1_full_model) # embryo_dry_weight_g_mean ~ stage + brood_ID
lab_cross_fit <- lm(embryo_dry_weight_g ~ species + stage + brood_ID, data=roof_chic_coac_f1_subset)
# Stats: compare means with ANOVA and Tukey
res <- summary(pairs(emmeans(lab_cross_fit, "species")))
res
res$p.value
# contrast                     estimate          SE df t.ratio p.value
# birxbir - birxmal_F1     0.000411 0.000671 141   0.613  0.9278
# birxbir - malxbir_F1    -0.001626 0.000216 141  -7.529  <.0001
# birxbir - malxmal       -0.002578 0.000349 141  -7.376  <.0001
# birxmal_F1 - malxbir_F1 -0.002037 0.000709 141  -2.872  0.0241
# birxmal_F1 - malxmal    -0.002989 0.000495 141  -6.038  <.0001
# malxbir_F1 - malxmal    -0.000952 0.000369 141  -2.579  0.0527

# save the stage partial residuals
lab_cross_partial_residuals <- visreg(lab_cross_fit, "species",plot=F)$res
#write.csv(lab_cross_partial_residuals,"Data/lab_cross_xmal-xbir-F1_species_partial_residuals.csv")

color_list <- c(bircol, birxmal_f1col, malxbir_f1col, malcol)
lab_cross_embryo_size <- ggplot(lab_cross_partial_residuals, aes(x=species, y=visregRes, fill=species, color=species)) +
  scale_color_manual(values=color_list) +
  scale_fill_manual(values=color_list) +
  geom_jitter(alpha=0.6,width=0.1, cex=0.7) +
  stat_summary(fun.data = mean_cl_normal, width=0.4, size=0.8, geom = "errorbar", fill=NA) +
  stat_summary(fun.y = mean, geom = "point", pch="-", size=10) +
  theme_minimal() +
  theme(
    legend.position="none",
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  scale_x_discrete(labels = c(expression('bir'%*%'bir'),expression('bir'%*%'mal F'[1]),expression('mal'%*%'bir F'[1]),expression('mal'%*%'mal'))) +
  xlab("lab cross") +
  ylab("partial residuals of embryo dry weight (g)")
lab_cross_embryo_size
ggsave(lab_cross_embryo_size,filename='Figures/Fig2C_lab-cross_xmal-xbir-f1_embryo-size_partial-residuals_black-text.pdf',height=8.5,width=8.5,bg = "transparent")


## Figure 2D: Histology


## Figure 3A: PCA plot of ovary RNAseq samples
library(PCAtools)
library(DESeq2)
samples <- read.table("Data/embryo_ovary_dge_combined_2023/ovary_embryo_kallisto_output_FebAug23_combined/ovary_embryo_rnaseq_samples_FebAug23_combined.txt", header = TRUE)
samples <- samples[samples$tissue == "ovary",] # subset tissue
samples <- samples[samples$sample_ID != "1M2O" & samples$sample_ID != "1B2O",] # remove outliers
vst <- readRDS("Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/ovary-xmacID-combined2023_dge_vst.rds")
x <- assay(vst)
rownames(samples) <- samples$sample
# create another variable
samples$stage_category_origin <- paste(samples$stage_category,samples$origin,sep="_")
p <- pca(x, metadata = samples, removeVar = 0.1)
# species and sample origin (lab v wild) cluster as expected
pdf(file='Figures/Fig3A_ovary_dge_xmacID_PCA.pdf',height=8.5,width=12,bg = "transparent")
biplot(p, x="PC1", y="PC2",
       colby = 'species', colkey=c('Xmalinche'=malcol,'Xbirchmanni'=bircol,'Xcortezi'=corcol),
       legendPosition = 'right', legendLabSize = 14, legendIconSize = 7,
       shape = 'stage_category_origin', shapekey = c('early_lab'=18, 'late_lab'=16, 'late_wild'=8), pointSize=4, labSize = 4,
       gridlines.major=FALSE, gridlines.minor=FALSE,
       drawConnectors = T,
       encircle=TRUE, encircleFill=TRUE,
       axisLabSize = 20,
       caption = '')
dev.off()


## Figure 3B: Prolactin receptor expression boxplots, as well as supplementary WGCNA cluster boxplots (Supp Fig S8 & S9)
library(plotrix)
# load ovary RNAseq data pseudoaligned against X. birchmanni pseudoreference
ovary_data <- read.csv("Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/ovary-xmacID-combined2023_dge_lfc-shr_all.csv_with-annots.csv",header=T)

# specify tissue
tissue <- "ovary"

# define error bar function
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

# define plotting function
dge_boxplot <- function(gene, tissue, ymin, ymax, legend_pos_x, legend_pos_y, save_plot) {
  # pull expression data
  if(tissue=="ovary") {
    n<-subset(ovary_data,Gene==gene)
    mal_early_lab <-c(n$X1M2O,n$X2M2O,n$X2M14O,n$X2M26O,n$X2M38O)
    mal_late_lab <- c(n$X1M6O,n$X1M4O,n$X1M3O)
    mal_late_wild <- c(n$X2M50OW,n$X2M58OW)
    bir_early_lab <-c(n$X1B2O,n$X2B30O,n$X2B42O)
    bir_late_lab <- c(n$X1B5O, n$X1B8O,n$X1B9O,n$X1B6O,n$X2B6O,n$X2B18O)
    bir_late_wild <- c(n$X2B54BOW,n$X2B70OW)
    cor_early_lab <-c(n$X1C2O,n$X1C4O,n$X2C10O,n$X2C22O,n$X2C34O,n$X2C46O)
    cor_late_lab <- c(n$X1C5O,n$X1C1O)
    ymin <- min(mal_early_lab,mal_late_lab,mal_late_wild,bir_early_lab,bir_late_lab,bir_late_wild,cor_early_lab,cor_late_lab)
    ymax <- max(mal_early_lab,mal_late_lab,mal_late_wild,bir_early_lab,bir_late_lab,bir_late_wild,cor_early_lab,cor_late_lab)
  }
  if(save_plot) {
    outfile = paste0(data2plot,"_ovary-xbir-gtf-exp-plot_",tissue,"_",n$Gene,".pdf")
    pdf(outfile,width=4,height=6)
  }
  ## visual settings
  malcol_early=rgb(0/255, 143/255, 147/255)
  malcol_late=rgb(0/255,191/255,196/255)
  bircol_early=rgb(186/255, 89/255, 82/255)
  bircol_late=rgb(248/255,118/255,109/255)
  corcol_early=rgb(230/255,172/255,0/255)
  corcol_late=rgb(255/255,220/255,116/255)
  malcol_wild=rgb(0/255,191/255,196/255)
  bircol_wild=rgb(248/255,118/255,109/255)

  BE<-c(bir_early_lab); ME<-c(mal_early_lab); CE<-c(cor_early_lab); BL<-c(bir_late_lab,bir_late_wild); ML<-c(mal_late_lab,mal_late_wild); CL<-c(cor_late_lab);

  plot(c(0.8,2.8,4.8),c(mean(ME),mean(BE),mean(CE)),col=c(malcol_early,bircol_early,corcol_early),ylim=c(ymin,ymax),ann=FALSE, xlab="",ylab="Normalized counts",xaxt="n",pch=18,cex=2,cex.axis=1.5,cex.label=1.5,cex.main=1.5, cex.sub=1.5, cex.names=1.5, xlim=c(0.4,6.6))
  error.bar(c(0.8,2.8,4.8),c(mean(ME),mean(BE),mean(CE)),c(2*std.error(ME),2*std.error(BE),2*std.error(CE)), col=c(malcol_early,bircol_early,corcol_early),lwd=2)

  points(c(1.8,3.8,5.8),c(mean(ML),mean(BL),mean(CL)),col=c(malcol_late,bircol_late,corcol_late),pch=20,cex=2)
  error.bar(c(1.8,3.8,5.8),c(mean(ML),mean(BL),mean(CL)),c(2*std.error(ML),2*std.error(BL),2*std.error(CL)),col=c(malcol_late,bircol_late,corcol_late),lwd=2)

  malcol_early=rgb(0/255,153/255,157/255,alpha=0.2)
  malcol_late=rgb(0/255,191/255,196/255,alpha=0.2)
  bircol_early=rgb(186/255, 89/255, 82/255,alpha=0.2)
  bircol_late=rgb(248/255,118/255,109/255,alpha=0.2)
  corcol_early=rgb(230/255,172/255,0/255,alpha=0.2)
  corcol_late=rgb(255/255,220/255,116/255,alpha=0.2)
  malcol_wild=rgb(0/255,191/255,196/255,alpha=0.2)
  bircol_wild=rgb(248/255,118/255,109/255,alpha=0.2)

  # add transparent points
  noise<-runif(length(mal_early_lab),0.2,0.3)
  points(rep(1,length(mal_early_lab))+noise,mal_early_lab,pch=18,cex=1.5,col=malcol_early)
  noise<-runif(length(mal_late_lab),0.2,0.3)
  points(rep(2,length(mal_late_lab))+noise,mal_late_lab,pch=20,cex=1.5,col=malcol_late)
  noise<-runif(length(mal_late_wild),0.2,0.3)
  points(rep(2,length(mal_late_wild))+noise,mal_late_wild,pch=8,cex=1.3,col=malcol_wild)

  noise<-runif(length(bir_early_lab),0.2,0.3)
  points(rep(3,length(bir_early_lab))+noise,bir_early_lab,pch=18,cex=1.5,col=bircol_early)
  noise<-runif(length(bir_late_lab),0.2,0.3)
  points(rep(4,length(bir_late_lab))+noise,bir_late_lab,pch=20,cex=1.5,col=bircol_late)
  noise<-runif(length(bir_late_wild),0.2,0.3)
  points(rep(4,length(bir_late_wild))+noise,bir_late_wild,pch=8,cex=1.3,col=bircol_wild)

  noise<-runif(length(cor_early_lab),0.2,0.3)
  points(rep(5,length(cor_early_lab))+noise,cor_early_lab,pch=18,cex=1.5,col=corcol_early)
  noise<-runif(length(cor_late_lab),0.2,0.3)
  points(rep(6,length(cor_late_lab))+noise,cor_late_lab,pch=20,cex=1.5,col=corcol_late)

  axis(1, at = c(1.4,3.4,5.4), labels = c("Xmal","Xbir","Xcor"), cex.axis = 1.5, tick=FALSE)
  title(ylab = "normalized counts", cex.lab = 1.5, line = 4)

  title(n$annot)

  if(save_plot) { dev.off() }
}

matrotrophy_genes <- c("ENSXMAG00000008852","ENSXMAG00000004299") # prolactin releaser (prlrh2a),  prlra
acylCoA_genes <- c("ENSXMAG00000017128","ENSXMAG00000010204","ENSXMAG00000015573","ENSXMAG00000008007","ENSXMAG00000010958","ENSXMAG00000014638","ENSXMAG00000000852","ENSXMAG00000013172") # acot15, aclya, elovl7a, mvk, dlat, dlst, acss2l, elovl1b
artery_dev_genes <- c("ENSXMAG00000001096","ENSXMAG00000005830","ENSXMAG00000016605","ENSXMAG00000012207") #cdkn1a, pth1ra, amotl2a, rbpja

# chnage which list of genes to plot with this variable
gene_list <- artery_dev_genes

# plot
dev.off()
par(mfrow=c(2,2), las=1,bg=NA,mar = c(5,6,4,2))
for (gene in gene_list) {

  ymin=0
  ymax=10000

  dge_boxplot(gene, "ovary", ymin, ymax, 0.1, 515, FALSE)
}
dev.copy(pdf,"Figures/FigX-artery-dev-boxplots.pdf",bg="transparent",width=8,height=8)
#dev.copy(pdf,"Figures/FigX-acyl-CoA-boxplots.pdf",bg="transparent",width=12,height=8)
#dev.copy(pdf,"Figures/Fig3B-prolactin-prlhr-prlr-boxplots.pdf",bg="transparent",width=8.5,height=5)
dev.off()


## Figure 3C: Heatmap of enriched biological pathways in WGCNA 'darkorange2' cluster
library("DESeq2")
library("genefilter")
library("pheatmap")
library("RColorBrewer")
library("reshape")
library("ggplot2")

# load sample metadata
samples <- read.table("Data/embryo_ovary_dge_combined_2023/ovary_embryo_kallisto_output_FebAug23_combined/ovary_embryo_rnaseq_samples_FebAug23_combined.txt", header = TRUE)
samples <- subset(samples, tissue == "ovary" & (samples$sample_ID != "1M2O" & samples$sample_ID != "1B2O"))
samples
# load DGE data
deseq_lfc <- read.csv("Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/ovary-xmacID-combined2023_dge_lfc-shr_all.csv_with-annots.csv",header=T)
# load variance stabilized DESeq2 normalized count data
vst_data <- readRDS("Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/ovary-xmacID-combined2023_dge_vst.rds")

# convert vst to dataframe
vst <- assay(vst_data)
vst <- as.data.frame(vst)
vst$Gene <- rownames(vst)

# define gene lists to map
acylCoA_genes <- c("ENSXMAG00000017128","ENSXMAG00000010204","ENSXMAG00000015573","ENSXMAG00000008007","ENSXMAG00000010958","ENSXMAG00000014638","ENSXMAG00000000852","ENSXMAG00000013172") # acot15, aclya, elovl7a, mvk, dlat, dlst, acss2l, elovl1b
artery_dev_genes <- c("ENSXMAG00000001096","ENSXMAG00000005830","ENSXMAG00000016605","ENSXMAG00000012207") #cdkn1a, pth1ra, amotl2a, rbpja

# XX change this variable depending on which list you want to map
gene_list <- acylCoA_genes

# subset genes
vst.sub <- vst[vst$Gene %in% gene_list,]
vst.sub <- subset(vst.sub, select = -Gene)

# calculate deviation of gene counts for each sample from mean
vst.sub <- vst.sub - rowMeans(vst.sub)
vst.sub$annot<-deseq_lfc$annot[match(rownames(vst.sub),deseq_lfc$Gene)]
vst.sub <- na.omit(vst.sub)
rownames(vst.sub)<-vst.sub$annot
vst.sub <- subset(vst.sub, select = -annot)

# make new dataframe
df <- data.frame(species_stage=vst_data$stage_group,stage_category=vst_data$stage_category,sample_ID=colnames(vst_data))

# sort metadata and expression data by sample species_stage category so that the
# species_stages cluster in the heatmap
df_sorted <- df[order(df$species_stage),]
rownames(df_sorted) <- df_sorted$sample_ID
df_sorted <- subset(df_sorted, select = -sample_ID)
vst.sub_sorted <- vst.sub[, row.names(df_sorted)]

# set colors
annotation_colors = list(species_stage=c(bir_early=bircol,bir_late=bircol,cor_early=corcol,cor_late=corcol,mal_early=malcol,mal_late=malcol),stage_category=c(early="grey",late="black"))

# set scale
range <- max(abs(vst.sub))

# plot
gene_list_heatmap <- pheatmap(vst.sub_sorted,
                              color = colorRampPalette(c("black","darkslategrey","darkcyan","white","goldenrod1","darkgoldenrod1","darkgoldenrod"))(100),
                              breaks = seq(-range, range, length.out = 100),
                              cluster_cols = F,
                              treeheight_row = 0,
                              annotation_col=df_sorted, annotation_colors = annotation_colors,
                              border_color=NA, show_rownames = TRUE, annotation_names_row = TRUE,
                              fontsize=14)
gene_list_heatmap
ggsave(gene_list_heatmap,filename='Figures/FigX_acyl-CoA_heatmap.pdf',height=6,width=9,bg = "transparent")
ggsave(gene_list_heatmap,filename='Figures/FigX_vasc-paths_heatmap.pdf',height=6,width=9,bg = "transparent")


## Figure 3D: Immunostains of prolactin in late-pregnancy ovaries


## Figure 4A: Fat content comparison between CHIC and COAC in February and September collections
library(ggpubr)
nonpreg_females <- read.csv("Data/II-2023_IX-2023_nonpreg-females_lipid-extraction.csv",header=T)
fat_content_data <- subset(nonpreg_females,FC_percent >= 0 & population!="PTHC" & population!="SELV")
fat_content_data$population <- as.factor(fat_content_data$population)

# stats - two-sided T-test
# subset data by species and season
fc_coac_feb <- subset(fat_content_data,population=="COAC" & collection=="II-2023")
fc_chic_feb <- subset(fat_content_data,population=="CHIC" & collection=="II-2023")
fc_coac_sep <- subset(fat_content_data,population=="COAC" & collection=="IX-2022")
fc_chic_sep <- subset(fat_content_data,population=="CHIC" & collection=="IX-2022")

# compare February 2023 Xmal versus Xbir
t.test(fc_chic_feb$FC_percent, fc_coac_feb$FC_percent, alternative = "two.sided", var.equal = FALSE) # t = 1.64, p-value 0.122
# compare September 2022 Xmal versus Xbir
t.test(fc_chic_sep$FC_percent, fc_coac_sep$FC_percent, alternative = "two.sided", var.equal = FALSE) # t = -12.45, p-value 1.297e-06

# plot
chiccol <- "#00BFC4"
coaccol <- "#F8766D"
color_list <- c(chiccol,coaccol)
group_order <- c("CHIC","COAC")
text_col="black"
feb_fc_sub <- subset(fat_content_data,collection=="II-2023")
sep_fc_sub <- subset(fat_content_data,collection=="IX-2022")
fat_content_plot_sep <-
  ggboxplot(sep_fc_sub, x="population", y="FC_percent",color="population", palette=color_list,fill=NA,lwd=1.5) +
  geom_jitter(aes(col=population),alpha=0.4,width=0.2) +
  theme_minimal() +
  theme(
    legend.position="none",
    legend.title=element_blank(),
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  scale_x_discrete(labels = c('Xmal','Xbir')) +
  xlab("September 2022") +
  ylab("female body fat (%)")

fat_content_plot_feb
fat_content_plot_sep
female_fc_plot <- ggarrange(fat_content_plot_feb,fat_content_plot_sep + rremove("y.title"))
ggsave(female_fc_plot,filename='Figures/Fig4C_population-fat-content.pdf',bg = "transparent",width=8.5,height=8.5)


## Figure 4B: Food deprivation fry experiment - dry mass (g)
fry_fat_content <- read.csv("Data/X-23_fry_starvation_fat_content.csv")
# remove rows with fry content < 0; remove trial run Xbir and Xmal broods (1m & 1b); remove NAs and negative FC percentages
fry_fat_content <- subset(fry_fat_content, FC_percent >= 0 & !is.na(FC_percent) & brood_no != "1m" & brood_no != "1b")
fry_fat_content$condition <- paste(fry_fat_content$species,fry_fat_content$treatment, sep="_")

# choose model
step(lm(dry_mass_1 ~ condition+species+treatment+brood_no,fry_fat_content)) # dry_mass_1 ~ condition + brood_no

# stats: Compare raw means, accounting for covariates
# ANOVA - posthoc Tukey
# no difference between Xmal v Xbir
# no difference between Xmal control v starved
# sig difference (decrease) between Xbir control v starved
# $`pairwise differences of condition`
#                           estimate     SE df t.ratio p.value
#Xbir_control - Xbir_starved  0.001043 0.000190 50   5.504  <.0001
#Xbir_control - Xmal_control -0.000256 0.000393 50  -0.651  0.9148
#Xbir_control - Xmal_starved  0.000552 0.000391 50   1.410  0.4992
#Xbir_starved - Xmal_control -0.001300 0.000390 50  -3.335  0.0085
#Xbir_starved - Xmal_starved -0.000492 0.000388 50  -1.269  0.5869
#Xmal_control - Xmal_starved  0.000808 0.000160 50   5.054  <.0001
dry_mass_fit <- lme(dry_mass_1 ~ condition, random=~1|brood_no,data=fry_fat_content)
emmeans(dry_mass_fit, list(pairwise ~ condition), adjust = "tukey")

# calculate partial residuals
dry_mass_fit <- lmer(standard_length_cm ~ condition+(1|brood_no),data=fry_fat_content)
dry_mass_partial_residuals <- visreg(dry_mass_fit, "condition",plot=F)$res

# plot partial residuals of dry mass
text_col <- "black"
color_list <- c(coaccol,coaccol,chiccol,chiccol)
condition_order <- c("Xmal_control","Xmal_starved","Xbir_control","Xbir_starved")
dm_fry_food_dep_plot <- ggplot(dry_mass_partial_residuals, aes(x=factor(condition,levels=condition_order), y=visregRes, fill=factor(condition,levels=condition_order))) +
  geom_boxplot(aes(col=condition),position="identity",alpha=0,lwd=1.5) +
  geom_jitter(aes(col=condition),alpha=0.4,width=0.2) +
  scale_color_manual(values=color_list) +
  theme_minimal() +
  theme(
    legend.position="none",
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  scale_x_discrete(labels = c("food", "no food", "food", "no food")) +
  xlab("fry treatment (3 days)") +
  ylab("partial residuals of dry mass (g)")
dm_fry_food_dep_plot
ggsave(dm_fry_food_dep_plot,filename='Figures/Fig4B_food-deprivation_dry-mass-boxplot_black-text.pdf',height=8.5,width=8.5,bg = "transparent")



##### SUPPLEMENTARY FIGURES

## Supp Fig S1: Image of measurements for standard length and head width


## Supp Fig S2A: 5 Xipho species newborn fry head width
fry_size_data <- read.table("Data/newborn_fry_size_data.csv",header=T,sep=',')
xipho_sp_size <- subset(fry_size_data, species_site %in% c("XbirCOAC","XmalCHIC","XcorPTHC","XpygPTHC","XvarCOAC"))
xipho_sp_hw_summed <- as.data.frame(
  xipho_sp_size %>%
    group_by(brood_ID) %>%
    reframe(hw_mm_mean = mean(hw_mm,na.rm=T),species=first(species)))
species_order <- c("Xmal","Xbir","Xcor","Xpyg","Xvar")
color_list <- c(malcol,bircol,corcol,"#FBAC87","#B3E561")
text_col <- "black"
xipho_sp_fry_size_plot_hw <-
  ggplot(xipho_sp_hw_summed, aes(x=factor(species,levels=species_order), y=hw_mm_mean, fill=factor(species,levels=species_order), color=factor(species,levels=species_order))) +
  scale_color_manual(values=color_list) +
  scale_fill_manual(values=color_list) +
  geom_jitter(alpha=0.6,width=0.1, cex=0.7) +
  stat_summary(fun.data = mean_se, width=0.4, size=0.8, geom = "crossbar", fill=NA) +
  geom_violin(color=NA, alpha=0.2, trim=FALSE) +
  theme_minimal() +
  theme(
    legend.position="none",
    legend.title=element_blank(),
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.1),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  ylim(1.3,3.3) +
  xlab("Xiphophorus species") +
  ylab("head width (mm)")
xipho_sp_fry_size_plot_hw
ggsave(xipho_sp_fry_size_plot_hw,filename='Figures/SuppFigS1_Xipho-5species-fry-head-width.pdf',height=8.5,width=8.5,bg = "transparent")

# get summary of mean and stdev for each species
xipho_sp_size_summed %>%
  group_by(species) %>%
  summarise_at(vars(hw_mm_mean),
               list(hw_mean = mean, hw_std = sd))

# pairwise comparisons using Wilcoxon rank sum exact test
pairwise.wilcox.test(xipho_sp_size_summed$hw_mm_mean, xipho_sp_size_summed$species,p.adjust.method = "BH")


## Supp. Fig. S2B: Newborn fry head width boxplot
fry_size_data <- read.table("Data/newborn_fry_size_data.csv",header=T,sep=',')
#mal_bir_hyb_size <- subset(fry_size_data, species_site %in% c("XbirCOAC","XmalCHIC","Xmal_xbir_F1","Xmal_xbir_F2","Xmal_xbirCALL","Xmal_xbirCAPS","Xmal_xbirPLNK"))
mal_bir_hyb_size <- subset(fry_size_data, species_site %in% c("XbirCOAC","XmalCHIC","Xmal_xbir_F1","Xmal_xbir_F2","Xmal_xbirCALL"))
mal_bir_hyb_hw_summed <- as.data.frame(
  mal_bir_hyb_size %>%
    group_by(brood_ID) %>%
    reframe(hw_mm_mean = mean(hw_mm,na.rm=T),species_site=first(species_site)))
group_order <- c('XmalCHIC','XbirCOAC','Xmal_xbir_F1','Xmal_xbir_F2','Xmal_xbirCALL')
color_list <- c(malcol, bircol, malxbir_f1col, malxbir_f2col, hybrid_col)
text_col <- "black"
hw <- ggplot(mal_bir_hyb_hw_summed, aes(x=factor(species_site,levels=group_order), y=hw_mm_mean, fill=factor(species_site,levels=group_order), color=factor(species_site,levels=group_order))) +
  scale_color_manual(values=color_list) +
  scale_fill_manual(values=color_list) +
  geom_jitter(alpha=0.6,width=0.1, cex=0.7) +
  stat_summary(fun.data = mean_se, width=0.4, size=0.8, geom = "crossbar", fill=NA) +
  geom_violin(color=NA, alpha=0.2, trim=FALSE) +
  theme_minimal() +
  theme(
    legend.position="none",
    legend.title=element_blank(),
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.1),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  ylim(1.3,3.3) +
  scale_x_discrete(labels = c('Xmal','Xbir',expression('F'[1]),expression('F'[2]),'CALL')) +
  xlab("Xiphophorus species") +
  ylab("head width (mm)")
hw
ggsave(hw,filename='Figures/SuppFig_Xipho-Xmal-Xbirch-hybrid-fry_head-width.pdf',bg = "transparent",width=8.5,height=8.5)

# Pairwise comparisons using Wilcoxon rank sum exact test
pairwise.wilcox.test(mal_bir_hyb_size_summed$hw_mm_mean, mal_bir_hyb_size_summed$species,p.adjust.method = "BH")


## Supp Fig S3: Brood size comparison between species


## Supp Fig S4: Adult lab-reared female size comparison between species


## Supp Fig S5: Embryo developmental profiles with raw data
## Code can be found under Figure 2A


## Supp Fig S6: Relationship of mother ancestry proportion and mitotype
# with embryo and ovary dry weights in hybrids from natural hybrid population Calnali Low
combined_embryo_data <- read.table("Data/all-pops_combined_embryo_datasets.csv",header=T,sep=',')
CALL_embryo_data <- subset(combined_embryo_data,population=="CALL")
CALL_mother_hi <- read.table("Data/CALL_mother_hybrid-index_mitotype.csv",header=T,sep=',')
CALL_embryo_data$mother_hi <- CALL_mother_hi$mom_index[match(CALL_embryo_data$brood_ID, CALL_mother_hi$brood_ID)]
CALL_embryo_data$mitotype <- CALL_mother_hi$mitotype[match(CALL_embryo_data$brood_ID, CALL_mother_hi$brood_ID)]
# remove early stages with few samples
CALL_embryo_data <- subset(CALL_embryo_data,stage!="-" & stage!="0" & stage!="10" & stage!="15")
# remove pure parent data
CALL_embryo_data_nomal <- subset(CALL_embryo_data,mother_hi<0.9 & mother_hi>0.1 )

## Supp Fig S6A - Relationship between embryo dry weight and mother hybrid index
# model fit
model_all<-lm(embryo_dry_weight_g~mother_hi*stage+brood_size,data=CALL_embryo_data_nomal)
selectedMod <- step(model_all) # embryo_dry_weight_g ~ mother_hi + stage + mother_hi:stage + brood_ID
summary(model_all)

# make new dataframe and remove NA
CALL_embryo_data_nomal_sub <- data.frame(brood_ID=CALL_embryo_data_nomal$brood_ID,stage=CALL_embryo_data_nomal$stage,brood_size=CALL_embryo_data_nomal$brood_size,mother_hi=CALL_embryo_data_nomal$mother_hi,mitotype=CALL_embryo_data_nomal$mitotype,embryo_dry_weight_g=CALL_embryo_data_nomal$embryo_dry_weight_g)
CALL_embryo_data_nomal_sub <- na.omit(CALL_embryo_data_nomal_sub)

# evaluate significance with likelihood ratio test (anova) on full and reduced model
CALL_hi_mlm_fit<-lme(embryo_dry_weight_g~mother_hi+stage+mother_hi:stage+brood_size,random=~1|brood_ID,data=CALL_embryo_data_nomal_sub)
CALL_hi_mlm_fit_reduced<-lme(embryo_dry_weight_g~stage,random=~1|brood_ID,data=CALL_embryo_data_nomal_sub)
anova(CALL_hi_mlm_fit,CALL_hi_mlm_fit_reduced)
# Model df       AIC       BIC   logLik   Test  L.Ratio p-value
# CALL_hi_mlm_fit             1 13 -10892.05 -10829.55 5459.026
# CALL_hi_mlm_fit_reduced     2  7 -10965.66 -10931.96 5489.829 1 vs 2 61.60707  <.0001

# use lmer to fit mixed linear model for visreg
CALL_hi_mlm_fit<-lmer(embryo_dry_weight_g~mother_hi+stage+mother_hi:stage+brood_size+(1|brood_ID),data=CALL_embryo_data_nomal_sub)

# plot mother hi and embryo dry weight
CALL_hi_mlm_vr <- visreg(CALL_hi_mlm_fit,"mother_hi")
plot(visregRes ~ mother_hi,CALL_hi_mlm_vr$res,xlab="mother ancestry proportion",ylab="partial residuals of embryo dry weight (g)", col=hybrid_col,pch=20,cex.lab=1.8,cex.axis=1.5)
abline(lm(CALL_hi_mlm_vr$res$visregRes~CALL_hi_mlm_vr$res$mother_hi), col="purple",lwd=4)
dev.copy(pdf,"Figures/FigSX_CALL_embryo-size_mother-hybrid-index.pdf",bg="transparent",width=8,height=8)
dev.off()
cor(CALL_hi_mlm_vr$res$visregRes,CALL_hi_mlm_vr$res$mother_hi,method="pearson") # 0.9877804

## Supp Fig S6B - Relationship between embryo dry weight and mother mitotype
# model fit
model_all<-lm(embryo_dry_weight_g~mitotype*stage+brood_size,data=CALL_embryo_data_nomal)
selectedMod <- step(model_all) # embryo_dry_weight_g ~ mitotype + stage + mitotype:stage + brood_size
summary(model_all)

# evaluate significance with likelihood ratio test (anova) on full and reduced model
CALL_mito_mlm_fit <- lme(embryo_dry_weight_g ~ mitotype + stage + brood_size,random=~1|brood_ID,data=CALL_embryo_data_nomal_sub)
CALL_mito_mlm_fit_reduced <- lme(embryo_dry_weight_g ~ stage + brood_size,random=~1|brood_ID,data=CALL_embryo_data_nomal_sub)
anova(CALL_mito_mlm_fit,CALL_mito_mlm_fit_reduced)
# CALL_mito_mlm_fit_reduced     2  8 -10943.09 -10904.58 5479.546 1 vs 2 11.81378   6e-04

# use lmer to fit mixed linear model for visreg
CALL_mito_mlm_fit<-lmer(embryo_dry_weight_g~mitotype+stage+brood_size+(1|brood_ID),data=CALL_embryo_data_nomal_sub)

# plot mitotype and embryo dry weight
CALL_mito_mlm_vr <- visreg(CALL_mito_mlm_fit,"mitotype")
plot(visregRes ~ mitotype,CALL_mito_mlm_vr$res,xlab="mitotype",ylab="embryo dry weight (g)", col=hybrid_col,pch=20)

CALL_mito_plot <-
  ggplot(CALL_mito_mlm_vr$res, aes(x=mitotype, y=visregRes, fill=mitotype, color=mitotype)) +
  scale_color_manual(values=c(hybrid_col,hybrid_col)) +
  scale_fill_manual(values=c(hybrid_col,hybrid_col)) +
  geom_jitter(alpha=0.6,width=0.1, cex=0.7) +
  stat_summary(fun.data = mean_se, width=0.4, size=0.8, geom = "crossbar", fill=NA) +
  theme_minimal() +
  theme(
    legend.position="none",
    legend.title=element_blank(),
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.1),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  xlab("mitotype") +
  ylab("partial residuals of embryo dry weight (g)")
CALL_mito_plot
ggsave(CALL_mito_plot,filename='Figures/FigSX_CALL_embryo-size_mitotype.pdf',height=8,width=8,bg = "transparent")


## Supp Fig S7: Allele-specific expression of X. malinche ovaries from F1 cross
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

indiv1<-read.csv(file="~/Data/ASE_overlap_AI_MxB1ov.R1_val_1.fq.gz.par1.bam.vcf_mod_rmdup_ASE_counts",sep="\t",head=FALSE)
indiv1<-subset(indiv1,(indiv1$V5+indiv1$V6) >=20)
indiv2<-read.csv(file="~/Data/ASE_overlap_AI_MxB2ov.R1_val_1.fq.gz.par1.bam.vcf_mod_rmdup_ASE_counts",sep="\t",head=FALSE)
indiv2<-subset(indiv2,(indiv2$V5+indiv2$V6) >=20)

prop_xmal_indiv1<- indiv1$V6/(indiv1$V6+indiv1$V5)
prop_xmal_indiv2<- indiv2$V6/(indiv2$V6+indiv2$V5)

plot(1:2,c(mean(prop_xmal_indiv1),mean(prop_xmal_indiv2)),col="cornflowerblue",xlab="Individual", ylab="Proportion of malinche alleles",pch=20,cex=2,ylim=c(0,1.02),xlim=c(0.5,2.5),xaxt="n")

noise=runif(length(prop_xmal_indiv1),-0.3,0.3)
points(rep(1,length(prop_xmal_indiv1))+noise,prop_xmal_indiv1,col=rgb(154/255,155/255,156/255,alpha=0.01),pch=20,cex=2)

noise=runif(length(prop_xmal_indiv2),-0.3,0.3)
points(rep(2,length(prop_xmal_indiv2))+noise,prop_xmal_indiv2,col=rgb(154/255,155/255,156/255,alpha=0.01),pch=20,cex=2)

points(1:2,c(mean(prop_xmal_indiv1),mean(prop_xmal_indiv2)),col="cornflowerblue",xlab="Individual", ylab="Proportion of malinche alleles",pch=20,cex=2,ylim=c(0,1),xlim=c(0.5,2.5),xaxt="n")

error.bar(c(1:2),c(mean(prop_xmal_indiv1),mean(prop_xmal_indiv2)),c(2*std.error(prop_xmal_indiv1),2*std.error(prop_xmal_indiv2)),lwd=2,col="cornflowerblue")


## Supp Fig S8: Artery development genes from WGCNA cluster associated with late-pregnancy X. malinche ovaries
## Code can be found under Figure 3B


## Supp Fig S9: Acyl- and acetyl-CoA genes from WGCNA cluster associated with late-pregnancy X. malinche ovaries
## Code can be found under Figure 3B


## Supp Fig S10: Heatmap of WGCNA module-trait correlations matrix
## Generated in WGCNA analysis, see ovary_WGCNA_xmac-IDs_combined2023.R


## Supp Fig S11: CHIC and ACUA temperature data
## Adapted from Payne et al 2022 temperature plot
library(ggplot2)
hobo_data <- read.csv("Data/HOBO_CHIC_ACUA_for_Payne-et-al-2021.csv")
datetime <- as.POSIXct(strptime(hobo_data$Datetime_GMT.06.00_CHIC_01, "%m/%d/%y %H:%M", tz = "CST6CDT"))
malcol <- "#00BFC4"
bircol <- "#F8766D"
text_col <- "black"
temp_plot <- ggplot() + theme_minimal() +
  theme(
    axis.text=element_text(color=text_col,size=22),
    text=element_text(color=text_col,size=22)
  ) +
  geom_point(data = hobo_data, aes(x = datetime, y = CHIC_01_C), colour=adjustcolor(malcol,alpha.f=0.2), size = 0.9) +
  geom_smooth(data = hobo_data, aes(x = datetime, y = CHIC_01_C), fill=malcol,colour=malcol,linewidth=1) +
  geom_point(data = hobo_data, aes(x = datetime, y = ACUA_01.16.06.17), colour=adjustcolor(bircol,alpha.f=0.2), size = 0.9) +
  geom_smooth(data = hobo_data, aes(x = datetime, y = ACUA_01.16.06.17), fill=bircol,colour=bircol,linewidth=1) +
  scale_x_datetime(date_breaks = "2 months", date_minor_breaks = "1 month", date_labels="%b") +
  scale_y_continuous(name="temperature (\u00B0C)", limits=c(10, 30)) +
  xlab("year-month-day-time (GMT-6)")
temp_plot
ggsave(temp_plot,filename="Figures/HOBO_CHIC-2020_ACUA-2016_combined-plot.pdf",width=12,height=8.5)


## Supp Fig S12 - Pregnancy rates in wild X. malinche and X. birchmanni populations
# load pregnancy rate data
pregnancy_data <- read.csv("Data/pregnancy_rate_CHIC_COAC_collections.csv",header=T)

# set colors
preg_col <- "saddle brown"
nonpreg_yolk_col    <- "goldenrod"
nonpreg_nonyolk_col <- "wheat"
color_list <- c(preg_col,nonpreg_yolk_col,nonpreg_nonyolk_col)

# reorder pregnancy status
status_order <- c("pregnant","nonpregnant_yolked","nonpregnant_nonyolked")

## split into two barplots, one for X. malinche and one for X. birchmanni
## unfortunately, was not able to stitch the two plots together with
## facet or ggarrange, something weird with creating two barplots at the same time?
## stitched manually
# Xbir plot
pregnancy_data_xbir <- subset(pregnancy_data, species=="Xbirchmanni")
collection_order <- c("COAC_II-23", "COAC_V-22", "COAC_VIII-20", "COAC_IX-22")
Xbir_preg_rate_plot <-
  ggplot(pregnancy_data_xbir, aes(x = factor(collection,levels=collection_order), y = number_of_females, fill = factor(pregnancy_status, levels=status_order))) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values=color_list) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  scale_x_discrete(labels = c("February", "May", "August", "November")) +
  xlab("X. birchmanni") +
  ylab("percentage of females collected")
Xbir_preg_rate_plot
ggsave(Xbir_preg_rate_plot,filename='Figures/SuppFigS12_Xbir_pregnancy-rates.pdf',bg = "transparent",width=11,height=8.5)

# Xmal plot
pregnancy_data_xmal <- subset(pregnancy_data, species=="Xmalinche")
collection_order <- c("CHIC_II-23", "CHIC_V-22", "CHIC_VIII-20", "CHIC_IX-22")
Xmal_preg_rate_plot <-
  ggplot(pregnancy_data_xmal, aes(x = factor(collection,levels=collection_order), y = number_of_females, fill = factor(pregnancy_status, levels=status_order))) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values=color_list) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  scale_x_discrete(labels = c("February", "May", "August", "November")) +
  xlab("X. malinche") +
  ylab("percentage of females collected")
Xmal_preg_rate_plot
ggsave(Xmal_preg_rate_plot,filename='Figures/SuppFigS12_Xmal_pregnancy-rates.pdf',bg = "transparent",width=11,height=8.5)


## Supp Fig S13A - Food deprivation fry experiment - standard length (cm)
## chose to plot partial residuals of standard lengths, instead of lengths normalized
## by average initial lengths, because it better shows that Xmal and Xbir fry achieve
## roughly the same size after 3 days, despite starting much smaller, and shows that
## Xbir gain less in no food conditions
fry_initial_std_lengths <- read.csv("Data/X-23_fry_starvation_initial_standard_lengths.csv")
# get average pre-trial standard length per brood, add to dataframe
fry_initial_std_lengths_avg <- as.data.frame(
  fry_initial_std_lengths %>%
    group_by(brood_no) %>%
    reframe(avg_initial_std_length_mm = mean(standard_length_cm,na.rm=T)*10))
fry_fat_content <- merge(fry_fat_content,fry_initial_std_lengths_avg,by="brood_no")
# convert std length cm to mm
fry_fat_content$standard_length_mm <- fry_fat_content$standard_length_cm*10

# choose model
step(lm(standard_length_mm ~ condition+species+avg_initial_std_length_mm+brood_no,fry_fat_content))
#standard_length_mm ~ condition + brood_no

# Stats: Compare raw means, accounting for covariates
# ANOVA - posthoc Tukey
#$`pairwise differences of condition`
#1                           estimate     SE df t.ratio p.value
#Xbir_control - Xbir_starved   0.1130 0.0325 50   3.481  0.0056
#Xbir_control - Xmal_control  -0.0393 0.0397 50  -0.989  0.7562
#Xbir_control - Xmal_starved   0.0280 0.0393 50   0.712  0.8920
#Xbir_starved - Xmal_control  -0.1522 0.0384 50  -3.963  0.0013
#Xbir_starved - Xmal_starved  -0.0850 0.0380 50  -2.237  0.1274
#Xmal_control - Xmal_starved   0.0672 0.0273 50   2.462  0.0786
std_length_fit <- lme(standard_length_mm ~ condition, random=~1|brood_no,data=fry_fat_content)
emmeans(std_length_fit, list(pairwise ~ condition), adjust = "tukey")

# calculate partial residuals of standard length
std_length_fit <- lmer(standard_length_mm ~ condition+(1|brood_no),data=fry_fat_content)
std_length_partial_residuals <- visreg(std_length_fit, "condition",plot=F)$res

# plot
text_col <- "black"
color_list <- c(coaccol,coaccol,chiccol,chiccol)
condition_order <- c("Xmal_control","Xmal_starved","Xbir_control","Xbir_starved")
sl_fry_food_dep_plot <- ggplot(std_length_partial_residuals, aes(x=factor(condition,levels=condition_order), y=visregRes, fill=condition)) +
  geom_boxplot(aes(col=condition),position="identity",alpha=0,lwd=1.5) +
  geom_jitter(aes(col=condition),alpha=0.4,width=0.2) +
  scale_color_manual(values=color_list) +
  theme_minimal() +
  theme(
    legend.position="none",
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  scale_x_discrete(labels = c("food", "no food", "food", "no food")) +
  ylim(c(8,13)) +
  xlab("fry treatment (3 days)") +
  ylab("partial residuals of standard length (mm)")
sl_fry_food_dep_plot
ggsave(sl_fry_food_dep_plot,filename='Figures/FigSX_food-deprivation_standard-length-boxplot_black-text.pdf',height=8.5,width=8.5,bg = "transparent")

## Supp Fig S13A: plot partial residuals of initial standard lengths
fry_initial_std_lengths <- subset(fry_initial_std_lengths,brood_no!="1m" & brood_no!="1b")
fry_initial_std_lengths$standard_length_mm <- fry_initial_std_lengths$standard_length_cm*10
initial_length_fit <- lmer(standard_length_mm ~ species+(1|brood_no),data=fry_initial_std_lengths)
# calculate partial residuals
initial_length_partial_residuals <- visreg(initial_length_fit, "species",plot=F)$res
# plot
text_col <- "black"
color_list <- c(coaccol,chiccol)
species_order <- c("Xmal","Xbir")
initial_sl_fry_food_dep_plot <- ggplot(initial_length_partial_residuals, aes(x=factor(species,levels=species_order), y=visregRes, fill=factor(species,levels=species_order))) +
  geom_boxplot(aes(col=species),position="identity",alpha=0,lwd=1.5) +
  geom_jitter(aes(col=species),alpha=0.4,width=0.2) +
  scale_color_manual(values=color_list) +
  theme_minimal() +
  theme(
    legend.position="none",
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  ylim(c(8,13)) +
  xlab("species") +
  ylab("partial residuals of standard length (mm)")
initial_sl_fry_food_dep_plot
ggsave(initial_sl_fry_food_dep_plot,filename='Figures/FigSX_food-deprivation_initial-standard-length-boxplot_black-text.pdf',height=8.5,width=5,bg = "transparent")


## Supp Fig S13B: Food deprivation fry experiment - fat content (total %)
fry_fat_content <- read.csv("Data/X-23_fry_starvation_fat_content.csv")
# remove rows with fry content < 0; remove trial run Xbir and Xmal broods (1m & 1b); remove NAs and negative FC percentages
fry_fat_content <- subset(fry_fat_content, FC_percent >= 0 & !is.na(FC_percent) & brood_no != "1m" & brood_no != "1b")
fry_fat_content$condition <- paste(fry_fat_content$species,fry_fat_content$treatment, sep="_")

## stats: ANOVA and Tukey
library(nlme)
library(emmeans)
fry_fat_fit <- lme(FC_percent ~ condition, random=~1|brood_no,data=fry_fat_content)
summary(fry_fat_fit)
emmeans(fry_fat_fit, list(pairwise ~ condition), adjust = "tukey")
# Stats: Compare raw means, accounting for covariates
# ANOVA - posthoc Tukey
# $`pairwise differences of condition`
#                           estimate   SE df t.ratio p.value
#Xbir_control - Xbir_starved    6.492 1.39 50   4.658  0.0001
#Xbir_control - Xmal_control   -0.898 4.51 50  -0.199  0.9972
#Xbir_control - Xmal_starved    8.006 4.50 50   1.779  0.2953
#Xbir_starved - Xmal_control   -7.390 4.49 50  -1.644  0.3637
#Xbir_starved - Xmal_starved    1.514 4.48 50   0.338  0.9866
#Xmal_control - Xmal_starved    8.903 1.18 50   7.573  <.0001

# choose model
step(lm(FC_percent ~ condition+species+treatment+brood_no,data=fry_fat_content)) # FC_percent ~ condition + brood_no

# calculate partial residuals
fry_fat_fit <- lmer(FC_percent ~ condition+(1|brood_no),data=fry_fat_content)
fry_fat_partial_residuals <- visreg(fry_fat_fit, "condition",plot=F)$res

# plot
text_col <- "black"
color_list <- c(coaccol,coaccol,chiccol,chiccol)
condition_order <- c("Xmal_control","Xmal_starved","Xbir_control","Xbir_starved")
fc_fry_food_dep_plot <- ggplot(fry_fat_partial_residuals, aes(x=factor(condition,levels=condition_order), y=visregRes, fill=factor(condition,levels=condition_order))) +
  geom_boxplot(aes(col=condition),position="identity",alpha=0,lwd=1.5) +
  geom_jitter(aes(col=condition),alpha=0.4,width=0.2) +
  scale_color_manual(values=color_list) +
  theme_minimal() +
  theme(
    legend.position="none",
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  scale_x_discrete(labels = c("food", "no food", "food", "no food")) +
  ylim(-3,30) +
  xlab("fry treatment (3 days)") +
  ylab("partial residuals of % fat content")
fc_fry_food_dep_plot
ggsave(fc_fry_food_dep_plot,filename='Figures/FigSX_food-deprivation_fat-content-boxplot_black-text.pdf',height=8.5,width=8.5,bg = "transparent")


## Supp Fig S14: Photo of premature birth phenotype


## Supp Fig S15: Within and between X. malinche lab population crosses [TETIxTETI (wild), CHICxCHIC (lab), TETIxCHIC (lab)]
tetixchic_data <- read.csv("Data/III-2023_TETI2_TETIxCHIC_CHICxCHIC_embryo-dry-weights.csv",header=T,sep=",")
# we only have stage 25 and stage 35 TETIxCHIC embryos, will need to use those stages
table(subset(tetixchic_data,population=="TETIxCHIC")$stage)
tetixchic_data <- subset(tetixchic_data, stage=="25" | stage=="35")
tetixchic_data$embryo_dry_weight_g <- as.numeric(tetixchic_data$embryo_dry_weight_g)
# average within brood/stage
tetixchic_data_avg <- as.data.frame(
  tetixchic_data %>%
    group_by(brood_ID,stage,population,mother_origin,brood_size,mother_std_length) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T)))
# model selection
tetixchic_full <- lm(embryo_dry_weight_g_mean~population+stage+brood_size+mother_std_length+mother_origin, data=tetixchic_data_avg)
step(tetixchic_full) # embryo_dry_weight_g ~ population + stage + mother_std_length
## stats: ANOVA and Tukey
tetixchic_fit <- lm(embryo_dry_weight_g_mean ~ population + stage + mother_std_length,data=tetixchic_data_avg)
summary(tetixchic_fit)
emmeans(tetixchic_fit, list(pairwise ~ population), adjust = "tukey")
#CHICxCHIC - TETI2      0.001120 0.000452 11   2.477  0.0729
#CHICxCHIC - TETIxCHIC  0.000825 0.000585 11   1.409  0.3700
#TETI2 - TETIxCHIC     -0.000295 0.000378 11  -0.782  0.7215

# pull partial residuals, accounting for covariates
tetixchic_partial_residuals <- visreg(tetixchic_fit, "population",plot=F)$res

text_col <- "black"
color_list <- c("#797EF6","#45b6fe",chiccol)
group_order <- c("TETI2","TETIxCHIC","CHICxCHIC")
tetixchic_embryo_size_plot <- ggplot(tetixchic_partial_residuals, aes(x=factor(population,levels=group_order), y=visregRes, fill=factor(population,levels=group_order), color=factor(population,levels=group_order))) +
  geom_jitter(aes(col=factor(population,levels=group_order)),alpha=0.4,width=0.2) +
  scale_color_manual(values=color_list) +
  scale_color_manual(values=color_list) +

  stat_summary(fun.data = mean_se, width=0.4, size=0.8, geom = "crossbar", fill=NA) +
  theme_minimal() +
  theme(
    legend.position="none",
    axis.ticks = element_line(colour = text_col, size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour=text_col, size=22),
    text=element_text(colour=text_col, size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = text_col, fill=NA)
  ) +
  scale_x_discrete(labels = c("TETIxTETI", "TETIxCHIC", "CHICxCHIC")) +
  xlab("population") +
  ylab("partial residuals of embryo dry mass (g)")
tetixchic_embryo_size_plot
ggsave(tetixchic_embryo_size_plot,filename='Figures/FigSX_TETIxCHIC-Xmal-intraspecific-embryo-sizes_black-text.pdf',height=8.5,width=8.5,bg = "transparent")


## Supp Fig S16: PCA plot of embryo RNAseq samples
library(PCAtools)
samples <- read.table("Data/embryo_ovary_dge_combined_2023/ovary_embryo_kallisto_output_FebAug23_combined/ovary_embryo_rnaseq_samples_FebAug23_combined.txt", header = TRUE)
samples <- samples[samples$tissue == "embryo",] # subset tissue
vst <- readRDS("Data/embryo_ovary_dge_combined_2023/embryo_xmac-gtf_dge/embryo-xmacID-combined2023_dge_vst.rds")
x <- assay(vst)
rownames(samples) <- samples$sample
p <- pca(x, metadata = samples, removeVar = 0.1)
# species and sample origin (lab v wild) cluster as expected
pdf(file='Figures/FigSX_embryo_dge_xmacID_PCA.pdf',height=8.5,width=12,bg = "transparent")
biplot(p, x="PC1", y="PC2",
       colby = 'species', colkey=c('Xmalinche'=malcol,'Xbirchmanni'=bircol,'Xcortezi'=corcol),
       legendPosition = 'right', legendLabSize = 14, legendIconSize = 7,
       shape = 'origin', shapekey = c('lab'=16, 'wild'=8), labSize = 4,
       gridlines.major=FALSE, gridlines.minor=FALSE,
       drawConnectors = T,
       encircle=TRUE, encircleFill=TRUE,
       axisLabSize = 20,
       caption = '')
dev.off()


## Supp Fig S17: Stage 0 egg size comparison between X. malinche, X. birchmanni, and X. cortezi
## Code can be found under Figure 2B


## Supp Fig S18: Fry CTmin comparison
ctmin_data <- read.csv("Data/CTmin_Xmal_Xbirch_newborn_trial_data.csv",header=T)
ctmin_data$CTmin_trial <- as.factor(ctmin_data$CTmin_trial)
# choose model
ctmin_species_full <- lm(CTmin_C~species+CTmin_trial+start_temp_C+born_on_date,ctmin_data)
step(ctmin_species_full)
ctmin_species_fit <- lm(CTmin_C ~ species + start_temp_C + born_on_date,ctmin_data) # CTmin_C ~ species + start_temp_C + born_on_date
summary(ctmin_species_fit)
#              Estimate Std. Error t value Pr(>|t|)
# speciesXmal              1.1992     0.3574   3.355  0.00216 **
# start_temp_C            -8.7418     6.2718  -1.394  0.17361
# born_on_date19-XIII-22   2.0563     1.6749   1.228  0.22910
# born_on_date25-XIII-22   2.4202     1.3300   1.820  0.07879 .

# calculate partial residuals
ctmin_partial_residuals <- visreg(ctmin_species_fit, "species",plot=F)$res

# plot
malcol <- "#00BFC4"
bircol <- "#F8766D"
color_list <- c(bircol,malcol)
ctmin_species <- ggplot(ctmin_data, aes(x=species, y=CTmin_C, fill=species)) +
  geom_boxplot(position="identity",color=color_list,alpha=0,lwd=1.5) +
  geom_jitter(aes(col=species),alpha=0.4,width=0.2) +
  theme_minimal() +
  theme(
    legend.position="none",
    axis.ticks = element_line(colour = "black", size = 0.2),
    panel.grid.major = element_line(colour = "grey", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text=element_text(colour="black", size=22),
    text=element_text(colour="black", size=26),
    plot.background = element_rect(fill = "transparent", color = NA),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = "black", fill=NA)
  ) +
  scale_x_discrete(labels = c('Xbir','Xmal')) +
  xlab("species") +
  ylab("CTmin (\u00B0C)")
ctmin_species
ggsave(ctmin_species,filename='Figures/FigSX_CTmin_newborn-Xmal-Xbirch.pdf',bg = "transparent",width=8.5,height=8.5)
