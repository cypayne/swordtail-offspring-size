##### SUPPLEMENTARY TABLES
##### Payne et al
##### 12/2023

## Supp Table S1 - Raw fry size data
fry_size_data <- read.table("Data/newborn_fry_size_data.csv",header=T,sep=',')


## Supp Table S2 - Compare standard length and head width means for
## A. Five Xiphophorus species
fry_size_data <- read.table("Data/newborn_fry_size_data.csv",header=T,sep=',')
xipho_sp_size <- subset(fry_size_data, species_site %in% c("XbirCOAC","XmalCHIC","XcorPTHC","XpygPTHC","XvarCOAC"))
# average fry size within a brood
xipho_sp_size_summed <- as.data.frame(
  xipho_sp_size %>%
    group_by(brood_ID) %>%
    reframe(sl_mm_mean = mean(sl_mm,na.rm=T),species=first(species)))

# summarize mean and standard deviation for each species
xipho_sp_size_summed %>%
  group_by(species) %>%
  summarise_at(vars(sl_mm_mean),
               list(sl_mean = mean, sl_std = sd))
#   species sl_mean sl_std
# 1 Xbir       9.90  0.785
# 2 Xcor       8.37  0.551
# 3 Xmal      10.8   1.13
# 4 Xpyg       8.31  0.569
# 5 Xvar       7.35  0.549

# Pairwise comparisons using Wilcoxon rank sum exact test
pairwise.wilcox.test(xipho_sp_size_summed$sl_mm_mean, xipho_sp_size_summed$species,p.adjust.method = "BH")
#      Xbir    Xcor    Xmal   Xpyg
# Xcor 2.1e-07 -       -       -
# Xmal 0.00043 2.1e-07 -       -
# Xpyg 0.00078 0.94505 0.00043 -
# Xvar 0.00466 0.07576 0.00466 0.29630
# P value adjustment method: BH

## B. Xmal, Xbir, F1, F2, and natural hybrids
fry_size_data <- read.table("Data/newborn_fry_size_data.csv",header=T,sep=',')
mal_bir_hyb_size <- subset(fry_size_data, species_site %in% c("XbirCOAC","XmalCHIC","Xmal_xbir_F1","Xmal_xbir_F2","Xmal_xbirCALL"))
# average fry size within a brood
mal_bir_hyb_size_summed <- as.data.frame(
  mal_bir_hyb_size %>%
    group_by(brood_ID) %>%
    reframe(sl_mm_mean = mean(sl_mm,na.rm=T),species_site=first(species_site)))
# Pairwise comparisons using Wilcoxon rank sum exact test
pairwise.wilcox.test(mal_bir_hyb_size_summed$sl_mm_mean, mal_bir_hyb_size_summed$species,p.adjust.method = "BH")
#              XbirCOAC Xmal_xbir_F1 Xmal_xbir_F2 Xmal_xbirCALL
# Xmal_xbir_F1  0.0109   -            -            -
# Xmal_xbir_F2  0.0281   0.3071       -            -
# Xmal_xbirCALL 0.1909   0.1909       0.3393       -
# XmalCHIC      0.0014   0.3071       0.3480       0.3393
# P value adjustment method: BH


## Supp Table S3 - Morphological stage descriptions


## Supp Table S4 - Raw embryo size data from wild collections
combined_embryo_data <- read.table("Data/all-pops_combined_embryo_datasets.csv",header=T,sep=',')


##  Supp Table S5 - Statistical comparisons of mean embryo size between species at each stage
combined_embryo_data <- read.table("Data/all-pops_combined_embryo_datasets.csv",header=T,sep=',')
combined_subset <- subset(combined_embryo_data,(popcoll=="CHIC_V-2022" | popcoll=="CHIC_VIII_20" | popcoll=="COAC_VIII_20" | popcoll=="COAC_II-2023" | popcoll=="COAC_IX-2022" | popcoll=="PTHC_II-2023") & stage!="0" & stage!="2" & stage!="5" & stage!="?" & stage!="-")
# average embryo size for embryos of the same stage within a brood
combined_subset_summed <- as.data.frame(
  combined_subset %>%
    group_by(brood_ID,species,population,collection,stage) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T),brood_size=mean(brood_size,na.rm=T),mother_std_length=first(mother_std_length)))## pairwise comparisons with linear model + Tukey
# choose model
combined_full_model <- lm(embryo_dry_weight_g_mean~species+stage+collection+brood_size+mother_std_length,data=combined_subset_summed)
step(combined_full_model) # embryo_dry_weight_g ~ species + stage + collection + mother_std_length

# compare means between species by stage
for(stg in c("10","15","20","25","30","35","40","45","50")){
  print(stg)
  coac_chic_subset <- subset(combined_subset_summed,(species=="Xmalinche" | species=="Xbirchmanni") & !is.na(stage) & stage==stg)
  # chosen model: embryo_dry_weight_g_mean ~ species + collection + mother_std_length
  coac_chic_stage_lm <- lm(embryo_dry_weight_g_mean ~ species + collection + mother_std_length,data=coac_chic_subset)
  # emmeans pairwise comparison with Tukey
  print(emmeans(coac_chic_stage_lm,pairwise~species)$contrasts)
}


##  Supp Table S6 - Raw stage 0 egg size data from wild collections
stage0_data <- read.table("Data/CHIC_COAC_PTHC_fully-yolked-stage0_dry-weights_mother-length.csv",header=T,sep=',')


## Supp Table S7 - Natural hybrid (CALL) mother hybrid index (X. malinche ancestry proportion) and mitotype data
combined_embryo_data <- read.table("Data/all-pops_combined_embryo_datasets.csv",header=T,sep=',')
CALL_embryo_data <- subset(combined_embryo_data,population=="CALL")
CALL_mother_hi <- read.table("Data/CALL_mother_hybrid-index_mitotype.csv",header=T,sep=',')
CALL_embryo_data$mother_hi <- CALL_mother_hi$mom_index[match(CALL_embryo_data$brood_ID, CALL_mother_hi$brood_ID)]
CALL_embryo_data$mitotype <- CALL_mother_hi$mitotype[match(CALL_embryo_data$brood_ID, CALL_mother_hi$brood_ID)]
# remove early stages with few samples
CALL_embryo_data <- subset(CALL_embryo_data,stage!="-" & stage!="0" & stage!="10" & stage!="15")
# remove pure parent data
CALL_embryo_data_nomal <- subset(CALL_embryo_data,mother_hi<0.9 & mother_hi>0.1 )


## Supp Table S8 - Matrotrophy Index partial residuals (see Supp. Materials S3 section below for calculations and bootstraps)
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


## Supp Table S9 - Embryo dry weights of multifactorial lab crosses between X. malinche and X. birchmanni
roof_chic_coac_f1 <- read.csv("Data/IV-2023_roof-tank_CHIC-COAC-F1_embryo-weights.csv",header=T,sep=",")


## Supp Table S10 - Compare mean embryo size of multifactorial lab crosses between X. malinche and X. birchmanni
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


## Supp Table S11 - Embryo size data from within and between population crosses within X. malinche
tetixchic_data <- read.csv("Data/III-2023_TETI2_TETIxCHIC_CHICxCHIC_embryo-dry-weights.csv",header=T,sep=",")


## Supp Table S12 - Compare mean embryo size between within and between population crosses within X. malinche
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
# choose model
tetixchic_full <- lm(embryo_dry_weight_g_mean~population+stage+brood_size+mother_std_length+mother_origin, data=tetixchic_data_avg)
step(tetixchic_full) # embryo_dry_weight_g ~ population + stage + mother_std_length
## stats: ANOVA and Tukey
tetixchic_fit <- lm(embryo_dry_weight_g_mean ~ population + stage + mother_std_length,data=tetixchic_data_avg)
summary(tetixchic_fit)
emmeans(tetixchic_fit, list(pairwise ~ population), adjust = "tukey")
#CHICxCHIC - TETI2      0.001120 0.000452 11   2.477  0.0729
#CHICxCHIC - TETIxCHIC  0.000825 0.000585 11   1.409  0.3700
#TETI2 - TETIxCHIC     -0.000295 0.000378 11  -0.782  0.7215


## Supp Table S13 - RNAseq metadata


## Supp Table S14 - Embryo differential gene expression results
## See embryo-DGE_DESeq2_xmac-IDs_2023.R to run DGE analysis
embryo_dge <- read.csv("Data/embryo_ovary_dge_combined_2023/embryo_xmac-gtf_dge/embryo-xmacID-combined2023_dge_lfc-shr_all.csv_with-annots.csv",header=T)


## Supp Table S15 - Ovary differential gene expression results
## See ovary-DGE_DESeq2_xmac-IDs_2023.R to run DGE analysis
ovary_dge <- read.csv("Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/ovary-xmacID-combined2023_dge_lfc-shr_all.csv_with-annots.csv",header=T)


## Supp Table S16 - Gene Ontology biological pathways enriched in differentially expressed genes between X. malinche and X. birchmanni ovaries

## Supp Table S17 - Significance of relationship between WGCNA co-expressed gene cluster ovary expression profile and species/stage
## See ovary_WGCNA_xmac-IDs_combined2023.R to run WGCNA analysis
ovary_wgcna_module_pvals <- read.csv("Data/embryo_ovary_dge_combined_2023/wgcna/ovary-xmac_combined-FebAug23_MEtraitpvals.csv",header=T)



## Table X - Heritability of standard length calculations
## prepare partial residuals to account for covariates in trait values
nb_fry_data <- read.csv("Data/newborn_fry_size_data.csv")
nb_fry_data$born_year <- as.character(nb_fry_data$born_year)

subset_nb_fry_data <- subset(nb_fry_data, species_site=="XbirCOAC" | site=="CHIC" | species=="Xmal_xbir_F2" | species=="Xmal_xbir_F1")
model_all<-lm(sl_mm~species+collection+born_month+born_year,data=subset_nb_fry_data)
selectedMod <- step(model_all) # sl_mm~species+collection
summary(model_all)

subset_nb_size_fit <- lmer(sl_mm ~ species+collection+(1|brood_ID),subset_nb_fry_data)
nb_size_residuals <- visreg(subset_nb_size_fit, "species",plot=F)$res
write.csv(nb_size_residuals,"Data/newborn_standard-length_Xbir-Xmal-F1-F2_partial-residuals_aug2023.csv")

## subset data by species
coac_fry <- subset(nb_size_residuals, species=="Xbir") # 479
chic_fry <- subset(nb_size_residuals, species=="Xmal") # 157
f2_fry <- subset(nb_size_residuals, species=="Xmal_xbir_F2") # 30
f1_fry <- subset(nb_size_residuals, species=="Xmal_xbir_F1") # 83

## calculate between-population variance in sl of two parental species fry
# calculate mean across parent populations
overall_sl_mean <- mean(c(coac_fry$visregRes,chic_fry$visregRes))
# calculate mean for each parent population
coac_sl_mean <- mean(coac_fry$visregRes)
chic_sl_mean <- mean(chic_fry$visregRes)
# get number of data points per parent population
coac_n <- length(coac_fry$visregRes)
chic_n <- length(chic_fry$visregRes)

## calculate variance between the means of the populations
# between parent population variance
between_var_sl_par <- var(c(coac_sl_mean,chic_sl_mean)) # 0.379264

## calculate within-population variance in size of parent, F1, and F2 fry
var_sl_coac <- var(coac_fry$visregRes) # 0.2047061
var_sl_chic <- var(chic_fry$visregRes) # 0.2222235
var_sl_f1 <- var(f1_fry$visregRes)     # 0.6866626
var_sl_f2 <- var(f2_fry$visregRes)     # 0.6179587
## within-F2 variance = 0.6179587

## calculate broad-sense heritability H^2 with F1 fry size data
# calculate environmental variance with Wright's weighted average of within-strain variance
var_sl_env <- var_sl_coac/4 + var_sl_chic/4 + var_sl_f1/2
## environmental variance = 0.4500637

# Calculate H^2 = (var_F2 - var_enviroment)/var_F2
h2 <- (var_sl_f2 - var_sl_env) / var_sl_f2
h2
## h^2 = 0.2716928


## Relationship between mother ancestry proportion and mitotype with embryo and ovary dry weights in hybrids from
## Natural hybrid population Calnali Low
combined_embryo_data <- read.table("Data/all-pops_combined_embryo_datasets.csv",header=T,sep=',')
CALL_embryo_data <- subset(combined_embryo_data,population=="CALL")
CALL_mother_hi <- read.table("Data/CALL_mother_hybrid-index_mitotype.csv",header=T,sep=',')
CALL_embryo_data$mother_hi <- CALL_mother_hi$mom_index[match(CALL_embryo_data$brood_ID, CALL_mother_hi$brood_ID)]
CALL_embryo_data$mitotype <- CALL_mother_hi$mitotype[match(CALL_embryo_data$brood_ID, CALL_mother_hi$brood_ID)]
# remove early stages with few samples
CALL_embryo_data <- subset(CALL_embryo_data,stage!="-" & stage!="0" & stage!="10" & stage!="15")
# remove pure parent data
CALL_embryo_data_nomal <- subset(CALL_embryo_data,mother_hi<0.9 & mother_hi>0.1 )

## embryo dry weight and mother hybrid index
# model fit
model_all<-lm(embryo_dry_weight_g~mother_hi*stage+brood_size+brood_ID,data=CALL_embryo_data_nomal)
selectedMod <- step(model_all) # embryo_dry_weight_g ~ mother_hi + stage + brood_ID + mother_hi:stage
summary(model_all)
# 2    879 0.00025928 -4 -4.6212e-06 3.9696 0.003368 **

# evaluate significance with likelihood ratio test (anova) on full and reduced model
CALL_hi_mlm_fit<-lm(embryo_dry_weight_g~mother_hi+stage+mother_hi:stage+brood_ID,data=CALL_embryo_data_nomal)
CALL_hi_mlm_fit_reduced<-lm(embryo_dry_weight_g~stage+brood_ID,data=CALL_embryo_data_nomal)
anova(CALL_hi_mlm_fit,CALL_hi_mlm_fit_reduced)

summary(CALL_hi_mlm_fit)

# plot mother hi and embryo dry weight
CALL_hi_mlm_vr <- visreg(CALL_hi_mlm_fit,"mother_hi")
plot(visregRes ~ mother_hi,CALL_hi_mlm_vr$res,xlab="mother ancestry proportion",ylab="partial residuals of embryo dry weight (g)", col=hybrid_col,pch=20,cex.lab=1.8,cex.axis=1.5)
abline(lm(CALL_hi_mlm_vr$res$visregRes~CALL_hi_mlm_vr$res$mother_hi), col="purple",lwd=4)
dev.copy(pdf,"Figures/FigSX_CALL_embryo-size_mother-hybrid-index.pdf",bg="transparent",width=8,height=8)
dev.off()
cor(CALL_hi_mlm_vr$res$visregRes,CALL_hi_mlm_vr$res$mother_hi,method="pearson") # 0.9877804

## embryo dry weight and mother mitotype
# model fit
model_all<-lm(embryo_dry_weight_g~mitotype*stage+brood_size+brood_ID,data=CALL_embryo_data_nomal)
selectedMod <- step(model_all) # embryo_dry_weight_g ~ mitotype + stage + brood_ID + mitotype:stage
summary(model_all)

# evaluate significance with likelihood ratio test (anova) on full and reduced model
CALL_mito_mlm_fit <- lm(embryo_dry_weight_g ~ mitotype + stage + mitotype:stage + brood_ID,CALL_embryo_data_nomal)
CALL_mito_mlm_fit_reduced <- lm(embryo_dry_weight_g ~ stage + brood_ID, CALL_embryo_data_nomal)
anova(CALL_mito_mlm_fit,CALL_mito_mlm_fit_reduced)
# 2    879 0.00025928 -3 -6.0832e-06 7.0155 0.0001153 ***

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

## ovary dry weight and mother hybrid index
# model fit
CALL_ovary_weight_unique <- as.data.frame(
  CALL_embryo_data_nomal %>%
    group_by(brood_ID,stage) %>%
    reframe(ovarian_tissue_dry_weight_g = first(ovarian_tissue_dry_weight_g),brood_size=first(brood_size),mother_hi=first(mother_hi),mitotype=first(mitotype)))

# model fit: mother hybrid index
CALL_hi_ovary_lm_fit <- lm(ovarian_tissue_dry_weight_g ~ mother_hi*stage + brood_size, CALL_ovary_weight_unique)
selectedMod <- step(CALL_hi_ovary_lm_fit) # ovarian_tissue_dry_weight_g ~ mother_hi + stage + mother_hi:stage + brood_size

# evaluate significance with likelihood ratio test (anova) on full and reduced model
CALL_hi_ovary_lm_fit <- lm(ovarian_tissue_dry_weight_g ~ mother_hi + stage + mother_hi:stage + brood_size, CALL_ovary_weight_unique)
summary(CALL_hi_ovary_lm_fit)
CALL_hi_ovary_lm_fit_reduced <- lm(ovarian_tissue_dry_weight_g ~ stage + brood_size, CALL_ovary_weight_unique)
anova(CALL_hi_ovary_lm_fit,CALL_hi_ovary_lm_fit_reduced)
# 2     75 0.0016711 -5 -0.00028619 2.8931 0.01971 *

# plot mother hi and embryo dry weight
CALL_hi_ovary_lm_vr <- visreg(CALL_hi_ovary_lm_fit,"mother_hi")
plot(visregRes ~ mother_hi,CALL_hi_ovary_lm_vr$res,xlab="mother ancestry proportion",ylab="partial residuals of ovary dry weight (g)", col=hybrid_col,pch=20,cex.lab=1.8,cex.axis=1.5)
abline(lm(CALL_hi_ovary_lm_vr$res$visregRes~CALL_hi_ovary_lm_vr$res$mother_hi), col="purple",lwd=4)

## ovary dry weight and mother mitotype
# model fit
CALL_mito_ovary_lm_fit <- lm(ovarian_tissue_dry_weight_g ~ mitotype*stage + brood_size, CALL_ovary_weight_unique)
selectedMod <- step(CALL_mito_ovary_lm_fit) # ovarian_tissue_dry_weight_g ~ stage + brood_size
summary(CALL_mito_ovary_lm_fit) # mitotypemalinche          1.210e-03  2.620e-03   0.462  0.64561
CALL_mito_ovary_lm_fit <- lm(ovarian_tissue_dry_weight_g ~ mitotype + stage + brood_size, CALL_ovary_weight_unique)
CALL_mito_ovary_lm_fit_reduced <- lm(ovarian_tissue_dry_weight_g ~ stage + brood_size, CALL_ovary_weight_unique)
anova(CALL_mito_ovary_lm_fit,CALL_mito_ovary_lm_fit_reduced)
# 2     75 0.0016711 -1 -6.2071e-07 0.0275 0.8687



## Supp. Materials S3: Matrotrophy Index bootstraps
## X. malinche dev profile plot
combined_embryo_data <- read.table("Data/all-pops_combined_embryo_datasets.csv",header=T,sep=',')
# average embryos at the same stage in the same brood, and filter out unused stages (0-5)
combined_subset <- subset(combined_embryo_data,(popcoll=="CHIC_V-2022" | popcoll=="CHIC_VIII_20" | popcoll=="COAC_VIII_20" | popcoll=="COAC_II-2023" | popcoll=="COAC_IX-2022") & stage!="0" & stage!="2" & stage!="5" & stage!="?" & stage!="-" & !is.na(embryo_dry_weight_g))
combined_subset_avg <- as.data.frame(
  combined_subset %>%
    group_by(brood_ID,species,stage) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T), popcoll=first(popcoll)))
combined_subset_avg

# Bootstrap MI for CHIC Xmal population, by post-2020 collections
chic_subset_avg <- subset(combined_subset_avg,popcoll=="CHIC_V-2022" | popcoll=="CHIC_VIII_20")
table(chic_subset_avg$stage)
chic_avg_stage10 <- subset(chic_subset_avg,stage==10)
chic_avg_stage45 <- subset(chic_subset_avg,stage==45)
chic_mi_bootstrap <- c()
for(i in 1:100) {
  resample_stage10 <- sample(chic_avg_stage10$embryo_dry_weight_g_mean, replace = TRUE)
  mean_stage10 <- mean(resample_stage10)
  resample_stage45 <- sample(chic_avg_stage45$embryo_dry_weight_g_mean, replace = TRUE)
  mean_stage45 <- mean(resample_stage45)
  chic_mi <- mean_stage45/mean_stage10
  chic_mi_bootstrap <- c(chic_mi_bootstrap,chic_mi)
}
mean(chic_mi_bootstrap) # 0.8581601
t.test(chic_mi_bootstrap,conf.level = 0.95) #  0.8419940 0.8743262

chic_avg_stage10 <- subset(chic_subset_avg,stage==10)
chic_avg_stage50 <- subset(chic_subset_avg,stage==50)
chic_mi_bootstrap <- c()
for(i in 1:100) {
  resample_stage10 <- sample(chic_avg_stage10$embryo_dry_weight_g_mean, replace = TRUE)
  mean_stage10 <- mean(resample_stage10)
  resample_stage50 <- sample(chic_avg_stage50$embryo_dry_weight_g_mean, replace = TRUE)
  mean_stage50 <- mean(resample_stage50)
  chic_mi <- mean_stage50/mean_stage10
  chic_mi_bootstrap <- c(chic_mi_bootstrap,chic_mi)
}
mean(chic_mi_bootstrap) # 0.947186
t.test(chic_mi_bootstrap,conf.level = 0.95) #  0.9293548 0.9650173

# Bootstrap MI for COAC Xbir population, by post-2020 collections
coac_subset_avg <- subset(combined_subset_avg,popcoll=="COAC_VIII_20" | popcoll=="COAC_II-2023" | popcoll=="COAC_IX-2022")
table(coac_subset_avg$stage)
coac_avg_stage10 <- subset(coac_subset_avg,stage==10)
coac_avg_stage45 <- subset(coac_subset_avg,stage==45)
coac_mi_bootstrap <- c()
for(i in 1:100) {
  resample_stage10 <- sample(coac_avg_stage10$embryo_dry_weight_g_mean, replace = TRUE)
  mean_stage10 <- mean(resample_stage10)
  resample_stage45 <- sample(coac_avg_stage45$embryo_dry_weight_g_mean, replace = TRUE)
  mean_stage45 <- mean(resample_stage45)
  coac_mi <- mean_stage45/mean_stage10
  coac_mi_bootstrap <- c(coac_mi_bootstrap,coac_mi)
}
mean(coac_mi_bootstrap) # 0.7735786
t.test(coac_mi_bootstrap,conf.level = 0.95) # 0.7668770 0.7802803



## Table S21: Food deprivation experiment results

## Fat content (total %)
fry_fat_content <- read.csv("Data/X-23_fry_starvation_fat_content.csv")
# remove rows with fry content < 0; remove trial run Xbir and Xmal broods (1m & 1b); remove NAs and negative FC percentages
fry_fat_content <- subset(fry_fat_content, FC_percent >= 0 & !is.na(FC_percent) & brood_no != "1m" & brood_no != "1b")
fry_fat_content$condition <- paste(fry_fat_content$species,fry_fat_content$treatment, sep="_")

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

# select model
step(lm(FC_percent ~ condition+species+treatment+brood_no,data=fry_fat_content))
# FC_percent ~ condition + brood_no

# Plot data, accounting for covariates: calculate partial residuals
fry_fat_fit <- lmer(FC_percent ~ condition+(1|brood_no),data=fry_fat_content)
fry_fat_partial_residuals <- visreg(fry_fat_fit, "condition",plot=F)$res
# Paper plot
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

## Standard length (cm)
## chose to plot partial residuals of standard lengths, instead of lengths normalized
## by average initial lengths, because it better shows that Xmal and Xbir fry achieve
## roughly the same size after 3 days, despite starting much smaller, and shows that
## xbir gain less in no food conditions
fry_initial_std_lengths <- read.csv("Data/X-23_fry_starvation_initial_standard_lengths.csv")
# get average pre-trial standard length per brood, add to dataframe
fry_initial_std_lengths_avg <- as.data.frame(
  fry_initial_std_lengths %>%
    group_by(brood_no) %>%
    reframe(avg_initial_std_length_mm = mean(standard_length_cm,na.rm=T)*10))
fry_fat_content <- merge(fry_fat_content,fry_initial_std_lengths_avg,by="brood_no")
# convert std length cm to mm
fry_fat_content$standard_length_mm <- fry_fat_content$standard_length_cm*10

# normalize standard lengths by avg length pre-trial
fry_fat_content$std_length_normalized <- fry_fat_content$standard_length_mm / fry_fat_content$avg_initial_std_length_mm
# calculate growth rate: (post trial standard length - pre-trial average lengths) / 3 days
fry_fat_content$std_length_growth_rate <- (fry_fat_content$standard_length_mm - fry_fat_content$avg_initial_std_length_mm) / 3
std_length_growth_rates_avg <- as.data.frame(
  fry_fat_content %>%
    group_by(condition) %>%
    reframe(std_length_growth_rate = mean(std_length_growth_rate,na.rm=T)))
std_length_growth_rates_avg

step(lm(standard_length_mm ~ condition+species+avg_initial_std_length_mm+brood_no,fry_fat_content))
#standard_length_mm ~ condition + brood_no

# Stats: Compare raw means, accounting for covariates
# ANOVA - posthoc Tukey
# no difference between Xmal v Xbir
# no difference between Xmal control v starved
# sig difference (decrease) between Xbir control v starved
#$`pairwise differences of condition`
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

# Plot data, accounting for covariates: calculate partial residuals
std_length_fit <- lmer(standard_length_mm ~ condition+(1|brood_no),data=fry_fat_content)
std_length_partial_residuals <- visreg(std_length_fit, "condition",plot=F)$res
# Paper plot: partial residuals of standard length
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

## plot partial residuals of initial standard lengths
fry_initial_std_lengths <- subset(fry_initial_std_lengths,brood_no!="1m" & brood_no!="1b")
fry_initial_std_lengths$standard_length_mm <- fry_initial_std_lengths$standard_length_cm*10
initial_length_fit <- lmer(standard_length_mm ~ species+(1|brood_no),data=fry_initial_std_lengths)
initial_length_partial_residuals <- visreg(initial_length_fit, "species",plot=F)$res
# Paper plot: partial residuals of standard length
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

## Dry mass (g)
step(lm(dry_mass_1 ~ condition+species+treatment+brood_no,fry_fat_content))
# dry_mass_1 ~ condition + brood_no

# Stats: Compare raw means, accounting for covariates
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
# Plot data, accounting for covariates: calculate partial residuals
dry_mass_fit <- lmer(standard_length_cm ~ condition+(1|brood_no),data=fry_fat_content)
dry_mass_partial_residuals <- visreg(dry_mass_fit, "condition",plot=F)$res
# Paper plot: partial residuals of standard length
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

x <- data.frame("no"=fry_fat_content$brood_no,"condition"=fry_fat_content$condition,"mass"=fry_fat_content$dry_mass_1)
subset(x,condition=="Xmal_starved")
subset(x,condition=="Xmal_starved" & mass<0.0033)
subset(x,condition=="Xmal_starved" & mass>0.0033)

write.csv(fry_fat_content,"Data/TableSX_food-deprivation-expt_data.csv")





