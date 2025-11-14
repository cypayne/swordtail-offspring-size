##### SUPPLEMENTARY TABLES (& MATERIALS)
##### Recent evolution of large offspring size and post-fertilization nutrient provisioning in swordtails
##### Payne et al 2025
##### 07/2025

library(lme4)
library(emmeans)

## Supp Table S1 - Raw fry size data
fry_size_data <- read.table("Data/newborn_fry_size_data.csv",header=T,sep=',')


## Supp Table S2 - Compare standard length and head width means for
## A. Five Xiphophorus species
fry_size_data <- read.table("Data/newborn_fry_size_data.csv",header=T,sep=',')
# subset species of interest and remove samples with missing data
xipho_sp_size <- subset(fry_size_data, species_site %in% c("XbirCOAC","XmalCHIC","XcorPTHC","XpygPTHC","XvarCOAC"))
xipho_sp_size <- subset(xipho_sp_size,!is.na(sl_mm) & !is.na(hw_mm) & !is.na(species) & !is.na(brood_ID))

# summarize mean and standard deviation for each species
xipho_sp_size %>%
  group_by(species) %>%
  summarise_at(vars(sl_mm, hw_mm),
               list(mean = mean, std = sd))

# lmm model for each morphometric, including brood ID as a random effect
# then make pairwise post hoc multiple comparisons between estimated marginal means
# with paired t test and tukey adjustment
# standard length
lmm.sl <- lme4::lmer(sl_mm ~ species + (1 | brood_ID), data = xipho_sp_size, REML=T )
emmeans(lmm.sl, list(pairwise ~ species), adjust = "tukey")
# head width
lmm.hw <- lme4::lmer(hw_mm ~ species + (1 | brood_ID), data = xipho_sp_size, REML=T )
emmeans(lmm.hw, list(pairwise ~ species), adjust = "tukey")

## B. Xmal, Xbir, F1, F2, and natural hybrids
fry_size_data <- read.table("Data/newborn_fry_size_data.csv",header=T,sep=',')
mal_bir_hyb_size <- subset(fry_size_data, species_site %in% c("XbirCOAC","XmalCHIC","Xmal_xbir_F1","Xmal_xbir_F2","Xmal_xbirCALL"))
mal_bir_hyb_size <- subset(mal_bir_hyb_size, !is.na(sl_mm) & !is.na(hw_mm) & !is.na(species) & !is.na(brood_ID))

# summarize mean and standard deviation for each species
mal_bir_hyb_size %>%
  group_by(species) %>%
  summarise_at(vars(sl_mm, hw_mm),
               list(mean = mean, std = sd))

# lmm model for each morphometric, including brood ID as a random effect
# then make pairwise post hoc multiple comparisons between estimated marginal means
# with paired t test and tukey adjustment
# standard length
lmm.sl <- lme4::lmer(sl_mm ~ species + (1 | brood_ID), data = mal_bir_hyb_size, REML=T )
emmeans(lmm.sl, list(pairwise ~ species), adjust = "tukey")
# head width
lmm.hw <- lme4::lmer(hw_mm ~ species + (1 | brood_ID), data = mal_bir_hyb_size, REML=T )
emmeans(lmm.hw, list(pairwise ~ species), adjust = "tukey")


## Supp Table S3 - Morphological stage descriptions


## Supp Table S4 - Raw embryo size data from wild collections
combined_embryo_data <- read.table("Data/all-pops_combined_embryo_datasets.csv",header=T,sep=',')


##  Supp Table S5 - Statistical comparisons of mean embryo size between species at each stage
combined_embryo_data <- read.table("Data/all-pops_combined_embryo_datasets.csv",header=T,sep=',')
combined_subset <- subset(combined_embryo_data,(popcoll=="CHIC_V-2022" | popcoll=="CHIC_VIII_20" | popcoll=="COAC_VIII_20" | popcoll=="COAC_II-2023" | popcoll=="COAC_IX-2022" | popcoll=="COAC_X-2023") & stage!="0" & stage!="2" & stage!="5" & stage!="?" & stage!="-")

# select best model by LRT
# include mother std length, season, collection date, and brood ID as a random effect
combined_subset <- subset(combined_subset,!is.na(mother_std_length))
fit_full <- lme4::lmer(embryo_dry_weight_g ~ species + stage + species:stage + mother_std_length + species:mother_std_length + season + collection + brood_size + (1 | brood_ID), data = combined_subset, REML=F)
# removed one variable at a time
fit_reduced <- lme4::lmer(embryo_dry_weight_g ~ species + stage + species:stage + mother_std_length + species:mother_std_length + season + collection + brood_size + (1 | brood_ID), data = combined_subset, REML=F)
anova(fit_full,fit_reduced)
# variable                  Pr(>Chisq)
# brood_size                0.81      [dropped]
# collection                0.27      [dropped]
# mother_std_length         0.04
# species:mother_std_length 0.85      [dropped]

# chosen fit
combined_fit <- lme4::lmer(embryo_dry_weight_g ~ species + stage + species:stage + mother_std_length + season + (1|brood_ID),data=combined_subset, REML=T)
# emmeans paired t-test and Tukey post-hoc adjustment  to compare species by stage
emmeans(combined_fit, list(pairwise ~ species | stage), adjust = "tukey")


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


## Supp Table S8 - Matrotrophy Index partial residuals
## See Supp. Materials S3 section below for partial residual, MI, and bootstrap calculations)


## Supp Table S9 - Embryo dry weights of multifactorial lab crosses between X. malinche and X. birchmanni
roof_chic_coac_f1 <- read.csv("Data/IV-2023_roof-tank_CHIC-COAC-F1_embryo-weights.csv",header=T,sep=",")


## Supp Table S10 - Compare mean embryo size of multifactorial lab crosses between X. malinche and X. birchmanni
roof_chic_coac_f1 <- read.csv("Data/IV-2023_roof-tank_CHIC-COAC-F1_embryo-weights.csv",header=T,sep=",")
roof_chic_coac_f1$stage <- as.factor(roof_chic_coac_f1$stage)
# drop early stages (have <=2 samples)
roof_chic_coac_f1_subset <- subset(roof_chic_coac_f1,stage!="0" & stage!="5" & stage!="10" & stage!="15" & stage!="20" & stage != "25")
# remove one outlier - likely erroneous data recording
roof_chic_coac_f1_subset <- subset(roof_chic_coac_f1_subset,embryo_dry_weight_g < 0.008)

# select best model by LRT
# include mother std length, brood size, and brood ID as a random effect
# first evaluate whether brood_ID should be included as a random effect
fit_full <- lme4::lmer(embryo_dry_weight_g ~ species + stage + species:stage + mother_std_length + species:mother_std_length + brood_size + (1 | brood_ID), data = roof_chic_coac_f1_subset, REML=F)
fit_reduced <- lm(embryo_dry_weight_g ~ species + stage + species:stage + mother_std_length + species:mother_std_length + brood_size, data = roof_chic_coac_f1_subset)
anova(fit_full,fit_reduced)
# variable                  Pr(>Chisq)
# (1|brood_ID)              1 [dropped]
# since random effect is dropped, switching to linear model
fit_full <- lm(embryo_dry_weight_g ~ species + stage + species:stage + mother_std_length + brood_size, data = roof_chic_coac_f1_subset)
# removed one variable at a time
fit_reduced <- lm(embryo_dry_weight_g ~ species + stage + species:stage + mother_std_length + brood_size, data = roof_chic_coac_f1_subset)
anova(fit_full,fit_reduced)
# variable                  Pr(>Chisq)
# stage                     0.06 [kept]
# species:stage             0.17 [dropped]
# mother_std_length         0.0002
# brood_size                1.2*10-7

# chosen fit
lab_cross_fit <- lm(embryo_dry_weight_g~species+stage+mother_std_length+brood_size,data=roof_chic_coac_f1_subset)

# Stats: compare means with ANOVA and Tukey
emmeans(lab_cross_fit, pairwise ~ species, adjust = "tukey")
# contrast                     estimate          SE df t.ratio p.value
# birxbir - birxmal_F1     0.000667 0.000712 112   0.936  0.7855
# birxbir - malxbir_F1    -0.000847 0.000639 112  -1.325  0.5491
# birxbir - malxmal       -0.001723 0.000426 112  -4.048  0.0005
# birxmal_F1 - malxbir_F1 -0.001514 0.001151 112  -1.315  0.5553
# birxmal_F1 - malxmal    -0.002390 0.000890 112  -2.684  0.0411
# malxbir_F1 - malxmal    -0.000876 0.000366 112  -2.393  0.0844


## Supp Table S11 - Embryo size data from within and between population crosses within X. malinche
tetixchic_data <- read.csv("Data/III-2023_TETI2_TETIxCHIC_CHICxCHIC_embryo-dry-weights.csv",header=T,sep=",")


## Supp Table S12 - Compare mean embryo size between within and between population crosses within X. malinche
tetixchic_data <- read.csv("Data/III-2023_TETI2_TETIxCHIC_CHICxCHIC_embryo-dry-weights.csv",header=T,sep=",")
# we only have stage 25 and stage 35 TETIxCHIC embryos, will need to use those stages
table(subset(tetixchic_data,population=="TETIxCHIC")$stage)
tetixchic_data <- subset(tetixchic_data, stage=="25" | stage=="30" | stage=="35")
tetixchic_data$embryo_dry_weight_g <- as.numeric(tetixchic_data$embryo_dry_weight_g)

# select best model by LRT
# include mother std length, mother origin brood size, and brood ID as a random effect
# first evaluate whether brood_ID should be included as a random effect
fit_full <- lme4::lmer(embryo_dry_weight_g ~ population + stage + population:stage + mother_std_length + mother_origin + brood_size + (1 | brood_ID), data = tetixchic_data, REML=F)
fit_reduced <- lm(embryo_dry_weight_g ~ population + stage + population:stage + mother_std_length + mother_origin + brood_size, data = tetixchic_data)
anova(fit_full,fit_reduced)
# variable                  Pr(>Chisq)
# (1|brood_ID)              0.04 [kept]
# keeping mixed linear model with brood_ID as random effect
fit_full <- lme4::lmer(embryo_dry_weight_g ~ population + stage + population:stage + mother_std_length + mother_origin + brood_size + (1 | brood_ID), data = tetixchic_data, REML=F)
fit_reduced <- lme4::lmer(embryo_dry_weight_g ~ population + stage + population:stage + mother_std_length + mother_origin +  (1 | brood_ID), data = tetixchic_data, REML=F)
anova(fit_full,fit_reduced)
# variable                  Pr(>Chisq)
# stage                     0.21  [dropped]
# population:stage          0.09  [dropped]
# mother_std_length         0.25  [dropped]
# mother_origin             NA    [dropped]
# brood_size                0.97  [dropped]

# chosen fit
tetixchic_fit <- lmer(embryo_dry_weight_g ~ population + (1|brood_ID),data=tetixchic_data, REML=T)

# Stats: compare means with ANOVA and Tukey
emmeans(tetixchic_fit, pairwise ~ population, adjust = "tukey")
# contrast               estimate       SE   df t.ratio p.value
# CHICxCHIC - TETI2      0.001770 0.000385 11.4   4.598  0.0018
# CHICxCHIC - TETIxCHIC  0.001637 0.000535 13.8   3.059  0.0220
# TETI2 - TETIxCHIC     -0.000133 0.000427 16.1  -0.312  0.9480

# get within cross averages
tetixchic_data_summed <- as.data.frame(
  tetixchic_data %>%
    group_by(brood_ID,population) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T)))
mean(subset(tetixchic_data_summed,population=="CHICxCHIC")$embryo_dry_weight_g) # 0.0060
mean(subset(tetixchic_data_summed,population=="TETIxCHIC")$embryo_dry_weight_g) # 0.0043
mean(subset(tetixchic_data_summed,population=="TETI2")$embryo_dry_weight_g) # 0.0042


## Supp Table S13 - RNAseq metadata


## Supp Table S14 - Embryo differential gene expression results
## See embryo-DGE_DESeq2_xmac-IDs_2023.R to run DGE analysis
embryo_dge <- read.csv("Data/embryo_ovary_dge_combined_2023/embryo_xmac-gtf_dge/embryo-xmacID-combined2023_dge_lfc-shr_all.csv_with-annots.csv",header=T)


## Supp Table S15 - Ovary differential gene expression results
## See ovary-DGE_DESeq2_xmac-IDs_2023.R to run DGE analysis
ovary_dge <- read.csv("Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/ovary-xmacID-combined2023_dge_lfc-shr_all.csv_with-annots.csv",header=T)


## Supp Table S16 - Gene Ontology biological pathways enriched in differentially expressed genes between X. malinche and X. birchmanni ovaries
## See ovary-GO-analysis_xmac-IDs_2023.R

## Supp Table S17 - Significance of relationship between WGCNA co-expressed gene cluster ovary expression profile and species/stage
## See ovary_WGCNA_xmac-IDs_combined2023.R to run WGCNA analysis
ovary_wgcna_module_pvals <- read.csv("Data/embryo_ovary_dge_combined_2023/wgcna/ovary-xmac_combined-FebAug23_MEtraitpvals.csv",header=T)


## Supp Table S18 - Enriched GO biological pathways in the WGCNA cluster 'darkorange2'
## See ovary-GO-analysis-WGCNA_xmac-IDs_2023.R


## Supp Table S19 - Critical thermal minimum (CTmin) of newborn fry
ctmin_data <- read.csv("Data/CTmin_Xmal_Xbirch_newborn_trial_data.csv",header=T)


## Supp Table S20 - Newborn fry food deprivation experiment data
fry_food_dep_data <- read.csv("Data/TableS21_food-deprivation-expt_data.csv",header=T)


## Table S21: Food deprivation experiment results
## A. Dry mass (g)
fry_fat_content <- read.csv("Data/X-23_fry_starvation_fat_content.csv")
# remove rows with fry content < 0; remove trial run Xbir and Xmal broods (1m & 1b); remove NAs and negative FC percentages
fry_fat_content <- subset(fry_fat_content, FC_percent >= 0 & !is.na(FC_percent) & brood_no != "1m" & brood_no != "1b")
fry_fat_content$condition <- paste(fry_fat_content$species,fry_fat_content$treatment, sep="_")

# select best model by LRT
# evaluate whether brood_ID should be included as a random effect
fit_full <- lme4::lmer(dry_mass_1 ~ condition+(1|brood_no), data = fry_fat_content, REML=F)
fit_reduced <- lm(dry_mass_1 ~ condition,data = fry_fat_content)
anova(fit_full,fit_reduced)
# variable                  Pr(>Chisq)
# (1|brood_ID)              1.4*10^-5 [kept]

# Stats: Compare raw means, accounting for covariates
dry_mass_fit <- lme4::lmer(dry_mass_1 ~ condition+(1|brood_no), data = fry_fat_content, REML=T)
emmeans(dry_mass_fit, list(pairwise ~ condition), adjust = "tukey")
# 1                            estimate       SE    df t.ratio p.value
# Xbir_control - Xbir_starved  0.001043 0.000190 51.16   5.500  <.0001
# Xbir_control - Xmal_control -0.000256 0.000394  6.27  -0.650  0.9121
# Xbir_control - Xmal_starved  0.000552 0.000392  6.17   1.409  0.5370
# Xbir_starved - Xmal_control -0.001300 0.000391  6.03  -3.328  0.0581
# Xbir_starved - Xmal_starved -0.000492 0.000388  5.93  -1.267  0.6131
# Xmal_control - Xmal_starved  0.000808 0.000160 51.35   5.047  <.0001

## old:
# ANOVA - posthoc Tukey
# no difference between Xmal v Xbir in either condition
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

## B. Standard length (cm)
fry_fat_content$standard_length_mm <- fry_fat_content$standard_length_cm*10 # convert post-trial std length cm to mm
std_length_fit <- lme4::lmer(standard_length_mm ~ condition+(1|brood_no), data = fry_fat_content, REML=T)

# Stats: Compare raw means, accounting for covariates
# ANOVA - posthoc Tukey
emmeans(std_length_fit, list(pairwise ~ condition), adjust = "tukey")
# $`pairwise differences of condition`
# 1                           estimate     SE df t.ratio p.value
# Xbir_control - Xbir_starved    1.130 0.325 51.78   3.471  0.0057
# Xbir_control - Xmal_control   -0.393 0.400  9.98  -0.980  0.7634
# Xbir_control - Xmal_starved    0.280 0.396  9.86   0.707  0.8921
# Xbir_starved - Xmal_control   -1.522 0.390  8.52  -3.901  0.0171
# Xbir_starved - Xmal_starved   -0.850 0.385  8.39  -2.206  0.1978
# Xmal_control - Xmal_starved    0.672 0.274 52.19   2.451  0.0800

## additional information
# add average pre-trial standard length per brood to dataframe
fry_initial_std_lengths <- read.csv("Data/X-23_fry_starvation_initial_standard_lengths.csv")
fry_initial_std_lengths_avg <- as.data.frame(
  fry_initial_std_lengths %>%
    group_by(brood_no) %>%
    reframe(avg_initial_std_length_mm = mean(standard_length_cm,na.rm=T)*10))
fry_fat_content <- merge(fry_fat_content,fry_initial_std_lengths_avg,by="brood_no")

# normalize standard lengths by avg length pre-trial
fry_fat_content$std_length_normalized <- fry_fat_content$standard_length_mm / fry_fat_content$avg_initial_std_length_mm
# calculate growth rate: (post trial standard length - pre-trial average lengths) / 3 days
fry_fat_content$std_length_growth_rate <- (fry_fat_content$standard_length_mm - fry_fat_content$avg_initial_std_length_mm) / 3
std_length_growth_rates_avg <- as.data.frame(
  fry_fat_content %>%
    group_by(condition) %>%
    reframe(std_length_growth_rate = mean(std_length_growth_rate,na.rm=T)))
std_length_growth_rates_avg
#    condition std_length_growth_rate
# Xbir_control             0.31943182
# Xbir_starved            -0.10678571
# Xmal_control             0.08593491
# Xmal_starved            -0.19475373

## C. Fat content (total %)
fry_fat_fit <- lme4::lmer(FC_percent ~ condition+(1|brood_no), data = fry_fat_content, REML=T)

# stats: ANOVA and Tukey
emmeans(fry_fat_fit, list(pairwise ~ condition), adjust = "tukey")
# Stats: Compare raw means, accounting for covariates
# ANOVA - posthoc Tukey
# $`pairwise differences of condition`
# 1                           estimate   SE    df t.ratio p.value
# Xbir_control - Xbir_starved    6.492 1.39 51.06   4.657  0.0001
# Xbir_control - Xmal_control   -0.898 4.51  5.48  -0.199  0.9969
# Xbir_control - Xmal_starved    8.006 4.50  5.44   1.779  0.3732
# Xbir_starved - Xmal_control   -7.390 4.50  5.40  -1.644  0.4303
# Xbir_starved - Xmal_starved    1.514 4.48  5.35   0.338  0.9853
# Xmal_control - Xmal_starved    8.903 1.18 51.14   7.568  <.0001

# write out complete table (Supp Table 20)
write.csv(fry_fat_content,"Data/TableS21_food-deprivation-expt_data.csv")


## Supp Table S22 - Pregnancy rates across seasons
pregnancy_data <- read.csv("Data/pregnancy_rate_CHIC_COAC_collections.csv",header=T)


## Supp Table S23 - Mortality data on X. birchmanni mother x X. malinche father F1 crosses


## REMOVED - Partial residuals of standard length used to calculate broad-sense heritability
# prepare partial residuals to account for covariates in trait values
nb_fry_data <- read.csv("Data/newborn_fry_size_data.csv")
nb_fry_data$born_year <- as.character(nb_fry_data$born_year)
subset_nb_fry_data <- subset(nb_fry_data, species_site=="XbirCOAC" | site=="CHIC" | species=="Xmal_xbir_F2" | species=="Xmal_xbir_F1")
# add season as variable
subset_nb_fry_data$season <- "warm" # IV-X
subset_nb_fry_data$season[subset_nb_fry_data$born_month=="XI"] <- "cold" # XI-III
subset_nb_fry_data$season[subset_nb_fry_data$born_month=="XII"] <- "cold" # XI-III
subset_nb_fry_data$season[subset_nb_fry_data$born_month=="I"] <- "cold" # XI-III
subset_nb_fry_data$season[subset_nb_fry_data$born_month=="II"] <- "cold" # XI-III
subset_nb_fry_data$season[subset_nb_fry_data$born_month=="III"] <- "cold" # XI-III

# choose model
## XXX
## # model selection by LRT with mother ID as random effect, mother length, season/collection date, etc
## fit0 <- lme4::glmer(sl ~ species + (1 | motherID), family = binomial, data = data)
## fit1 <- lme4::glmer(sl ~ species, family = binomial, data = data)
## lrtest(fit0,fit1)

model_all<-lm(sl_mm~species+season,data=subset_nb_fry_data)
selectedMod <- step(model_all) # sl_mm~species+season
summary(model_all)
subset_nb_size_fit <- lmer(sl_mm ~ species+season+(1|brood_ID),subset_nb_fry_data)

# calculate partial residuals
nb_size_residuals <- visreg(subset_nb_size_fit, "species",plot=F)$res
write.csv(nb_size_residuals,"Data/newborn_standard-length_Xbir-Xmal-F1-F2_partial-residuals_april2024.csv")



#### Additional Supp Materials calculations and statistics

## Supp Materials S1: Life history traits
combined_embryo_data <- read.table("Data/all-pops_combined_embryo_datasets.csv",header=T,sep=',')
combined_subset <- subset(combined_embryo_data,(popcoll=="CHIC_V-2022" | popcoll=="CHIC_VIII_20" | popcoll=="COAC_VIII_20" | popcoll=="COAC_II-2023" | popcoll=="COAC_IX-2022" | popcoll=="COAC_X-2023") & stage!="0" & stage!="2" & stage!="5" & stage!="?" & stage!="-")
# average embryo size within brood
combined_subset_summed <- as.data.frame(
  combined_subset %>%
    group_by(brood_ID,species,population,season,collection) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T),brood_size=mean(brood_size,na.rm=T),mother_std_length=first(mother_std_length)))

# effect of mother standard length on average embryo weight within brood
step(lm(embryo_dry_weight_g_mean~species+season+brood_size+mother_std_length,data=combined_subset_summed)) # embryo_dry_weight_g_mean ~ species + season + mother_std_length
anova(lm(embryo_dry_weight_g_mean ~ species + season + mother_std_length,combined_subset_summed))
#                    Df     Sum Sq    Mean Sq F value    Pr(>F)
# species            1 1.4169e-05 1.4169e-05 43.9352 2.269e-09 ***
# season             1 4.4797e-06 4.4797e-06 13.8908 0.0003343 ***
# mother_std_length  1 2.1559e-06 2.1559e-06  6.6852 0.0112925 *

# check whether there is a significant interaction between species and mother standard length
# there is not (p=0.65), and species:mother_std_length is dropped from the model
step(lm(embryo_dry_weight_g_mean~species+season+brood_size+mother_std_length+species:mother_std_length,data=combined_subset_summed)) # embryo_dry_weight_g_mean ~ species + season + mother_std_length
anova(lm(embryo_dry_weight_g_mean ~ species + season + mother_std_length + species:mother_std_length,combined_subset_summed))

# relationship between brood size and mother standard length
step(lm(brood_size~mother_std_length+species+season,data=combined_subset_summed)) # brood_size ~ mother_std_length + season
anova(lm(brood_size ~ mother_std_length + species + season,combined_subset_summed))
#                    Df Sum Sq Mean Sq  F value  Pr(>F)
# mother_std_length  1 3347.6  3347.6 95.3841 6.381e-16 ***
# species            1  163.2   163.2  4.6510   0.03362 *
# season             1   90.2    90.2  2.5698   0.11231

# mean brood sizes by species
mean(subset(combined_subset_summed,species=="Xbirchmanni")$brood_size) # 13.15
mean(subset(combined_subset_summed,species=="Xmalinche")$brood_size) # 17.18

# range of mother standard lengths, by population
mother_lengths <- as.data.frame(
  combined_subset %>%
    group_by(population) %>%
    reframe(min_mother_length = min(mother_std_length,na.rm=T),max_mother_length = max(mother_std_length,na.rm=T),min_brood_size=min(brood_size,na.rm=T),max_brood_size=max(brood_size,na.rm=T)))
mother_lengths

# number of pregnant mothers per population and per collection
dim(subset(combined_subset_summed,species=="Xbirchmanni" & !is.na(embryo_dry_weight_g_mean))) # 59
table(subset(combined_subset_summed,species=="Xbirchmanni" & !is.na(embryo_dry_weight_g_mean))$collection)
dim(subset(combined_subset_summed,species=="Xmalinche")) # 38
table(subset(combined_subset_summed,species=="Xmalinche" & !is.na(embryo_dry_weight_g_mean))$collection)

# number of pregnant + nonpregnant mothers per population and per collection
combined_embryo_data <- read.table("Data/all-pops_combined_embryo_datasets.csv",header=T,sep=',')
combined_subset <- subset(combined_embryo_data,(popcoll=="CHIC_V-2022" | popcoll=="CHIC_VIII_20" | popcoll=="COAC_VIII_20" | popcoll=="COAC_II-2023" | popcoll=="COAC_IX-2022" | popcoll=="COAC_X-2023" | popcoll=="PTHC_II-2023"))
combined_subset_summed <- as.data.frame(
  combined_subset %>%
    group_by(brood_ID,species,population,season,collection) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T),brood_size=mean(brood_size,na.rm=T),mother_std_length=first(mother_std_length)))
dim(subset(combined_subset_summed,species=="Xbirchmanni")) # 110
table(subset(combined_subset_summed,species=="Xbirchmanni")$collection)
dim(subset(combined_subset_summed,species=="Xmalinche")) # 49
table(subset(combined_subset_summed,species=="Xmalinche")$collection)
dim(subset(combined_subset_summed,species=="Xcortezi")) # 19 (only pregnant)



## REMOVED: Heritability of standard length calculation
# prepare partial residuals to account for covariates in trait values
nb_fry_data <- read.csv("Data/newborn_fry_size_data.csv")
nb_fry_data$born_year <- as.character(nb_fry_data$born_year)
subset_nb_fry_data <- subset(nb_fry_data, species_site=="XbirCOAC" | site=="CHIC" | species=="Xmal_xbir_F2" | species=="Xmal_xbir_F1")
# add season as variable
subset_nb_fry_data$season <- "warm" # IV-X
subset_nb_fry_data$season[subset_nb_fry_data$born_month=="XI"] <- "cold" # XI-III
subset_nb_fry_data$season[subset_nb_fry_data$born_month=="XII"] <- "cold" # XI-III
subset_nb_fry_data$season[subset_nb_fry_data$born_month=="I"] <- "cold" # XI-III
subset_nb_fry_data$season[subset_nb_fry_data$born_month=="II"] <- "cold" # XI-III
subset_nb_fry_data$season[subset_nb_fry_data$born_month=="III"] <- "cold" # XI-III

# choose model
model_all<-lm(sl_mm~species+season,data=subset_nb_fry_data)
selectedMod <- step(model_all) # sl_mm~species+season
summary(model_all)
subset_nb_size_fit <- lmer(sl_mm ~ species+season+(1|brood_ID),subset_nb_fry_data)

# calculate partial residuals
nb_size_residuals <- visreg(subset_nb_size_fit, "species",plot=F)$res
write.csv(nb_size_residuals,"Data/newborn_standard-length_Xbir-Xmal-F1-F2_partial-residuals_april2024.csv")

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
between_var_sl_par <- var(c(coac_sl_mean,chic_sl_mean)) # 0.3604126

## calculate within-population variance in size of parent, F1, and F2 fry
var_sl_coac <- var(coac_fry$visregRes) # 0.2035312
var_sl_chic <- var(chic_fry$visregRes) # 0.2255541
var_sl_f1 <- var(f1_fry$visregRes)     # 0.6902981
var_sl_f2 <- var(f2_fry$visregRes)     # 0.6065563
## within-F2 variance = 0.6065563

## calculate broad-sense heritability H^2 with F1 fry size data
# calculate environmental variance with Wright's weighted average of within-strain variance
var_sl_env <- var_sl_coac/4 + var_sl_chic/4 + var_sl_f1/2
## environmental variance = 0.4524204

# Calculate H^2 = (var_F2 - var_enviroment)/var_F2
h2 <- (var_sl_f2 - var_sl_env) / var_sl_f2
h2
## h^2 = 0.2541165


## Supp Materials S2: Matrotrophy Index calculations and bootstraps
## calculated MI based on average embryo size by stage within brood
combined_embryo_data <- read.table("Data/all-pops_combined_embryo_datasets.csv",header=T,sep=',')
combined_subset <- subset(combined_embryo_data,(popcoll=="CHIC_V-2022" | popcoll=="CHIC_VIII_20" | popcoll=="COAC_VIII_20" | popcoll=="COAC_II-2023" | popcoll=="COAC_IX-2022" | popcoll=="COAC_X-2023" | popcoll=="PTHC_II-2023") & stage!="0" & stage!="2" & stage!="5" & stage!="?" & stage!="-")

## Xbirchmanni MI
coac_subset <- subset(combined_subset, population=="COAC")

# summarize all COAC mothers
coac_uniq_IDs <- na.omit(as.data.frame(coac_subset %>%group_by(brood_ID) %>% reframe(collection = first(collection))))
table(coac_uniq_IDs$collection)

# average embryo size by stage within brood
coac_stage_size_avged <- as.data.frame(
  coac_subset %>%
    group_by(brood_ID,stage,season) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T),mother_std_length=mean(mother_std_length,na.rm=T)))
coac_stage_size_avged <- na.omit(coac_stage_size_avged)
coac_stage_size_avged

# Xbir MI with raw values, with embryo dry weight averaged within broods
coac_stage10 <- subset(coac_stage_size_avged,stage=="10") # 4 broods
coac_stage10_mean <- mean(coac_stage10$embryo_dry_weight_g_mean) # 0.0037875
coac_stage50 <- subset(coac_stage_size_avged, stage=="50") # 3 broods
coac_stage50_mean <- mean(coac_stage50$embryo_dry_weight_g_mean) # 0.002487908
coac_stage50_mean/coac_stage10_mean # 0.6568735
# output the raw averages used for this calculation
coac_raw_table <- subset(coac_stage_size_avged,stage=="10" | stage=="50")
write.csv(coac_raw_table,"Data/MI-calculation_COAC_raw-avg-values_stage-10-50.csv")

# Xbir MI with partial residuals, with embryo dry weight averaged within broods *
# calculate partial residuals
coac_fit <- lm(embryo_dry_weight_g_mean~stage+season+mother_std_length,data=coac_stage_size_avged)
coac_avg_pr <- visreg(coac_fit, "stage",plot=F)$res
coac_stage10 <- subset(coac_avg_pr,stage=="10")
coac_stage10_mean <- mean(coac_stage10$visregRes)
coac_stage50 <- subset(coac_avg_pr, stage=="50")
coac_stage50_mean <- mean(coac_stage50$visregRes)
coac_stage50_mean/coac_stage10_mean # 0.6588837
# output table
coac_pr_table <- subset(coac_avg_pr,stage=="10" | stage=="50")
coac_pr_table$species <- "Xbirchmanni"
write.csv(coac_pr_table,"Data/MI-calculation_COAC_partial-residuals-avg_stage-10-50.csv")

## Bootstrap MI for COAC Xbir population, by post-2020 collections
# average embryos at the same stage in the same brood, and filter out unused stages (0-5)
combined_subset <- subset(combined_embryo_data,(popcoll=="CHIC_V-2022" | popcoll=="CHIC_VIII_20" | popcoll=="COAC_VIII_20" | popcoll=="COAC_II-2023" | popcoll=="COAC_IX-2022" | popcoll=="COAC_X-2023") & stage!="0" & stage!="2" & stage!="5" & stage!="?" & stage!="-")
combined_subset_avg <- as.data.frame(
  combined_subset %>%
    group_by(brood_ID,species,population,stage) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T), popcoll=first(popcoll)))
combined_subset_avg

coac_subset_avg <- subset(combined_subset_avg,population=="COAC")
table(coac_subset_avg$stage)
coac_avg_stage10 <- subset(coac_subset_avg,stage==10)
coac_avg_stage50 <- subset(coac_subset_avg,stage==50)
coac_mi_bootstrap <- c()
for(i in 1:1000) {
  resample_stage10 <- sample(coac_avg_stage10$embryo_dry_weight_g_mean, replace = TRUE)
  mean_stage10 <- mean(resample_stage10)
  resample_stage50 <- sample(coac_avg_stage50$embryo_dry_weight_g_mean, replace = TRUE)
  mean_stage50 <- mean(resample_stage50)
  coac_mi <- mean_stage50/mean_stage10
  coac_mi_bootstrap <- c(coac_mi_bootstrap,coac_mi)
}
mean(coac_mi_bootstrap) # 0.658127
t.test(coac_mi_bootstrap,conf.level = 0.95) # 0.6565315 0.6597226


## Xmalinche MI
chic_subset <- subset(combined_subset, population=="CHIC")

# summarize all CHIC mothers
chic_uniq_IDs <- na.omit(as.data.frame(chic_subset %>%group_by(brood_ID) %>% reframe(collection = first(collection))))
table(chic_uniq_IDs$collection)

# average embryo size by stage within brood
chic_stage_size_avged <- as.data.frame(
  chic_subset %>%
    group_by(brood_ID,stage,season) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T),mother_std_length=mean(mother_std_length,na.rm=T)))
chic_stage_size_avged <- na.omit(chic_stage_size_avged)
chic_stage_size_avged
# note: only one season for X. malinche collections (warm)

# Xmal MI with raw values, with embryo dry weight averaged within broods
chic_stage10 <- subset(chic_stage_size_avged,stage=="10") # 8 broods
chic_stage10_mean <- mean(chic_stage10$embryo_dry_weight_g_mean) # 0.004010417
chic_stage50 <- subset(chic_stage_size_avged,stage=="50") # 7 broods
chic_stage50_mean <- mean(chic_stage50$embryo_dry_weight_g_mean) # 0.003810819
chic_stage50_mean/chic_stage10_mean # 0.9502302
# output table
chic_raw_table <- subset(chic_stage_size_avged,stage=="10" | stage=="50")
write.csv(chic_raw_table,"Data/MI-calculation_CHIC_raw-avg-values_stage-10-50.csv")

# Xmal MI with partial residuals, with embryo dry weight averaged within broods *
# calculate partial residuals
chic_fit <- lm(embryo_dry_weight_g_mean~stage+mother_std_length,data=chic_stage_size_avged)
chic_avg_pr <- visreg(chic_fit, "stage",plot=F)$res
# calculate MI
chic_stage10 <- subset(chic_avg_pr,stage=="10")
chic_stage10_mean <- mean(chic_stage10$visregRes)
chic_stage50 <- subset(chic_avg_pr,stage=="50")
chic_stage50_mean <- mean(chic_stage50$visregRes)
chic_stage50_mean/chic_stage10_mean # 0.9836957
# output table
chic_pr_table <- subset(chic_stage_partial_residuals,stage=="10" | stage=="50")
chic_pr_table$species <- "Xmalinche"
write.csv(chic_pr_table,"Data/MI-calculation_CHIC_partial-residuals-avg_stage-10-50.csv")

## Bootstrap MI for CHIC Xmal population
# average embryos at the same stage in the same brood, and filter out unused stages (0-5)
combined_subset <- subset(combined_embryo_data,(popcoll=="CHIC_V-2022" | popcoll=="CHIC_VIII_20" | popcoll=="COAC_VIII_20" | popcoll=="COAC_II-2023" | popcoll=="COAC_IX-2022" | popcoll=="COAC_X-2023") & stage!="0" & stage!="2" & stage!="5" & stage!="?" & stage!="-")
combined_subset_avg <- as.data.frame(
  combined_subset %>%
    group_by(brood_ID,species,population,stage) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T), popcoll=first(popcoll)))
combined_subset_avg

chic_subset_avg <- subset(combined_subset_avg,population=="CHIC")
table(chic_subset_avg$stage)
chic_avg_stage10 <- subset(chic_subset_avg,stage==10)
chic_avg_stage50 <- subset(chic_subset_avg,stage==50)
chic_mi_bootstrap <- c()
for(i in 1:1000) {
  resample_stage10 <- sample(chic_avg_stage10$embryo_dry_weight_g_mean, replace = TRUE)
  mean_stage10 <- mean(resample_stage10)
  resample_stage50 <- sample(chic_avg_stage50$embryo_dry_weight_g_mean, replace = TRUE)
  mean_stage50 <- mean(resample_stage50)
  chic_mi <- mean_stage50/mean_stage10
  chic_mi_bootstrap <- c(chic_mi_bootstrap,chic_mi)
}
mean(chic_mi_bootstrap) # 0.9599878
t.test(chic_mi_bootstrap,conf.level = 0.95) #  0.9547505 0.9652250


## Xcortezi MI
pthc_subset <- subset(combined_subset, population=="PTHC")

# Xcor MI with partial residuals
pthc_fit <- lm(embryo_dry_weight_g_mean~stage+mother_std_length,data=pthc_stage_size_avged)
pthc_avg_pr <- visreg(pthc_fit, "stage",plot=F)$res
# calculate MI
pthc_stage10 <- subset(pthc_avg_pr,stage=="10") # 3 broods
pthc_stage10_mean <- mean(pthc_stage10$visregRes)
pthc_stage50 <- subset(pthc_avg_pr, stage=="50") # 1 brood
pthc_stage50_mean <- mean(pthc_stage50$visregRes)
pthc_stage50_mean/pthc_stage10_mean # 0.6995778
pthc_stage50_45 <- subset(pthc_avg_pr,stage=="45" | stage=="50") # 5 broods
pthc_stage50_45_mean <- mean(pthc_stage50_45$visregRes)
pthc_stage50_45_mean/pthc_stage10_mean # 0.9404641
# output table
pthc_pr_table <- subset(pthc_stage_partial_residuals,stage=="10" | stage=="45" | stage=="50")
write.csv(pthc_pr_table,"Data/MI-calculation_PTHC_partial-residuals-avg_stage-10-45-50.csv")

## calculated MI based on average embryo size by stage within brood

# average embryo size by stage within brood
pthc_stage_size_avged <- as.data.frame(
  pthc_subset %>%
    group_by(brood_ID,stage,season) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T),mother_std_length=mean(mother_std_length,na.rm=T)))
pthc_stage_size_avged

# Xcor MI with raw values, with embryo dry weight averaged within broods
pthc_stage10 <- subset(pthc_stage_size_avged,stage=="10") # 3 broods
pthc_stage10_mean <- mean(pthc_stage10$embryo_dry_weight_g_mean)
pthc_stage50 <- subset(pthc_stage_size_avged, stage=="50") # 1 brood
pthc_stage50_mean <- mean(pthc_stage50$embryo_dry_weight_g_mean)
pthc_stage50_mean/pthc_stage10_mean # 0.7242604
pthc_stage50_45 <- subset(pthc_stage_size_avged,stage=="45" | stage=="50") # 5 broods
pthc_stage50_45_mean <- mean(pthc_stage50_45$embryo_dry_weight_g_mean)
pthc_stage50_45_mean/pthc_stage10_mean # 0.9484674
# output table
pthc_raw_table <- subset(pthc_stage_size_avged,stage=="10" | stage=="45" | stage=="50")
write.csv(pthc_raw_table,"Data/MI-calculation_PTHC_raw-avg-values_stage-10-45-50.csv")


## Calculate MI with partial residuals

## Calculate MI for Xmal with partial residuals instead of averages, with brood ID as random effect
# calculate partial residuals, using the model selected for Fig 2A
chic_fit <- lmer(embryo_dry_weight_g~stage+mother_std_length+(1|brood_ID),data=chic_subset, REML=T)
chic_pr <- visreg(chic_fit, "stage", plot=F)$res
# calculate MI
chic_stage10 <- subset(chic_pr,stage=="10")
chic_stage10_mean <- mean(chic_stage10$visregRes)
chic_stage50 <- subset(chic_pr,stage=="50")
chic_stage50_mean <- mean(chic_stage50$visregRes)
chic_stage50_mean/chic_stage10_mean # 0.9439981

# Bootstrap MI for CHIC Xmal population

table(chic_pr$stage)
chic_stage10 <- subset(chic_pr,stage==10)
chic_stage50 <- subset(chic_pr,stage==50)
chic_mi_bootstrap <- c()
for(i in 1:1000) {
  resample_stage10 <- sample(chic_stage10$visregRes, replace = TRUE)
  mean_stage10 <- mean(resample_stage10)
  resample_stage50 <- sample(chic_stage50$visregRes, replace = TRUE)
  mean_stage50 <- mean(resample_stage50)
  chic_mi <- mean_stage50/mean_stage10
  chic_mi_bootstrap <- c(chic_mi_bootstrap,chic_mi)
}
mean(chic_mi_bootstrap) # 0.9451785
t.test(chic_mi_bootstrap,conf.level = 0.95) #  0.9437161 0.9466408


## Calculate MI for Xbir with partial residuals instead of averages, with brood ID as random effect
# calculate partial residuals, using the model selected for Fig 2A
coac_fit <- lme4::lmer(embryo_dry_weight_g~stage+season+mother_std_length+(1|brood_ID),data=coac_subset, REML=T)
coac_pr <- visreg(coac_fit, "stage", plot=F)$res
# calculate MI
coac_stage10 <- subset(coac_pr,stage=="10")
coac_stage10_mean <- mean(coac_stage10$visregRes) # 0.003508688
coac_stage50 <- subset(coac_pr, stage=="50")
coac_stage50_mean <- mean(coac_stage50$visregRes) # 0.002540325
coac_stage50_mean/coac_stage10_mean # 0.7240099

# Bootstrap MI for COAC Xbir population, by post-2020 collections
table(coac_pr$stage)
coac_stage10 <- subset(coac_pr,stage==10)
coac_stage50 <- subset(coac_pr,stage==50)
coac_mi_bootstrap <- c()
for(i in 1:1000) {
  resample_stage10 <- sample(coac_stage10$visregRes, replace = TRUE)
  mean_stage10 <- mean(resample_stage10)
  resample_stage50 <- sample(coac_stage50$visregRes, replace = TRUE)
  mean_stage50 <- mean(resample_stage50)
  coac_mi <- mean_stage50/mean_stage10
  coac_mi_bootstrap <- c(coac_mi_bootstrap,coac_mi)
}
mean(coac_mi_bootstrap) # 0.7251131
t.test(coac_mi_bootstrap,conf.level = 0.95) # 0.7229252 0.7273011



## Supp Materials S3: Relationship between ovary dry weight and mother hybrid index / mitotype
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

## Relationship between ovary dry weight and mother hybrid index
# average embryo dry weight by stage within brood
CALL_ovary_weight_unique <- as.data.frame(
  CALL_embryo_data_nomal %>%
    group_by(brood_ID,stage) %>%
    reframe(ovarian_tissue_dry_weight_g = first(ovarian_tissue_dry_weight_g),brood_size=first(brood_size),mother_hi=first(mother_hi),mitotype=first(mitotype)))
# model fit with LRT: mother hybrid index
CALL_hi_ovary_lm_full <- lm(ovarian_tissue_dry_weight_g ~ mother_hi*stage + brood_size, CALL_ovary_weight_unique)
CALL_hi_ovary_lm_reduced <- lm(ovarian_tissue_dry_weight_g ~ stage + mother_hi + brood_size, CALL_ovary_weight_unique)
anova(CALL_hi_ovary_lm_full,CALL_hi_ovary_lm_reduced)
# brood_size                0.001242 [kept]
# stage                     0.004342 [kept]
# mother_hi+mother_hi:stage 0.01971 [kept]
# mother_hi:stage           0.01221 [kept]
# ovarian_tissue_dry_weight_g ~ mother_hi + stage + mother_hi:stage + brood_size

# evaluate significance with likelihood ratio test (anova) on full and reduced model
CALL_hi_ovary_lm_fit <- lm(ovarian_tissue_dry_weight_g ~ mother_hi + stage + mother_hi:stage + brood_size, CALL_ovary_weight_unique)
summary(CALL_hi_ovary_lm_fit)
CALL_hi_ovary_lm_fit_reduced <- lm(ovarian_tissue_dry_weight_g ~ stage + brood_size, CALL_ovary_weight_unique)
anova(CALL_hi_ovary_lm_fit,CALL_hi_ovary_lm_fit_reduced)
#   Res.Df       RSS Df   Sum of Sq      F  Pr(>F)
# 2     75 0.0016711 -5 -0.00028619 2.8931 0.01971 *

# plot mother hi and ovary dry weight
CALL_hi_ovary_lm_vr <- visreg(CALL_hi_ovary_lm_fit,"mother_hi")
plot(visregRes ~ mother_hi,CALL_hi_ovary_lm_vr$res,xlab="mother ancestry proportion",ylab="partial residuals of ovary dry weight (g)", col=hybrid_col,pch=20,cex.lab=1.8,cex.axis=1.5)
abline(lm(CALL_hi_ovary_lm_vr$res$visregRes~CALL_hi_ovary_lm_vr$res$mother_hi), col="purple",lwd=4)

## Relationship between ovary dry weight and mother mitotype
# model fit with LRT: mitotype
CALL_hi_ovary_lm_full <- lm(ovarian_tissue_dry_weight_g ~ mitotype + stage + mitotype:stage + brood_size, CALL_ovary_weight_unique)
CALL_hi_ovary_lm_reduced <- lm(ovarian_tissue_dry_weight_g ~ mitotype + stage + mitotype:stage, CALL_ovary_weight_unique)
anova(CALL_hi_ovary_lm_full,CALL_hi_ovary_lm_reduced)
# mitotype:stage        0.4323 [dropped]
# stage+mitotype:stage  0.01671 [kept]
# brood_size            0.001541 [kept]
# ovarian_tissue_dry_weight_g ~ mitotype + stage + brood_size

# evaluate significance of mitotype with likelihood ratio test
summary(CALL_mito_ovary_lm_fit) # mitotypemalinche          1.210e-03  2.620e-03   0.462  0.64561
CALL_mito_ovary_lm_fit <- lm(ovarian_tissue_dry_weight_g ~ mitotype + stage + brood_size, CALL_ovary_weight_unique)
CALL_mito_ovary_lm_fit_reduced <- lm(ovarian_tissue_dry_weight_g ~ stage + brood_size, CALL_ovary_weight_unique)
anova(CALL_mito_ovary_lm_fit,CALL_mito_ovary_lm_fit_reduced)
#   Res.Df       RSS Df   Sum of Sq      F Pr(>F)
# 2     75 0.0016711 -1 -6.2071e-07 0.0275 0.8687

## correlation between mother hybrid index and mitotype
temp <- CALL_embryo_data
temp$mito_num <- ifelse(temp$mitotype=="birchmanni",0,1)
temp <- subset(temp,!is.na(mother_hi) & !is.na(mito_num))
cor(temp$mother_hi,temp$mito_num,method="spearman") # spearman rho=0.61

