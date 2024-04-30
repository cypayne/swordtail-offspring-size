##### SUPPLEMENTARY TABLES (& MATERIALS)
##### Recent evolution of large offspring size and post-fertilization nutrient provisioning in swordtails
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
combined_subset <- subset(combined_embryo_data,(popcoll=="CHIC_V-2022" | popcoll=="CHIC_VIII_20" | popcoll=="COAC_VIII_20" | popcoll=="COAC_II-2023" | popcoll=="COAC_IX-2022" | popcoll=="COAC_X-2023") & stage!="0" & stage!="2" & stage!="5" & stage!="?" & stage!="-")
# chosen fit
combined_fit <- lmer(embryo_dry_weight_g~species+stage+species:stage+season+mother_std_length+(1|brood_ID),data=combined_subset, REML=T)
# ANOVA and Tukey-test to compare species by stage
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
roof_chic_coac_f1_subset <- subset(roof_chic_coac_f1,stage!="0" & stage!="5" & stage!="10" & stage!="15" & stage!="20" & stage!="25")
# remove one outlier - likely erroneous data recording
roof_chic_coac_f1_subset <- subset(roof_chic_coac_f1_subset,embryo_dry_weight_g < 0.008)

# choose model
step(lm(embryo_dry_weight_g~species*stage+brood_size+mother_std_length,roof_chic_coac_f1_subset)) # embryo_dry_weight_g ~ species + stage + brood_size + mother_std_length
lab_cross_fit <- lmer(embryo_dry_weight_g~species+stage+brood_size+mother_std_length+(1|brood_ID),data=roof_chic_coac_f1_subset, REML=T)
# the above model is singular because brood_ID and maternal standard length are perfectly correlated (R^2=1)
summary(lm(mother_std_length~brood_ID,roof_chic_coac_f1_subset))
# importantly, the visreg partial residuals are exactly the same whether we use a mixed linear model with brood_ID or
# a simple linear model without it
# therefore, we drop brood_ID as a random effect and use a simple linear model instead
lab_cross_fit <- lm(embryo_dry_weight_g~species+stage+brood_size+mother_std_length,data=roof_chic_coac_f1_subset)

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

# choose model
tetixchic_full <- lm(embryo_dry_weight_g~population+stage+brood_size+mother_std_length+mother_origin, data=tetixchic_data)
step(tetixchic_full) # embryo_dry_weight_g ~ population + stage + brood_size + mother_std_length

# stats: ANOVA and Tukey
tetixchic_fit <- lmer(embryo_dry_weight_g ~ population + stage + brood_size + mother_std_length + (1|brood_ID),data=tetixchic_data, REML=T)
summary(tetixchic_fit)
emmeans(tetixchic_fit, list(pairwise ~ population | stage), adjust = "tukey")
# 1                      estimate       SE   df t.ratio p.value
# CHICxCHIC - TETI2      0.001680 0.000647 10.5   2.599  0.0612
# CHICxCHIC - TETIxCHIC  0.001103 0.000799 11.4   1.380  0.3827
# TETI2 - TETIxCHIC     -0.000577 0.000595 12.1  -0.970  0.6085

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

# choose model
step(lm(dry_mass_1 ~ condition+species+treatment+brood_no,fry_fat_content)) # dry_mass_1 ~ condition + brood_no

# Stats: Compare raw means, accounting for covariates
dry_mass_fit <- lme(dry_mass_1 ~ condition, random=~1|brood_no,data=fry_fat_content)
emmeans(dry_mass_fit, list(pairwise ~ condition), adjust = "tukey")
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

## B. Standard length (cm)
## chose to plot partial residuals of standard lengths, instead of lengths normalized
## by average initial lengths, because it better shows that Xmal and Xbir fry achieve
## roughly the same size after 3 days, despite starting much smaller, and shows that
## xbir gain less in no food conditions
# add average pre-trial standard length per brood to dataframe
fry_initial_std_lengths <- read.csv("Data/X-23_fry_starvation_initial_standard_lengths.csv")
fry_initial_std_lengths_avg <- as.data.frame(
  fry_initial_std_lengths %>%
    group_by(brood_no) %>%
    reframe(avg_initial_std_length_mm = mean(standard_length_cm,na.rm=T)*10))
fry_fat_content <- merge(fry_fat_content,fry_initial_std_lengths_avg,by="brood_no")
# convert post-trial std length cm to mm
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

# choose model
step(lm(standard_length_mm ~ condition+species+avg_initial_std_length_mm+brood_no,fry_fat_content)) #standard_length_mm ~ condition + brood_no

# Stats: Compare raw means, accounting for covariates
# ANOVA - posthoc Tukey
std_length_fit <- lme(standard_length_mm ~ condition, random=~1|brood_no,data=fry_fat_content)
emmeans(std_length_fit, list(pairwise ~ condition), adjust = "tukey")
#$`pairwise differences of condition`
#1                           estimate     SE df t.ratio p.value
#Xbir_control - Xbir_starved   0.1130 0.0325 50   3.481  0.0056
#Xbir_control - Xmal_control  -0.0393 0.0397 50  -0.989  0.7562
#Xbir_control - Xmal_starved   0.0280 0.0393 50   0.712  0.8920
#Xbir_starved - Xmal_control  -0.1522 0.0384 50  -3.963  0.0013
#Xbir_starved - Xmal_starved  -0.0850 0.0380 50  -2.237  0.1274
#Xmal_control - Xmal_starved   0.0672 0.0273 50   2.462  0.0786

## C. Fat content (total %)
# choose model
step(lm(FC_percent ~ condition+species+treatment+brood_no,data=fry_fat_content)) # FC_percent ~ condition + brood_no

# stats: ANOVA and Tukey
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

# write out complete table (Supp Table 20)
write.csv(fry_fat_content,"Data/TableS21_food-deprivation-expt_data.csv")


## Supp Table S22 - Pregnancy rates across seasons
pregnancy_data <- read.csv("Data/pregnancy_rate_CHIC_COAC_collections.csv",header=T)


## Supp Table S23 - Mortality data on X. birchmanni mother x X. malinche father F1 crosses


## Supp Table S24 - Partial residuals of standard length used to calculate broad-sense heritability
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



## Supp Materials S2: Heritability of standard length calculation
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


## Supp Materials S3: Matrotrophy Index calculations and bootstraps
## calculate matrotrophy index
combined_embryo_data <- read.table("Data/all-pops_combined_embryo_datasets.csv",header=T,sep=',')

## Xbirchmanni MI
coac_subset <- subset(combined_embryo_data,stage!="0" & population=="COAC" & (popcoll=="COAC_II-2023" | popcoll=="COAC_IX-2022" | popcoll=="COAC_VIII_20"))
# average embryo size by stage within brood
coac_stage_size_avged <- as.data.frame(
  coac_subset %>%
    group_by(brood_ID,stage,collection) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T),mother_std_length=mean(mother_std_length,na.rm=T)))
coac_stage_size_avged <- na.omit(coac_stage_size_avged)

# summarize all COAC mothers
coac_uniq_IDs <- na.omit(as.data.frame(subset(combined_embryo_data,(popcoll=="COAC_II-2023" | popcoll=="COAC_IX-2022" | popcoll=="COAC_VIII_20")) %>%group_by(brood_ID) %>% reframe(collection = first(collection))))
table(coac_uniq_IDs$collection)

# Xbir MI with raw values
coac_fit <- lm(embryo_dry_weight_g_mean~collection+stage+mother_std_length,data=coac_stage_size_avged)
step(coac_fit) # embryo_dry_weight_g_mean~collection+stage+mother_std_length
coac_stage10 <- subset(coac_stage_size_avged,stage=="10") # 3
coac_stage10_mean <- mean(coac_stage10$embryo_dry_weight_g_mean)
coac_stage45 <- subset(coac_stage_size_avged,stage=="45") # 10
coac_stage45_mean <- mean(coac_stage45$embryo_dry_weight_g_mean)
coac_stage50 <- subset(coac_stage_size_avged, stage=="50") # only 1
coac_stage50_mean <- mean(coac_stage50$embryo_dry_weight_g_mean)
coac_stage50_45 <- subset(coac_stage_size_avged,stage=="45" | stage=="50")
coac_stage50_45_mean <- mean(coac_stage50_45$embryo_dry_weight_g_mean)
coac_stage45_mean/coac_stage10_mean # 0.7727921
coac_stage50_mean/coac_stage10_mean # 0.593617
coac_stage50_45_mean/coac_stage10_mean # 0.7548746
coac_raw_table <- subset(coac_stage_size_avged,stage=="10" | stage=="45" | stage=="50")
write.csv(coac_raw_table,"Data/MI-calculation_COAC_raw-values_stage-10-45-50.csv")

# Xbir MI with partial residuals
coac_fit <- lm(embryo_dry_weight_g_mean~collection+stage+mother_std_length,data=coac_stage_size_avged)
step(coac_fit) # embryo_dry_weight_g_mean~collection+stage+mother_std_length
coac_stage_partial_residuals <- visreg(coac_fit, "stage",plot=F)$res
coac_stage10 <- subset(coac_stage_partial_residuals,stage=="10")
coac_stage10_mean <- mean(coac_stage10$visregRes)
coac_stage45 <- subset(coac_stage_partial_residuals,stage=="45")
coac_stage45_mean <- mean(coac_stage45$visregRes)
coac_stage50 <- subset(coac_stage_partial_residuals, stage=="50")
coac_stage50_mean <- mean(coac_stage50$visregRes) # only 1
coac_stage50_45 <- subset(coac_stage_partial_residuals,stage=="45" | stage=="50")
coac_stage50_45_mean <- mean(coac_stage50_45$visregRes)
coac_stage45_mean/coac_stage10_mean # 0.689198
coac_stage50_mean/coac_stage10_mean # 0.6510765
coac_stage50_45_mean/coac_stage10_mean # 0.6853859
coac_pr_table <- subset(coac_stage_partial_residuals,stage=="10" | stage=="45" | stage=="50")
write.csv(coac_pr_table,"Data/MI-calculation_COAC_partial-residuals_stage-10-45-50.csv")

## Xmalinche MI
chic_subset <- subset(combined_embryo_data,stage!="0" & population=="CHIC" & (popcoll=="CHIC_V-2022" | popcoll=="CHIC_VIII_20"))
# average embryo size by stage within brood
chic_stage_size_avged <- as.data.frame(
  chic_subset %>%
    group_by(brood_ID,stage,collection) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T),mother_std_length=mean(mother_std_length,na.rm=T)))
chic_stage_size_avged <- na.omit(chic_stage_size_avged)
chic_stage_size_avged

# summarize all CHIC mothers
chic_uniq_IDs <- na.omit(as.data.frame(subset(combined_embryo_data,(popcoll=="CHIC_V-2022" | popcoll=="CHIC_VIII_20")) %>%group_by(brood_ID) %>% reframe(collection = first(collection))))
table(chic_uniq_IDs$collection)

# Xmal MI with raw values
chic_stage10 <- subset(chic_stage_size_avged,stage=="10") # 8
chic_stage10_mean <- mean(chic_stage10$embryo_dry_weight_g_mean)
chic_stage45 <- subset(chic_stage_size_avged,stage=="45") # 6
chic_stage45_mean <- mean(chic_stage45$embryo_dry_weight_g_mean)
chic_stage50 <- subset(chic_stage_size_avged,stage=="50") # 7
chic_stage50_mean <- mean(chic_stage50$embryo_dry_weight_g_mean)
chic_stage50_45 <- subset(chic_stage_size_avged,stage=="45" | stage=="50")
chic_stage50_45_mean <- mean(chic_stage50_45$embryo_dry_weight_g_mean)
chic_stage45_mean/chic_stage10_mean # 0.8494234
chic_stage50_mean/chic_stage10_mean # 0.9502302
chic_stage50_45_mean/chic_stage10_mean # 0.903704
chic_raw_table <- subset(chic_stage_size_avged,stage=="10" | stage=="45" | stage=="50")
write.csv(chic_raw_table,"Data/MI-calculation_CHIC_raw-values_stage-10-45-50.csv")

# Xmal MI with partial residuals
chic_fit <- lm(embryo_dry_weight_g_mean~collection+stage+mother_std_length,data=chic_stage_size_avged)
step(chic_fit) # embryo_dry_weight_g_mean ~ collection + mother_std_length, kept stage for downstream steps
chic_stage_partial_residuals <- visreg(chic_fit, "stage",plot=F)$res
chic_stage10 <- subset(chic_stage_partial_residuals,stage=="10")
chic_stage10_mean <- mean(chic_stage10$visregRes)
chic_stage45 <- subset(chic_stage_partial_residuals,stage=="45")
chic_stage45_mean <- mean(chic_stage45$visregRes)
chic_stage50 <- subset(chic_stage_partial_residuals,stage=="50")
chic_stage50_mean <- mean(chic_stage50$visregRes)
chic_stage50_45 <- subset(chic_stage_partial_residuals,stage=="45" | stage=="50")
chic_stage50_45_mean <- mean(chic_stage50_45$visregRes)
chic_stage45_mean/chic_stage10_mean # 0.8945757
chic_stage50_mean/chic_stage10_mean # 0.9902006
chic_stage50_45_mean/chic_stage10_mean # 0.946066
chic_pr_table <- subset(chic_stage_partial_residuals,stage=="10" | stage=="45" | stage=="50")
write.csv(chic_pr_table,"Data/MI-calculation_CHIC_partial-residuals_stage-10-45-50.csv")

## Xcortezi MI
pthc_subset <- subset(combined_embryo_data,stage!="0" & population=="PTHC" & popcoll=="PTHC_II-2023")
# average embryo size by stage within brood
pthc_stage_size_avged <- as.data.frame(
  pthc_subset %>%
    group_by(brood_ID,stage,collection) %>%
    reframe(embryo_dry_weight_g_mean = mean(embryo_dry_weight_g,na.rm=T),mother_std_length=mean(mother_std_length,na.rm=T)))
pthc_stage_size_avged

# Xcor MI with raw values
pthc_stage10 <- subset(pthc_stage_size_avged,stage=="10") # 3
pthc_stage10_mean <- mean(pthc_stage10$embryo_dry_weight_g_mean)
pthc_stage45 <- subset(pthc_stage_size_avged,stage=="45") # 4
pthc_stage45_mean <- mean(pthc_stage45$embryo_dry_weight_g_mean)
pthc_stage50 <- subset(pthc_stage_size_avged, stage=="50") # 1
pthc_stage50_mean <- mean(pthc_stage50$embryo_dry_weight_g_mean)
pthc_stage50_45 <- subset(pthc_stage_size_avged,stage=="45" | stage=="50")
pthc_stage50_45_mean <- mean(pthc_stage50_45$embryo_dry_weight_g_mean)
pthc_stage45_mean/pthc_stage10_mean # 1.000686
pthc_stage50_mean/pthc_stage10_mean # 0.7242604
pthc_stage50_45_mean/pthc_stage10_mean # 0.9484674
pthc_raw_table <- subset(pthc_stage_size_avged,stage=="10" | stage=="45" | stage=="50")
write.csv(pthc_raw_table,"Data/MI-calculation_PTHC_raw-values_stage-10-45-50.csv")

# Xcor MI with partial residuals
pthc_fit <- lm(embryo_dry_weight_g_mean~stage+mother_std_length,data=pthc_stage_size_avged)
step(pthc_fit) # embryo_dry_weight_g_mean ~ mother_std_length, kept stage for downstream steps
pthc_stage_partial_residuals <- visreg(pthc_fit, "stage",plot=F)$res
pthc_stage10 <- subset(pthc_stage_partial_residuals,stage=="10")
pthc_stage10_mean <- mean(pthc_stage10$visregRes)
pthc_stage45 <- subset(pthc_stage_partial_residuals,stage=="45")
pthc_stage45_mean <- mean(pthc_stage45$visregRes)
pthc_stage50 <- subset(pthc_stage_partial_residuals, stage=="50")
pthc_stage50_mean <- mean(pthc_stage50$visregRes)
pthc_stage50_45 <- subset(pthc_stage_partial_residuals,stage=="45" | stage=="50")
pthc_stage50_45_mean <- mean(pthc_stage50_45$visregRes)
pthc_stage45_mean/pthc_stage10_mean # 1.000686
pthc_stage50_mean/pthc_stage10_mean # 0.6995778
pthc_stage50_45_mean/pthc_stage10_mean # 0.9404641
pthc_pr_table <- subset(pthc_stage_partial_residuals,stage=="10" | stage=="45" | stage=="50")
write.csv(pthc_pr_table,"Data/MI-calculation_PTHC_partial-residuals_stage-10-45-50.csv")

## Boostrap MI using raw embryo sizes
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


## Supp Materials S4: Relationship between ovary dry weight and mother hybrid index / mitotype
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
# model fit: mother hybrid index
CALL_hi_ovary_lm_fit <- lm(ovarian_tissue_dry_weight_g ~ mother_hi*stage + brood_size, CALL_ovary_weight_unique)
selectedMod <- step(CALL_hi_ovary_lm_fit) # ovarian_tissue_dry_weight_g ~ mother_hi + stage + mother_hi:stage + brood_size
# evaluate significance with likelihood ratio test (anova) on full and reduced model
CALL_hi_ovary_lm_fit <- lm(ovarian_tissue_dry_weight_g ~ mother_hi + stage + mother_hi:stage + brood_size, CALL_ovary_weight_unique)
summary(CALL_hi_ovary_lm_fit)
CALL_hi_ovary_lm_fit_reduced <- lm(ovarian_tissue_dry_weight_g ~ stage + brood_size, CALL_ovary_weight_unique)
anova(CALL_hi_ovary_lm_fit,CALL_hi_ovary_lm_fit_reduced)
# 2     75 0.0016711 -5 -0.00028619 2.8931 0.01971 *

# plot mother hi and ovary dry weight
CALL_hi_ovary_lm_vr <- visreg(CALL_hi_ovary_lm_fit,"mother_hi")
plot(visregRes ~ mother_hi,CALL_hi_ovary_lm_vr$res,xlab="mother ancestry proportion",ylab="partial residuals of ovary dry weight (g)", col=hybrid_col,pch=20,cex.lab=1.8,cex.axis=1.5)
abline(lm(CALL_hi_ovary_lm_vr$res$visregRes~CALL_hi_ovary_lm_vr$res$mother_hi), col="purple",lwd=4)

## Relationship between ovary dry weight and mother mitotype
# model fit
CALL_mito_ovary_lm_fit <- lm(ovarian_tissue_dry_weight_g ~ mitotype*stage + brood_size, CALL_ovary_weight_unique)
selectedMod <- step(CALL_mito_ovary_lm_fit) # ovarian_tissue_dry_weight_g ~ stage + brood_size
summary(CALL_mito_ovary_lm_fit) # mitotypemalinche          1.210e-03  2.620e-03   0.462  0.64561
CALL_mito_ovary_lm_fit <- lm(ovarian_tissue_dry_weight_g ~ mitotype + stage + brood_size, CALL_ovary_weight_unique)
CALL_mito_ovary_lm_fit_reduced <- lm(ovarian_tissue_dry_weight_g ~ stage + brood_size, CALL_ovary_weight_unique)
anova(CALL_mito_ovary_lm_fit,CALL_mito_ovary_lm_fit_reduced)
# 2     75 0.0016711 -1 -6.2071e-07 0.0275 0.8687

## correlation between mother hybrid index and mitotype
temp <- CALL_embryo_data
temp$mito_num <- ifelse(temp$mitotype=="birchmanni",0,1)
temp <- subset(temp,!is.na(mother_hi) & !is.na(mito_num))
cor(temp$mother_hi,temp$mito_num,method="spearman") # spearman rho=0.69

