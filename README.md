# Recent evolution of large offspring size and post-fertilization nutrient provisioning in swordtails
### Payne et al 2025. Proc B.

This repository contains all code and accompanying input data files that generated the results
presented in this study.

## Abstract

Organisms have evolved diverse reproductive strategies that impact the probability that their offspring survive to adulthood. Using morphological measurements in embryos and fry, gene expression analysis, and hybrid crosses, we describe divergence in reproductive strategy between two closely related species of swordtail fish (Xiphophorus), which have internal fertilization and give birth to free-swimming fry. We find that one species, X. malinche, which lives in high-elevation environments, has evolved larger offspring than its closest relative X. birchmanni and dwarfs the offspring size of other species in the genus. The fry of X. malinche are more resilient to starvation than X. birchmanni, hinting that the evolution of large offspring size could be an adaptation to the challenging environments in which X. malinche are born. We also find evidence that X. malinche mothers provision nutrients to their offspring during embryonic development, the first time this process has been documented in the Xiphophorus genus. Moreover, in the ovary, we observe differential regulation of genes associated with maternal nutrient provisioning in other groups that use this reproductive strategy. Finally, we generated hybrid crosses between X. malinche and X. birchmanni to explore the impact of genetics and maternal environment on offspring size, finding that offspring size is at least in part genetically determined. Intriguingly, we find a low rate of survival in one cross direction and investigate the links between reproductive strategy and this asymmetric hybrid incompatibility.

---

## Table of Contents
1. [File Inventory](#file-inventory)
2. [Software Requirements](#software-requirements)
3. [Data Files and Variable Definitions](#data-files-and-variable-definitions)
4. [Glossary of Terms](#glossary-of-terms)
5. [Contact Information](#contact-information)

---

## File Inventory

### Analysis Scripts (Root Directory)

**Main Analysis Scripts:**
- `Payne-et-al-2025_Figures.R` - Code to generate all main text and supplementary figures presented in the manuscript
- `Payne-et-al-2025_Tables.R` - Code to generate supplementary tables and perform statistical analyses referenced in the manuscript

**RNA-Seq Analysis Scripts:**
- `embryo-DGE_DESeq2_xmac-IDs_2023.R` - Differential gene expression (DGE) analysis for embryo RNA-seq data using DESeq2
- `ovary-DGE_DESeq2_xmac-IDs_2023.R` - Differential gene expression (DGE) analysis for ovary RNA-seq data using DESeq2
- `ovary-GO-analysis_xmac-IDs_2023.R` - Gene Ontology (GO) enrichment analysis for differentially expressed genes in ovary tissue
- `ovary-WGCNA_xmac-IDs_combined2023.R` - Weighted Gene Co-expression Network Analysis (WGCNA) to identify co-expressed gene modules in ovary tissue
- `ovary-GO-analysis-WGCNA_xmac-IDs_2023.R` - Gene Ontology enrichment analysis for gene modules identified by WGCNA

### Data Files

**Data/** (Main data directory containing all input files for analyses)

#### Morphological and Experimental Data Files:
- `newborn_fry_size_data.csv` - Morphological measurements of newborn fry across multiple species and populations
- `all-pops_combined_embryo_datasets.csv` - Combined embryo dry weight measurements across all populations and developmental stages
- `CHIC_COAC_PTHC_fully-yolked-stage0_dry-weights_mother-length.csv` - Dry weights of fully-yolked (stage 0) embryos with maternal length data
- `IV-2023_roof-tank_CHIC-COAC-F1_embryo-weights.csv` - Embryo dry weights from F1 hybrid crosses (X. malinche × X. birchmanni)
- `III-2023_TETI2_TETIxCHIC_CHICxCHIC_embryo-dry-weights.csv` - Embryo dry weights of intraspecific X. malinche crosses (TETIxTETI, TETIxCHIC, CHICxCHIC)
- `CALL_mother_hybrid-index_mitotype.csv` - Genetic hybrid index and mitochondrial type data for mothers from the natural hybrid zone (Calnali)
- `newborn_fry_size_data.csv` - Standard length and head width measurements for newborn fry
- `X-23_fry_starvation_fat_content.csv` - Fat content measurements from fry food deprivation experiments
- `X-23_fry_starvation_initial_standard_lengths.csv` - Initial size measurements for fry included in food deprivation experiments
- `TableS21_food-deprivation-expt_data.csv` - Formatted and combined data from food deprivation experiment
- `CTmin_Xmal_Xbirch_newborn_trial_data.csv` - Critical thermal minimum (CTmin) data for newborn X. malinche and X. birchmanni fry
- `pregnancy_rate_embryo_collections.csv` - Pregnancy rates across different collection sites, seasons, and species
- `II-2023_IX-2023_nonpreg-females_lipid-extraction.csv` - Lipid content measurements from non-pregnant females

#### RNA-Seq Data Directory:
**Data/embryo_ovary_dge_combined_2023/** (Contains all RNA-seq analysis results)

**Subdirectories:**

1. **embryo_xmac-gtf_dge/** - Embryo tissue differential gene expression results
   - `embryo-xmacID-combined2023_dge_dds.rds` - DESeq2 dataset object (R binary format)
   - `embryo-xmacID-combined2023_dge_vst.rds` - Variance stabilized transformed expression data (R binary format)
   - `embryo-xmacID-combined2023_dge_lfc-shr_all.csv_with-annots.csv` - Log2 fold changes and adjusted p-values for all pairwise comparisons with gene annotations

2. **ovary_xmac-gtf_dge/** - Ovary tissue differential gene expression results
   - `ovary-xmacID-combined2023_dge_dds.rds` - DESeq2 dataset object (R binary format)
   - `ovary-xmacID-combined2023_dge_vst.rds` - Variance stabilized transformed expression data (R binary format)
   - `ovary-xmacID-combined2023_dge_lfc-shr_all.csv_with-annots.csv` - Log2 fold changes and adjusted p-values for all pairwise comparisons with gene annotations

3. **ovary_embryo_kallisto_output_FebAug23_combined/** - RNA-seq sample metadata
   - `ovary_embryo_rnaseq_samples_FebAug23_combined.txt` - Sample information including species, tissue type, developmental stage, and sequencing batch

4. **wgcna/** - Weighted Gene Co-expression Network Analysis results (185 files total)
   - For each co-expression module (identified by color name), there are 5 files:
     - `ovary-xmac_combined-FebAug23_GO-WGCNA_[MODULE]_genes.csv` - List of genes in the module
     - `ovary-xmac_combined-FebAug23_GO-WGCNA_[MODULE]_GO-gene-annots.csv` - Gene Ontology annotations for genes in the module
     - `ovary-xmac_combined-FebAug23_GO-WGCNA_[MODULE]_pval0.05_BP.tsv` - Enriched Biological Process GO terms (p < 0.05)
     - `ovary-xmac_combined-FebAug23_GO-WGCNA_[MODULE]_pval0.05_CC.tsv` - Enriched Cellular Component GO terms (p < 0.05)
     - `ovary-xmac_combined-FebAug23_GO-WGCNA_[MODULE]_pval0.05_MF.tsv` - Enriched Molecular Function GO terms (p < 0.05)
   - `ovary-xmac_combined-FebAug23_MEtraitpvals.csv` - Module eigengene-trait correlation p-values
   - `ovary-xmac_combined-FebAug23_GO-WGCNA_sigMEs.RData` - Significant module eigengenes (R binary format)

   Module colors included: black, blue, brown, brown4, cyan, darkgreen, darkgrey, darkmagenta, darkolivegreen, darkorange, darkorange2, darkturquoise, floralwhite, green, greenyellow, ivory, lightcyan1, lightgreen, lightsteelblue1, lightyellow, magenta, mediumpurple3, midnightblue, orangered4, paleturquoise, pink, plum1, purple, red, royalblue, salmon, sienna3, skyblue3, steelblue, turquoise, violet, white, yellow, yellowgreen

---

## Software Requirements

### Required Software (Free and Open Source)

**For running R analysis scripts:**
- **R** (version 4.0 or higher) - Download from [https://www.r-project.org/](https://www.r-project.org/)
- **RStudio** (recommended IDE) - Download from [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/)

**Required R packages:**
- `lme4` - Linear mixed-effects models
- `emmeans` - Estimated marginal means for post-hoc comparisons
- `DESeq2` - Differential gene expression analysis
- `tximport` - Import transcript-level abundance estimates
- `WGCNA` - Weighted gene co-expression network analysis
- `GO.db` - Gene Ontology database
- `ggplot2` - Data visualization
- `dplyr` - Data manipulation
- `tidyr` - Data tidying

---

## Data Files and Variable Definitions

### 1. newborn_fry_size_data.csv
Morphological measurements of newborn fry (free-swimming offspring at birth) from multiple Xiphophorus species, populations, and crosses.

**Variables:**
- `born_on_date` - Date when fry were born (format: DAY-MONTH-YEAR)
- `born_month` - Month of birth (Roman numerals: I-XII)
- `born_year` - Year of birth (two-digit format)
- `collection` - Collection identifier combining month and year
- `species_site` - Combined species and collection site identifier
- `species` - Species name (abbreviated: Xmal = X. malinche, Xbir = X. birchmanni, Xcor = X. cortezi, Xpyg = X. pygmaeus, Xvar = X. variatus, hyb = hybrid)
- `site` - Geographic collection site or population code
- `parents` - Information about parental cross (if applicable)
- `tank_id` - Laboratory tank identifier
- `brood_no` - Brood number for a given female
- `brood_ID` - Unique identifier for each brood
- `individual` - Individual fry number within a brood
- `sl_mm` - Standard length in millimeters (mm) - measured from tip of snout to posterior end of caudal peduncle
- `hw_mm` - Head width in millimeters (mm) - measured at widest point of head
- `hw/sl` - Ratio of head width to standard length (dimensionless)
- `notes` - Any additional observations or comments
- `gen` - Generation (F1, F2, or UNK = unknown)

---

### 2. all-pops_combined_embryo_datasets.csv
Combined dataset of embryo dry weight measurements across all populations, species, and developmental stages.

**Variables:**
- `brood_ID` - Unique identifier for each brood
- `species` - Species name (Xmalinche, Xbirchmanni, Xcortezi)
- `population` - Population or site code (CHIC, COAC, PTHC, etc.)
- `collection` - Collection identifier (month-year format)
- `female_ID` - Female individual identifier within a collection
- `embryo_ID` - Embryo individual identifier within a brood
- `stage` - Developmental stage as defined in Reznick (1981)
- `brood_size` - Total number of embryos in the brood
- `embryo_dry_weight_g` - Dry weight of individual embryo in grams (g) - measured after drying to constant mass
- `ovarian_tissue_dry_weight_g` - Dry weight of attached ovarian tissue in grams (g) - NA if not included
- `mother_std_length` - Mother's standard length in centimeters (cm)
- `mother_total_wet_mass` - Mother's total wet mass in grams (g)
- `mother_no_ovary_wet_mass` - Mother's wet mass excluding ovary in grams (g)
- `popcoll` - Combined population and collection identifier
- `season` - Season of collection (warm, cold)

---

### 3. CHIC_COAC_PTHC_fully-yolked-stage0_dry-weights_mother-length.csv
Dry weights of fully-yolked eggs (developmental stage 0, unfertilized, even distribution of yolk) with maternal morphological data.

**Variables:**
- `brood_ID` - Unique identifier for each brood
- `species` - Species name (Xmalinche, Xbirchmanni, Xcortezi)
- `population` - Population code (CHIC-UP, COAC, PTHC)
- `collection` - Collection identifier
- `female_ID` - Female individual identifier
- `fresh_or_etoh.preserved` - Sample preservation method (fresh or ethanol-preserved)
- `brood_size` - Total number of embryos in the brood
- `sample_storage` - Storage conditions (e.g., "2ml etoh")
- `embryo_ID` - Embryo individual identifier within a brood
- `stage` - Developmental stage (all 0 = fully-yolked, unfertilized egg, as defined in Reznick (1981))
- `embryo_dry_weight_g` - Dry weight of individual embryo in grams (g)
- `ovarian_tissue_dry_weight_g` - Dry weight of attached ovarian tissue in grams (g)
- `Comments` - Descriptive notes about lipid distribution and yolk appearance
- `mother_std_length` - Mother's standard length in centimeters (cm)

---

### 4. IV-2023_roof-tank_CHIC-COAC-F1_embryo-weights.csv
Embryo dry weight measurements from lab-reared F1 hybrid crosses between X. malinche (CHIC population) and X. birchmanni (COAC population).

**Variables:**
- `brood_ID` - Unique identifier for each brood
- `species` - Cross type (malxbir_F1 = X. malinche × X. birchmanni F1 hybrid)
- `population` - Population identifier for hybrid cross
- `collection` - Collection identifier
- `female_ID` - Female individual identifier
- `fresh_or_etoh-preserved` - Sample preservation method
- `brood_size` - Total number of embryos in the brood
- `sample_storage` - Storage conditions
- `ovary_tissue_included` - Whether ovarian tissue was included in measurements (yes/no)
- `embryo_ID` - Embryo individual identifier within a brood
- `stage` - Developmental stage as defined in Reznick (1981)
- `embryo_dry_weight_g` - Dry weight of individual embryo in grams (g)
- `ovarian_tissue_dry_weight_g` - Dry weight of attached ovarian tissue in grams (g) - NA if not included
- `comments` - Descriptive notes about developmental stage
- `mother_std_length` - Mother's standard length in centimeters (cm)
- `mother_total_wet_mass` - Mother's total wet mass in grams (g)
- `mother_no_ovary_wet_mass` - Mother's wet mass excluding ovary in grams (g)

---

### 5. III-2023_TETI2_TETIxCHIC_CHICxCHIC_embryo-dry-weights.csv
Embryo dry weight measurements of embryos from intraspecific X. malinche crosses.

**Variables:**
- `brood_ID` - Unique identifier for each brood
- `species` - Species name (Xmalinche)
- `population` - Population code (TETI2 = TETIxTETI embryos from wild-caught mothers; TETIxCHIC = TETI×CHIC embryos from lab cross; CHICxCHIC = CHICxCHIC embryos from lab-reared mothers)
- `collection` - Collection identifier
- `female_ID` - Female individual identifier
- `fresh_or_etoh-preserved` - Sample preservation method
- `brood_size` - Total number of embryos in the brood
- `sample_storage` - Storage conditions
- `ovary_tissue_included` - Whether ovarian tissue was included in measurements (yes/no)
- `embryo_ID` - Embryo individual identifier within a brood
- `stage` - Developmental stage (0 or 35)
- `embryo_dry_weight_g` - Dry weight of individual embryo in grams (g)
- `ovarian_tissue_dry_weight_g` - Dry weight of attached ovarian tissue in grams (g)
- `comments` - Descriptive notes
- `mother_std_length` - Mother's standard length in centimeters (cm)
- `mother_total_wet_mass` - Mother's total wet mass in grams (g)
- `mother_no_ovary_wet_mass` - Mother's wet mass excluding ovary in grams (g)
- `mother_origin` - Origin of mother (lab or wild)

---

### 6. CALL_mother_hybrid-index_mitotype.csv
Genetic characterization of mothers from the natural hybrid zone in Calnali (CALL), including hybrid index (proportion of X. malinche ancestry) and mitochondrial genotype.

**Variables:**
- `mother_ID` - Unique identifier for each mother
- `brood_ID` - Unique identifier for each brood
- `mom_index` - Hybrid index (proportion of X. malinche ancestry; 0 = pure X. birchmanni, 1 = pure X. malinche)
- `ndufs5` - Genotype at ndufs5 mitochondrial-diagnostic nuclear marker (birchmanni, malinche, or heterozygous)
- `ndufa13` - Genotype at ndufa13 mitochondrial-diagnostic nuclear marker (birchmanni, malinche, heterozygous, or NA)
- `ScyDAA6-5984-HRSCAF-6694:3367085` - Genotype at additional diagnostic locus (birchmanni, malinche, or heterozygous)
- `mitotype` - Mitochondrial haplotype (birchmanni or malinche)

---

### 7. TableS21_food-deprivation-expt_data.csv
Data from experimental food deprivation trials measuring fry growth, survival, and fat content under starvation versus control (fed) conditions.

**Variables:**
- `brood_no` - Brood number identifier
- `dissection_date` - Date of final measurement and dissection (format: D-MONTH-YY)
- `species` - Species (Xbir = X. birchmanni, Xmal = X. malinche)
- `population` - Population code (COAC, CHIC)
- `number` - Individual fish number within treatment group
- `treatment` - Experimental treatment (starved or control/fed)
- `standard_length_cm` - Final standard length in centimeters (cm)
- `wet_mass` - Final wet mass in grams (g)
- `dry_mass_1` - First dry mass measurement in grams (g)
- `dry_mass_2` - Second dry mass measurement in grams (g)
- `FC_percent` - Fat content as percentage (%) calculated from difference between dry mass measurements
- `notes` - Additional observations
- `condition` - Combined species and treatment identifier
- `avg_initial_std_length_mm` - Average initial standard length in millimeters (mm) for the brood
- `standard_length_mm` - Final standard length in millimeters (mm) (converted from cm)
- `std_length_normalized` - Standard length normalized by initial average (dimensionless ratio)
- `std_length_growth_rate` - Growth rate in millimeters per day (mm/day)

---

### 8. X-23_fry_starvation_fat_content.csv
Fat content measurements from fry starvation experiments.

**Variables:**
- `dissection_date` - Date of measurement (format: D-MONTH-YY)
- `species` - Species (Xmal = X. malinche, Xbir = X. birchmanni)
- `population` - Population code (CHIC, COAC)
- `brood_no` - Brood number identifier
- `number` - Individual fish number
- `treatment` - Experimental treatment (starved or control)
- `standard_length_cm` - Standard length in centimeters (cm)
- `wet_mass` - Wet mass in grams (g)
- `dry_mass_1` - First dry mass measurement in grams (g) - before lipid extraction
- `dry_mass_2` - Second dry mass measurement in grams (g) - after lipid extraction
- `FC_percent` - Fat content as percentage (%) calculated as ((dry_mass_1 - dry_mass_2) / dry_mass_1) × 100
- `notes` - Additional observations

---

### 9. X-23_fry_starvation_initial_standard_lengths.csv
Initial size measurements for fry used in starvation experiments (measured at start of experiment).

**Variables:**
- `dissection_date` - Date experiment began
- `species` - Species (Xmal, Xbir)
- `population` - Population code (CHIC, COAC)
- `brood_no` - Brood number identifier
- `treatment` - Always "initial" (indicating initial measurement)
- `standard_length_cm` - Initial standard length in centimeters (cm)

---

### 10. CTmin_Xmal_Xbirch_newborn_trial_data.csv
Critical thermal minimum (CTmin) data for newborn X. malinche and X. birchmanni. CTmin is the temperature at which fish lose equilibrium and swimming ability when temperature is gradually decreased.

**Variables:**
- `newborn_ID` - Unique identifier for each newborn fish
- `species` - Species (Xmal = X. malinche, Xbir = X. birchmanni)
- `trial_date` - Date of CTmin trial
- `trial_number` - Trial replicate number
- `fish_number` - Fish number within trial
- `CTmin_C` - Critical thermal minimum temperature in degrees Celsius (°C)
- `CTmin_Time` - Time to reach CTmin in minutes (min)
- `start_temp_C` - Starting temperature of trial in degrees Celsius (°C)
- `born_on_date` - Date fish was born
- `CTmin_trial` - Trial identifier
- `individual` - Individual fish identifier
- `sl(mm)` - Standard length in millimeters (mm)
- `hw(mm)` - Head width in millimeters (mm)
- `hw/sl` - Ratio of head width to standard length (dimensionless)

---

### 11. pregnancy_rate_embryo_collections.csv
Pregnancy rates across different collection sites, seasons, and species. Pregnancy rate is the proportion of females carrying developing embryos at time of collection.

**Variables:**
- `collection` - Collection identifier (site_month-year)
- `month` - Month of collection (full name)
- `year` - Year of collection (four-digit)
- `site` - Geographic collection site code
- `species` - Species name (Xmalinche, Xbirchmanni, Xcortezi)
- `sample_size` - Total number of females sampled
- `pregnant_N` - Number of pregnant females (containing developing embryos)
- `yolking_yolked_nonpreg_N` - Number of non-pregnant females with yolking or yolked follicles
- `nonyolked_nonpreg_N` - Number of non-pregnant females without yolked follicles
- `pregnancy_rate` - Proportion of pregnant females (pregnant_N / sample_size)
- `data_origin` - Source of data
- `notes` - Additional notes about sampling design

---

### 12. II-2023_IX-2023_nonpreg-females_lipid-extraction.csv
Lipid content measurements from non-pregnant female fish collected across seasons.

**Variables:**
- `dissection_date` - Date of dissection
- `species` - Species name (X. birchmanni, X. malinche)
- `population` - Population code
- `collection` - Collection identifier
- `popcoll` - Combined population and collection identifier
- `female_ID` - Female individual identifier
- `female_std_length` - Female standard length in centimeters (cm)
- `female_total_wet_mass` - Female total wet mass in grams (g)
- `dry_mass_1` - First dry mass measurement in grams (g) - before lipid extraction
- `dry_mass_2` - Second dry mass measurement in grams (g) - after lipid extraction
- `FC_g` - Fat content in grams (g)
- `FC_percent` - Fat content as percentage (%)
- `sc` - Spawning condition indicator
- `notes` - Additional observations

---

### 13. RNA-Seq Data Files

#### ovary_embryo_rnaseq_samples_FebAug23_combined.txt
Sample metadata for RNA-sequencing experiments on ovary and embryo tissues.

**Variables:**
- `file_basename` - Base filename for sequencing data
- `sample_ID` - Short sample identifier
- `og_sample_id` - Original sample identifier
- `batch` - Sequencing batch identifier (1R1, 1R2, etc.)
- `species` - Species (Xbirchmanni, Xmalinche, Xcortezi)
- `population` - Population code (CHIC, COAC, PTHC)
- `tissue` - Tissue type (ovary or embryo)
- `embryo_stage` - Developmental stage as defined in Reznick (1981)
- `stage_category` - Developmental category (early or late)
- `stage_group` - Combined species and stage category (e.g., mal_late, bir_early)
- `origin` - Sample origin (lab)
- `collection` - Collection date identifier (2023Aug or 2023Feb)

---

#### ovary-xmacID-combined2023_dge_lfc-shr_all.csv_with-annots.csv
Differential gene expression results for ovary tissue comparing species and developmental stages.

**Variables:**
- `Unnamed: 0` - Row index
- `Gene` - Ensembl gene identifier (ENSXMAG########)
- `LFC_res.corvbirch_early` - Log2 fold change for X. cortezi vs. X. birchmanni at early stage
- `padj_res.corvbirch_early` - Adjusted p-value (Benjamini-Hochberg) for X. cortezi vs. X. birchmanni at early stage
- `LFC_res.corvbirch_late` - Log2 fold change for X. cortezi vs. X. birchmanni at late stage
- `padj_res.corvbirch_late` - Adjusted p-value for X. cortezi vs. X. birchmanni at late stage
- `LFC_res.latevearly_birch` - Log2 fold change for late vs. early stage in X. birchmanni
- `padj_res.latevearly_birch` - Adjusted p-value for late vs. early stage in X. birchmanni
- `LFC_res.latevearly_cor` - Log2 fold change for late vs. early stage in X. cortezi
- `padj_res.latevearly_cor` - Adjusted p-value for late vs. early stage in X. cortezi
- `LFC_res.latevearly_mal` - Log2 fold change for late vs. early stage in X. malinche
- `padj_res.latevearly_mal` - Adjusted p-value for late vs. early stage in X. malinche
- `LFC_res.malvbirch_early` - Log2 fold change for X. malinche vs. X. birchmanni at early stage
- `padj_res.malvbirch_early` - Adjusted p-value for X. malinche vs. X. birchmanni at early stage
- `LFC_res.malvbirch_late` - Log2 fold change for X. malinche vs. X. birchmanni at late stage
- `padj_res.malvbirch_late` - Adjusted p-value for X. malinche vs. X. birchmanni at late stage
- `LFC_res.malvcor_early` - Log2 fold change for X. malinche vs. X. cortezi at early stage
- `padj_res.malvcor_early` - Adjusted p-value for X. malinche vs. X. cortezi at early stage
- `LFC_res.malvcor_late` - Log2 fold change for X. malinche vs. X. cortezi at late stage
- `padj_res.malvcor_late` - Adjusted p-value for X. malinche vs. X. cortezi at late stage
- Sample columns (1C2O, 1M6O, etc.) - Variance-stabilized expression values for each sample
- `annot` - Gene name/symbol annotation

Note: Positive log2 fold changes indicate higher expression in the first listed species/stage; negative values indicate higher expression in the second listed species/stage. Adjusted p-values < 0.05 are typically considered statistically significant.

---

#### embryo-xmacID-combined2023_dge_lfc-shr_all.csv_with-annots.csv
Differential gene expression results for embryo tissue (same format as ovary file described above).

---

#### WGCNA Gene Module Files

**Gene List Files** (`ovary-xmac_combined-FebAug23_GO-WGCNA_[MODULE]_genes.csv`):
- `x` - Ensembl gene identifier for genes belonging to the co-expression module

**Gene Annotation Files** (`ovary-xmac_combined-FebAug23_GO-WGCNA_[MODULE]_GO-gene-annots.csv`):
- `frame.go_id` - Gene Ontology term identifier
- `frame.gene_name` - Gene symbol/name
- `frame.gene_id` - Gene identifier
- `frame.ensembl` - Ensembl gene identifier
- `frame.Evidence` - Evidence code for GO annotation

**GO Enrichment Files** (`ovary-xmac_combined-FebAug23_GO-WGCNA_[MODULE]_pval0.05_[BP/CC/MF].tsv`):
BP = Biological Process, CC = Cellular Component, MF = Molecular Function

- `GOBPID` (or GOCCID, GOMFID) - Gene Ontology term identifier
- `Pvalue` - Enrichment p-value (unadjusted)
- `OddsRatio` - Odds ratio for enrichment
- `ExpCount` - Expected number of genes in this GO term
- `Count` - Observed number of genes from module in this GO term
- `Size` - Total number of genes annotated to this GO term in reference
- `Term` - GO term description
- `GENEID_list` - Comma-separated list of Ensembl gene IDs
- `GENENAME_list` - Comma-separated list of gene names/symbols

---

#### Binary R Data Files (.rds and .RData)

These files contain R objects and require R software to open:

**DESeq2 Dataset Objects** (`*_dge_dds.rds`):
- Complete DESeq2 analysis object containing count data, size factors, dispersion estimates, and statistical results
- Load in R using: `dds <- readRDS("filename.rds")`

**Variance Stabilized Transformed Data** (`*_dge_vst.rds`):
- Expression values after variance stabilizing transformation, suitable for visualization and clustering
- Load in R using: `vst <- readRDS("filename.rds")`

**WGCNA Module Eigengenes** (`ovary-xmac_combined-FebAug23_GO-WGCNA_sigMEs.RData`):
- Module eigengenes (first principal component of each gene module) and their correlations with traits
- Load in R using: `load("filename.RData")`

---

## Glossary of Terms

### Species Abbreviations
- **Xmal** or **X. malinche** - *Xiphophorus malinche*, high-elevation swordtail species from Mexico
- **Xbir** or **X. birchmanni** - *Xiphophorus birchmanni*, lowland swordtail from Mexico, sister species to *X. malinche*
- **Xcor** or **X. cortezi** - *Xiphophorus cortezi*, swordtail from Mexico, sister species to the *X. malinche* and *X. birchmanni* clade
- **Xpyg** or **X. pygmaeus** - *Xiphophorus pygmaeus*, swordtail species from Mexico
- **Xvar** or **X. variatus** - *Xiphophorus variatus*, platyfish species from Mexico

### Population Codes
- **CHIC** - Chicayotla, high-elevation *X. malinche* population in a tributary of Río Xontla in Hidalgo, Mexico 
- **COAC** - Coacuilco, lowland *X. birchmanni* population on the Río Coacuilco in Hidalgo, Mexico 
- **CALL** - Calnali Low, natural hybrid population between *X. malinche* and *X. birchmanni* on the Río Calnali in Hidalgo, Mexico
- **PTHC** - Puente de Huichihuayán, *X. cortezi* population on the Río Santa Cruz in San Luis Potosí, Mexico
- **TETI** - Tetipanchalco, *X. malinche* population

### Morphology/Trait Terms
- **Standard length (SL)** - Distance from the tip of the snout to the posterior end of the caudal peduncle (base of tail fin), measured in mm or cm
- **Head width (HW)** - Maximum width of the head, measured in mm
- **Dry weight** - Weight after drying to constant mass measured in grams (g)
- **Wet mass** - Fresh weight before drying, measured in grams (g)
- **Fat content (FC)** - Lipid content measured as difference in dry mass before and after lipid extraction, expressed as percentage
- **CTmin (Critical Thermal Minimum)** - Temperature at which fish lose equilibrium and cannot maintain swimming position when cooled gradually

### Developmental Stages (as defined in Reznick (1981))
Reznick D. 1981 ‘Grandfather Effects’: The Genetics of Interpopulation Differences in Offspring Size in the Mosquito Fish. Evolution 35, 941–953. (doi:10.2307/2407865)

Embryonic development in *Xiphophorus* is divided into stages based on morphological landmarks:
- **0** -	no development, even lipid distribution
- **2** -	blastodisc formed
- **5** - nerula starts to form or is already present
- **10** - head and optic cups are visible (under high magnification)
- **15** - faint eye pigmentation (grey eyed)
- **20** - eye pigmentation (black) and caudal fin bud
- **25** - pectoral fin buds and body pigmentation
- **30** - caudal fin rays clearly visible and dorsal and anal fin buds
- **35** - pectoral fin rays clearly visible
- **40** - dorsal and anal fin rays visible
- **45** - oval shape, operculae conspicuous, scales visible
- **50** - pericardial cavity almost or completely closed

### Bioinformatics and Statistical Terms
- **DGE (Differential Gene Expression)** - Statistical comparison of gene expression levels between conditions
- **DESeq2** - R package for differential expression analysis of RNA-seq data using negative binomial models
- **Log2 fold change (LFC)** - Base-2 logarithm of expression ratio between conditions; LFC=1 means 2-fold higher, LFC=-1 means 2-fold lower
- **Adjusted p-value (padj)** - P-value corrected for multiple testing using Benjamini-Hochberg method; padj < 0.05 typically indicates statistical significance
- **VST (Variance Stabilizing Transformation)** - Transformation that normalizes variance across expression levels for visualization and clustering
- **WGCNA (Weighted Gene Co-expression Network Analysis)** - Method to identify clusters (modules) of genes with correlated expression patterns
- **Module** - Group of co-expressed genes identified by WGCNA, labeled by color names
- **Module eigengene** - First principal component of gene expression in a module, representing the module's overall expression pattern
- **GO (Gene Ontology)** - Standardized vocabulary describing gene functions in three categories: Biological Process (BP), Cellular Component (CC), and Molecular Function (MF)
- **GO enrichment** - Statistical test to identify GO terms that are over-represented in a gene set compared to background
- **Ensembl gene ID** - Unique identifier for genes in the Ensembl database (format: ENSXMAG########)


---

## Contact Information

Questions, comments, and suggestions are always welcome. Please feel free to contact the repository owner:

**Cheyenne Y. Payne**:
Email: cypayne [at] ucsc [dot] edu

For issues or bugs with code, please open an issue on the GitHub repository.

---

## Citation

If you use data or code from this repository, please cite:

Payne et al. (2025). Recent evolution of large offspring size and post-fertilization nutrient provisioning in swordtails. *Proceedings of the Royal Society B*.

---

