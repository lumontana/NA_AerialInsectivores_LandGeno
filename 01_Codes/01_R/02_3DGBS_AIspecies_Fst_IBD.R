# Info --------------------------------------------------------------------
# 
# Author: Luca Montana
# Affiliation: University of Lethbridge
# Group: T. Burg research group
# Location: Rimouski, QC
# Date: 2024-10-31
# 
# Overview: Genetic distance and its causes among four North American AI species -
#           Alder flycatchers, cliff swallows, purple martins, violet-green swallows
#

# Library -----------------------------------------------------------------

library(readxl)
library(tidyverse)
library(vcfR)
# library(adegenet)  # tab function
library(dartR)  # gl.fst.pop
library(QuickPop)
library(patchwork)
library(hierfstat)

# Internal functions
# not in
`%nin%` <- Negate(`%in%`)

# count proportion of NA loci per individual in a genlight object
count.ind.na.gl <- function(gl){
  res <- apply(tab(gl,  NA.method = c("asis")), MARGIN = 1, FUN = function(l){n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  return(res)
}

# PCA functions
plot_pca_eig <- list(geom_bar(stat = "identity"),
                     geom_line(),
                     scale_fill_viridis_c(),
                     labs(y = "% variance", title = NULL, x = "PC axis"),
                     theme_bw(),
                     theme(axis.text.x = element_blank(),
                           panel.grid = element_blank(),
                           axis.ticks.length.y = unit(0.1, "in"),
                           axis.ticks.length.x = unit(0, "in"),
                           legend.position = "none"))

# Fst functions
table.fst <- function(fst){
  res <-  fst$Bootstraps %>% dplyr::select(Population1, Population2, "Lower bound CI limit", "Upper bound CI limit", "p-value", "Fst")
  
  return(res)
  
}

heat.fst <- function(fst){
  res <- bind_rows(table.fst(fst),
                   table.fst(fst) %>% mutate(Pop1.int = Population2,
                                             Pop2.int = Population1,
                                             Population1 = Pop1.int,
                                             Population2 = Pop2.int) %>%
                     dplyr::select(-c(Pop1.int, Pop2.int))
  )
  
  return(res)
  
}


# Data --------------------------------------------------------------------

## Upload metadata --------------------------------------------------------

d_meta <- read.csv("./00_Data/00_fileinfo/AI_metadata.csv", stringsAsFactors = F)
# I made a manual modification to this file after its creation on script 00_DataPrep.R
# I addedd that PUMA from SK were sampled during the breeding season.


## Upload vcf files and convert to genlight formant -----------------------

# The following VCF files are issued from 3DGBS reads that have been aligned on species-specific reference genomes,
# loci have been assembled using stacks 2 (gstacks), and the fist VCFs were the outputs of stacks 'populations' unit (min-mac 3, random snp)
# VCFs have been filters for lmiss, imiss, and maxDP (outlier SNPs with very high coverage deptth removed only). 
# Outlier SNPs (under directional and balancing selection) have also been removed

### ALFL ------------------------------------------------------------------

# Populations: --min-mac 3 --write-random-snp
# 58 samples
# Removed 0 loci that did not pass sample/population constraints from 2469285 loci.
# Kept 2469285 loci, composed of 417099871 sites; 2025363 of those sites were filtered, 266864 variant sites remained.
# Mean genotyped sites per locus: 132.93bp (stderr 0.03).

# Missingness
# --max-missing 0.3: After filtering, kept 46995 out of a possible 266864 Sites
# imiss --keep inds_to_keep_70_90.txt: After filtering, kept 36 out of 58 Individuals

# MaxDP
# --max-meanDP 1.6e+03: After filtering, kept 46992 out of a possible 46995 Sites

# Neutral - bayescan
# What the code does is actually removing loci that should have been removed with the normal filtering prior to ideintifying 
# loci under selection...
# After filtering, kept 20673 out of a possible 46992 Sites

#### Neutral SNPs ---------------------------------------------------------

# from ALFL_ALL/ALFL_ALL_merged/01_reference/00_alfl/8_neutral/neutral/
vcf.alfl <- vcfR::read.vcfR("./00_Data/03_filtering/03_neutral/ALFL/alfl_7090_neutral_maxDP.recode.vcf")  # 20673 SNPs
head(vcf.alfl)

# Count how many RAD loci left from lmiss and imiss filters
SCAFFOLD.info <- vcf.alfl@fix %>% as.data.frame() %>%  # 
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, "_"), `[`,2), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )
length(unique(SCAFFOLD.info$RADloc))  # 20673 loci for 20673 SNPs - 1 SNP per locus, excellent

# Check names of vcf fields
vcf_field_names(vcf.alfl, tag = "FORMAT")

# Conversion to genlight and genind R formats - necessary to run glPCA from adegenet package
gl.alfl <- vcfR::vcfR2genlight(vcf.alfl)
# gi.data.alfl  <- vcfR::vcfR2genind(vcf.alfl) 

# Rename samples to uniformize suffix in genlight, genind and d_meta (also to find out which one have been merged...)
c.names.gl <- indNames(gl.alfl)
n.names.gl <- gsub("\\.sorted$", "", c.names.gl)  # remove the ".sorted" suffix
indNames(gl.alfl) <- n.names.gl  # assign the new names back to the genlight object
# c.names.gi <- indNames(gi.alfl)
# n.names.gi <- gsub("\\.sorted$", "", c.names.gi)  # remove the ".sorted" suffix
# indNames(gi.alfl) <- n.names.gi  # assign the new names back to the genlight object

# Identify merged samples
merged_alfl <- grep("_merged", n.names.gl, value = TRUE)
merged_alfl_IDs <- gsub("_merged$", "", merged_alfl)  # remove the '_merged' suffix

# Remove duplicated samples in d_meta and specify they have been merged
d_meta <- d_meta %>% 
  # Step 1: update the `ID_dra` column for the specified IDs
  mutate(ID_dra = ifelse(ID %in% merged_alfl_IDs, paste0(ID, "_merged"), ID_dra),
         Plate_ID = ifelse(ID %in% merged_alfl_IDs, "Merged", Plate_ID),
         Sent = ifelse(ID %in% merged_alfl_IDs, "Merged", Sent),
         Plate = ifelse(ID %in% merged_alfl_IDs, "Merged", Plate),
         Well = ifelse(ID %in% merged_alfl_IDs, "Merged", Well))

# Remove duplicate rows based on original ID: do it for each species one at a time
d_meta_alfl <- d_meta %>%
  filter(Species %in% "ALFL") %>% 
  group_by(ID_dra) %>%  # Use ID_dra here instead of ID because ID_dra has the duplicated names I need to remove
  slice(1) %>%
  ungroup() %>% 
  mutate(ID_dra = gsub("YT_", "YT", ID_dra))
# Two samples have mismatching names between the overall ALFL analysis and the neutral-SNP-only one:
# ALFL_YT_007 and ALFL_YT_012_seq2 --> ALFL_YT007 and ALFL_YT012_seq2
# Taken care here

# Extract raw info from vcf file
gt.tidy.alfl <- extract_gt_tidy(vcf.alfl, format_types = NULL)
gt.tidy.alfl <- gt.tidy.alfl %>% mutate(gt_DP = as.numeric(as.character(gt_DP)))  # convert gt_DP into numeric values
head(gt.tidy.alfl)

# Verify if any ID is present multiple times (do we have duplicated samples in the dataset?)
gt.tidy.alfl %>% count(Indiv, name = "N") %>% print(n = 100)  # all is good, same number of duplicate for each sample

# Add metadata to raw info of vcf file
gt.meta.tidy.alfl <- gt.tidy.alfl %>%
  mutate(Indiv = gsub(".sorted", "", Indiv)) %>% 
  mutate(Batch = ifelse(str_detect(Indiv, "_seq2"), 2, 1),
         ID = gsub("_seq2", "", Indiv),  # I want to maintain the old name as well since it's the key to append data to gl files
         Merged = ifelse(str_detect(ID, "_merged"), 1, 0),
         ID = str_replace(ID, "_merged", "")) %>%
  group_by(ID) %>% summarise(Indiv = first(Indiv),  # Maintain the original Indiv column
                             Nsnps = length(gt_GT[!is.na(gt_GT)]), 
                             DP = mean(gt_DP, na.rm=T),
                             Batch = max(Batch),  #  this is ok since we don't have duplicates - see L148
                             Merged = max(Merged)) %>%  # since it's binary it doesn't really matter if I keep the max or the min
  left_join(d_meta_alfl %>% select(3:18), by = c("Indiv"="ID_dra"))
head(gt.meta.tidy.alfl)
gt.meta.tidy.alfl %>% group_by(Region) %>% summarise(N = n())
# Region     N
# AB         8
# GU         2
# HO         1
# LA         3
# MT         1
# NBC       11
# NL         1
# SK         7
# YT         2

# Data frame with NA content per ID (see 4th graph below)
na.info.alfl <- data.frame(ID = indNames(gl.alfl),
                           NNA = count.ind.na.gl(gl.alfl))
na.info.alfl <- na.info.alfl %>% mutate(ID = gsub(".sorted","", ID))

# Data frame with Ho per ID
ho.alfl <- gl.report.heterozygosity(gl.alfl, method='ind')  # Verify what gl.report.heterozigosity does
ho.alfl <- ho.alfl %>% mutate(ind.name = gsub(".sorted","",ind.name))
colnames(ho.alfl)[2:4] <- c("Ho","f.hom.ref","f.hom.alt")

# Merge NA info on gt.meta dataset
gt.meta.tidy.alfl <- gt.meta.tidy.alfl %>% 
  left_join(na.info.alfl, by = c("Indiv" = "ID")) %>% 
  left_join(ho.alfl, by = c("Indiv" = "ind.name")) %>% 
  select(1:4,22:25,5:21)


### CLSW ------------------------------------------------------------------

# Populations: -R 0.5 --min-mac 3 --write-random-snp
# 101 samples
# Removed 2127559 loci that did not pass sample/population constraints from 2160011 loci.
# Kept 32452 loci, composed of 5705100 sites; 259623 of those sites were filtered, 21950 variant sites remained.
# Mean genotyped sites per locus: 162.13bp (stderr 0.27).

# Missingness
# --max-missing 0.5: After filtering, kept 21950 out of a possible 21950 Sites
# imiss --keep inds_to_keep_s50_i70.txt: After filtering, kept 60 out of 101 Individuals

# MaxDP
# --max-meanDP 1e+03: After filtering, kept 21946 out of a possible 21950 Sites

# Neutral - bayescan
# What the code does is actually removing loci that should have been removed with the normal filtering prior to ideintifying 
# loci under selection...
# After filtering, kept 9507 out of a possible 21946 Sites

#### Neutral SNPs ---------------------------------------------------------

# from CLSW_21_24_oct/5_DP_maxonly/
vcf.clsw <- vcfR::read.vcfR("./00_Data/03_filtering/03_neutral/CLSW/re_filtered_CLSW.recode.vcf")  # 9507 SNPs
head(vcf.clsw)

# Count how many RAD loci left from lmiss and imiss filters
SCAFFOLD.info <- vcf.clsw@fix %>% as.data.frame() %>%  # 
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, "_"), `[`,2), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )
length(unique(SCAFFOLD.info$RADloc))  # 9507 loci for 9507 SNPs - 1 SNP per locus, excellent

# Check names of vcf fields
vcf_field_names(vcf.clsw, tag = "FORMAT")

# Conversion to genlight and genind R formats - necessary to run glPCA from adegenet package
gl.clsw <- vcfR::vcfR2genlight(vcf.clsw)
# gi.data.clsw  <- vcfR::vcfR2genind(vcf.clsw) 

# Rename samples to uniformize suffix in genlight, genind and d_meta (also to find out which one have been merged...)
c.names.gl <- indNames(gl.clsw)
n.names.gl <- gsub("\\.sorted$", "", c.names.gl)  # remove the ".sorted" suffix
indNames(gl.clsw) <- n.names.gl  # assign the new names back to the genlight object

# Identify merged samples
merged_clsw <- grep("_merged", n.names.gl, value = TRUE)
merged_clsw_IDs <- gsub("_merged$", "", merged_clsw)  # remove the '_merged' suffix

# Remove duplicated samples in d_meta and specify they have been merged
d_meta <- d_meta %>% 
  # Step 1: update the `ID_dra` column for the specified IDs
  mutate(ID_dra = ifelse(ID %in% merged_clsw_IDs, paste0(ID, "_merged"), ID_dra),
         Plate_ID = ifelse(ID %in% merged_clsw_IDs, "Merged", Plate_ID),
         Sent = ifelse(ID %in% merged_clsw_IDs, "Merged", Sent),
         Plate = ifelse(ID %in% merged_clsw_IDs, "Merged", Plate),
         Well = ifelse(ID %in% merged_clsw_IDs, "Merged", Well))

# Remove duplicate rows based on original ID
d_meta_clsw <- d_meta %>%
  filter(Species %in% "CLSW") %>% 
  group_by(ID_dra) %>%  # Use ID_dra here instead of ID because ID_dra has the duplicated names I need to remove
  slice(1) %>%
  ungroup()

# Extract raw info from vcf file
gt.tidy.clsw <- extract_gt_tidy(vcf.clsw, format_types = NULL)
gt.tidy.clsw <- gt.tidy.clsw %>% mutate(gt_DP = as.numeric(as.character(gt_DP)))  # convert gt_DP into numeric values
head(gt.tidy.clsw)

# Verify if any ID is present multiple times (do we have duplicated samples in the dataset?)
gt.tidy.clsw %>% count(Indiv, name = "N") %>% print(n = 100)  # all is good, same number of duplicate for each sample

# Add metadata to raw info of vcf file
gt.meta.tidy.clsw <- gt.tidy.clsw %>%
  mutate(Batch = ifelse(str_detect(Indiv, "_seq2"), 2, 1),
         ID = gsub("_seq2", "", Indiv),  # I want to maintain the old name as well since it's the key to append data to gl files
         Merged = ifelse(str_detect(ID, "_merged"), 1, 0),
         ID = str_replace(ID, "_merged", "")) %>% 
  group_by(ID) %>% summarise(Indiv = first(Indiv),  # Maintain the original Indiv column
                             Nsnps = length(gt_GT[!is.na(gt_GT)]), 
                             DP = mean(gt_DP, na.rm=T),
                             Batch = max(Batch),  #  this is ok since we don't have duplicates - see L148
                             Merged = max(Merged)) %>%  # since it's binary it doesn't really matter if I keep the max or the min
  left_join(d_meta_clsw %>% select(3:18), by = c("Indiv"="ID_dra"))
head(gt.meta.tidy.clsw)
gt.meta.tidy.clsw %>% group_by(Region) %>% summarise(N = n())
# Region     N
# AZ         5
# BC         5
# CO         8
# MB         6
# MS         4
# MX         5
# ON         7
# SK        12
# WA         8

# Data frame with NA content per ID (see 4th graph below)
na.info.clsw <- data.frame(ID = indNames(gl.clsw),
                           NNA = count.ind.na.gl(gl.clsw))

# Data frame with Ho per ID
ho.clsw <- gl.report.heterozygosity(gl.clsw, method='ind')  # Verify what gl.report.heterozigosity does
colnames(ho.clsw)[2:4] <- c("Ho","f.hom.ref","f.hom.alt")

# Merge NA info on gt.meta dataset
gt.meta.tidy.clsw <- gt.meta.tidy.clsw %>% 
  left_join(na.info.clsw, by = c("Indiv" = "ID")) %>% 
  left_join(ho.clsw, by = c("Indiv" = "ind.name")) %>% 
  select(1:4,22:25,5:21)


### PUMA ------------------------------------------------------------------

# Populations: -R 0.5 --min-mac 3 --write-random-snp
# 50 samples
# Removed 1939511 loci that did not pass sample/population constraints from 2141609 loci.
# Kept 202098 loci, composed of 40131018 sites; 236450 of those sites were filtered, 108698 variant sites remained.
# Mean genotyped sites per locus: 165.58bp (stderr 0.10).

# Missingness
# --max-missing 0.3: After filtering, kept 108698 out of a possible 108698 Sites
# imiss --keep inds_to_keep_7070.txt: After filtering, kept 43 out of 50 Individuals

# MaxDP
# --max-meanDP 1.6e+03: After filtering, kept 108694 out of a possible 108698 Sites

# Neutral - bayescan
# What the code does is actually removing loci that should have been removed with the normal filtering prior to ideintifying 
# loci under selection...
# After filtering, kept 40927 out of a possible 108694 Sites

#### Neutral SNPs ---------------------------------------------------------

# From PUMA_ALL/PUMA_allsamples_preanalysis/8_neutral/
vcf.puma <- vcfR::read.vcfR("./00_Data/03_filtering/03_neutral/PUMA/re_filtered_PUMA_7070.recode.vcf")  # 40927 SNPs
head(vcf.puma)

# Count how many RAD loci left from lmiss and imiss filters
SCAFFOLD.info <- vcf.puma@fix %>% as.data.frame() %>%  # 
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, "_"), `[`,2), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )
length(unique(SCAFFOLD.info$RADloc))  # 40927 loci for 40927 SNPs - 1 SNP per locus, excellent

# Check names of vcf fields
vcf_field_names(vcf.puma, tag = "FORMAT")

# Conversion to genlight and genind R formats - necessary to run glPCA from adegenet package
gl.puma <- vcfR::vcfR2genlight(vcf.puma)
# gi.data.puma  <- vcfR::vcfR2genind(vcf.puma) 

# Rename samples to uniformize suffix in genlight, genind and d_meta (also to find out which one have been merged...)
c.names.gl <- indNames(gl.puma)
n.names.gl <- gsub("\\.sorted$", "", c.names.gl)  # remove the ".sorted" suffix
indNames(gl.puma) <- n.names.gl  # assign the new names back to the genlight object
# c.names.gi <- indNames(gi.puma)
# n.names.gi <- gsub("\\.sorted$", "", c.names.gi)  # remove the ".sorted" suffix
# indNames(gi.puma) <- n.names.gi  # assign the new names back to the genlight object

# Identify merged samples
merged_puma <- grep("_merged", n.names.gl, value = TRUE)
merged_puma_IDs <- gsub("_merged$", "", merged_puma)  # remove the '_merged' suffix

# Remove duplicated samples in d_meta and specify they have been merged
d_meta <- d_meta %>% 
  # Step 1: update the `ID_dra` column for the specified IDs
  mutate(ID_dra = ifelse(ID %in% merged_puma_IDs, paste0(ID, "_merged"), ID_dra),
         Plate_ID = ifelse(ID %in% merged_puma_IDs, "Merged", Plate_ID),
         Sent = ifelse(ID %in% merged_puma_IDs, "Merged", Sent),
         Plate = ifelse(ID %in% merged_puma_IDs, "Merged", Plate),
         Well = ifelse(ID %in% merged_puma_IDs, "Merged", Well))

# Remove duplicate rows based on original ID
d_meta_puma <- d_meta %>%
  filter(Species %in% "PUMA") %>% 
  group_by(ID_dra) %>%  # Use ID_dra here instead of ID because ID_dra has the duplicated names I need to remove
  slice(1) %>%
  ungroup()

# Extract raw info from vcf file
gt.tidy.puma <- extract_gt_tidy(vcf.puma, format_types = NULL)
gt.tidy.puma <- gt.tidy.puma %>% mutate(gt_DP = as.numeric(as.character(gt_DP)))  # convert gt_DP into numeric values
head(gt.tidy.puma)

# Verify if any ID is present multiple times (do we have duplicated samples in the dataset?)
gt.tidy.puma %>% count(Indiv, name = "N") %>% print(n = 100)  # all is good, same number of duplicate for each sample

# Add metadata to raw info of vcf file
gt.meta.tidy.puma <- gt.tidy.puma %>%
  mutate(Indiv = gsub(".sorted", "", Indiv)) %>% 
  mutate(Batch = ifelse(str_detect(Indiv, "_seq2"), 2, 1),
         ID = gsub("_seq2", "", Indiv),  # I want to maintain the old name as well since it's the key to append data to gl files
         Merged = ifelse(str_detect(ID, "_merged"), 1, 0),
         ID = str_replace(ID, "_merged", "")) %>% 
  group_by(ID) %>% summarise(Indiv = first(Indiv),  # Maintain the original Indiv column
                             Nsnps = length(gt_GT[!is.na(gt_GT)]), 
                             DP = mean(gt_DP, na.rm=T),
                             Batch = max(Batch),  #  this is ok since we don't have duplicates - see L148
                             Merged = max(Merged)) %>%  # since it's binary it doesn't really matter if I keep the max or the min
  left_join(d_meta_puma %>% select(3:18), by = c("Indiv"="ID_dra"))
head(gt.meta.tidy.puma)
gt.meta.tidy.puma %>% group_by(Region) %>% summarise(N = n())
# Region     N
# NC         4
# NV         2
# PI_BC      8
# SBC        1
# SK        12
# SWON      10
# VI         1
# WA         5

# Data frame with NA content per ID (see 4th graph below)
na.info.puma <- data.frame(ID = indNames(gl.puma),
                           NNA = count.ind.na.gl(gl.puma))
na.info.puma <- na.info.puma %>% mutate(ID = gsub(".sorted","", ID))

# Data frame with Ho per ID
ho.puma <- gl.report.heterozygosity(gl.puma, method='ind')  # Verify what gl.report.heterozigosity does
ho.puma <- ho.puma %>% mutate(ind.name = gsub(".sorted","",ind.name))
colnames(ho.puma)[2:4] <- c("Ho","f.hom.ref","f.hom.alt")

# Merge NA info on gt.meta dataset
gt.meta.tidy.puma <- gt.meta.tidy.puma %>% 
  left_join(na.info.puma, by = c("Indiv" = "ID")) %>% 
  left_join(ho.puma, by = c("Indiv" = "ind.name")) %>% 
  select(1:4,22:25,5:21)


### VGSW ------------------------------------------------------------------

# Populations: -R 0.5 --min-mac 3 --write-random-snp
# 51 samples
# populations.log file incomplete (110197 variant sites)

# Missingness
# --max-missing 0.6: After filtering, kept 110197 out of a possible 110197 Sites
# imiss --keep inds_to_keep_s60_i70.txt: After filtering, kept 41 out of 51 Individuals

# MaxDP
# --max-meanDP 900: After filtering, kept 110189 out of a possible 110197 Sites

# Neutral - bayescan
# What the code does is actually removing loci that should have been removed with the normal filtering prior to ideintifying 
# loci under selection...
# After filtering, kept 44211 out of a possible 110189 Sites

#### Neutral SNPs ---------------------------------------------------------

# From VGSW_21_24_oct/8_neutral/neutral/
vcf.vgsw <- vcfR::read.vcfR("./00_Data/03_filtering/03_neutral/VGSW/vgsw_s60_i70_neutral.vcf.recode.vcf")  # 44211 SNPs
head(vcf.vgsw)

# Count how many RAD loci left from lmiss and imiss filters
SCAFFOLD.info <- vcf.vgsw@fix %>% as.data.frame() %>%  # 
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, "_"), `[`,2), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )
length(unique(SCAFFOLD.info$RADloc))  # 44211 loci for 44211 SNPs - 1 SNP per locus, excellent

# Check names of vcf fields
vcf_field_names(vcf.vgsw, tag = "FORMAT")

# Conversion to genlight and genind R formats - necessary to run glPCA from adegenet package
gl.vgsw <- vcfR::vcfR2genlight(vcf.vgsw)
# gi.data.vgsw  <- vcfR::vcfR2genind(vcf.vgsw) 

# Rename samples to uniformize suffix in genlight, genind and d_meta (also to find out which one have been merged...)
c.names.gl <- indNames(gl.vgsw)
n.names.gl <- gsub("\\.sorted$", "", c.names.gl)  # remove the ".sorted" suffix
indNames(gl.vgsw) <- n.names.gl  # assign the new names back to the genlight object
# c.names.gi <- indNames(gi.vgsw)
# n.names.gi <- gsub("\\.sorted$", "", c.names.gi)  # remove the ".sorted" suffix
# indNames(gi.vgsw) <- n.names.gi  # assign the new names back to the genlight object

# Identify merged samples
merged_vgsw <- grep("_merged", n.names.gl, value = TRUE)
merged_vgsw_IDs <- gsub("_merged$", "", merged_vgsw)  # remove the '_merged' suffix

# Remove duplicated samples in d_meta and specify they have been merged
d_meta <- d_meta %>% 
  # Step 1: update the `ID_dra` column for the specified IDs
  mutate(ID_dra = ifelse(ID %in% merged_vgsw_IDs, paste0(ID, "_merged"), ID_dra),
         Plate_ID = ifelse(ID %in% merged_vgsw_IDs, "Merged", Plate_ID),
         Sent = ifelse(ID %in% merged_vgsw_IDs, "Merged", Sent),
         Plate = ifelse(ID %in% merged_vgsw_IDs, "Merged", Plate),
         Well = ifelse(ID %in% merged_vgsw_IDs, "Merged", Well))

# Remove duplicate rows based on original ID
d_meta_vgsw <- d_meta %>%
  filter(Species %in% "VGSW") %>% 
  group_by(ID_dra) %>%  # Use ID_dra here instead of ID because ID_dra has the duplicated names I need to remove
  slice(1) %>%
  ungroup()

# Extract raw info from vcf file
gt.tidy.vgsw <- extract_gt_tidy(vcf.vgsw, format_types = NULL)
gt.tidy.vgsw <- gt.tidy.vgsw %>% mutate(gt_DP = as.numeric(as.character(gt_DP)))  # convert gt_DP into numeric values
head(gt.tidy.vgsw)

# Verify if any ID is present multiple times (do we have duplicated samples in the dataset?)
gt.tidy.vgsw %>% count(Indiv, name = "N") %>% print(n = 100)  # all is good, same number of duplicate for each sample

# Add metadata to raw info of vcf file
gt.meta.tidy.vgsw <- gt.tidy.vgsw %>%
  mutate(Batch = ifelse(str_detect(Indiv, "_seq2"), 2, 1),
         ID = gsub("_seq2", "", Indiv),  # I want to maintain the old name as well since it's the key to append data to gl files
         Merged = ifelse(str_detect(ID, "_merged"), 1, 0),
         ID = str_replace(ID, "_merged", "")) %>% 
  group_by(ID) %>% summarise(Indiv = first(Indiv),  # Maintain the original Indiv column
                             Nsnps = length(gt_GT[!is.na(gt_GT)]), 
                             DP = mean(gt_DP, na.rm=T),
                             Batch = max(Batch),  #  this is ok since we don't have duplicates - see L148
                             Merged = max(Merged)) %>%  # since it's binary it doesn't really matter if I keep the max or the min
  left_join(d_meta_vgsw %>% select(3:18), by = c("Indiv"="ID_dra"))
head(gt.meta.tidy.vgsw)
gt.meta.tidy.vgsw %>% group_by(Region) %>% summarise(N = n())
# Region     N
# CO         8
# NM         2
# NV        10
# PI_BC      7
# SBC        5
# WA         9

# Data frame with NA content per ID (see 4th graph below)
na.info.vgsw <- data.frame(ID = indNames(gl.vgsw),
                           NNA = count.ind.na.gl(gl.vgsw))

# Data frame with Ho per ID
ho.vgsw <- gl.report.heterozygosity(gl.vgsw, method='ind')  # Verify what gl.report.heterozigosity does
colnames(ho.vgsw)[2:4] <- c("Ho","f.hom.ref","f.hom.alt")

# Merge NA info on gt.meta dataset
gt.meta.tidy.vgsw <- gt.meta.tidy.vgsw %>% 
  left_join(na.info.vgsw, by = c("Indiv" = "ID")) %>% 
  left_join(ho.vgsw, by = c("Indiv" = "ind.name")) %>% 
  select(1:4,22:25,5:21)


### Merge species-specific gt.meta files together -------------------------

gt.meta.tidy <- bind_rows(gt.meta.tidy.alfl,
                          gt.meta.tidy.clsw,
                          gt.meta.tidy.puma,
                          gt.meta.tidy.vgsw)
gt.meta.tidy %>% group_by(Region, Species) %>% summarise(N= n()) %>% print(n = 50)
# write.csv(gt.meta.tidy, file = "./02_Results/05_filtering/01_reference/05_summary_stats/AI_meta_filtering_tidy.csv", row.names = F)
# gt.meta.tidy <- read.csv("./02_Results/05_filtering/01_reference/05_summary_stats/AI_meta_filtering_tidy.csv", stringsAsFactors = F)

# Plots looking at SNPs and coverage
gSNP.by.region <- gt.meta.tidy %>% ggplot(aes(x = Region, y = Nsnps)) +
  # geom_point(alpha = 3/5, size = 3, aes(col = Region)) +
  geom_jitter(size = 3, height = 0, alpha = 3/5, aes(col = Region)) +
  geom_boxplot(alpha = 0) +
  labs(y = "N snps", x = NULL) +
  facet_wrap(~ Species, scales = "free_y") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
gSNP.by.region
# ggsave(filename = "./02_Results/05_filtering/01_reference/05_summary_stats/Summarise_SNPs.png",
#        plot = gSNP.by.region,
#        width = 20, height = 12, units = "in")

gNA.by.region <- gt.meta.tidy %>% ggplot(aes(x = Region, y = NNA)) +
  # geom_point(alpha = 3/5, size = 3, aes(col = Region)) +
  geom_jitter(size = 3, height = 0, alpha = 3/5, aes(col = Region)) +
  geom_boxplot(alpha = 0) +
  labs(y = "Proportion of missing data", x = NULL) +
  facet_wrap(~ Species, scales = "free_y") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
gNA.by.region
# ggsave(filename = "./02_Results/05_filtering/01_reference/05_summary_stats/Summarise_NNAs.png",
#        plot = gNA.by.region,
#        width = 20, height = 12, units = "in")

gSNP.by.cov <- gt.meta.tidy %>% ggplot(aes(x = DP, y = Nsnps, col = Region)) +
  geom_point(size = 3, alpha = 0.6)+
  labs(y = "N snps", x = "Coverage") +
  facet_wrap(~ Species, scales = "free_y") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
gSNP.by.cov
# ggsave(filename = "./02_Results/05_filtering/01_reference/05_summary_stats/Summarise_SNPsbyDP.png",
#        plot = gSNP.by.cov,
#        width = 20, height = 12, units = "in")

gNA.by.cov <- gt.meta.tidy %>% ggplot(aes(x = DP, y = NNA, col = Species)) +
  geom_point(size = 3, alpha = 0.7)+
  labs(y = "Proportion of missing data", x = "Coverage") +
  scale_x_continuous(breaks = seq(0, max(gt.meta.tidy$DP, na.rm = TRUE), by = 10)) +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
gNA.by.cov
# ggsave(filename = "./02_Results/05_filtering/01_reference/05_summary_stats/Summarise_NAbyDP.png",
#        plot = gNA.by.cov,
#        width = 20, height = 12, units = "in")

gHo.by.cov <- gt.meta.tidy %>% ggplot(aes(x = DP, y = Ho, col = Species)) +
  geom_point(size = 3, alpha = 0.7)+
  labs(y = "Ho", x = "Coverage") +
  scale_x_continuous(breaks = seq(0, max(gt.meta.tidy$DP, na.rm = TRUE), by = 10)) +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
gHo.by.cov
# ggsave(filename = "./02_Results/05_filtering/01_reference/05_summary_stats/Summarise_HobyDP.png",
#        plot = gHo.by.cov,
#        width = 20, height = 12, units = "in")

gHo.by.NNA <- gt.meta.tidy %>% ggplot(aes(x = NNA, y = Ho, col = Species)) +
  geom_point(size = 3, alpha = 0.7)+
  labs(y = "Ho", x = "Proportion of missing data") +
  scale_x_continuous(breaks = seq(0, max(gt.meta.tidy$NNA, na.rm = TRUE), by = 0.10)) +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
gHo.by.NNA
# ggsave(filename = "./02_Results/05_filtering/01_reference/05_summary_stats/Summarise_HobyNNA.png",
#        plot = gHo.by.NNA,
#        width = 20, height = 12, units = "in")

gHo.by.region <- gt.meta.tidy %>% ggplot(aes(x = Region, y = Ho)) +
  geom_jitter(size = 3, height = 0, alpha = 3/5, aes(col = DP)) +
  geom_boxplot(alpha = 0) +
  scale_colour_viridis_c(option = "turbo") +
  labs(y = "Ho", x = NULL) +
  facet_wrap(~ Species, scales = "free_y") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
gHo.by.region
# ggsave(filename = "./02_Results/05_filtering/01_reference/05_summary_stats/Summarise_HobyDP_perRegion.png",
#        plot = gHo.by.region,
#        width = 20, height = 12, units = "in")



# PCA ---------------------------------------------------------------------

# Make a list with gl files to run PCA with lapply
# gl.list <- list(gl.alfl, gl.clsw, gl.puma, gl.vgsw)
# names(gl.list) <- c("ALFL","CLSW","PUMA","VGSW")
# 
# # Run PCA
# pca.list <- setNames(lapply(seq_along(gl.list), function(i) {
#   gl.file <- gl.list[[i]]
#   gl.file.name <- names(gl.list)[i]  # Get the name of the current element
#   print(paste0("Estimating PCA for ", gl.file.name))
# 
#   # Perform PCA
#   pca.res <- glPca(gl.file, center = TRUE, scale = FALSE, parallel = TRUE, n.core = 8, nf = 1000)
#   print(paste0("PCA estimated for ", gl.file.name))
# 
#   return(pca.res)  # Store result in the list
# }), names(gl.list))  # Assign names to pca.list based on gl.list
# 
# # Extract PCA eigen score
# pca.mat.list <- lapply(pca.list, function(y){
#   pca.mat <- data.frame(y$scores)  # score matrices
#   print(y)
#   pca.mat
# })
# 
# save(list = c("pca.list", "pca.mat.list"), file = "./02_Results/06_PCA/PCA_alfl_clsw_puma_vgsw.Rdata")
load("./02_Results/06_PCA/PCA_alfl_clsw_puma_vgsw.Rdata")


## Eig var ----------------------------------------------------------------

pca.var.alfl <- pca_var(pca.list[["ALFL"]], nInd(gl.alfl)-1) %>%
  mutate(Species = "ALFL")
pca.var.clsw <- pca_var(pca.list[["CLSW"]], nInd(gl.clsw)-1) %>%
  mutate(Species = "CLSW")
pca.var.puma <- pca_var(pca.list[["PUMA"]], nInd(gl.puma)-1) %>%
  mutate(Species = "PUMA")
pca.var.vgsw <- pca_var(pca.list[["VGSW"]], nInd(gl.vgsw)-1) %>%
  mutate(Species = "VGSW")
# pca.var.puma <- pca_var(pca.list[["PUMA"]], length(pca.list[["PUMA"]]$eig)) %>%  # originally using nInd(gl.clsw)-1, issue here in the number of eig estimated...
#   mutate(Species = "PUMA")
# pca.var.vgsw <- pca_var(pca.list[["VGSW"]], length(pca.list[["VGSW"]]$eig)) %>%  # originally using nInd(gl.clsw)-1, issue here in the number of eig estimated...
#   mutate(Species = "VGSW")
pca.var <- bind_rows(pca.var.alfl,
                     pca.var.clsw,
                     pca.var.puma,
                     pca.var.vgsw)

gPCA.var <- pca.var %>%
  # filter(SNPs %in% "all") %>% 
  ggplot(aes(x = axis, y = p.eig * 100, fill = axis)) +
  plot_pca_eig +
  # scale_y_continuous(limits=c(0,1.2),breaks=c(0,0.3,0.6,0.9,1.2)) +
  # geom_hline(yintercept = 0, col = "grey20") +
  facet_wrap(~ Species, scales = "free") +
  guides(fill = "none") +
  theme_bw(base_size = 20) +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
gPCA.var
# ggsave(filename = "./02_Results/06_PCA/PCA.var.species.png",
#        plot = gPCA.var,
#        width = 15, height = 12, units = "in")


## PCA figures ------------------------------------------------------------

### Data ------------------------------------------------------------------

pca.alfl <- pca.list[["ALFL"]] %>% QuickPop::pca_scoretable(naxe = 6) %>%
  mutate(ID = gsub(".sorted", "", ID),
         SNPs = "all") %>%
  left_join(gt.meta.tidy, by = c("ID" = "Indiv")) %>%
  droplevels()  # to remove unused Pop levels

pca.clsw <- pca.list[["CLSW"]] %>% QuickPop::pca_scoretable(naxe = 6) %>%
  mutate(SNPs = "all") %>% 
  left_join(gt.meta.tidy, by = c("ID" = "Indiv")) %>%
  droplevels()  # to remove unused Pop levels

pca.puma <- pca.list[["PUMA"]] %>% QuickPop::pca_scoretable(naxe = 6) %>%
  mutate(ID = gsub(".sorted", "", ID),
         SNPs = "all") %>%
  left_join(gt.meta.tidy, by = c("ID" = "Indiv")) %>%
  droplevels()  # to remove unused Pop levels

pca.vgsw <- pca.list[["VGSW"]] %>% QuickPop::pca_scoretable(naxe = 6) %>%
  mutate(SNPs = "all") %>% 
  left_join(gt.meta.tidy, by = c("ID" = "Indiv")) %>%
  droplevels()  # to remove unused Pop levels

d.pca <- bind_rows(pca.alfl,
                   pca.clsw,
                   pca.puma,
                   pca.vgsw)

# Add Season defined by Thilini
d.pca.tt <- read.csv("./00_Data/00_fileinfo/PCA.data_TT_March_25.csv")
# which(d.pca$ID %nin% d.pca.tt$ID)
# which(d.pca.tt$ID %nin% d.pca$ID)
table(duplicated(d.pca.tt$ID))  # F 190 T 84 --> PUMA samples duplicated, removing the duplicated ones
d.pca.tt <- d.pca.tt %>% 
  group_by(ID) %>%
  slice(1) %>%
  ungroup()

d.pca <- d.pca %>% 
  left_join(d.pca.tt[, c("ID","Season..TT.")])
colnames(d.pca)[33] <- "Season.TT"
d.pca %>% group_by(Season.TT) %>% summarise(N = n())
# Season.TT     N
# Breeding    140
# Migrating    40
# write.csv(d.pca, file = "./02_Results/06_PCA/PCA.data.csv", row.names = F)


## Figures ---------------------------------------------------------------

d.pca <- d.pca %>%
  mutate(Region = factor(Region, levels = c("YT","NBC","BC","SBC","VI","PI_BC","WA","MT","NV","AZ",  
                                            "AB","MB","SK","CO","NM","MX","GU","HO",  
                                            "ON","SWON","NL","NC","MS","LA")),
         Season = factor(Season.TT, levels = c("Breeding","Migrating")))

region_colors <- c("YT" = "#A2E8F0", "NBC" = "#0000FF", "BC" = "#0000FF", "SBC" = "#0000FF", 
                   "VI" = "#0086FF", "PI_BC" = "#0086FF", "WA" = "#8AB3FF", "MT" = "#1380B6",
                   "NV" = "#00327E", "AZ" = "#18579A", "AB" = "#E32DD7", "MB" = "#CD6EF7",
                   "SK" = "#942FEF", "CO" = "#8556AB", "NM" = "#D998E8", "MX" = "#ff84ff",
                   "GU" = "#D998E8", "HO" = "#8B4CC1", "ON" = "#C40000", "SWON" = "#C40000",
                   "NL" = "#FF2D15", "NC" = "#F96D6D", "MS" = "#8B0000", "LA" = "#FF7C00")
season_shapes <- c(19,17)

# Create separate plots
l_gPCA.1v2 <- list()
l_species <- unique(d.pca$Species)

for(sp in l_species){
  species_data <- d.pca %>% filter(Species == sp)
  pcavar <- pca.var %>% filter(Species == sp)
  
  # Subset colors for only the regions present in this species
  species_colors <- region_colors[names(region_colors) %in% unique(species_data$Region)]

  l_gPCA.1v2[[sp]] <- ggplot(species_data, aes(x = score.PC1, y = score.PC2, col = Region, shape = Season.TT)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(stroke = 1, size = 6, alpha = 0.5) +
    scale_color_manual(values = species_colors, guide = guide_legend(ncol = 1), name = "Region",) +
    scale_shape_manual(values = season_shapes, name = "Season",) +
    guides(color = guide_legend(order = 1),  # Fix legend order 
           shape = guide_legend(order = 2)) +
    labs(x = paste0("PC1 (", round(pcavar$p.eig[pcavar$axis %in% 1 & pcavar$Species %in% sp], 4)*100, "%)"), 
         y = paste0("PC2 (", round(pcavar$p.eig[pcavar$axis %in% 2 & pcavar$Species %in% sp], 4)*100, "%)"), 
         title = sp) +
    theme_bw(base_size = 20) +
    theme(axis.text = element_text(size = 18, colour = "black"),
          axis.title.y.right = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
}

# Combine plots with patchwork
gPCA.1v2 <- wrap_plots(l_gPCA.1v2) 
gPCA.1v2
# ggsave(filename = "./02_Results/06_PCA/PCA.1v2.species.png",
#        plot = gPCA.1v2,
#        width = 20, height = 15, units = "in")


l_gPCA.3v4 <- list()
for(sp in l_species){
  species_data <- d.pca %>% filter(Species == sp)
  pcavar <- pca.var %>% filter(Species == sp)
  
  # Subset colors for only the regions present in this species
  species_colors <- region_colors[names(region_colors) %in% unique(species_data$Region)]
  
  l_gPCA.3v4[[sp]] <- ggplot(species_data, aes(x = score.PC3, y = score.PC4, col = Region, shape = Season.TT)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(stroke = 1, size = 6, alpha = 0.5) +
    scale_color_manual(values = species_colors, guide = guide_legend(ncol = 1), name = "Region",) +
    scale_shape_manual(values = season_shapes, name = "Season",) +
    guides(color = guide_legend(order = 1),  # Fix legend order 
           shape = guide_legend(order = 2)) +
    labs(x = paste0("PC3 (", round(pcavar$p.eig[pcavar$axis %in% 3 & pcavar$Species %in% sp], 4)*100, "%)"), 
         y = paste0("PC4 (", round(pcavar$p.eig[pcavar$axis %in% 4 & pcavar$Species %in% sp], 4)*100, "%)"), 
         title = sp) +
    theme_bw(base_size = 20) +
    theme(axis.text = element_text(size = 18, colour = "black"),
          axis.title.y.right = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
}

# Combine plots with patchwork
gPCA.3v4 <- wrap_plots(l_gPCA.3v4) 
gPCA.3v4
# ggsave(filename = "./02_Results/06_PCA/PCA.3v4.species.png",
#        plot = gPCA.3v4,
#        width = 20, height = 15, units = "in")
