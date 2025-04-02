# Info --------------------------------------------------------------------
# 
# Author: Luca Montana
# Affiliation: University of Lethbridge
# Group: T. Burg research group
# Location: Rimouski, QC
# Date: 2025-03-12
# 
# Overview: Identification of putative neutral SNPs with BayeScan for
#           Alder flycatchers, cliff swallows, purple martins, violet-green swallows
#

# Library -----------------------------------------------------------------

library(tidyverse)
library(vcfR)
library(adegenet)
library(coda)  #  evaluation convergence MCMC
library(OutFLANK)
library(bigstatsr)
library(bigsnpr)   # package for SNP trimming


# Internal functions
# not in
`%nin%` <- Negate(`%in%`)

# From mfoll/BayeScan github page
plot_bayescan<-function(res,FDR=0.05,size=1,pos=0.35,highlight=NULL,name_highlighted=F,add_text=T)
{
  if (is.character(res))
    res=read.table(res)
  
  colfstat=5
  colq=colfstat-2
  
  highlight_rows=which(is.element(as.numeric(row.names(res)),highlight))
  non_highlight_rows=setdiff(1:nrow(res),highlight_rows)
  
  outliers=as.integer(row.names(res[res[,colq]<=FDR,]))
  
  ok_outliers=TRUE
  if (sum(res[,colq]<=FDR)==0)
    ok_outliers=FALSE;
  
  res[res[,colq]<=0.0001,colq]=0.0001
  
  # plot
  plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),xlab="log10(q value)",ylab=names(res[colfstat]),type="n")
  points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,colfstat],pch=19,cex=size)
  
  if (name_highlighted) {
    if (length(highlight_rows)>0) {
      text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
    }
  }
  else {
    points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],col="red",pch=19,cex=size)
    # add names of loci over p and vertical line
    if (ok_outliers & add_text) {
      text(log10(res[res[,colq]<=FDR,][,colq])+pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),cex=size)
    }
  }
  lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=2)
  
  return(list("outliers"=outliers,"nb_outliers"=length(outliers)))
}



# Data --------------------------------------------------------------------

## Upload metadata --------------------------------------------------------

d_meta <- read.csv("./00_Data/00_fileinfo/AI_metadata.csv", stringsAsFactors = F)


### VCF files -------------------------------------------------------------

### ALFL ------------------------------------------------------------------

# from ALFL_ALL/ALFL_ALL_merged/01_reference/00_alfl/alfl_population_final/
vcf.alfl <- vcfR::read.vcfR("./00_Data/03_filtering/02_DPmax/ALFL/alfl_7090_maxDP.recode.vcf")  # 46992 SNPs
head(vcf.alfl)

# Count how many RAD loci left from lmiss and imiss filters
SCAFFOLD.info <- vcf.alfl@fix %>% as.data.frame() %>%  # 
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, "_"), `[`,2), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )
length(unique(SCAFFOLD.info$RADloc))  # 46992 loci for 46992 SNPs - 1 SNP per locus, excellent

# Conversion to genlight and genind R formats - necessary to run glPCA from adegenet package
gl.alfl <- vcfR::vcfR2genlight(vcf.alfl)

# Rename samples to uniformize suffix in genlight, genind and d_meta (also to find out which one have been merged...)
c.names.gl <- indNames(gl.alfl)
n.names.gl <- gsub("\\.sorted$", "", c.names.gl)  # remove the ".sorted" suffix
indNames(gl.alfl) <- n.names.gl  # assign the new names back to the genlight object

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
  ungroup()

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
                             Nsnps.all = length(gt_GT[!is.na(gt_GT)]), 
                             DP.all = mean(gt_DP, na.rm=T),
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

# Two samples have mismatching names between the overall ALFL analysis and the neutral-SNP-only one:
# ALFL_YT_007 and ALFL_YT_012_seq2 --> ALFL_YT007 and ALFL_YT012_seq2
# Take care of this in next section


### CLSW ------------------------------------------------------------------

# from CLSW_21_24_oct/5_DP_maxonly/
vcf.clsw <- vcfR::read.vcfR("./00_Data/03_filtering/02_DPmax/CLSW/clsw.s50_i70.21946snps.60ids.recode.vcf")  # 21946 SNPs
head(vcf.clsw)

# Count how many RAD loci left from lmiss and imiss filters
SCAFFOLD.info <- vcf.clsw@fix %>% as.data.frame() %>%  # 
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, "_"), `[`,2), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )
length(unique(SCAFFOLD.info$RADloc))  # 21946 loci for 21946 SNPs - 1 SNP per locus, excellent

# Conversion to genlight and genind R formats - necessary to run glPCA from adegenet package
gl.clsw <- vcfR::vcfR2genlight(vcf.clsw)

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
                             Nsnps.all = length(gt_GT[!is.na(gt_GT)]), 
                             DP.all = mean(gt_DP, na.rm=T),
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


### PUMA ------------------------------------------------------------------

# from PUMA_ALL/PUMA_allsamples_preanalysis/5_DP_maxonly/
vcf.puma <- vcfR::read.vcfR("./00_Data/03_filtering/02_DPmax/PUMA/puma_7070_maxDP.recode.vcf")  # 108694 SNPs
head(vcf.puma)

# Count how many RAD loci left from lmiss and imiss filters
SCAFFOLD.info <- vcf.puma@fix %>% as.data.frame() %>%  # 
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, "_"), `[`,2), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )
length(unique(SCAFFOLD.info$RADloc))  # 108694 loci for 108694 SNPs - 1 SNP per locus, excellent

# Conversion to genlight and genind R formats - necessary to run glPCA from adegenet package
gl.puma <- vcfR::vcfR2genlight(vcf.puma)

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
                             Nsnps.all = length(gt_GT[!is.na(gt_GT)]), 
                             DP.all = mean(gt_DP, na.rm=T),
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


### VGSW ------------------------------------------------------------------

# from VGSW_21_24_oct/5_DP_maxonly/
vcf.vgsw <- vcfR::read.vcfR("./00_Data/03_filtering/02_DPmax/VGSW/vgsw.s60_i70.110189snps.41ids.recode.vcf")  # 110189 SNPs
head(vcf.vgsw)

# Count how many RAD loci left from lmiss and imiss filters
SCAFFOLD.info <- vcf.vgsw@fix %>% as.data.frame() %>%  # 
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, "_"), `[`,2), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )
length(unique(SCAFFOLD.info$RADloc))  # 110189 loci for 110189 SNPs - 1 SNP per locus, excellent

# Conversion to genlight and genind R formats - necessary to run glPCA from adegenet package
gl.vgsw <- vcfR::vcfR2genlight(vcf.vgsw)

# Rename samples to uniformize suffix in genlight, genind and d_meta (also to find out which one have been merged...)
c.names.gl <- indNames(gl.vgsw)
n.names.gl <- gsub("\\.sorted$", "", c.names.gl)  # remove the ".sorted" suffix
indNames(gl.vgsw) <- n.names.gl  # assign the new names back to the genlight object

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
                             Nsnps.all = length(gt_GT[!is.na(gt_GT)]), 
                             DP.all = mean(gt_DP, na.rm=T),
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



# Bayescan ----------------------------------------------------------------

## ALFL -------------------------------------------------------------------

loc.alfl <- data.frame(vcf.alfl@fix[, c("CHROM","POS","ID")])

out.bayescan.alfl <- read.table("./00_Data/03_filtering/03_neutral/ALFL/ALFL_outlier_fst.txt",
                                header = T)

# Evaluate MCMC chain convergence
# https://evomics.org/learning/population-and-speciation-genomics/2016-population-and-speciation-genomics/bayescan-exercise/
chain.alfl <- read.table("./00_Data/03_filtering/03_neutral/ALFL/ALFL_outlier.sel",
                         header = T)
chain.alfl <- mcmc(chain.alfl, thin = 10)
# plot(chain.alfl)  # check trace plots and posterior distribution
summary(chain.alfl)
autocorr.diag(chain.alfl)  # check correlation between sampled parameter values (is sample from the MCMC large enough)
effectiveSize(chain.alfl)  # if there is correlation, the effective sample size will be smaller than the sample size (number of iterations) used as input (100000)
geweke.diag(chain.alfl, frac1=0.1, frac2=0.5)  # Z-scores (reject H0 (equality of means => convergence) if z < -1.96 or z > +1.96, when α = 0.05)
heidel.diag(chain.alfl, eps=0.1, pvalue=0.05)  # Heidelberg and Welch’s convergence diagnostic, which is also based on a single chain
# We run only 1 chain with BayeScan. TO compare multiple chains, a popular test is the Gelman and Rubin’s diagnostic (but see: https://stats.stackexchange.com/questions/507/what-is-the-best-method-for-checking-convergence-in-mcmc).
# combined = mcmc.list(chain1,chain2)
# plot(combined)
# gelman.diag(combined)
# gelman.plot(combined,ask)
# gelman.diag is based on a comparison of between and within chain variances. If the chains converged, these two variances should be equal.
# The output of gelman.diag are the scale reduction factors for each parameter. A factor of 1 means that between- and within-chain variance are the same; 
# larger values mean that they are fairly different, indicating convergence problems. 
# The rule of thumb is that values below 1.1 or so are OK but, being a rule of thumb, this is not an extremely rigorous criterion.

# Combine list of loci and bayescan output
d.bayescan.alfl <- cbind(loc.alfl, out.bayescan.alfl)

# Assign selection categories based on alpha and qval
plot_bayescan(out.bayescan.alfl, FDR = 0.05)
d.bayescan.alfl <- d.bayescan.alfl %>% 
  mutate(selection = factor(ifelse(qval <= 0.05 & alpha > 0, "diversifying",
                                   ifelse(qval <= 0.05 & alpha < 0, "balancing",
                                          "neutral")), levels = c("neutral","balancing","diversifying")))
d.bayescan.alfl %>% group_by(selection) %>% summarise(N = n())
# selection     N
# neutral   46992

# Outlier loci plot
gg.bayescan.alfl <- d.bayescan.alfl %>% 
  mutate(log10qval = -log10(qval)) %>% 
  ggplot(aes(x = log10qval, y = fst)) + 
  geom_point(aes(fill = selection), pch = 21, size = 2) + 
  scale_fill_manual(name = "Selection", values = c("white","orange","red")) +
  labs(x = "Log(q-value)", y = "Fst") +
  ggtitle('Bayescan Outlier Analysis (n = 46992, FDR = 0.05)') +
  xlim(0,4) +
  geom_vline(aes(xintercept = 1.30103), color = "black") +
  geom_vline(aes(xintercept = 4), colour = "blue", linetype = "dashed") +
  theme_classic()
gg.bayescan.alfl
# ggsave(filename = "./02_Results/05_filtering/01_reference/04_neutral/ALFL/Bayescan.alfl.outliers.png",
#        plot = gg.bayescan.alfl,
#        width = 7, height = 6, units = "in")

# List of neutral loci to keep for vcf file: no outlier SNPs detected
bayescan.outliers.alfl <- d.bayescan.alfl %>%
  filter(selection %nin% "neutral") %>% select(ID)
# write.table(bayescan.outliers.alfl,
#             "./00_Data/03_filtering/03_neutral/ALFL/outlier.bayescan.alfl", 
#             quote = F, row.names = F, col.names = F, sep = "\t")


## CLSW -------------------------------------------------------------------

loc.clsw <- data.frame(vcf.clsw@fix[, c("CHROM","POS","ID")])

out.bayescan.clsw <- read.table("./00_Data/03_filtering/03_neutral/CLSW/CLSW_outlier_fst_fst.txt",
                                          header = T)

# Evaluate MCMC chain convergence
# https://evomics.org/learning/population-and-speciation-genomics/2016-population-and-speciation-genomics/bayescan-exercise/
chain.clsw <- read.table("./00_Data/03_filtering/03_neutral/CLSW/CLSW_outlier_fst.sel",
                          header = T)
chain.clsw <- mcmc(chain.clsw, thin = 10)
# plot(chain.clsw)  # check trace plots and posterior distribution
summary(chain.clsw)
autocorr.diag(chain.clsw)  # check correlation between sampled parameter values (is sample from the MCMC large enough)
effectiveSize(chain.clsw)  # if there is correlation, the effective sample size will be smaller than the sample size used as input (100000)
geweke.diag(chain.clsw, frac1=0.1, frac2=0.5)  # Z-scores (reject H0 (equality of means => convergence) if z < -1.96 or z > +1.96, when α = 0.05)
heidel.diag(chain.clsw, eps=0.1, pvalue=0.05)  # Heidelberg and Welch’s convergence diagnostic, which is also based on a single chain

# Combine list of loci and bayescan output
d.bayescan.clsw <- cbind(loc.clsw, out.bayescan.clsw)

# Assign selection categories based on alpha and qval
plot_bayescan(out.bayescan.clsw, FDR = 0.05)
d.bayescan.clsw <- d.bayescan.clsw %>% 
  mutate(selection = factor(ifelse(qval <= 0.05 & alpha > 0, "diversifying",
                                   ifelse(qval <= 0.05 & alpha < 0, "balancing",
                                          "neutral")), levels = c("neutral","balancing","diversifying")))
d.bayescan.clsw %>% group_by(selection) %>% summarise(N = n())
# selection     N
# neutral   21946

# Outlier loci plot
gg.bayescan.clsw <- d.bayescan.clsw %>% 
  mutate(log10qval = -log10(qval)) %>% 
  ggplot(aes(x = log10qval, y = fst)) + 
  geom_point(aes(fill = selection), pch = 21, size = 2) + 
  scale_fill_manual(name = "Selection", values = c("white","orange","red")) +
  labs(x = "Log(q-value)", y = "Fst") +
  ggtitle('Bayescan Outlier Analysis (n = 21946, FDR = 0.05)') +
  xlim(0,4) +
  geom_vline(aes(xintercept = 1.30103), color = "black") +
  geom_vline(aes(xintercept = 4), colour = "blue", linetype = "dashed") +
  theme_classic()
gg.bayescan.clsw
# ggsave(filename = "./02_Results/05_filtering/01_reference/04_neutral/CLSW/Bayescan.clsw.outliers.png",
#        plot = gg.bayescan.clsw,
#        width = 7, height = 6, units = "in")

# List of neutral loci to keep for vcf file: no outlier SNPs detected
bayescan.outliers.clsw <- d.bayescan.clsw %>%
  filter(selection %nin% "neutral") %>% select(ID)
# write.table(bayescan.outliers.clsw,
#             "./00_Data/03_filtering/03_neutral/CLSW/outlier.bayescan.clsw", 
#             quote = F, row.names = F, col.names = F, sep = "\t")


## PUMA -------------------------------------------------------------------

loc.puma <- data.frame(vcf.puma@fix[, c("CHROM","POS","ID")])

out.bayescan.puma <- read.table("./00_Data/03_filtering/03_neutral/PUMA/PUMA_outlier_fst_fst.txt",
                                header = T)

# Evaluate MCMC chain convergence
# https://evomics.org/learning/population-and-speciation-genomics/2016-population-and-speciation-genomics/bayescan-exercise/
chain.puma <- read.table("./00_Data/03_filtering/03_neutral/PUMA/PUMA_outlier_fst.sel",
                         header = T)
chain.puma <- mcmc(chain.puma, thin = 10)
# plot(chain.puma)  # check trace plots and posterior distribution
summary(chain.puma)
autocorr.diag(chain.puma)  # check correlation between sampled parameter values (is sample from the MCMC large enough)
effectiveSize(chain.puma)  # if there is correlation, the effective sample size will be smaller than the sample size used as input (100000)
geweke.diag(chain.puma, frac1=0.1, frac2=0.5)  # Z-scores (reject H0 (equality of means => convergence) if z < -1.96 or z > +1.96, when α = 0.05)
heidel.diag(chain.puma, eps=0.1, pvalue=0.05)  # Heidelberg and Welch’s convergence diagnostic, which is also based on a single chain

# Combine list of loci and bayescan output
d.bayescan.puma <- cbind(loc.puma, out.bayescan.puma)

# Assign selection categories based on alpha and qval
plot_bayescan(out.bayescan.puma, FDR = 0.05)
d.bayescan.puma[d.bayescan.puma$qval %in% 0, "qval"] <- 0.0000001  #assign this value to those SNPs with qval = 0
d.bayescan.puma <- d.bayescan.puma %>% 
  mutate(selection = factor(ifelse(qval <= 0.05 & alpha > 0, "diversifying",
                                   ifelse(qval <= 0.05 & alpha < 0, "balancing",
                                          "neutral")), levels = c("neutral","balancing","diversifying")))
d.bayescan.puma %>% group_by(selection) %>% summarise(N = n())
# selection         N
# neutral      108506
# diversifying    188

# Outlier loci plot
gg.bayescan.puma <- d.bayescan.puma %>% 
  mutate(log10qval = -log10(qval)) %>% 
  ggplot(aes(x = log10qval, y = fst)) + 
  geom_point(aes(fill = selection), pch = 21, size = 2) + 
  scale_fill_manual(name = "Selection", values = c("white","orange","red")) +
  labs(x = "Log(q-value)", y = "Fst") +
  ggtitle('Bayescan Outlier Analysis (n = 108694, FDR = 0.05)') +
  xlim(0,7) +
  geom_vline(aes(xintercept = 1.30103), color = "black") +
  geom_vline(aes(xintercept = 4), colour = "blue", linetype = "dashed") +
  theme_classic()
gg.bayescan.puma
# ggsave(filename = "./02_Results/05_filtering/01_reference/04_neutral/PUMA/Bayescan.puma.outliers.png",
#        plot = gg.bayescan.puma,
#        width = 7, height = 6, units = "in")

# List of neutral loci to keep for vcf file
bayescan.outliers.puma <- d.bayescan.puma %>%
  filter(selection %nin% "neutral") %>% select(ID)
# write.table(bayescan.outliers.puma,
#             "./00_Data/03_filtering/03_neutral/PUMA/outlier.bayescan.puma",
#             quote = F, row.names = F, col.names = F, sep = "\t")


## VGSW -------------------------------------------------------------------

loc.vgsw <- data.frame(vcf.vgsw@fix[, c("CHROM","POS","ID")])

out.bayescan.vgsw <- read.table("./00_Data/03_filtering/03_neutral/VGSW/VGSW_outlier_fst.txt",
                                header = T)

# Evaluate MCMC chain convergence
# https://evomics.org/learning/population-and-speciation-genomics/2016-population-and-speciation-genomics/bayescan-exercise/
chain.vgsw <- read.table("./00_Data/03_filtering/03_neutral/VGSW/VGSW_outlier.sel",
                         header = T)
chain.vgsw <- mcmc(chain.vgsw, thin = 10)
# plot(chain.vgsw)  # check trace plots and posterior distribution
summary(chain.vgsw)
autocorr.diag(chain.vgsw)  # check correlation between sampled parameter values (is sample from the MCMC large enough)
effectiveSize(chain.vgsw)  # if there is correlation, the effective sample size will be smaller than the sample size (number of iterations) used as input (100000)
geweke.diag(chain.vgsw, frac1=0.1, frac2=0.5)  # Z-scores (reject H0 (equality of means => convergence) if z < -1.96 or z > +1.96, when α = 0.05)
heidel.diag(chain.vgsw, eps=0.1, pvalue=0.05)  # Heidelberg and Welch’s convergence diagnostic, which is also based on a single chain

# Combine list of loci and bayescan output
d.bayescan.vgsw <- cbind(loc.vgsw, out.bayescan.vgsw)

# Assign selection categories based on alpha and qval
plot_bayescan(out.bayescan.vgsw, FDR = 0.05)
d.bayescan.vgsw <- d.bayescan.vgsw %>% 
  mutate(selection = factor(ifelse(qval <= 0.05 & alpha > 0, "diversifying",
                                   ifelse(qval <= 0.05 & alpha < 0, "balancing",
                                          "neutral")), levels = c("neutral","balancing","diversifying")))
d.bayescan.vgsw %>% group_by(selection) %>% summarise(N = n())
# selection         N
# neutral      110184
# diversifying      5

# Outlier loci plot
gg.bayescan.vgsw <- d.bayescan.vgsw %>% 
  mutate(log10qval = -log10(qval)) %>% 
  ggplot(aes(x = log10qval, y = fst)) + 
  geom_point(aes(fill = selection), pch = 21, size = 2) + 
  scale_fill_manual(name = "Selection", values = c("white","orange","red")) +
  labs(x = "Log(q-value)", y = "Fst") +
  ggtitle('Bayescan Outlier Analysis (n = 110189, FDR = 0.05)') +
  xlim(0,4) +
  geom_vline(aes(xintercept = 1.30103), color = "black") +
  geom_vline(aes(xintercept = 4), colour = "blue", linetype = "dashed") +
  theme_classic()
gg.bayescan.vgsw
# ggsave(filename = "./02_Results/05_filtering/01_reference/04_neutral/VGSW/Bayescan.vgsw.outliers.png",
#        plot = gg.bayescan.vgsw,
#        width = 7, height = 6, units = "in")

# List of neutral loci to keep for vcf file: no outlier SNPs detected
bayescan.outliers.vgsw <- d.bayescan.vgsw %>%
  filter(selection %nin% "neutral") %>% select(ID)
# write.table(bayescan.outliers.vgsw,
#             "./00_Data/03_filtering/03_neutral/VGSW/outlier.bayescan.vgsw",
#             quote = F, row.names = F, col.names = F, sep = "\t")

