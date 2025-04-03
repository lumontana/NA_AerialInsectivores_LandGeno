# Info --------------------------------------------------------------------
# 
# Author: Luca Montana
# Affiliation: University of Lethbridge
# Group: T. Burg research group
# Location: Rimouski, QC
# Date: 2024-04-04
# 
# Overview: Creation of array files for arrayjobs on bash
# 
#

# Library -----------------------------------------------------------------

library(readxl)
library(tidyverse)

# Internal functions
for(i in 1:length(list.files("./01_Codes/00_Functions"))){
  source(file.path("./01_Codes/00_Functions",  list.files("./01_Codes/00_Functions")[i]))  
}
`%nin%` <- Negate(`%in%`)



# Data --------------------------------------------------------------------

d_meta <- read_excel("./00_Data/00_fileinfo/NGS 3D samples master sheet July 2024 TB_LM241129.xlsx", sheet = "all in one sheet")  # metadata

d_meta <- d_meta %>% filter(spp %in% c("ALFL","CLSW","PUMA","VGSW")) %>% 
  select(1:3,5:8)
colnames(d_meta) <- c("Sort","Sent","Plate_ID","ID","Species","Band_no","Location")

table(d_meta$Species)
# ALFL CLSW PUMA VGSW 
# 103  139   66   74

d_meta <- d_meta %>% 
  mutate(Plate = paste0("Plate_",sub("^.*Plate[_]?([0-9]+).*", "\\1", Plate_ID)),  # Extract Plate number
         Well = sub("^.*/([0-9]+[a-zA-Z]+)$", "\\1", Plate_ID)) %>%  # Extract Well info 
  mutate(Region = gsub("^[A-Z]{4}_|_[0-9]+$|[0-9]+$", "", ID)) %>%  # I manually changed some things on excel file that were causing troubles:
  # space after the species name for some ALFL individuals from AB. Also, some ALFL YT samples had just no _ between species name, location, 
  # and numeric identifier
  mutate(Region = ifelse(Region %in% c("CAB_Edm", "AB"), "AB",
                         ifelse(Region %in% c("CBC_Pri","BC_Smc"), "BC",
                                ifelse(Region %in% c("MT_Tal"), "MT",
                                       ifelse(Region %in% "NB_Que", "NB",
                                              ifelse(Region %in% c("NL","NL_Pas","NL_Ter"), "NL",
                                                     ifelse(Region %in% c("SK","SK_Pan","SK_Bls","SK_Ind","SK_xxx"), "SK",
                                                            ifelse(Region %in% c("SWON","SWON_xxx"), "SWON",
                                                                   ifelse(Region %in% c("VI", "VI_BC"), "VI",
                                                                          Region))))))))) #%>% 
  #filter(!duplicated(ID))  # keeping metadata to be as clean as possible. Removing duplicated ID because some were sent to genotype multiple times
table(d_meta$Species)
# ALFL CLSW PUMA VGSW 
# 103  139   66   74
table(d_meta$Sent)
# 20210222 20231102 20240711 
#       52      192      138



## Identify duplicated and merged samples ---------------------------------

# How I would have thought it was done...
# dup.id <- d_meta$ID[duplicated(d_meta$ID)]
# d_meta.dup <- d_meta[d_meta$ID %in% dup.id,] %>% 
#   arrange(ID)
# d_meta.uni <- d_meta[d_meta$ID %nin% dup.id,] %>% 
#   arrange(ID)
# # Rename seq2 for duplicated samples
# d_meta.dup <- d_meta.dup %>% 
#   group_by(ID) %>%  # group by original ID without suffix
#   mutate(Dup = row_number()) %>%  # Number the duplicates among original IDs
#   ungroup() %>% 
#   mutate(ID_dra = ifelse(Dup %in% 1, ID, 
#                          paste0(ID, "_seq2")))
# 
# Add Dup column in uni dataset as well
# d_meta.uni <- d_meta.uni %>% 
#   group_by(ID) %>%  # group by original ID without suffix
#   mutate(Dup = row_number()) %>%  # Number the duplicates among original IDs
#   ungroup() %>% 
#   mutate(ID_dra = ID)
# 
# d_meta <- bind_rows(d_meta.uni, d_meta.dup) %>% 
#   arrange(ID_dra) %>% 
#   select(1,4,12,11,5,10,7,6,3,2,8,9)

# How it is actually done...
d_meta <- d_meta %>% 
  group_by(ID) %>%  # group by original ID without suffix
  mutate(Dup = row_number()) %>%  # Number the duplicates among original IDs
  ungroup() %>% 
  mutate(ID_dra = ifelse(Sent %in% "20240711", paste0(ID, "_seq2"),
                         ID)) %>% 
  arrange(ID_dra) %>%
  select(1,4,12,11,5,10,7,6,3,2,8,9)

d_meta %>% filter(Dup %in% 1) %>% group_by(Species, Region) %>% summarise(N = n()) %>% print(n = 50)
# Species Region     N
# ALFL    AB        17
# ALFL    BC         1
# ALFL    GU         2
# ALFL    HO         1
# ALFL    LA         4
# ALFL    MT         1
# ALFL    NB         1
# ALFL    NBC       14
# ALFL    NL         4
# ALFL    QC         4
# ALFL    SK        20
# ALFL    YT        19
# CLSW    AZ         5
# CLSW    BC        13
# CLSW    CO        10
# CLSW    MB        14
# CLSW    MS         8
# CLSW    MX         5
# CLSW    NBC        4
# CLSW    ON        14
# CLSW    SK        19
# CLSW    TX         4
# CLSW    VI         1
# CLSW    WA         8
# PUMA    NC         4
# PUMA    NV         3
# PUMA    PI_BC      8
# PUMA    SBC        3
# PUMA    SK        14
# PUMA    SWON      12
# PUMA    VI         1
# PUMA    WA         6
# VGSW    CO         9
# VGSW    NM         2
# VGSW    NV        10
# VGSW    PI_BC     10
# VGSW    SBC       10
# VGSW    WA        10


## Dates ------------------------------------------------------------------

d.date <- read_excel("./00_Data/00_fileinfo/AI data with dates.xlsx", sheet = "combined", col_types = c("text","date",rep("text",18)))  # metadata
d.date <- d.date %>% mutate(species = toupper(species)) %>% 
  filter(!is.na(`band number`))
# dup <- d.date$`band number`[duplicated(d.date$`band number`)]
# d.date.dup <- d.date[d.date$`band number` %in% dup,]
d.date <- d.date[!duplicated(d.date$`band number`),]

# bandno_dup <- d_meta$Band_no[duplicated(d_meta$Band_no)]  #  2 duplicated band no.
# d.dup_bandno <- d_meta[d_meta$Band_no %in% bandno_dup,]

d_meta <- d_meta %>% left_join(d.date %>% select(2:4), by = c("Species"="species", "Band_no"="band number"))

# ALFL: BS 21 Jun-12 Jul, Win 15 Nov-29 Mar
# CLSW: BS 7 Jun-12 Jul, Win 29 Nov-11 Jan
# PUMA: BS 17 May-21 Jun, Win 1 Nov-27 Dec
# VGSW: BS 17 May-12 Jul (Northern Mex + US + CA; year-round Mex + southern regions), Win 22 Nov-18 Jan
d_meta <- d_meta %>%
  mutate(Month = as.integer(format(date, "%m")),
         Day = as.integer(format(date, "%d")),
         Year = as.integer(format(date, "%Y")),
         JD = as.integer(format(date, "%j"))) #%>% 
  # mutate(Season = ifelse(is.na(date), "Unknown",
  #                        ifelse(Month %in% c(6:7), "Breeding",
  #                               ifelse(Month %in% c(11,12,1,2,3), "Wintering",
  #                                      "Migrating"))))


# Save dataset ------------------------------------------------------------

write.csv(d_meta, file = "./00_Data/00_fileinfo/AI_metadata.csv", row.names = F)
