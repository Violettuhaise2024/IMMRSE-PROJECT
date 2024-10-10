###################################################
###          Install required packages          ###
###################################################
library(readxl)
library(janitor)
library(tidyverse)

###################################################
###                 Import Data                 ###
###################################################
rm(list=ls())
detach(package:plotly, unload=TRUE)
setwd("D:/Computer/IMMRSE/Data-from-pipeline/Drug_resistance/Round3/")

## import variant table
var_tab <- read.csv("C:/Users/Thomas Katairo/Box/IMMRSE-U Study/Melissa Drug Resistance Code/variant_table.csv")
var_labs <- subset(var_tab, !is.na(mutation_name), select = -c(no_alleles, AA))
var_no_alleles <- subset(var_tab, !is.na(no_alleles), select = c(name, mutation_name, no_alleles))
##  paragon DR genotype data
IM_R101 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-05-04-IM-RUN-101/resistance_marker_module/resmarker_table.txt")
IM_R102 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-06-07-IM-RUN-102/resistance_marker_module/resmarker_table.txt")
IM_R103 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-06-11-IM-RUN-103/resistance_marker_module/resmarker_table.txt")
IM_R104 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-06-19-IM-RUN-104/resistance_marker_module/resmarker_table.txt")
IM_R105 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-06-30-IM-RUN-105/resistance_marker_module/resmarker_table.txt")
IM_R106 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-07-01-IM-RUN-106/resistance_marker_module/resmarker_table.txt")
IM_R107 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-07-16-IM-RUN-107/resistance_marker_module/resmarker_table.txt")
IM_R108 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-07-18-IM-RUN-108/resistance_marker_module/resmarker_table.txt")
IM_R109 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-07-25-IM-RUN-109/resistance_marker_module/resmarker_table.txt")
IM_R110 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-08-07-IM-RUN-110/resistance_marker_module/resmarker_table.txt")
IM_R111 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-08-09-IM-RUN-111/resistance_marker_module/resmarker_table.txt")
IM_R112 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-08-12-IM-RUN-112/resistance_marker_module/resmarker_table.txt")
IM_R113 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-08-14-IM-RUN-113/resistance_marker_module/resmarker_table.txt")
IM_R114 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-08-20-IM-RUN-114/resistance_marker_module/resmarker_table.txt")
IM_R116 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-09-02-IM-RUN-116/resistance_marker_module/resmarker_table.txt")
IM_R117 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-09-09-IM-RUN-117/resistance_marker_module/resmarker_table.txt")
IM_R118 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-09-13-IM-RUN-118/resistance_marker_module/resmarker_table.txt")
IM_R119 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-09-14-IM-RUN-119/resistance_marker_module/resmarker_table.txt")
IM_R120 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-09-16-IM-RUN-120/resistance_marker_module/resmarker_table.txt")
IM_R121 <- read.delim("D:/Computer/IMMRSE/Data-from-pipeline/Results-2024-09-17-IM-RUN-121/resistance_marker_module/resmarker_table.txt")

paragon_df_R3 <- bind_rows(IM_R101, IM_R102, IM_R103, IM_R104, IM_R105, IM_R106, IM_R107, IM_R108, IM_R109, IM_R110, IM_R111, IM_R112, 
                           IM_R113, IM_R114, IM_R116, IM_R117, IM_R118, IM_R119, IM_R120, IM_R121)

metadata <- read.csv("C:/Users/Thomas Katairo/Box/IMMRSE-U Study/UMSP and Travel Survey Databases/qPCR data/R3 Preliminary/qpcr_r3_immrse.csv")
metadata <- metadata %>% select(StudySubject, Health.facility)
metadata <- metadata %>% rename(barcode= StudySubject) 
metadata <- metadata %>% rename(Sample.Site= Health.facility)

###################################################
###         Format Paragon Data                ###
###################################################
#paragon_df_R3 <- paragon_df_R3 %>% left_join(metadata)

paragon_df_R3 <- 
  paragon_df_R3 %>%
  separate(SampleID, sep = "_", into = c("barcode", NA, NA))

paragon_df <- paragon_df_R3 %>% left_join(metadata)

paragon_df$barcode <- gsub("-(20|15)$", "", paragon_df$barcode)

paragon_grouped <- paragon_df %>% group_by(barcode, Sample.Site, Gene, CodonID, AA) %>% summarize(reads = sum(Reads))

#create name variable
paragon_grouped$name <- paste(paragon_grouped$Gene, paragon_grouped$CodonID, sep = "_")

##############################
###   Link read count      ###
##############################
paragon_grouped_m <- merge(paragon_grouped, var_tab, by = c("Gene", "CodonID", "name", "AA"), all.x = T)

## see if any variants need to be added ot the var_tab table
pending <- subset(paragon_grouped_m, is.na(type) | type=="")
View(pending)
#write.csv(pending, "pending.csv")

paragon_grouped_m_ss <- subset(paragon_grouped_m, select = c("barcode", "Sample.Site", "Gene", "CodonID", "name", "reads", "type"))
wide <- pivot_wider(paragon_grouped_m_ss, names_from = "type", names_prefix = "bc_", values_from = "reads", values_fill = 0)

## check max number of alt alleles; 
table(paragon_grouped_m$type)

## calculate coverage
wide$coverage <- wide$bc_ref + wide$bc_alt + wide$bc_alt2 + wide$bc_alt3 + wide$bc_alt4

###################################################
###                     Filter                  ###
###################################################
geno_filter <- subset(wide, coverage > 10) # filters out all negative controls with reads 
geno_filter <- geno_filter %>% mutate(ref = case_when(bc_ref   > 10 & geno_filter$bc_ref/geno_filter$coverage > 0.05 ~ 1, .default = 0),
                                      alt = case_when(bc_alt   > 10 & geno_filter$bc_alt/geno_filter$coverage > 0.05 ~ 1, .default = 0),
                                      alt2 = case_when(bc_alt2 > 10 & geno_filter$bc_alt2/geno_filter$coverage > 0.05 ~ 1, .default = 0),
                                      alt3 = case_when(bc_alt3 > 10 & geno_filter$bc_alt3/geno_filter$coverage > 0.05 ~ 1, .default = 0),
                                      alt4 = case_when(bc_alt4 > 10 & geno_filter$bc_alt4/geno_filter$coverage > 0.05 ~ 1, .default = 0)) %>% 
                          mutate(genotype = case_when(ref==1 & alt == 0 & alt2==0 & alt3==0 & alt4==0 ~ 0,
                                                      ref==1 & (alt == 1 | alt2==1 | alt3==1 | alt4==1) ~ 1,
                                                      ref==0 & (alt == 1 | alt2==1 | alt3==1 | alt4==1) ~ 2),
                                 mixed_alt = case_when(alt + alt2 +alt3 +alt4 > 1 ~ 1, .default=0))

geno_filter_lab <- merge(geno_filter, var_no_alleles, by = "name")

geno_filter_bi <- subset(geno_filter_lab, no_alleles==2)
geno_filter_3 <- subset(geno_filter_lab, no_alleles==3)
geno_filter_4 <- subset(geno_filter_lab, no_alleles==4)
geno_filter_5 <- subset(geno_filter_lab, no_alleles==5)

###################################################
###       Assign Genotype and Mutation Name     ###
###################################################
### Dealing with bialleleic loci:
geno_filter_bi_ss <- subset(geno_filter_bi, select = c("barcode", "Sample.Site", "name", "bc_ref", "bc_alt", "bc_alt2", "bc_alt3", "bc_alt4", "genotype", "mutation_name"))

##############
### Dealing with triallelic loci:
all_3 <- unique(subset(geno_filter_3, select = name))$name
## make a data frame for the loci that have 3 alleles in our data.
geno_filter_3_ss <- geno_filter_3 %>% subset(., name %in% all_3, select = -c(genotype, mixed_alt, mutation_name)) %>% 
                                   mutate(genotype_alt = case_when(ref==1 & alt==0 ~ 0,
                                                                     ref==0 & alt==1 ~ 2,
                                                                     ref==1 & alt==1 ~ 1,
                                                                     alt==1 & alt2==1 ~ 3),
                                            genotype_alt2 = case_when(ref==1 & alt2==0 ~ 0,
                                                                   ref==0 & alt2==1 ~ 2,
                                                                   ref==1 & alt2==1 ~ 1,
                                                                   alt==1 & alt2==1 ~ 3))
geno_filter_3_long <- pivot_longer(geno_filter_3_ss, cols = genotype_alt:genotype_alt2, names_to = "variant", values_to = "genotype")






geno_filter_3_long <- merge(geno_filter_3_long, var_labs, by = c("Gene", "CodonID", "name", "variant"))## Here Here
geno_filter_3_long <- subset(geno_filter_3_long, select = c("barcode", "Sample.Site", "name", "bc_ref", "bc_alt", "bc_alt2", "bc_alt3", "bc_alt4", "genotype", "mutation_name"))

##############
### Dealing with 4 allele loci:
all_4 <- unique(subset(geno_filter_4, select = name))$name
## make a data frame for the loci that have >2 alleles in our data.
geno_filter_4_ss <- geno_filter_4 %>% subset(., name %in% all_4, select = -c(genotype, mixed_alt, mutation_name)) %>% 
  mutate(genotype_alt = case_when(ref==1 & alt==0 ~ 0,
                                  ref==0 & alt==1 ~ 2,
                                  ref==1 & alt==1 ~ 1,
                                  alt==1 & alt2==1 ~ 3,
                                  alt==1 & alt3==1 ~ 3,
                                  alt2==1 & alt3==1 ~3),
         genotype_alt2 = case_when(ref==1 & alt2==0 ~ 0,
                                   ref==0 & alt2==1 ~ 2,
                                   ref==1 & alt2==1 ~ 1,
                                   alt==1 & alt2==1 ~ 3,
                                   alt==1 & alt3==1 ~ 3,
                                   alt2==1 & alt3==1 ~3),
         genotype_alt3 = case_when(ref==1 & alt3==0 ~ 0,
                                   ref==0 & alt3==1 ~ 2,
                                   ref==1 & alt3==1 ~ 1,
                                   alt==1 & alt2==1 ~ 3,
                                   alt==1 & alt3==1 ~ 3,
                                   alt2==1 & alt3==1 ~3))
geno_filter_4_long <- pivot_longer(geno_filter_4_ss, cols = genotype_alt:genotype_alt3, names_to = "variant", values_to = "genotype")
geno_filter_4_long <- merge(geno_filter_4_long, var_labs, by = c("Gene", "CodonID", "name", "variant"))
geno_filter_4_long <- subset(geno_filter_4_long, select = c("barcode", "Sample.Site", "name", "bc_ref", "bc_alt", "bc_alt2", "bc_alt3", "bc_alt4", "genotype", "mutation_name"))

### Dealing with 5 allele loci:
all_5 <- unique(subset(geno_filter_5, select = name))$name

## make a data frame for the loci that have >2 alleles in our data.
geno_filter_5_ss <- geno_filter_5 %>% subset(., name %in% all_5, select = -c(genotype, mixed_alt, mutation_name)) %>% 
  mutate(genotype_alt = case_when(ref==1 & alt==0 ~ 0,
                                  ref==0 & alt==1 ~ 2,
                                  ref==1 & alt==1 ~ 1,
                                  alt==1 & alt2==1 ~ 3,
                                  alt==1 & alt3==1 ~ 3,
                                  alt2==1 & alt3==1 ~3),
         genotype_alt2 = case_when(ref==1 & alt2==0 ~ 0,
                                   ref==0 & alt2==1 ~ 2,
                                   ref==1 & alt2==1 ~ 1,
                                   alt==1 & alt2==1 ~ 3,
                                   alt==1 & alt3==1 ~ 3,
                                   alt2==1 & alt3==1 ~3),
         genotype_alt3 = case_when(ref==1 & alt3==0 ~ 0,
                                   ref==0 & alt3==1 ~ 2,
                                   ref==1 & alt3==1 ~ 1,
                                   alt==1 & alt2==1 ~ 3,
                                   alt==1 & alt3==1 ~ 3,
                                   alt2==1 & alt3==1 ~3),
         genotype_alt4 = case_when(ref==1 & alt4==0 ~ 0,
                                   ref==0 & alt4==1 ~ 2,
                                   ref==1 & alt4==1 ~ 1,
                                   alt==1 & alt4==1 ~ 3,
                                   alt==1 & alt4==1 ~ 3,
                                   alt2==1 & alt4==1 ~3))
geno_filter_5_long <- pivot_longer(geno_filter_5_ss, cols = genotype_alt:genotype_alt4, names_to = "variant", values_to = "genotype")
geno_filter_5_long <- merge(geno_filter_5_long, var_labs, by = c("Gene", "CodonID", "name", "variant"))
geno_filter_5_long <- subset(geno_filter_5_long, select = c("barcode", "Sample.Site", "name", "bc_ref", "bc_alt", "bc_alt2", "bc_alt3", "bc_alt4", "genotype", "mutation_name"))


###################################################
###               Bind data frames              ###
###################################################
data_pgn <- rbind(geno_filter_bi_ss, geno_filter_3_long, geno_filter_4_long, geno_filter_5_long)

## make sure all alleles got assigned a mutation_name. This should be an empty data frame
View(subset(data_pgn, mutation_name==""))

rm(list = setdiff(ls(), c("data_pgn", "metadata", "sample_metadata", "pending")))

###########
# IMMRSE-U EXPERIMENTAL
###########
data_exp <- subset(data_pgn, select = c("barcode", "Sample.Site", "genotype", "mutation_name"))
data_exp_wide <- pivot_wider(data_exp, names_from = mutation_name, values_from = genotype)

nosite <- data_exp_wide %>% dplyr::filter(is.na(Sample.Site))
#Ignore the chahafi samples because they were in another round

data_exp_wide <- data_exp_wide %>%
  mutate(Sample.Site = case_when(
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-011", barcode) ~ "Aduku",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-016", barcode) ~ "Kasambya",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-019", barcode) ~ "Nagongera",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-030", barcode) ~ "Amolatar",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-034", barcode) ~ "Orum",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-036", barcode) ~ "Aboke",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-041", barcode) ~ "Lalogi",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-042", barcode) ~ "Patongo",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-047", barcode) ~ "Opia",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-054", barcode) ~ "Kigorobya",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-056", barcode) ~ "Karambi",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-060", barcode) ~ "Lokolia",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-063", barcode) ~ "Lobule",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-064", barcode) ~ "Nadunget",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-087", barcode) ~ "Kyatiri",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-089", barcode) ~ "Nawaikoke",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-090", barcode) ~ "Budondo",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-096", barcode) ~ "Kigandalo",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-110", barcode) ~ "Busitema",
    is.na(Sample.Site) & grepl("IM-[0-9]{2}-120", barcode) ~ "Muko",
    is.na(Sample.Site) & grepl("Z2", barcode) ~ "Kamwezi",
    TRUE ~ as.character(Sample.Site)
  ))

nosite <- data_exp_wide %>% dplyr::filter(is.na(Sample.Site))

write.csv(data_exp_wide, "databases/sample_genotypes_R3.csv", row.names = FALSE)
