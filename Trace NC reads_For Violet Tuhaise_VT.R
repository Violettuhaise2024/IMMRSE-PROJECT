# Fri Jan 12 11:34:24 2024 ------------------------------
#Clearing the environment and console
rm(list = ls())
cat("\14")
#load the necessary library
library(tidyverse)
library(dplyr)

#need to pull in the sample sheet for run 52 so we have well positions
sample_sheet_Run52 <- read_csv("C:/Users/USER/OneDrive/Desktop/Bienvenu Files/IMMRSE-U/Phase2/QC/Samplesheet_Run_52_Samplesheet_Run_52.csv.csv")
sample_sheet_Run52$sampleID <- sample_sheet_Run52$SampleID
sample_sheet_Run52 <- sample_sheet_Run52 %>% select(-sampleID)
start_row <- 97
sample_sheet_Run52 <- sample_sheet_Run52[-(start_row:nrow(sample_sheet_Run52)), ]

#pull in run52 data 
run52 <- read_delim("C:/Users/USER/OneDrive/Desktop/Bienvenu Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-01-03-IM-RUN-52/Results-2024-01-03-IM-RUN-52/allele_data.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE) %>% 
  mutate(sampleID = word(sampleID, 1, sep = "_")) 

#join in sample sheet info so that you have plate position 
Run52_samplesheet_merged <- run52 %>% left_join(sample_sheet_Run52, by = "sampleID")

#if no plate position, then it was the other pool on the run, so filter it out
Run52_samplesheet_merged <- Run52_samplesheet_merged %>% filter(!is.na(Position))

#check number of samples
check <- Run52_samplesheet_merged %>% select(sampleID) %>% unique()

#Extract the allele data for that contaminated NC
bad_nc_allele <- Run52_samplesheet_merged %>% filter(sampleID == "Negative-Control-5") %>% pull(allele) 

# Mark the alleles presented in that contaminated NC
Run52_samplesheet_merged <- Run52_samplesheet_merged %>% 
  mutate(nc_allele = ifelse(allele %in% bad_nc_allele, 1, 0))

# Proportion of alleles that are present in contaminated NC 
Plate52_norm <- Run52_samplesheet_merged %>% group_by(sampleID, Position, nc_allele) %>% summarize(count = n())
Plate52_norm <- Plate52_norm %>% group_by(sampleID) %>% mutate(total_allele = sum(count))
Plate52_norm <- Plate52_norm %>% mutate(prop_nc_allele = count/total_allele)
Plate52_norm %<>% 
  separate(col = Position, into = c("Row", "Col"), sep = 1) %>%  #separate columns and rows 
  mutate(Col = as.numeric(Col)) %>% 
  mutate(group = case_when(
    grepl("Negative", sampleID) ~ "Negative",
    grepl("Positive", sampleID) ~ "Positive",
    TRUE ~ "Sample"))

Plate52_norm$Row <- as.factor(Plate52_norm$Row)
Plate52_norm$Col <- as.factor(Plate52_norm$Col)

Plate52_norm <- Plate52_norm %>% 
  mutate(group = factor(group, levels = c("Negative", "Positive", "Sample"))) #level the factor

my_colors <- c("Negative" = "red", "Positive" = "blue", "Sample" = "white")

heat_lab <- seq(0, 1, 0.1)
heat_break <- seq(0, 1, 0.1)

#I don't know why the blue outline of the positive control is not appearing
Plate52_norm %>% 
  ggplot(aes(x = Col, y = Row, text = paste(sampleID))) +
  geom_tile(aes(fill = prop_nc_allele, color = group), size = 2) +
  scale_colour_manual(values = my_colors) + 
  scale_size_manual("group", values = c(0, 1, 1)) +
  scale_fill_gradient(low = "blue", high = "yellow", trans = "log",
                      breaks = heat_break, labels = heat_lab) +
  ylim(rev(levels(Plate52_norm$Row))) + 
 # theme_bw() +
  theme(axis.title = element_blank()) 
 

table(Plate37_norm$group)


# No. of negative control reads on the plate
Plate37_nc_reads <- Plate37 %>% 
  group_by(sampleID, nc_allele) %>% 
  summarize(mean_reads = round(mean(reads)), reads = sum(reads)) %>%
  mutate(total_reads = sum(reads)) %>% mutate(prop_nc_reads = reads/total_reads) %>%
  mutate(group = case_when(
    grepl("Negative", sampleID) ~ "Negative",
    grepl("Positive", sampleID) ~ "Positive",
    TRUE ~ "Sample"))

Plate37_nc_reads %>%
  filter(nc_allele == 1) %>%
  ggplot(aes(x=sampleID, y=log10(mean_reads), color = group)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = 6)) +
  geom_hline(aes(yintercept = 2), color = "red", linetype='dashed') 

# Sun Feb  4 16:28:35 2024 ------------------------------
#Clearing the environment and console
rm(list = ls())
cat("\14")
#load the necessary library
library(tidyverse)
library(dplyr)
library(ggplot2)

#need to pull in the sample sheet for run 55 so we have well positions
sample_sheet_Run55 <- read_csv("C:/Users/USER/OneDrive/Desktop/Bienvenu Files/IMMRSE-U/Phase2/QC/Samplesheet-Run-55.csv")
sample_sheet_Run55$sampleID <- sample_sheet_Run55$Sample_ID
sample_sheet_Run55 <- sample_sheet_Run55 %>% select(-Sample_ID)
start_row <- 97
sample_sheet_Run55 <- sample_sheet_Run55[(start_row:nrow(sample_sheet_Run55)), ]

#pull in run52 data 
run55 <- read_delim("C:/Users/USER/OneDrive/Desktop/Bienvenu Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-02-02-IM-RUN-55/Results-2024-02-02-IM-RUN-55/allele_data.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE) %>% 
  mutate(sampleID = word(sampleID, 1, sep = "_")) 

#join in sample sheet info so that you have plate position 
Run55_samplesheet_merged <- run55 %>% left_join(sample_sheet_Run55, by = "sampleID")

#if no plate position, then it was the other pool on the run, so filter it out
Run55_samplesheet_merged <- Run55_samplesheet_merged %>% filter(!is.na(Position))

#check number of samples
check <- Run55_samplesheet_merged %>% select(sampleID) %>% unique()

#Extract the allele data for that contaminated NC
bad_nc_allele <- Run55_samplesheet_merged %>% filter(sampleID == "Negative-Control-8") %>% pull(allele) 

# Mark the alleles presented in that contaminated NC
Run55_samplesheet_merged <- Run55_samplesheet_merged %>% 
  mutate(nc_allele = ifelse(allele %in% bad_nc_allele, 1, 0))

# Proportion of alleles that are present in contaminated NC 
Plate55_norm <- Run55_samplesheet_merged %>% group_by(sampleID, Position, nc_allele) %>% summarize(count = n())
Plate55_norm <- Plate55_norm %>% group_by(sampleID) %>% mutate(total_allele = sum(count))
Plate55_norm <- Plate55_norm %>% mutate(prop_nc_allele = count/total_allele)
Plate55_norm %<>% 
  separate(col = Position, into = c("Row", "Col"), sep = 1) %>%  #separate columns and rows 
  mutate(Col = as.numeric(Col)) %>% 
  mutate(group = case_when(
    grepl("Negative", sampleID) ~ "Negative",
    grepl("Positive", sampleID) ~ "Positive",
    TRUE ~ "Sample"))

Plate55_norm$Row <- as.factor(Plate55_norm$Row)
Plate55_norm$Col <- as.factor(Plate55_norm$Col)

Plate55_norm <- Plate55_norm %>% 
  mutate(group = factor(group, levels = c("Negative", "Positive", "Sample"))) #level the factor

my_colors <- c("Negative" = "red", "Positive" = "black", "Sample" = "white")

heat_lab <- seq(0, 1, 0.1)
heat_break <- seq(0, 1, 0.1)

#I don't know why the black outline of the positive control is not appearing
Plate55_norm %>% 
  ggplot(aes(x = Col, y = Row, text = paste(sampleID))) +
  geom_tile(aes(fill = prop_nc_allele, color = group), size = 2) +
  scale_colour_manual(values = my_colors) + 
  scale_size_manual("group", values = c(0, 1, 1)) +
  scale_fill_gradient(low = "blue", high = "yellow", trans = "log",
                      breaks = heat_break, labels = heat_lab) +
  ylim(rev(levels(Plate55_norm$Row))) + 
  # theme_bw() +
  theme(axis.title = element_blank()) 


# table(Plate55_norm$group)
# 
# 
# # No. of negative control reads on the plate
# Plate55_nc_reads <- Run55_samplesheet_merged %>% 
#   group_by(sampleID, nc_allele) %>% 
#   summarize(mean_reads = round(mean(reads)), reads = sum(reads)) %>%
#   mutate(total_reads = sum(reads)) %>% mutate(prop_nc_reads = reads/total_reads) %>%
#   mutate(group = case_when(
#     grepl("Negative", sampleID) ~ "Negative",
#     grepl("Positive", sampleID) ~ "Positive",
#     TRUE ~ "Sample"))
# 
# Plate55_nc_reads %>%
#   filter(nc_allele == 1) %>%
#   ggplot(aes(x=sampleID, y=log10(mean_reads), color = group)) + 
#   geom_point() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = 6)) +
#   geom_hline(aes(yintercept = 2), color = "red", linetype='dashed') 

# Thu Apr  4 11:48:59 2024 ------------------------------
#Clearing the environment and console
rm(list = ls())
cat("\14")
#load the necessary library
library(tidyverse)
library(dplyr)
library(ggplot2)

#need to pull in the sample sheet for run 55 so we have well positions
sample_sheet_Run17 <- read_csv("C:/Users/USER/OneDrive/Desktop/Bienvenu Files/IMMRSE-U/R_codes/QC/Run17_SampleSheet.csv")
sample_sheet_Run17$sampleID <- sample_sheet_Run17$`Sample-ID`
sample_sheet_Run17 <- sample_sheet_Run17 %>% select(-`Sample-ID`)

#pull in run52 data 
run17 <- read_delim("C:/Users/USER/OneDrive/Desktop/Bienvenu Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2023-07-27-AG-IM-RUN-17/Results-2023-07-27-AG-IM-RUN-17/allele_data.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE) %>% 
  mutate(sampleID = word(sampleID, 1, sep = "_")) 

#join in sample sheet info so that you have plate position 
Run17_samplesheet_merged <- run17 %>% left_join(sample_sheet_Run17, by = "sampleID")

#if no plate position, then it was the other pool on the run, so filter it out
Run17_samplesheet_merged <- Run17_samplesheet_merged %>% filter(!is.na(Position))

#check number of samples
check <- Run17_samplesheet_merged %>% select(Position) %>% unique()

#Extract the allele data for that contaminated NC
bad_nc_allele <- Run17_samplesheet_merged %>% filter(sampleID == "Negative-control-8") %>% pull(allele) 

# Mark the alleles presented in that contaminated NC
Run17_samplesheet_merged <- Run17_samplesheet_merged %>% 
  mutate(nc_allele = ifelse(allele %in% bad_nc_allele, 1, 0))

# Proportion of alleles that are present in contaminated NC 
Plate17_norm <- Run17_samplesheet_merged %>% group_by(sampleID, Position, nc_allele) %>% summarize(count = n())
Plate17_norm <- Plate17_norm %>% group_by(sampleID) %>% mutate(total_allele = sum(count))
Plate17_norm <- Plate17_norm %>% mutate(prop_nc_allele = count/total_allele)
Plate17_norm %<>% 
  separate(col = Position, into = c("Row", "Col"), sep = 1) %>%  #separate columns and rows 
  mutate(Col = as.numeric(Col)) %>% 
  mutate(group = case_when(
    grepl("Negative", sampleID) ~ "Negative",
    grepl("Positive", sampleID) ~ "Positive",
    TRUE ~ "Sample"))

Plate17_norm$Row <- as.factor(Plate17_norm$Row)
Plate17_norm$Col <- as.factor(Plate17_norm$Col)

Plate17_norm <- Plate17_norm %>% 
  mutate(group = factor(group, levels = c("Negative", "Positive", "Sample"))) #level the factor

my_colors <- c("Negative" = "red", "Positive" = "black", "Sample" = "white")

heat_lab <- seq(0, 1, 0.1)
heat_break <- seq(0, 1, 0.1)

#I don't know why the black outline of the positive control is not appearing
Plate17_norm %>% 
  ggplot(aes(x = Col, y = Row, text = paste(sampleID))) +
  geom_tile(aes(fill = prop_nc_allele, color = group), size = 2) +
  scale_colour_manual(values = my_colors) + 
  scale_size_manual("group", values = c(0, 1, 1)) +
  scale_fill_gradient(low = "blue", high = "yellow", trans = "log",
                      breaks = heat_break, labels = heat_lab) +
  ylim(rev(levels(Plate17_norm$Row))) + 
  # theme_bw() +
  theme(axis.title = element_blank()) 

