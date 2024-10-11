# Mon Apr 15 07:45:40 2024 ------------------------------

#################################################################################
################################       #########################################
############################### Run01 #########################################
##############################       #########################################
#############################################################################
# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run1 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2023-05-16-IM-RUN-1/Results-2023-05-16-IM-RUN-1/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run1_wide <- run1 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run1_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run1_wide$Input, na.rm = TRUE) - sum(run1_wide$`No Dimers`, na.rm = TRUE)) / sum(run1_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run1_wide_Per_DADA2_QC <- run1_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run1_wide$OutputDada2, na.rm = TRUE)/sum(run1_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide file and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run1_no_controls <- run1_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run1_no_controls$`No Dimers`, na.rm = TRUE)
sd(run1_no_controls$`No Dimers`, na.rm = TRUE)


#################################################################################
################################       #########################################
############################### Run47 #########################################
##############################       #########################################
#############################################################################

# Clean the environmemt
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run47 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2023-12-26-IM-RUN-47/Results-2023-12-26-IM-RUN-47/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run47_wide <- run47%>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  )

# Remove controls
run47_no_controls <- run47_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID))

# Make valid data and sum up the reads
valid_data <- subset(run47_wide, run47_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)

# Calculate median and standard deviation of reads
median(run47_no_controls$`No Dimers`, na.rm = TRUE)
sd(run47_no_controls$`No Dimers`, na.rm = TRUE)


# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run47_wide, run47_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

#################################################################################
################################       #########################################
############################### Run48 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run48 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2023-12-28-IM-RUN-48/Results-2023-12-28-IM-RUN-48/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run48_wide <- run48%>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  )

# Remove controls
run48_no_controls <- run48_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID))

# Make valid data and sum up the reads
valid_data <- subset(run48_wide, run48_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)

# Calculate median and standard deviation of reads
median(run48_no_controls$`No Dimers`, na.rm = TRUE)
sd(run48_no_controls$`No Dimers`, na.rm = TRUE)


# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run48_wide, run48_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Wed Jan  3 11:19:29 2024 ------------------------------

#################################################################################
################################       #########################################
############################### Run49 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run49 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2023-12-30-IM-RUN-49/Results-2023-12-30-IM-RUN-49/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run49_wide <- run49%>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  )

# Remove controls
run49_no_controls <- run49_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID))

# Make valid data and sum up the reads
valid_data <- subset(run49_wide, run49_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)

# Calculate median and standard deviation of reads
median(run49_no_controls$`No Dimers`, na.rm = TRUE)
sd(run49_no_controls$`No Dimers`, na.rm = TRUE)


# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run49_wide, run49_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Wed Jan  3 11:45:34 2024 ------------------------------

#################################################################################
################################       #########################################
############################### Run50 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run50 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2023-12-31-IM-RUN-50/Results-2023-12-31-IM-RUN-50/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run50_wide <- run50%>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  )

# Remove controls
run50_no_controls <- run50_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID))

# Make valid data and sum up the reads
valid_data <- subset(run50_wide, run50_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)

# Calculate median and standard deviation of reads
run50_no_controls_as_numeric <- as.numeric(run50_no_controls$`No Dimers`) #(This was introduced because there is 'null' in "No Dimers" column) 
median(run50_no_controls_as_numeric, na.rm = TRUE)
sd(run50_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run50_wide, run50_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Wed Jan  3 13:25:25 2024 ------------------------------

#################################################################################
################################       #########################################
############################### Run51 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run51 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-01-01-IM-RUN-51/Results-2024-01-01-IM-RUN-51/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run51_wide <- run51%>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  )

# Remove controls
run51_no_controls <- run51_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("Control", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run51_wide, run51_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)

# Calculate median and standard deviation of reads
median(as.numeric(run51_no_controls$`No Dimers`), na.rm = TRUE)
sd(run51_no_controls$`No Dimers`, na.rm = TRUE)


# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run51_wide, run51_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Fri Jan  5 22:37:07 2024 ------------------------------

#################################################################################
################################       #########################################
############################### Run52 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run52 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-01-03-IM-RUN-52/Results-2024-01-03-IM-RUN-52/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run52_wide <- run52%>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  )

# Remove controls
run52_no_controls <- run52_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run52_wide, run52_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)

# Calculate median and standard deviation of reads
median(as.numeric(run52_no_controls$`No Dimers`), na.rm = TRUE)
sd(run52_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run52_wide, run52_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Thu Mar  7 07:10:25 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run57 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run57 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-02-12-IM-RUN-57/Results-2024-02-12-IM-RUN-57/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run57_wide <- run57 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  )

# Remove controls
run57_no_controls <- run57_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run57_wide, run57_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run61_wide$OutputDada2, na.rm = TRUE)/sum(run61_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run57_no_controls$`No Dimers`), na.rm = TRUE)
sd(run57_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run57_wide, run57_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)



# Mon Mar  4 13:18:36 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run60 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run60 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-02-22-IM-RUN-60/Results-2024-02-22-IM-RUN-60/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run60_wide <- run60%>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% mutate_at(vars(2:6), as.numeric)

# Remove controls
run60_no_controls <- run60_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run60_wide, run60_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run60_wide$OutputDada2, na.rm = TRUE)/sum(run60_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run61_no_controls$`No Dimers`), na.rm = TRUE)
sd(run61_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run61_wide, run61_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Mon Mar  4 13:49:20 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run61 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run61 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-03-01-IM-RUN-61/Results-2024-03-01-IM-RUN-61/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run61_wide <- run61%>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% mutate_at(vars(2:6), as.numeric)

# Remove controls
run61_no_controls <- run61_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run61_wide, run61_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run61_wide$OutputDada2, na.rm = TRUE)/sum(run61_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run61_no_controls$`No Dimers`), na.rm = TRUE)
sd(run61_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run61_wide, run61_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Mon Mar  4 13:57:59 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run62 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run62 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-02-28-IM-RUN-62/Results-2024-02-28-IM-RUN-62/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run62_wide <- run62%>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% mutate_at(vars(2:6), as.numeric)

# Remove controls
run62_no_controls <- run62_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run62_wide, run62_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run62_wide$OutputDada2, na.rm = TRUE)/sum(run62_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run62_no_controls$`No Dimers`), na.rm = TRUE)
sd(run61_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run62_wide, run62_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Mon Mar  4 14:00:08 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run63 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run63 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-03-03-IM-RUN-63/Results-2024-03-03-IM-RUN-63/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run63_wide <- run63%>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% mutate_at(vars(2:6), as.numeric)

# Remove controls
run63_no_controls <- run63_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run63_wide, run63_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run63_wide$OutputDada2, na.rm = TRUE)/sum(run63_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run63_no_controls$`No Dimers`), na.rm = TRUE)
sd(run63_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run63_wide, run63_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Thu Mar  7 06:59:00 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run64 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run64 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-03-04-IM-RUN-64/Results-2024-03-04-IM-RUN-64/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run64_wide <- run64 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% mutate_at(vars(2:6), as.numeric)

# Remove controls
run64_no_controls <- run64_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run64_wide, run64_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run64_wide$OutputDada2, na.rm = TRUE)/sum(run64_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run64_no_controls$`No Dimers`), na.rm = TRUE)
sd(run64_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run64_wide, run64_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Sat Mar 16 22:39:19 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run65 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run65 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-03-11-IM-RUN-65//Results-2024-03-11-IM-RUN-65/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run65_wide <- run65 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% mutate_at(vars(2:6), as.numeric)

# Remove controls
run65_no_controls <- run65_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run65_wide, run65_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run65_wide$OutputDada2, na.rm = TRUE)/sum(run65_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run65_no_controls$`No Dimers`), na.rm = TRUE)
sd(run65_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run65_wide, run65_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


# Sat Mar 16 21:48:29 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run66 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run66 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-03-15-IM-RUN-66/Results-2024-03-15-IM-RUN-66/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run66_wide <- run66 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% mutate_at(vars(2:6), as.numeric)

# Remove controls
run66_no_controls <- run66_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run66_wide, run66_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run66_wide$OutputDada2, na.rm = TRUE)/sum(run66_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)


# Calculate median and standard deviation of reads
median(as.numeric(run66_no_controls$`No Dimers`), na.rm = TRUE)
sd(run66_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run66_wide, run66_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Mon Mar 18 19:31:48 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run67 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run67 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-03-16-IM-RUN-67/Results-2024-03-16-IM-RUN-67/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run67_wide <- run67 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% mutate_at(vars(2:6), as.numeric)

# Remove controls
run67_no_controls <- run67_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run67_wide, run67_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run67_wide$OutputDada2, na.rm = TRUE)/sum(run67_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run67_no_controls$`No Dimers`), na.rm = TRUE)
sd(run67_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run67_wide, run67_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Fri Mar 22 07:48:06 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run68 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run68 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-03-18-IM-RUN-68/Results-2024-03-18-IM-RUN-68/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run68_wide <- run68 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% mutate_at(vars(2:6), as.numeric)

# Remove controls
run68_no_controls <- run68_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run68_wide, run68_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run68_wide$OutputDada2, na.rm = TRUE)/sum(run68_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run68_no_controls$`No Dimers`), na.rm = TRUE)
sd(run68_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run68_wide, run68_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Fri Mar 22 07:57:52 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run69 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run69 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-03-19-IM-RUN-69/Results-2024-03-19-IM-RUN-69/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run69_wide <- run69 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

run69_wide_Per_DADA2_QC <- run69_wide %>% 
  mutate(`%age_passing_DADA2_QC` = run69_wide$OutputPostprocessing/run69_wide$`No Dimers`*100)

write.csv(run69_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_miseq_M08585_Run69.csv", row.names = FALSE)

# Remove controls
run69_no_controls <- run69_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run69_wide, run69_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run69_wide$OutputDada2, na.rm = TRUE)/sum(run69_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run69_no_controls$`No Dimers`), na.rm = TRUE)
sd(run69_no_controls$`No Dimers`, na.rm = TRUE)
sum(as.numeric(run69_wide$OutputDada2))/sum(as.numeric(run69_wide$`No Dimers`))*100

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run69_wide, run69_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Sat Mar 30 21:18:30 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run70 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file  
run70 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-03-15-IM-RUN-70/Results-2024-03-15-IM-RUN-70/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run70_wide <- run70 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  )

# Remove controls
run70_no_controls <- run70_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run70_wide, run70_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)


# Calculate median and standard deviation of reads
median(as.numeric(run70_no_controls$`No Dimers`), na.rm = TRUE)
sd(run70_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run70_wide, run70_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Sat Mar 30 21:27:50 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run71 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file  
run71 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-03-27-IM-RUN-71/Results-2024-03-27-IM-RUN-71/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run71_wide <- run71 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% mutate_at(vars(2:6), as.numeric)

# Remove controls
run71_no_controls <- run71_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run71_wide, run71_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run71_wide$OutputDada2, na.rm = TRUE)/sum(run71_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)


# Calculate median and standard deviation of reads
median(as.numeric(run71_no_controls$`No Dimers`), na.rm = TRUE)
sd(run71_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run71_wide, run71_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Thu Apr  4 12:50:56 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run17 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file  
run17 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2023-07-27-AG-IM-RUN-17/Results-2023-07-27-AG-IM-RUN-17/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run17_wide <- run17 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% mutate_at(vars(2:6), as.numeric)

# Remove controls
run17_no_controls <- run17_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run17_wide, run17_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run17_wide$OutputDada2, na.rm = TRUE)/sum(run17_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run71_no_controls$`No Dimers`), na.rm = TRUE)
sd(run71_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run71_wide, run71_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


# Thu Apr  4 15:58:58 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run73 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file  
run73 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-04-02-IM-RUN-73/Results-2024-04-02-IM-RUN-73/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run73_wide <- run73 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% mutate_at(vars(2:6), as.numeric)

# Remove controls
run73_no_controls <- run73_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run73_wide, run73_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run73_wide$OutputDada2, na.rm = TRUE)/sum(run73_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run73_no_controls$`No Dimers`), na.rm = TRUE)
sd(run73_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run73_wide, run73_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Sun Apr  7 23:55:18 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run74 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file  
run74 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-04-05-IM-RUN-74/Results-2024-04-05-IM-RUN-74/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run74_wide <- run74 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate(Per_passing = (OutputDada2/`No Dimers`)*100)

# Remove controls
run74_no_controls <- run74_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run74_wide, run74_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run74_wide$OutputDada2, na.rm = TRUE)/sum(run74_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)


# Calculate median and standard deviation of reads
median(as.numeric(run74_no_controls$`No Dimers`), na.rm = TRUE)
sd(run74_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run74_wide, run74_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Thu Apr 11 14:23:37 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run75 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run75 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-04-05-IM-RUN-75/Results-2024-04-05-IM-RUN-75/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run75_wide <- run75 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

run75_wide_Per_DADA2_QC <- run75_wide %>% 
  mutate(`%age_passing_DADA2_QC` = run75_wide$OutputPostprocessing/run75_wide$`No Dimers`*100)

write.csv(run69_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_miseq_M08585_Run69.csv", row.names = FALSE)

# Remove controls
run75_no_controls <- run75_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run75_wide, run75_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run75_wide$OutputDada2, na.rm = TRUE)/sum(run75_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run75_no_controls$`No Dimers`), na.rm = TRUE)
sd(run75_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run75_wide, run75_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Thu Apr 11 14:38:04 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run76 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run76 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-04-05-IM-RUN-76/Results-2024-04-05-IM-RUN-76/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run76_wide <- run76 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

run76_wide_Per_DADA2_QC <- run76_wide %>% 
  mutate(`%age_passing_DADA2_QC` = run76_wide$OutputPostprocessing/run76_wide$`No Dimers`*100)

# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_miseq_M08585_Run69.csv", row.names = FALSE)

# Remove controls
run76_no_controls <- run76_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run76_wide, run76_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run76_wide$OutputDada2, na.rm = TRUE)/sum(run76_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run76_no_controls$`No Dimers`), na.rm = TRUE)
sd(run76_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run75_wide, run75_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Thu Apr 11 14:48:50 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run78 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run78 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-04-08-IM-RUN-78/Results-2024-04-08-IM-RUN-78/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run78_wide <- run78 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

run78_wide_Per_DADA2_QC <- run78_wide %>% 
  mutate(`%age_passing_DADA2_QC` = run78_wide$OutputPostprocessing/run78_wide$`No Dimers`*100)

# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_miseq_M08585_Run69.csv", row.names = FALSE)

# Remove controls
run78_no_controls <- run78_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run78_wide, run78_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run78_wide$OutputDada2, na.rm = TRUE)/sum(run78_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run78_no_controls$`No Dimers`), na.rm = TRUE)
sd(run78_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run78_wide, run78_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Thu Apr 11 15:18:25 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run79 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run79 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-04-07-IM-RUN-79/Results-2024-04-07-IM-RUN-79/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run79_wide <- run79 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

run79_wide_Per_DADA2_QC <- run79_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)

# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_miseq_M08585_Run69.csv", row.names = FALSE)

# Remove controls
run79_no_controls <- run79_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Make valid data and sum up the reads
valid_data <- subset(run79_wide, run79_wide$`No Dimers` != "null")
sum(as.numeric(valid_data$`No Dimers`), na.rm = TRUE)
`%age_passing_DADA2_QC` <- sum(run79_wide$OutputDada2, na.rm = TRUE)/sum(run79_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Calculate median and standard deviation of reads
median(as.numeric(run79_no_controls$`No Dimers`), na.rm = TRUE)
sd(run79_no_controls$`No Dimers`, na.rm = TRUE)

# Filter out rows with "null" in the `No Dimers` column
valid_data <- subset(run79_wide, run79_wide$`No Dimers` != "null")

# Convert columns to numeric
valid_data$Input <- as.numeric(valid_data$Input)
valid_data$`No Dimers` <- as.numeric(valid_data$`No Dimers`)

# Calculate percentage_dimmers
percentage_dimmers <- ((sum(valid_data$Input, na.rm = TRUE) - sum(valid_data$`No Dimers`, na.rm = TRUE)) / sum(valid_data$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)

# Mon Apr 15 07:46:08 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run80 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run80 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-04-13-IM-RUN-80/Results-2024-04-13-IM-RUN-80/sample_coverage.txt", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)

# Change sample coverage file orientation to wide
run80_wide <- run80 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run80_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run80_wide$Input, na.rm = TRUE) - sum(run80_wide$`No Dimers`, na.rm = TRUE)) / sum(run80_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run80_wide_Per_DADA2_QC <- run80_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run80_wide$OutputDada2, na.rm = TRUE)/sum(run80_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run80_no_controls <- run80_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run80_no_controls$`No Dimers`, na.rm = TRUE)
sd(run80_no_controls$`No Dimers`, na.rm = TRUE)

# Wed Apr 17 14:35:52 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run81 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run81 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-04-16-IM-RUN-81/Results-2024-04-16-IM-RUN-81/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run81_wide <- run81 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run81_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run81_wide$Input, na.rm = TRUE) - sum(run81_wide$`No Dimers`, na.rm = TRUE)) / sum(run81_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run81_wide_Per_DADA2_QC <- run81_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run81_wide$OutputDada2, na.rm = TRUE)/sum(run81_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run81_no_controls <- run81_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run81_no_controls$`No Dimers`, na.rm = TRUE)
sd(run81_no_controls$`No Dimers`, na.rm = TRUE)

# Fri Apr 19 00:03:05 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run82 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run82 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-04-17-IM-RUN-82/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run82_wide <- run82 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run82_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run82_wide$Input, na.rm = TRUE) - sum(run82_wide$`No Dimers`, na.rm = TRUE)) / sum(run82_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run82_wide_Per_DADA2_QC <- run82_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run82_wide$OutputDada2, na.rm = TRUE)/sum(run82_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run82_no_controls <- run82_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run82_no_controls$`No Dimers`, na.rm = TRUE)
sd(run82_no_controls$`No Dimers`, na.rm = TRUE)

# Tue Apr 23 19:04:50 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run83 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run83 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-04-18-IM-RUN-83/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run83_wide <- run83 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run83_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run83_wide$Input, na.rm = TRUE) - sum(run83_wide$`No Dimers`, na.rm = TRUE)) / sum(run83_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run83_wide_Per_DADA2_QC <- run83_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run83_wide$OutputDada2, na.rm = TRUE)/sum(run83_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run83_no_controls <- run83_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run83_no_controls$`No Dimers`, na.rm = TRUE)
sd(run83_no_controls$`No Dimers`, na.rm = TRUE)

# Tue Apr 23 20:03:38 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run84 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run84 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-04-19-IM-RUN-84/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run84_wide <- run84 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run84_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run84_wide$Input, na.rm = TRUE) - sum(run84_wide$`No Dimers`, na.rm = TRUE)) / sum(run84_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run84_wide_Per_DADA2_QC <- run84_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run84_wide$OutputDada2, na.rm = TRUE)/sum(run84_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run84_no_controls <- run84_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run84_no_controls$`No Dimers`, na.rm = TRUE)
sd(run84_no_controls$`No Dimers`, na.rm = TRUE)

# Tue Apr 23 21:23:40 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run85 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run85 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-04-22-IM-RUN-85/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run85_wide <- run85 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run85_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run85_wide$Input, na.rm = TRUE) - sum(run85_wide$`No Dimers`, na.rm = TRUE)) / sum(run85_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run85_wide_Per_DADA2_QC <- run85_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run85_wide$OutputDada2, na.rm = TRUE)/sum(run85_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run85_no_controls <- run85_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run85_no_controls$`No Dimers`, na.rm = TRUE)
sd(run85_no_controls$`No Dimers`, na.rm = TRUE)

# Tue Apr 30 13:58:53 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run85 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run86 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-04-23-IM-RUN-86/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run86_wide <- run86 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run86_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run86_wide$Input, na.rm = TRUE) - sum(run86_wide$`No Dimers`, na.rm = TRUE)) / sum(run86_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run86_wide_Per_DADA2_QC <- run86_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run86_wide$OutputDada2, na.rm = TRUE)/sum(run86_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run86_no_controls <- run86_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run86_no_controls$`No Dimers`, na.rm = TRUE)
sd(run86_no_controls$`No Dimers`, na.rm = TRUE)

# Mon May  6 11:26:55 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run101 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run101 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-05-04-IM-RUN-101/sample_coverage.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Change sample coverage file orientation to wide
run101_wide <- run101 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run101_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run101_wide$Input, na.rm = TRUE) - sum(run101_wide$`No Dimers`, na.rm = TRUE)) / sum(run101_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run101_wide_Per_DADA2_QC <- run101_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run101_wide$OutputDada2, na.rm = TRUE)/sum(run101_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run101_no_controls <- run101_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run101_no_controls$`No Dimers`, na.rm = TRUE)
sd(run101_no_controls$`No Dimers`, na.rm = TRUE)

# Mon Jun 17 10:22:34 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run102 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run102 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-06-07-IM-RUN-102/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run102_wide <- run102 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run102_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run102_wide$Input, na.rm = TRUE) - sum(run102_wide$`No Dimers`, na.rm = TRUE)) / sum(run10_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run102_wide_Per_DADA2_QC <- run102_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run102_wide$OutputDada2, na.rm = TRUE)/sum(run102_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run102_no_controls <- run102_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run102_no_controls$`No Dimers`, na.rm = TRUE)
sd(run102_no_controls$`No Dimers`, na.rm = TRUE)

# Wed Jun 26 20:02:49 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run103 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
  run103 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-06-11-IM-RUN-103//sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run103_wide <- run103 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run103_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run103_wide$Input, na.rm = TRUE) - sum(run103_wide$`No Dimers`, na.rm = TRUE)) / sum(run103_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run103_wide_Per_DADA2_QC <- run103_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run103_wide$OutputDada2, na.rm = TRUE)/sum(run103_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run103_no_controls <- run103_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run103_no_controls$`No Dimers`, na.rm = TRUE)
sd(run103_no_controls$`No Dimers`, na.rm = TRUE)

# Wed Jun 26 20:06:51 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run104 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run104 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-06-19-IM-RUN-104/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run104_wide <- run104 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run104_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run104_wide$Input, na.rm = TRUE) - sum(run104_wide$`No Dimers`, na.rm = TRUE)) / sum(run104_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run104_wide_Per_DADA2_QC <- run104_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run104_wide$OutputDada2, na.rm = TRUE)/sum(run104_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run102_no_controls <- run102_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run102_no_controls$`No Dimers`, na.rm = TRUE)
sd(run102_no_controls$`No Dimers`, na.rm = TRUE)

# Mon Jul  1 14:38:21 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run105 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run105 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-06-30-IM-RUN-105/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run105_wide <- run105 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run105_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run105_wide$Input, na.rm = TRUE) - sum(run105_wide$`No Dimers`, na.rm = TRUE)) / sum(run105_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run105_wide_Per_DADA2_QC <- run105_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run105_wide$OutputDada2, na.rm = TRUE)/sum(run105_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run105_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-06-30-IM-RUN-105/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >=16)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run105_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE))
print(filtered_Negative_Controls)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run105_no_controls <- run105_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run105_no_controls$`No Dimers`, na.rm = TRUE)
sd(run105_no_controls$`No Dimers`, na.rm = TRUE)

# Wed Jul 24 00:06:53 2024 ------------------------------
#################################################################################
################################       #########################################
############################### Run106 #########################################
##############################       #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run106 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-07-01-IM-RUN-106/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run106_wide <- run106 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run106_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run106_wide$Input, na.rm = TRUE) - sum(run106_wide$`No Dimers`, na.rm = TRUE)) / sum(run106_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run106_wide_Per_DADA2_QC <- run106_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run106_wide$OutputDada2, na.rm = TRUE)/sum(run106_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run106_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-07-01-IM-RUN-106/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >1)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run106_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE))
print(filtered_Negative_Controls)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run106_no_controls <- run106_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run106_no_controls$`No Dimers`, na.rm = TRUE)
sd(run106_no_controls$`No Dimers`, na.rm = TRUE)

# Wed Jul 24 00:45:56 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## Run107 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run107 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-07-16-IM-RUN-107/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run107_wide <- run107 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run107_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run107_wide$Input, na.rm = TRUE) - sum(run107_wide$`No Dimers`, na.rm = TRUE)) / sum(run107_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run107_wide_Per_DADA2_QC <- run107_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run107_wide$OutputDada2, na.rm = TRUE)/sum(run107_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run107_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-07-16-IM-RUN-107/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >18)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run107_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE))
print(filtered_Negative_Controls)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run107_no_controls <- run107_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run107_no_controls$`No Dimers`, na.rm = TRUE)
sd(run107_no_controls$`No Dimers`, na.rm = TRUE)

# Wed Jul 24 01:16:11 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## Run108 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run108 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-07-18-IM-RUN-108/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run108_wide <- run108 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run108_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run108_wide$Input, na.rm = TRUE) - sum(run108_wide$`No Dimers`, na.rm = TRUE)) / sum(run108_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run108_wide_Per_DADA2_QC <- run108_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run108_wide$OutputDada2, na.rm = TRUE)/sum(run108_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run108_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-07-18-IM-RUN-108/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >22)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run108_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE))
print(filtered_Negative_Controls)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run108_no_controls <- run108_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run108_no_controls$`No Dimers`, na.rm = TRUE)
sd(run108_no_controls$`No Dimers`, na.rm = TRUE)

# Wed Jul 31 16:10:11 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## DPSP_001 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
DPSP_001 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/DPSP/DPSP-rawdata/Results-2024-07-27-DPSP-001/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
DPSP_001_wide <- DPSP_001 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(DPSP_001_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(DPSP_001_wide$Input, na.rm = TRUE) - sum(DPSP_001_wide$`No Dimers`, na.rm = TRUE)) / sum(DPSP_001_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
DPSP_001_wide_Per_DADA2_QC <- DPSP_001_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(DPSP_001_wide$OutputDada2, na.rm = TRUE)/sum(DPSP_001_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- DPSP_001_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/DPSP/DPSP-rawdata/Results-2024-07-27-DPSP-001/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "DPSP_001", sep = "-"), SampleID)) %>% 
  filter(Reads >31)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- DPSP_001_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE))
print(filtered_Negative_Controls)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
DPSP_001_no_controls <- DPSP_001_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(DPSP_001_no_controls$`No Dimers`, na.rm = TRUE)
sd(DPSP_001_no_controls$`No Dimers`, na.rm = TRUE)

# Thu Aug  8 17:22:27 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## DPSP_002 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
DPSP_002 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/DPSP/DPSP-rawdata/Results-2024-07-31-DPSP-002//sample_coverage.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

# Change sample coverage file orientation to wide
DPSP_002_wide <- DPSP_002 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(DPSP_002_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(DPSP_002_wide$Input, na.rm = TRUE) - sum(DPSP_002_wide$`No Dimers`, na.rm = TRUE)) / sum(DPSP_002_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
DPSP_002_wide_Per_DADA2_QC <- DPSP_002_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(DPSP_002_wide$OutputDada2, na.rm = TRUE)/sum(DPSP_002_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- DPSP_002_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/DPSP/DPSP-rawdata/Results-2024-07-31-DPSP-002//allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "DPSP_002", sep = "-"), SampleID)) %>% 
  filter(Reads >11)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- DPSP_002_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE))
print(filtered_Negative_Controls)


# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
DPSP_002_no_controls <- DPSP_002_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(DPSP_002_no_controls$`No Dimers`, na.rm = TRUE)
sd(DPSP_002_no_controls$`No Dimers`, na.rm = TRUE)

# Mon Aug 19 18:29:54 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## Run109 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run109 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-07-25-IM-RUN-109/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run109_wide <- run109 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run109_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run109_wide$Input, na.rm = TRUE) - sum(run109_wide$`No Dimers`, na.rm = TRUE)) / sum(run109_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run109_wide_Per_DADA2_QC <- run109_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run109_wide$OutputDada2, na.rm = TRUE)/sum(run109_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run109_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-08-07-IM-RUN-110/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >17)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run109_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE))
print(filtered_Negative_Controls)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run109_no_controls <- run109_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run109_no_controls$`No Dimers`, na.rm = TRUE)
sd(run109_no_controls$`No Dimers`, na.rm = TRUE)

# Mon Aug 19 18:55:19 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## Run110 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run110 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-08-07-IM-RUN-110/sample_coverage.txt", 
run110 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-08-07-IM-RUN-110/sample_coverage.txt", 
run110 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-08-07-IM-RUN-110/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run110_wide <- run110 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run110_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run110_wide$Input, na.rm = TRUE) - sum(run110_wide$`No Dimers`, na.rm = TRUE)) / sum(run110_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run110_wide_Per_DADA2_QC <- run110_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run110_wide$OutputDada2, na.rm = TRUE)/sum(run110_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run110_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-08-07-IM-RUN-110/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >17)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run110_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE))
print(filtered_Negative_Controls)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run110_no_controls <- run110_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run110_no_controls$`No Dimers`, na.rm = TRUE)
sd(run110_no_controls$`No Dimers`, na.rm = TRUE)

# Mon Aug 19 19:06:25 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## Run111 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run111 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-08-09-IM-RUN-111/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run111_wide <- run111 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run111_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run111_wide$Input, na.rm = TRUE) - sum(run111_wide$`No Dimers`, na.rm = TRUE)) / sum(run111_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run111_wide_Per_DADA2_QC <- run111_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run111_wide$OutputDada2, na.rm = TRUE)/sum(run111_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run111_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-08-09-IM-RUN-111/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >12)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run111_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE))
print(filtered_Negative_Controls)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run111_no_controls <- run111_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run111_no_controls$`No Dimers`, na.rm = TRUE)
sd(run111_no_controls$`No Dimers`, na.rm = TRUE)

# Mon Aug 19 19:21:20 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## Run112 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run112 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-08-12-IM-RUN-112/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run112_wide <- run112 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run112_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run112_wide$Input, na.rm = TRUE) - sum(run112_wide$`No Dimers`, na.rm = TRUE)) / sum(run112_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run112_wide_Per_DADA2_QC <- run112_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run112_wide$OutputDada2, na.rm = TRUE)/sum(run112_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run112_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-08-12-IM-RUN-112//allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >14)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run112_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE))
print(filtered_Negative_Controls)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run112_no_controls <- run112_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run112_no_controls$`No Dimers`, na.rm = TRUE)
sd(run112_no_controls$`No Dimers`, na.rm = TRUE)

# Mon Aug 19 19:29:34 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## Run113 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run113 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-08-14-IM-RUN-113/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run113_wide <- run113 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run113_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run113_wide$Input, na.rm = TRUE) - sum(run113_wide$`No Dimers`, na.rm = TRUE)) / sum(run113_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run113_wide_Per_DADA2_QC <- run113_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run113_wide$OutputDada2, na.rm = TRUE)/sum(run113_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run113_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-08-14-IM-RUN-113/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >26)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run113_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE))
print(filtered_Negative_Controls)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run113_no_controls <- run113_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run113_no_controls$`No Dimers`, na.rm = TRUE)
sd(run113_no_controls$`No Dimers`, na.rm = TRUE)

# Thu Aug 22 08:52:53 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## Run114 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run114 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-08-20-IM-RUN-114/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run114_wide <- run114 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run114_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run114_wide$Input, na.rm = TRUE) - sum(run114_wide$`No Dimers`, na.rm = TRUE)) / sum(run114_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run114_wide_Per_DADA2_QC <- run114_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run114_wide$OutputDada2, na.rm = TRUE)/sum(run114_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run114_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-08-20-IM-RUN-114/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >12)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run114_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE) | grepl("Neagtive", SampleID, ignore.case = TRUE))
head(filtered_Negative_Controls, n = 10)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run114_no_controls <- run114_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run114_no_controls$`No Dimers`, na.rm = TRUE)
sd(run114_no_controls$`No Dimers`, na.rm = TRUE)

# Mon Sep 30 13:06:02 2024 ------------------------------
# Thu Aug 22 08:52:53 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## Run116 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run116 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-09-02-IM-RUN-116/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run116_wide <- run116 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

run116_wide <- run116_wide %>%
  mutate(SampleID = case_when(
    grepl("1K-Positive-Control-8070942023-1_S35_L001", SampleID) ~ gsub("1K-Positive-Control-8070942023-1_S35_L001", "NTC1", SampleID),  # Replace "1K" with "NTC1"
    grepl("Negative-Control-1_S10_L001", SampleID) ~ gsub("Negative-Control-1_S10_L001", "1K-Positive-Control-8070942023-1_S35_L001", SampleID),  # Replace "NTC1" with "1K"
    TRUE ~ SampleID  # Leave other SampleIDs unchanged
  )) %>% 
  mutate(SampleID = gsub("NTC1", "Negative-Control-1_S10_L001", SampleID))  # Replace the placeholder with "NTC1"

# Sum up the reads without dimmers for the whole run
sum(run116_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run116_wide$Input, na.rm = TRUE) - sum(run116_wide$`No Dimers`, na.rm = TRUE)) / sum(run116_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run116_wide_Per_DADA2_QC <- run116_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run116_wide$OutputDada2, na.rm = TRUE)/sum(run116_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run116_wide %>%
  filter(grepl("Pos", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-09-02-IM-RUN-116/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >5)

df <- df %>%
  mutate(SampleID = case_when(
    grepl("1K-Positive-Control-8070942023-1", SampleID) ~ gsub("1K-Positive-Control-8070942023-1", "NTC1", SampleID),  # Replace "1K" with "NTC1"
    grepl("Negative-Control-1", SampleID) ~ gsub("Negative-Control-1", "1K-Positive-Control-8070942023-1", SampleID),  # Replace "NTC1" with "1K"
    TRUE ~ SampleID  # Leave other SampleIDs unchanged
  )) %>% 
  mutate(SampleID = gsub("NTC1", "Negative-Control-1", SampleID))  # Replace the placeholder with "NTC1"


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run116_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE) | grepl("Neagtive", SampleID, ignore.case = TRUE))
head(filtered_Negative_Controls, n = 10)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run116_no_controls <- run116_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run116_no_controls$`No Dimers`, na.rm = TRUE)
sd(run116_no_controls$`No Dimers`, na.rm = TRUE)

# Mon Sep 30 14:23:07 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## Run117 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run117 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-09-09-IM-RUN-117/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run117_wide <- run117 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run117_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run117_wide$Input, na.rm = TRUE) - sum(run117_wide$`No Dimers`, na.rm = TRUE)) / sum(run117_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run117_wide_Per_DADA2_QC <- run117_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run117_wide$OutputDada2, na.rm = TRUE)/sum(run117_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run117_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-09-09-IM-RUN-117/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >39)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run117_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE) | grepl("Neagtive", SampleID, ignore.case = TRUE))
head(filtered_Negative_Controls, n = 10)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run117_no_controls <- run117_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run117_no_controls$`No Dimers`, na.rm = TRUE)
sd(run117_no_controls$`No Dimers`, na.rm = TRUE)

#################################################################################
###############################        #########################################
############################## Run118 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run118 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-09-13-IM-RUN-118/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run118_wide <- run118 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run118_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run118_wide$Input, na.rm = TRUE) - sum(run118_wide$`No Dimers`, na.rm = TRUE)) / sum(run118_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run118_wide_Per_DADA2_QC <- run118_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run118_wide$OutputDada2, na.rm = TRUE)/sum(run118_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run118_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-09-13-IM-RUN-118/allele_data.txt",
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-09-13-IM-RUN-118/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >8)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run118_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE) | grepl("Neagtive", SampleID, ignore.case = TRUE))
head(filtered_Negative_Controls, n = 10)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run118_no_controls <- run118_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run118_no_controls$`No Dimers`, na.rm = TRUE)
sd(run118_no_controls$`No Dimers`, na.rm = TRUE)

# Mon Sep 30 14:59:59 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## run119 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run119 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-09-14-IM-RUN-119/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run119_wide <- run119 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run119_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run119_wide$Input, na.rm = TRUE) - sum(run119_wide$`No Dimers`, na.rm = TRUE)) / sum(run119_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run119_wide_Per_DADA2_QC <- run119_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run119_wide$OutputDada2, na.rm = TRUE)/sum(run119_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run119_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-09-14-IM-RUN-119/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >17)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run119_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE) | grepl("Neagtive", SampleID, ignore.case = TRUE))
head(filtered_Negative_Controls, n = 10)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run119_no_controls <- run119_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run119_no_controls$`No Dimers`, na.rm = TRUE)
sd(run119_no_controls$`No Dimers`, na.rm = TRUE)

# Mon Sep 30 15:06:14 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## run120 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run120 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-09-16-IM-RUN-120/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run120_wide <- run120 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run120_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run120_wide$Input, na.rm = TRUE) - sum(run120_wide$`No Dimers`, na.rm = TRUE)) / sum(run120_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run120_wide_Per_DADA2_QC <- run120_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run120_wide$OutputDada2, na.rm = TRUE)/sum(run120_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run120_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-09-16-IM-RUN-120/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >0)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run120_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE) | grepl("Neagtive", SampleID, ignore.case = TRUE))
head(filtered_Negative_Controls, n = 10)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run120_no_controls <- run120_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run120_no_controls$`No Dimers`, na.rm = TRUE)
sd(run120_no_controls$`No Dimers`, na.rm = TRUE)

# Mon Sep 30 15:14:18 2024 ------------------------------
#################################################################################
###############################        #########################################
############################## run121 #########################################
#############################        #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
run121 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-09-17-IM-RUN-121/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
run121_wide <- run121 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(run121_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(run121_wide$Input, na.rm = TRUE) - sum(run121_wide$`No Dimers`, na.rm = TRUE)) / sum(run121_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
run121_wide_Per_DADA2_QC <- run121_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(run121_wide$OutputDada2, na.rm = TRUE)/sum(run121_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- run121_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/IMMRSE-U raw data/Results-2024-09-17-IM-RUN-121/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >4)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- run121_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE) | grepl("Neagtive", SampleID, ignore.case = TRUE))
head(filtered_Negative_Controls, n = 10)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
run121_no_controls <- run121_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(run121_no_controls$`No Dimers`, na.rm = TRUE)
sd(run121_no_controls$`No Dimers`, na.rm = TRUE)

# Thu Oct  3 13:16:53 2024 ------------------------------
#################################################################################
############################           #########################################
########################### smaart001 #########################################
##########################           #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
smaart001 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/SMAART/SMAART_raw data/Results-2024-09-10-SMAART-RUN-001/sample_coverage.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Change sample coverage file orientation to wide
smaart001_wide <- smaart001 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(smaart001_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(smaart001_wide$Input, na.rm = TRUE) - sum(smaart001_wide$`No Dimers`, na.rm = TRUE)) / sum(smaart001_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
smaart001_wide_Per_DADA2_QC <- smaart001_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(smaart001_wide$OutputDada2, na.rm = TRUE)/sum(smaart001_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- smaart001_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/SMAART/SMAART_raw data/Results-2024-09-10-SMAART-RUN-001/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >13)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- smaart001_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE) | grepl("Neagtive", SampleID, ignore.case = TRUE))
head(filtered_Negative_Controls, n = 20)

# Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
smaart001_no_controls <- smaart001_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(smaart001_no_controls$`No Dimers`, na.rm = TRUE)
sd(smaart001_no_controls$`No Dimers`, na.rm = TRUE)

# Thu Oct  3 13:43:27 2024 ------------------------------
#################################################################################
############################           #########################################
########################### smaart002 #########################################
##########################           #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
smaart002 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/SMAART/SMAART_raw data/Results-2024-08-28-SMAART-RUN-002/sample_coverage.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)

# Change sample coverage file orientation to wide
smaart002_wide <- smaart002 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(smaart002_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(smaart002_wide$Input, na.rm = TRUE) - sum(smaart002_wide$`No Dimers`, na.rm = TRUE)) / sum(smaart002_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
smaart002_wide_Per_DADA2_QC <- smaart002_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(smaart002_wide$OutputDada2, na.rm = TRUE)/sum(smaart002_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- smaart002_wide %>%
  filter(grepl("1K", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/SMAART/SMAART_raw data/Results-2024-08-28-SMAART-RUN-002/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >11)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- smaart002_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE) | grepl("Neagtive", SampleID, ignore.case = TRUE))
head(filtered_Negative_Controls, n = 20)

#Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
smaart002_no_controls <- smaart002_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(smaart002_no_controls$`No Dimers`, na.rm = TRUE)
sd(smaart002_no_controls$`No Dimers`, na.rm = TRUE)

# Thu Oct  3 13:51:48 2024 ------------------------------
#################################################################################
############################           #########################################
########################### smaart003 #########################################
##########################           #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
smaart003 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/SMAART/SMAART_raw data/Results-2024-09-06-SMAART-RUN-003/sample_coverage.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)

# Change sample coverage file orientation to wide
smaart003_wide <- smaart003 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(smaart003_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(smaart003_wide$Input, na.rm = TRUE) - sum(smaart003_wide$`No Dimers`, na.rm = TRUE)) / sum(smaart003_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
smaart003_wide_Per_DADA2_QC <- smaart003_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(smaart003_wide$OutputDada2, na.rm = TRUE)/sum(smaart003_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- smaart003_wide %>%
  filter(grepl("pos", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/SMAART/SMAART_raw data/Results-2024-09-06-SMAART-RUN-003/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >1)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- smaart003_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE) | grepl("Neagtive", SampleID, ignore.case = TRUE))
head(filtered_Negative_Controls, n = 20)

#Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
smaart003_no_controls <- smaart003_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(smaart003_no_controls$`No Dimers`, na.rm = TRUE)
sd(smaart003_no_controls$`No Dimers`, na.rm = TRUE)

# Thu Oct  3 14:20:19 2024 ------------------------------
#################################################################################
############################           #########################################
########################### smaart006 #########################################
##########################           #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
smaart006 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/SMAART/SMAART_raw data/Results-2024-09-12-SMAART-RUN-006/sample_coverage.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)

# Change sample coverage file orientation to wide
smaart006_wide <- smaart006 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(smaart006_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(smaart006_wide$Input, na.rm = TRUE) - sum(smaart006_wide$`No Dimers`, na.rm = TRUE)) / sum(smaart006_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
smaart006_wide_Per_DADA2_QC <- smaart006_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(smaart006_wide$OutputDada2, na.rm = TRUE)/sum(smaart006_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- smaart006_wide %>%
  filter(grepl("pos", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/SMAART/SMAART_raw data/Results-2024-09-12-SMAART-RUN-006/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >13)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- smaart006_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE) | grepl("Neagtive", SampleID, ignore.case = TRUE))
head(filtered_Negative_Controls, n = 20)

#Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
smaart006_no_controls <- smaart006_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(smaart006_no_controls$`No Dimers`, na.rm = TRUE)
sd(smaart006_no_controls$`No Dimers`, na.rm = TRUE)

# Thu Oct  3 14:27:48 2024 ------------------------------
#################################################################################
############################           #########################################
########################### smaart007 #########################################
##########################           #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
smaart007 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/SMAART/SMAART_raw data/Results-2024-09-20-SMAART-RUN-007/sample_coverage.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)

# Change sample coverage file orientation to wide
smaart007_wide <- smaart007 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(smaart007_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(smaart007_wide$Input, na.rm = TRUE) - sum(smaart007_wide$`No Dimers`, na.rm = TRUE)) / sum(smaart007_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
smaart007_wide_Per_DADA2_QC <- smaart007_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(smaart007_wide$OutputDada2, na.rm = TRUE)/sum(smaart007_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- smaart007_wide %>%
  filter(grepl("POS", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/SMAART/SMAART_raw data/Results-2024-09-20-SMAART-RUN-007/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >0)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- smaart007_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE) | grepl("Neagtive", SampleID, ignore.case = TRUE))
head(filtered_Negative_Controls, n = 20)

#Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
smaart007_no_controls <- smaart007_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(smaart007_no_controls$`No Dimers`, na.rm = TRUE)
sd(smaart007_no_controls$`No Dimers`, na.rm = TRUE)

# Thu Oct  3 14:37:46 2024 ------------------------------
#################################################################################
############################           #########################################
########################### smaart008 #########################################
##########################           #########################################
#############################################################################

# Clean the environment
rm(list = ls())

# Clean the console
cat("\14")

# Load the necessary library
library(tidyverse)

# Read the sample coverage text file
smaart008 <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/SMAART/SMAART_raw data/Results-2024-09-21-SMAART-RUN-008/sample_coverage.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)

# Change sample coverage file orientation to wide
smaart008_wide <- smaart008 %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = "Stage",
    values_from = Reads
  ) %>% 
  mutate_at(vars(2:6), as.numeric)

# Sum up the reads without dimmers for the whole run
sum(smaart008_wide$`No Dimers`, na.rm = TRUE)

# Calculate percentage_dimmers on the run
percentage_dimmers <- ((sum(smaart008_wide$Input, na.rm = TRUE) - sum(smaart008_wide$`No Dimers`, na.rm = TRUE)) / sum(smaart008_wide$Input, na.rm = TRUE)) * 100
print(percentage_dimmers)


#Calculate Percentage reads without dimmers passing DADA2 QC and eyeball for every sample in the file
smaart008_wide_Per_DADA2_QC <- smaart008_wide %>% 
  mutate(`%age_passing_DADA2_QC` = OutputPostprocessing/`No Dimers`*100)
`%age_passing_DADA2_QC` <- sum(smaart008_wide$OutputDada2, na.rm = TRUE)/sum(smaart008_wide$`No Dimers`, na.rm = TRUE)*100
print(`%age_passing_DADA2_QC`)

#OutputPostProcessing reads in the 3D7 positive controls
filtered_3D7 <- smaart008_wide %>%
  filter(grepl("POS", SampleID, ignore.case = TRUE))
print(filtered_3D7)

#Checking for 3D7 positive control clonality
df <- read_delim("C:/Users/USER/OneDrive/Desktop/Violet Files/SMAART/SMAART_raw data/Results-2024-09-21-SMAART-RUN-008/allele_data.txt",
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE) %>% 
  mutate(SampleID = word(SampleID, 1, sep = "_")) %>% 
  mutate(SampleID = ifelse(grepl("control", SampleID) == TRUE, paste(SampleID, "run101", sep = "-"), SampleID)) %>% 
  filter(Reads >8)


clone <- df %>% 
  distinct(SampleID, Locus, ASV, .keep_all = TRUE) %>%
  group_by(SampleID, Locus) %>% 
  summarise(total_alleles = n()) %>%
  mutate(Clone = ifelse(total_alleles == 1, "Mono", "Poly"))
clone %<>% group_by(SampleID) %>% count(Clone)
clone_wide <- tidyr::spread(clone, Clone, n)
clone_wide[is.na(clone_wide)] <- 0
clone_wide %<>% mutate(Total_alleles = Mono + Poly) %>% 
  mutate(Mono_prop = Mono/Total_alleles) %>% 
  mutate(Poly_prop = Poly/Total_alleles)
print(sprintf("Mean of total locus = %s, Mean prop. of monoclonal loci = %s, Mean prop. of polyclonal loci = %s",
              round(mean(clone_wide$Total_alleles),2), 
              round(mean(clone_wide$Mono_prop),2), 
              round(mean(clone_wide$Poly_prop),2)))


plot_clones <- function(clone){
  g <- clone %>%
    ggplot(aes(x= SampleID, y= n)) +
    geom_col(aes(fill = Clone)) +
    theme(axis.text.y = element_text(size = 5)) +
    coord_flip()+
    scale_y_continuous(name = "Locus No.", breaks = seq(0,200, by = 50)) 
  
  g
}
library(plotly)
clone_plot <-  ggplotly(plot_clones(clone))
clone_plot

#OutputPostProcessing reads in the Negative controls
#OutputPostProcessing reads in the 3D7 positive controls
filtered_Negative_Controls <- smaart008_wide %>%
  filter(grepl("Negative", SampleID, ignore.case = TRUE) | grepl("Neagtive", SampleID, ignore.case = TRUE))
head(filtered_Negative_Controls, n = 20)

#Write the CSV including variable `%age_passing_DADA2_QC` if it seems unusual for discussion with colleagues
# write.csv(run76_wide_Per_DADA2_QC, "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/QC/Good_Miseq_M08585_Run69.csv", row.names = FALSE)

# Read through the run_wide dataframe and record in your post sequencing records log the OutputPostprocessing reads for each positive and negative controls

# Remove negative and positive controls to get the picture of coverage of samples only without controls
smaart008_no_controls <- smaart008_wide %>%
  filter(!grepl("Neg", SampleID) & !grepl("1K-Positive", SampleID) & !grepl("Undeter", SampleID) & !grepl("null", `No Dimers`))

# Calculate median and standard deviation of reads in samples only without controls
median(smaart008_no_controls$`No Dimers`, na.rm = TRUE)
sd(smaart008_no_controls$`No Dimers`, na.rm = TRUE)
