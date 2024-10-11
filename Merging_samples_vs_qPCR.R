setwd("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase1/Pooling/")

# library(dplyr)
library(dplyr)
library(writexl)
library(readxl)

#Pool 45
qPCR <- read.csv("Paragon Sample Preparation Master Spreadsheet - Reprep.csv")
pool45 <- read.csv("IM-23-135.csv")
merged <- pool47 %>% left_join(qPCR, by = c("Barcode", "Study.Subject"))




pool45 <- merged %>% 
  mutate(pool = case_when(
    Quantity < 10 ~ 30,
    Quantity >= 10 & Quantity < 100 ~ 20,  
    Quantity >= 100 & Quantity < 1000 ~ 15,
    Quantity >= 1000 & Quantity < 10000 ~ 10,
    Quantity >= 10000 & Quantity < 100000 ~ 6,
    Quantity >= 100000 ~ 3,
    TRUE ~ 6
  ))


write.csv(pool45, "Pool45.csv")


#Pool 46

qPCR <- read.csv("Paragon Sample Preparation Master Spreadsheet - Reprep.csv")
pool46 <- read.csv("IM-23-140.csv")
merged <- pool46 %>% left_join(qPCR, by = c("Barcode", "Study.Subject"))




pool46 <- merged %>% 
  mutate(pool = case_when(
    Quantity < 10 ~ 30,
    Quantity >= 10 & Quantity < 100 ~ 20,  
    Quantity >= 100 & Quantity < 1000 ~ 15,
    Quantity >= 1000 & Quantity < 10000 ~ 10,
    Quantity >= 10000 & Quantity < 100000 ~ 6,
    Quantity >= 100000 ~ 3,
    TRUE ~ 6
  ))


write.csv(pool46, "Pool46.csv")

#Pool 47
qPCR <- read.csv("Paragon Sample Preparation Master Spreadsheet - Reprep.csv")
pool47 <- read.csv("IM-23-141.csv")
merged <- pool47 %>% left_join(qPCR, by = c("Barcode", "Study.Subject"))




pool47 <- merged %>% 
  mutate(pool = case_when(
    Quantity < 10 ~ 30,
    Quantity >= 10 & Quantity < 100 ~ 20,  
    Quantity >= 100 & Quantity < 1000 ~ 15,
    Quantity >= 1000 & Quantity < 10000 ~ 10,
    Quantity >= 10000 & Quantity < 100000 ~ 6,
    Quantity >= 100000 ~ 3,
    TRUE ~ 6
    ))


write.csv(pool47, "Pool47.csv")

# Tue Jan  9 10:42:18 2024 ------------------------------
#Run-B7

#Clear the environment and the console
rm(list = ls())
cat("\14")

#Call the necessary libraries
library(dplyr)
library(writexl)
library(readxl)

#Read through the files to merge
qPCR <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/Paragon Sample Preparation Master Spreadsheet-Phase-2 - Samples to be prepped_ 2nd round.csv")
Run_B7 <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/IM-24-101.csv")

#Merge the two files, the plate sample list and qPCR file
Run_B7.qPCR_merged <- Run_B7 %>% 
  merge(qPCR, by = c("Barcode", "Study.Subject"), all.x = TRUE)

#Introducing the volumes to pool per sample according to parasite quantity.
Run_B7_Set_F_Merged <- Run_B7.qPCR_merged %>% 
  mutate(pool = case_when(
    Quantity < 10 ~ 30,
    Quantity >= 10 & Quantity < 100 ~ 20,  
    Quantity >= 100 & Quantity < 1000 ~ 15,
    Quantity >= 1000 & Quantity < 10000 ~ 10,
    Quantity >= 10000 & Quantity < 100000 ~ 6,
    Quantity >= 100000 ~ 3,
    TRUE ~ 6))

#Save the merged file with the volumes to pool
write.csv(Run_B7_Set_F_Merged, file = "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/2024-01-08-Run-B7-Set_F-Merged.csv")

# Fri Jan 12 10:20:43 2024 ------------------------------
#Clear the environment and the console
rm(list = ls())
cat("\14")

#Call the necessary libraries
library(dplyr)
library(writexl)
library(readxl)

#Read through the files to merge
qPCR <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/Paragon Sample Preparation Master Spreadsheet-Phase-2 - Samples to be prepped_ 2nd round.csv")
Run_B9 <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/IM-24-103.csv")

#Merge the two files, the plate sample list and qPCR file
Run_B9.qPCR_merged <- Run_B9 %>% 
  merge(qPCR, by = c("Barcode", "StudySubject"), all.x = TRUE)

#Introducing the volumes to pool per sample according to parasite quantity.
Run_B9_Set_D_Merged <- Run_B9.qPCR_merged %>% 
  mutate(pool = case_when(
    Quantity < 10 ~ 30,
    Quantity >= 10 & Quantity < 100 ~ 20,  
    Quantity >= 100 & Quantity < 1000 ~ 15,
    Quantity >= 1000 & Quantity < 10000 ~ 10,
    Quantity >= 10000 & Quantity < 100000 ~ 6,
    Quantity >= 100000 ~ 3,
    TRUE ~ 6))

#Save the merged file with the volumes to pool
write.csv(Run_B9_Set_D_Merged, file = "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/2024-01-11-Run-B9-Set_D-Merged.csv", na = "")

# Thu Jan 18 16:40:46 2024 ------------------------------

#Clear the environment and the console
rm(list = ls())
cat("\14")

#Call the necessary libraries
library(dplyr)
library(writexl)
library(readxl)

#Read through the files to merge
qPCR <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/Paragon Sample Preparation Master Spreadsheet-Phase-2 - Samples to be prepped_ 2nd round.csv")
Run_B12 <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/IM-24-105.csv")

#Merge the two files, the plate sample list and qPCR file
Run_B12.qPCR_merged <- Run_B12 %>% 
  left_join(qPCR, by = c("Barcode", "StudySubject"), all.x)

#Introducing the volumes to pool per sample according to parasite quantity.
Run_B12_Set_C_Merged <- Run_B12.qPCR_merged %>% 
  mutate(pool = case_when(
    Quantity < 10 ~ 30,
    Quantity >= 10 & Quantity < 100 ~ 20,  
    Quantity >= 100 & Quantity < 1000 ~ 15,
    Quantity >= 1000 & Quantity < 10000 ~ 10,
    Quantity >= 10000 & Quantity < 100000 ~ 6,
    Quantity >= 100000 ~ 3,
    TRUE ~ 6))

#Save the merged file with the volumes to pool
write.csv(Run_B12_Set_C_Merged, file = "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/2024-01-18-Run-B12-Set_C-Merged.csv", na = "")

# Mon Jan 22 09:14:57 2024 ------------------------------

#Clear the environment and the console
rm(list = ls())
cat("\14")

#Call the necessary libraries
library(dplyr)
library(writexl)
library(readxl)

#Read through the files to merge
qPCR <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/Paragon Sample Preparation Master Spreadsheet-Phase-2 - Samples to be prepped_ 2nd round.csv")
Run_B14 <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/IM-24-108.csv")

#Merge the two files, the plate sample list and qPCR file
Run_B14.qPCR_merged <- Run_B14 %>% 
  left_join(qPCR, by = c("Barcode", "StudySubject"), all.x)

#Introducing the volumes to pool per sample according to parasite quantity.
Run_B14_Set_E_Merged <- Run_B14.qPCR_merged %>% 
  mutate(pool = case_when(
    Quantity < 10 ~ 30,
    Quantity >= 10 & Quantity < 100 ~ 20,  
    Quantity >= 100 & Quantity < 1000 ~ 15,
    Quantity >= 1000 & Quantity < 10000 ~ 10,
    Quantity >= 10000 & Quantity < 100000 ~ 6,
    Quantity >= 100000 ~ 3,
    TRUE ~ 6))

#Save the merged file with the volumes to pool
write.csv(Run_B14_Set_E_Merged, file = "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/2024-01-20-Run-B14-Set_E-Merged.csv", na = "")

# Wed Mar  6 22:37:00 2024 ------------------------------
#Clear the environment and the console
rm(list = ls())
cat("\14")

#Call the necessary libraries
library(dplyr)
library(writexl)
library(readxl)

#Read through the files to merge
qPCR <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/Paragon Sample Preparation Master Spreadsheet-Phase-2 - Samples to be prepped_ 2nd round.csv") %>% 
  rename("StudySubject" = "Study.Subject") %>% select(-"X", -"site", -"Prepped.")
Run_B24 <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/IM-24-115.csv") %>% 
  select("Position", "Barcode", "StudySubject", "PlateName")

#Merge the two files, the plate sample list and qPCR file
Run_B24.qPCR_merged <- Run_B24 %>% 
  left_join(qPCR, by = c("Barcode", "StudySubject"), all.x)

#Introducing the volumes to pool per sample according to parasite quantity.
Run_B24_Set_C_Merged <- Run_B24.qPCR_merged %>% 
  mutate(pool = case_when(
    Quantity < 10 ~ 30,
    Quantity >= 10 & Quantity < 100 ~ 20,  
    Quantity >= 100 & Quantity < 1000 ~ 15,
    Quantity >= 1000 & Quantity < 10000 ~ 10,
    Quantity >= 10000 & Quantity < 100000 ~ 6,
    Quantity >= 100000 ~ 3,
    TRUE ~ 6))

#Save the merged file with the volumes to pool
write.csv(Run_B24_Set_C_Merged, file = "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/2024-03-01-Run-B24-Set_C-Merged.csv", na = "", row.names = TRUE)

# Thu Mar 14 15:27:31 2024 ------------------------------
#Clear the environment and the console
rm(list = ls())
cat("\14")

#Call the necessary libraries
library(dplyr)
library(writexl)
library(readxl)

#Read through the files to merge
qPCR <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/Paragon Sample Preparation Master Spreadsheet-Phase-2 - Samples to be prepped_ 2nd round.csv") %>% 
  rename("StudySubject" = "Study.Subject") %>% select(-"X", -"site", -"Prepped.")
Run_B28 <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/IM-24-119.csv") %>% 
  select("Position", "Barcode", "StudySubject", "PlateName")

#Merge the two files, the plate sample list and qPCR file
Run_B28.qPCR_merged <- Run_B28 %>% 
  left_join(qPCR, by = c("Barcode", "StudySubject"), all.x)

#Introducing the volumes to pool per sample according to parasite quantity.
Run_B28_Set_C_Merged <- Run_B28.qPCR_merged %>% 
  mutate(pool_vol = case_when(
    Quantity < 10 ~ 30,
    Quantity >= 10 & Quantity < 100 ~ 20,  
    Quantity >= 100 & Quantity < 1000 ~ 15,
    Quantity >= 1000 & Quantity < 10000 ~ 10,
    Quantity >= 10000 & Quantity < 100000 ~ 6,
    Quantity >= 100000 ~ 3,
    TRUE ~ 6))

#Save the merged file with the volumes to pool
write.csv(Run_B28_Set_C_Merged, file = "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/2024-03-14-Run-B28-Set_C-Merged.csv", na = "", row.names = TRUE)


# Tue Mar 26 09:13:01 2024 ------------------------------
#Clear the environment and the console
rm(list = ls())
cat("\14")

#Call the necessary libraries
library(dplyr)
library(writexl)
library(readxl)

#Read through the files to merge
qPCR <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/Paragon Sample Preparation Master Spreadsheet-Phase-2 - Samples to be prepped_ 2nd round.csv") %>% 
  rename("StudySubject" = "Study.Subject") %>% select(-"X", -"site", -"Prepped.")
Run_B33 <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/IM-24-124.csv") %>% 
  rename("Barcode" = "Tube.ID") %>% select("Position", "Barcode", "StudySubject", "PlateName")

#Merge the two files, the plate sample list and qPCR file
Run_B33.qPCR_merged <- Run_B33 %>% 
  left_join(qPCR, by = c("Barcode", "StudySubject"), all.x)

#Introducing the volumes to pool per sample according to parasite quantity.
Run_B33_Set_C_Merged <- Run_B33.qPCR_merged %>% 
  mutate(pool_vol = case_when(
    Quantity < 10 ~ 30,
    Quantity >= 10 & Quantity < 100 ~ 20,  
    Quantity >= 100 & Quantity < 1000 ~ 15,
    Quantity >= 1000 & Quantity < 10000 ~ 10,
    Quantity >= 10000 & Quantity < 100000 ~ 6,
    Quantity >= 100000 ~ 3,
    TRUE ~ 6))

#Save the merged file with the volumes to pool
write.csv(Run_B33_Set_C_Merged, file = "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/2024-03-14-Run-B33-Set_C-Merged.csv", na = "", row.names = TRUE)

# Fri May  3 08:28:27 2024 ------------------------------
#Clear the environment and the console
rm(list = ls())
cat("\14")

#Call the necessary libraries
library(dplyr)
library(writexl)
library(readxl)

#Read through the files to merge
qPCR <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase3/Library-prep/Paragon Sample Preparation Master Spreadsheet-Phase-3 - Samples to be prepped_ 3rd round.csv") %>% 
  rename("StudySubject" = "Study.Subject") %>% select(-"site", -"Prepped.")
Plate_C01 <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase3/Library-prep/IM-24-131.csv") %>% 
  rename("Barcode" = "Tube.ID") %>% select("Position", "Barcode", "StudySubject", "PlateName")

#Merge the two files, the plate sample list and qPCR file
Plate_CO1.qPCR_merged <- Plate_C01 %>% 
  left_join(qPCR, by = c("Barcode", "StudySubject"), all.x)

#Introducing the volumes to pool per sample according to parasite quantity.
Plate_C01_Set_C_Merged <- Plate_CO1.qPCR_merged %>% 
  mutate(pool_vol = case_when(
    Quantity < 10 ~ 30,
    Quantity >= 10 & Quantity < 100 ~ 20,  
    Quantity >= 100 & Quantity < 1000 ~ 15,
    Quantity >= 1000 & Quantity < 10000 ~ 10,
    Quantity >= 10000 & Quantity < 100000 ~ 6,
    Quantity >= 100000 ~ 3,
    TRUE ~ 6))

#Save the merged file with the volumes to pool
write.csv(Plate_C01_Set_C_Merged, file = "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase3/Library-prep/2024-05-02-Plate-C02-Set_C-Merged.csv", na = "", row.names = TRUE)

# Mon May 13 14:00:44 2024 ------------------------------
#Clear the environment and the console
rm(list = ls())
cat("\14")

#Call the necessary libraries
library(dplyr)
library(writexl)
library(readxl)

#Read through the files to merge
qPCR <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Reprep_samples_Jessica.csv") %>% 
  rename("StudySubject" = "sampleID", "Quantity" = "Parasitemia.by.qPCR") %>% select(-"Sample.Site")
qPCR$Barcode <- as.character(qPCR$Barcode)
Plate_B41 <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/2024-05-07-Run-B41-Set-E.csv") %>% 
  rename("Barcode" = "Tube.ID") %>% select("Position", "Barcode", "StudySubject", "PlateName")

#Merge the two files, the plate sample list and qPCR file
Plate_B41.qPCR_merged <- Plate_B41 %>% 
  left_join(qPCR, by = c("Barcode", "StudySubject"), all.x)

#Introducing the volumes to pool per sample according to parasite quantity.
Plate_B41_Set_E_Merged <- Plate_B41.qPCR_merged %>% 
  mutate(pool_vol = case_when(
    Quantity < 10 ~ 30,
    Quantity >= 10 & Quantity < 100 ~ 20,  
    Quantity >= 100 & Quantity < 1000 ~ 15,
    Quantity >= 1000 & Quantity < 10000 ~ 10,
    Quantity >= 10000 & Quantity < 100000 ~ 6,
    Quantity >= 100000 ~ 4,
    TRUE ~ 6))

#Save the merged file with the volumes to pool
write.csv(Plate_B41_Set_E_Merged, file = "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/2024-05-07-Plate-B41-Set_E-Merged.csv", na = "", row.names = TRUE)

# Fri May 17 10:06:57 2024 ------------------------------
#Clear the environment and the console
rm(list = ls())
cat("\14")

#Call the necessary libraries
library(dplyr)
library(writexl)
library(readxl)

#Read through the files to merge
qPCR <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/Pooling/Paragon Sample Preparation Master Spreadsheet-Phase-2 - Samples to be prepped_ 2nd round.csv") %>% 
  rename("StudySubject" = "Study.Subject") %>% select(-"X", -"site", -"Prepped.")
Plate_B42 <- read.csv("C:/Users/USER/OneDrive/Desktop/Bienvenu Files/IMMRSE-U/Phase2/Pooling/IM-24-134.csv") %>% 
  rename("Barcode" = "Tube.ID") %>% select("Position", "Barcode", "StudySubject", "PlateName")

#Merge the two files, the plate sample list and qPCR file
Plate_B42.qPCR_merged <- Plate_B42 %>% 
  left_join(qPCR, by = c("Barcode", "StudySubject"), all.x)

#Introducing the volumes to pool per sample according to parasite quantity.
Plate_B42_Set_F_Merged <- Plate_B42.qPCR_merged %>% 
  mutate(pool_vol = case_when(
    Quantity < 10 ~ 30, 
    Quantity >= 10 & Quantity < 100 ~ 20, 
    Quantity >= 100 & Quantity < 1000 ~ 17, 
    Quantity >= 1000 & Quantity < 10000 ~ 10, 
    Quantity >= 10000 & Quantity < 100000 ~ 6, 
    Quantity >= 100000 ~ 4, 
    TRUE ~ 6 ))

#Save the merged file with the volumes to pool
write.csv(Plate_B42_Set_F_Merged, file = "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase2/2024-05-15-Plate-B42-Set_F-Merged.csv", na = "", row.names = TRUE)

# Mon Aug  5 14:45:11 2024 ------------------------------
#Clear the environment and the console
rm(list = ls())
cat("\14")

#Call the necessary libraries
library(dplyr)
library(writexl)
library(readxl)

#Read through the files to merge
qPCR <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase3/Library-prep/Paragon Sample Preparation Master Spreadsheet-Phase-3 - Samples to be prepped_ 3rd round.csv") %>% 
  rename("StudySubject" = "Study.Subject") %>% select(-"site", -"Prepped.")
Plate_C15_C16Upper <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase3/Library-prep/Plate_C15_C16Upper.csv") %>% 
  select("Position", "Barcode", "StudySubject", "PlateName", "PoolID", "Set")

#Merge the two files, the plate sample list and qPCR file
Plate_C15_C16Upper.qPCR_merged <- Plate_C15_C16Upper %>% 
  left_join(qPCR, by = "Barcode", all.x)

#Introducing the volumes to pool per sample according to parasite quantity.
Plate_C15_C16Upper_Set_F_C_Merged <- Plate_C15_C16Upper.qPCR_merged %>% 
  mutate(pool_vol = case_when(
    Quantity < 10 ~ 30, 
    Quantity >= 10 & Quantity < 100 ~ 20, 
    Quantity >= 100 & Quantity < 1000 ~ 17, 
    Quantity >= 1000 & Quantity < 10000 ~ 10, 
    Quantity >= 10000 & Quantity < 100000 ~ 6, 
    Quantity >= 100000 ~ 4, 
    TRUE ~ 6 ))

#Save the merged file with the volumes to pool
write.csv(Plate_C15_C16Upper_Set_F_C_Merged, file = "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase3/Library-prep//2024-08-05-Plate_C15_C16Upper_Set_F_C_Merged.csv", na = "", row.names = TRUE)

# Mon Aug  5 15:06:54 2024 ------------------------------
#Clear the environment and the console
rm(list = ls())
cat("\14")

#Call the necessary libraries
library(dplyr)
library(writexl)
library(readxl)

#Read through the files to merge
qPCR <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase3/Library-prep/Paragon Sample Preparation Master Spreadsheet-Phase-3 - Samples to be prepped_ 3rd round.csv") %>% 
  rename("StudySubject" = "Study.Subject") %>% select(-"site", -"Prepped.")
Plate_C17_C16Lower <- read.csv("C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase3/Library-prep/Plate_C17_Lower_C16.csv") %>% 
  select("Position", "Barcode", "StudySubject", "PlateName", "PoolID", "Set")

#Merge the two files, the plate sample list and qPCR file
Plate_C17_C16Lower.qPCR_merged <- Plate_C17_C16Lower %>% 
  left_join(qPCR, by = "Barcode", all.x)

#Introducing the volumes to pool per sample according to parasite quantity.
Plate_C17_C16Lower_Set_C_D_Merged <- Plate_C17_C16Lower.qPCR_merged %>% 
  mutate(pool_vol = case_when(
    Quantity < 10 ~ 30, 
    Quantity >= 10 & Quantity < 100 ~ 20, 
    Quantity >= 100 & Quantity < 1000 ~ 17, 
    Quantity >= 1000 & Quantity < 10000 ~ 10, 
    Quantity >= 10000 & Quantity < 100000 ~ 6, 
    Quantity >= 100000 ~ 4, 
    TRUE ~ 6 ))

#Save the merged file with the volumes to pool
write.csv(Plate_C17_C16Lower_Set_C_D_Merged, file = "C:/Users/USER/OneDrive/Desktop/Violet Files/IMMRSE-U/Phase3/Library-prep/2024-08-05-Plate_C17_C16Lower_Set_C_D_Merged.csv", na = "", row.names = TRUE)
