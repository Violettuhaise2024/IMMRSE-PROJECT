###################################################
###          Install required packages          ###
###################################################
library(dplyr)
library(psych)
library(stringr)
library(naniar)
library(tidyr)


###################################################
###                 Import Data                 ###
###################################################
rm(list=ls())
setwd("C:/Users/conradm1/Desktop/Projects/STARTUP/molecular")

## metadata
kamwezi <- read.csv("metadata/2023_05_29_Kamwezi_samples_for_geno.csv")
bikurungu <- read.csv("metadata/Bikurungu_samples_for_Paragon.csv")
nagongera <- read.csv("metadata/Nagongera samples for Paragon.csv")
namokora <- read.csv("metadata/Namokora_for_genotyping.csv")
lobule <- read.csv("metadata/Paragon_Lobule_samples.csv")

##  paragon data
IM_R1 <- read.delim("raw_data/Results_2023-05-16-IM-RUN-1/resmarkers_summary.txt")
IM_R2 <- read.delim("raw_data/Results-2023-05-29-IM-RUN-2/resmarkers_summary.txt")
IM_R3 <- read.delim("raw_data/Results-2023-05-29-IM-RUN-3/resmarkers_summary.txt")
IM_R4 <- read.delim("raw_data/Results-2023-06-08-IM-RUN-4/resmarkers_summary.txt")

###################################################
###         Merge Metadata Files              ###
###################################################
kamwezi_ss <- subset(kamwezi, select = c("site", "barcode", "Sample"))
bikurungu_ss <- subset(bikurungu, select = c("site", "barcode", "Sample"))
nagongera_ss <- subset(nagongera, select = c("site", "barcode", "Sample"))
namokora_ss <- subset(namokora, select = c("site", "Study.Subject", "Barcode"))
names(namokora_ss) <- c("site", "barcode", "Sample")
lobule_ss <- subset(lobule, select = c("site", "barcode", "Sample"))

metadata <- rbind(kamwezi_ss, bikurungu_ss, nagongera_ss, namokora_ss, lobule_ss)

###################################################
###         Format Genotype Data              ###
###################################################
IM_R1_L <- pivot_longer(IM_R1, cols = dhfr_16:k13_596)
IM_R2_L <- pivot_longer(IM_R2, cols = dhfr_16:k13_596)
IM_R3_L <- pivot_longer(IM_R3, cols = dhfr_16:k13_596)
IM_R4_L <- pivot_longer(IM_R4, cols = dhfr_16:k13_596)

IM_R1_L<- 
  IM_R1_L %>%
    separate(value, sep = " ", into = c("genotype", "reads"))
IM_R2_L<- 
  IM_R2_L %>%
  separate(value, sep = " ", into = c("genotype", "reads"))
IM_R3_L<- 
  IM_R3_L %>%
  separate(value, sep = " ", into = c("genotype", "reads"))
IM_R4_L<- 
  IM_R4_L %>%
  separate(value, sep = " ", into = c("genotype", "reads"))

IM <- rbind(IM_R1_L, IM_R2_L, IM_R3_L, IM_R4_L)

IM_2<- 
  IM %>%
  separate(genotype, sep = "_", into = c("genotype1", "genotype2"))

IM_2<- 
  IM_2 %>%
  separate(reads, sep = "_", into = c("reads1", "reads2"))

IM_2a <- subset(IM_2, select = c("SampleName", "name", "genotype1", "reads1"))
names(IM_2a) <- c("SampleName", "name", "genotype", "reads")
IM_2b <- subset(IM_2, !is.na(genotype2), select = c("SampleName", "name", "genotype2", "reads2"))
names(IM_2b) <- c("SampleName", "name", "genotype", "reads")

IM_3 <- rbind(IM_2a, IM_2b)

IM_3$reads <- gsub("\\[|\\]", "", IM_3$reads) 

IM_3<- 
  IM_3 %>%
  separate(SampleName, sep = "_", into = c("SampleName", NA, NA))

###################################################
###             Link read count              ###
###################################################
IM_3$ref <- ifelse(IM_3$name=="dhfr_16" & IM_3$genotype=="A", IM_3$reads, NA)
IM_3$alt <- ifelse(IM_3$name=="dhfr_16" & IM_3$genotype=="V", IM_3$reads, NA)
IM_3$alt2 <- ifelse(IM_3$name=="dhfr_16" & IM_3$genotype=="X", IM_3$reads, NA)

IM_3$ref <- ifelse(IM_3$name=="dhfr_51" & IM_3$genotype=="N", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="dhfr_51" & IM_3$genotype=="I", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="dhfr_59" & IM_3$genotype=="C", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="dhfr_59" & IM_3$genotype=="R", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="dhfr_108" & IM_3$genotype=="S", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="dhfr_108" & IM_3$genotype=="N", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="dhfr_164" & IM_3$genotype=="I", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="dhfr_164" & IM_3$genotype=="L", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="dhfr_185" & IM_3$genotype=="T", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="dhfr_185" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="mdr1_86" & IM_3$genotype=="N", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="mdr1_86" & IM_3$genotype=="Y", IM_3$reads, IM_3$alt)
IM_3$alt2 <- ifelse(IM_3$name=="mdr1_86" & IM_3$genotype=="F", IM_3$reads, IM_3$alt2)

IM_3$ref <- ifelse(IM_3$name=="mdr1_184" & IM_3$genotype=="Y", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="mdr1_184" & IM_3$genotype=="F", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="mdr1_1034" & IM_3$genotype=="S", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="mdr1_1034" & IM_3$genotype=="C", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="mdr1_1042" & IM_3$genotype=="N", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="mdr1_1042" & IM_3$genotype=="D", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="mdr1_1246" & IM_3$genotype=="D", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="mdr1_1246" & IM_3$genotype=="Y", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="crt_72" & IM_3$genotype=="C", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="crt_72" & IM_3$genotype=="S", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="crt_73" & IM_3$genotype=="V", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="crt_73" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="crt_74" & IM_3$genotype=="M", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="crt_74" & IM_3$genotype=="I", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="crt_75" & IM_3$genotype=="N", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="crt_75" & IM_3$genotype=="E", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="crt_76" & IM_3$genotype=="K", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="crt_76" & IM_3$genotype=="T", IM_3$reads, IM_3$alt)
IM_3$alt2 <- ifelse(IM_3$name=="crt_76" & IM_3$genotype=="R", IM_3$reads, IM_3$alt2)

#IM_3$ref <- ifelse(IM_3$name=="crt_326" & IM_3$genotype=="", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="crt_326" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="dhps_436" & IM_3$genotype=="S", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="dhps_436" & IM_3$genotype=="H", IM_3$reads, IM_3$alt)
IM_3$alt2 <- ifelse(IM_3$name=="dhps_436" & IM_3$genotype=="F", IM_3$reads, IM_3$alt2)

IM_3$ref <- ifelse(IM_3$name=="dhps_437" & IM_3$genotype=="A", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="dhps_437" & IM_3$genotype=="G", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="dhps_540" & IM_3$genotype=="K", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="dhps_540" & IM_3$genotype=="E", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="dhps_581" & IM_3$genotype=="A", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="dhps_581" & IM_3$genotype=="G", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="dhps_613" & IM_3$genotype=="A", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="dhps_613" & IM_3$genotype=="S", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_675" & IM_3$genotype=="A", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="k13_675" & IM_3$genotype=="V", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_622" & IM_3$genotype=="R", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="k13_622" & IM_3$genotype=="I", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_580" & IM_3$genotype=="C", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="k13_580" & IM_3$genotype=="Y", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_574" & IM_3$genotype=="P", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_574" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_568" & IM_3$genotype=="V", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_568" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_561" & IM_3$genotype=="R", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="k13_561" & IM_3$genotype=="H", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_553" & IM_3$genotype=="P", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_553" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_543" & IM_3$genotype=="I", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="k13_543" & IM_3$genotype=="T", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_539" & IM_3$genotype=="R", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="k13_539" & IM_3$genotype=="T", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_538" & IM_3$genotype=="G", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="k13_538" & IM_3$genotype=="D", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_537" & IM_3$genotype=="N", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_537" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_527" & IM_3$genotype=="P", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="k13_527" & IM_3$genotype=="L", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_515" & IM_3$genotype=="R", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_515" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_493" & IM_3$genotype=="Y", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_493" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_481" & IM_3$genotype=="A", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_481" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_476" & IM_3$genotype=="M", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="k13_476" & IM_3$genotype=="I", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_469" & IM_3$genotype=="C", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="k13_469" & IM_3$genotype=="Y", IM_3$reads, IM_3$alt)
IM_3$alt2 <- ifelse(IM_3$name=="k13_469" & IM_3$genotype=="F", IM_3$reads, IM_3$alt2)

IM_3$ref <- ifelse(IM_3$name=="k13_458" & IM_3$genotype=="N", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_458" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_449" & IM_3$genotype=="G", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_449" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_446" & IM_3$genotype=="F", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_446" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_441" & IM_3$genotype=="P", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="k13_441" & IM_3$genotype=="L", IM_3$reads, IM_3$alt)
IM_3$alt2 <- ifelse(IM_3$name=="k13_441" & IM_3$genotype=="A", IM_3$reads, IM_3$alt2)

IM_3$ref <- ifelse(IM_3$name=="exo_415" & IM_3$genotype=="E", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="exo_415" & IM_3$genotype=="G", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="fd_193" & IM_3$genotype=="D", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="fd_193" & IM_3$genotype=="Y", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="mdr2_484" & IM_3$genotype=="T", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="mdr2_484" & IM_3$genotype=="I", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="arps10_127" & IM_3$genotype=="V", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="arps10_127" & IM_3$genotype=="M", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="dhps_431" & IM_3$genotype=="I", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="dhps_431" & IM_3$genotype=="V", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_662" & IM_3$genotype=="F", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_662" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_661" & IM_3$genotype=="Q", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_661" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_659" & IM_3$genotype=="R", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_659" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_656" & IM_3$genotype=="F", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_656" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_641" & IM_3$genotype=="D", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_641" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_631" & IM_3$genotype=="L", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_631" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_578" & IM_3$genotype=="A", IM_3$reads, IM_3$ref)
IM_3$alt <- ifelse(IM_3$name=="k13_578" & IM_3$genotype=="S", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_546" & IM_3$genotype=="Y", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_546" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_544" & IM_3$genotype=="G", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_544" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_520" & IM_3$genotype=="V", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_520" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_494" & IM_3$genotype=="V", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_494" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_485" & IM_3$genotype=="S", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_485" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_483" & IM_3$genotype=="F", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_483" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_480" & IM_3$genotype=="K", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_480" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_454" & IM_3$genotype=="F", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_454" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_607" & IM_3$genotype=="K", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_607" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_605" & IM_3$genotype=="E", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_605" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_599" & IM_3$genotype=="N", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_599" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_598" & IM_3$genotype=="L", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_598" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_597" & IM_3$genotype=="R", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_597" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

IM_3$ref <- ifelse(IM_3$name=="k13_596" & IM_3$genotype=="E", IM_3$reads, IM_3$ref)
#IM_3$alt <- ifelse(IM_3$name=="k13_596" & IM_3$genotype=="", IM_3$reads, IM_3$alt)

ref_data <- subset(IM_3, !is.na(ref))
alt_data <- subset(IM_3, !is.na(alt))
alt2_data <- subset(IM_3, !is.na(alt2))

pending <- subset(IM_3, is.na(ref) & is.na(alt) & is.na(alt2) & !is.na(genotype) & genotype!="")



ref_data_ss <- subset(ref_data, select = c("SampleName", "name", "ref"))
alt_data_ss <- subset(alt_data, select = c("SampleName", "name", "alt"))
alt2_data_ss <- subset(alt2_data, select = c("SampleName", "name", "alt2"))

IM_4 <- merge(ref_data_ss, alt_data_ss, by = c("SampleName", "name"), all = T)
IM_4b <- merge(IM_4, alt2_data_ss, by = c("SampleName", "name"), all = T)

IM_4b$ref <- as.numeric(IM_4b$ref)
IM_4b$alt <- as.numeric(IM_4b$alt)
IM_4b$alt2 <- as.numeric(IM_4b$alt2)

IM_4b$ref <- replace_na(IM_4b$ref, 0)
IM_4b$alt <- replace_na(IM_4b$alt, 0)
IM_4b$alt2 <- replace_na(IM_4b$alt2, 0)

IM_4b$coverage <- IM_4b$ref + IM_4b$alt + IM_4b$alt2


IM_5 <- IM_4b %>% group_by(SampleName, name) %>% summarise(merge_ref = sum(ref), merge_alt = sum(alt), merge_alt2 = sum(alt2), merge_coverage = sum(coverage))


###################################################
###                   Filter                   ###
###################################################
IM_6 <- subset(IM_5, merge_coverage > 10)
IM_6$ref <- 
  ifelse(IM_6$merge_ref > 10 &
           IM_6$merge_ref/IM_6$merge_coverage > 0.05, 1, 0)
IM_6$alt <- 
  ifelse(IM_6$merge_alt > 10 &
           IM_6$merge_alt/IM_6$merge_coverage > 0.05, 1, 0)

IM_6$alt2 <- 
  ifelse(IM_6$merge_alt2 > 10 &
           IM_6$merge_alt2/IM_6$merge_coverage > 0.05, 1, 0)

IM_6$genotype_bi <- ifelse((IM_6$alt == 1 | IM_6$alt2 == 1) & IM_6$ref==0, 2,
                            ifelse(IM_6$alt == 0 & IM_6$alt2 == 0 & IM_6$ref==1, 0,
                                   ifelse((IM_6$alt == 1 | IM_6$alt2 == 1) & IM_6$ref==1, 1, NA)))

IM_6$genotype_tri <- ifelse(IM_6$alt == 0 & IM_6$alt2 == 0 & IM_6$ref==1, 0,
                        ifelse(IM_6$alt == 1 & IM_6$alt2 == 0 & IM_6$ref==0, 1, 
                           ifelse(IM_6$alt == 0 & IM_6$alt2 == 1 & IM_6$ref==0, 2,
                                  
                              ifelse(IM_6$alt == 1 & IM_6$alt2 == 0 & IM_6$ref== 1, 3, 
                                 ifelse(IM_6$alt == 0 & IM_6$alt2 == 1 & IM_6$ref==1, 4, 
                                    ifelse(IM_6$alt == 1 & IM_6$alt2 == 1 & IM_6$ref==0, 5,
                                           
                                      ifelse(IM_6$alt == 1 & IM_6$alt2 == 1 & IM_6$ref==1, 6, NA)))))))




IM_ctl_bi <- subset(IM_6, str_detect(IM_6$SampleName, "ontrol") | str_detect(IM_6$SampleName, "egative") | str_detect(IM_6$SampleName, "ositive") , select = c("SampleName", "name", "genotype_bi"))
IM_ctl_bi_wide <- pivot_wider(IM_ctl_bi, names_from = "name", values_from = "genotype_bi")


IM_exp_bi <- subset(IM_6, !str_detect(IM_6$SampleName, "ontrol") & !str_detect(IM_6$SampleName, "egative") & !str_detect(IM_6$SampleName, "ositive"), select = c("SampleName", "name", "genotype_bi"))
IM_exp_bi <- subset(IM_exp_bi, select = c("SampleName", "name", "genotype_bi"))
IM_exp_bi_wide <- pivot_wider(IM_exp_bi, names_from = "name", values_from = "genotype_bi")


table(IM_exp_bi_wide$dhfr_16)
table(IM_exp_bi_wide$dhfr_51)
table(IM_exp_bi_wide$dhfr_59)
table(IM_exp_bi_wide$dhfr_108)
table(IM_exp_bi_wide$dhfr_164)
table(IM_exp_bi_wide$dhfr_185)
table(IM_exp_bi_wide$mdr1_86)
table(IM_exp_bi_wide$mdr1_184)
table(IM_exp_bi_wide$mdr1_1034)
table(IM_exp_bi_wide$mdr1_1042)
table(IM_exp_bi_wide$mdr1_1246)
table(IM_exp_bi_wide$crt_72)
table(IM_exp_bi_wide$crt_73)
table(IM_exp_bi_wide$crt_74)
table(IM_exp_bi_wide$crt_75)
table(IM_exp_bi_wide$crt_76)
#table(IM_exp_bi_wide$crt_326)
table(IM_exp_bi_wide$dhps_436)
table(IM_exp_bi_wide$dhps_437)
table(IM_exp_bi_wide$dhps_540)
table(IM_exp_bi_wide$dhps_581)
table(IM_exp_bi_wide$dhps_613)
table(IM_exp_bi_wide$k13_675)
table(IM_exp_bi_wide$k13_622)
table(IM_exp_bi_wide$k13_580)
table(IM_exp_bi_wide$k13_574)
table(IM_exp_bi_wide$k13_568)
table(IM_exp_bi_wide$k13_561)
table(IM_exp_bi_wide$k13_553)
table(IM_exp_bi_wide$k13_543)
table(IM_exp_bi_wide$k13_539)
table(IM_exp_bi_wide$k13_538)
table(IM_exp_bi_wide$k13_537)
table(IM_exp_bi_wide$k13_527)
table(IM_exp_bi_wide$k13_515)
table(IM_exp_bi_wide$k13_493)
table(IM_exp_bi_wide$k13_481)
table(IM_exp_bi_wide$k13_476)
table(IM_exp_bi_wide$k13_469)
table(IM_exp_bi_wide$k13_458)
table(IM_exp_bi_wide$k13_449)
table(IM_exp_bi_wide$k13_446)
table(IM_exp_bi_wide$k13_441)
table(IM_exp_bi_wide$exo_415)
table(IM_exp_bi_wide$fd_193)
table(IM_exp_bi_wide$mdr2_484)
table(IM_exp_bi_wide$arps10_127)
table(IM_exp_bi_wide$dhps_431)
table(IM_exp_bi_wide$k13_662)
table(IM_exp_bi_wide$k13_661)
table(IM_exp_bi_wide$k13_659)
table(IM_exp_bi_wide$k13_656)
table(IM_exp_bi_wide$k13_641)
table(IM_exp_bi_wide$k13_631)
table(IM_exp_bi_wide$k13_578)
table(IM_exp_bi_wide$k13_546)
table(IM_exp_bi_wide$k13_544)
table(IM_exp_bi_wide$k13_520)
table(IM_exp_bi_wide$k13_494)
table(IM_exp_bi_wide$k13_485)
table(IM_exp_bi_wide$k13_483)
table(IM_exp_bi_wide$k13_480)
table(IM_exp_bi_wide$k13_454)
table(IM_exp_bi_wide$k13_607)
table(IM_exp_bi_wide$k13_605)
table(IM_exp_bi_wide$k13_599)
table(IM_exp_bi_wide$k13_598)
table(IM_exp_bi_wide$k13_597)
table(IM_exp_bi_wide$k13_596)

IM_ctl_tri <- subset(IM_6, str_detect(IM_6$SampleName, "ontrol") | str_detect(IM_6$SampleName, "egative") | str_detect(IM_6$SampleName, "ositive") , select = c("SampleName", "name", "genotype_tri"))
IM_ctl_tri_wide <- pivot_wider(IM_ctl_tri, names_from = "name", values_from = "genotype_tri")


IM_exp_tri <- subset(IM_6, !str_detect(IM_6$SampleName, "ontrol") & !str_detect(IM_6$SampleName, "egative") & !str_detect(IM_6$SampleName, "ositive"), select = c("SampleName", "name", "genotype_tri"))
IM_exp_tri <- subset(IM_exp_tri, select = c("SampleName", "name", "genotype_tri"))
IM_exp_tri_wide <- pivot_wider(IM_exp_tri, names_from = "name", values_from = "genotype_tri")


table(IM_exp_tri_wide$dhfr_16)
table(IM_exp_tri_wide$dhfr_51)
table(IM_exp_tri_wide$dhfr_59)
table(IM_exp_tri_wide$dhfr_108)
table(IM_exp_tri_wide$dhfr_164)
table(IM_exp_tri_wide$dhfr_185)
table(IM_exp_tri_wide$mdr1_86)
table(IM_exp_tri_wide$mdr1_184)
table(IM_exp_tri_wide$mdr1_1034)
table(IM_exp_tri_wide$mdr1_1042)
table(IM_exp_tri_wide$mdr1_1246)
table(IM_exp_tri_wide$crt_72)
table(IM_exp_tri_wide$crt_73)
table(IM_exp_tri_wide$crt_74)
table(IM_exp_tri_wide$crt_75)
table(IM_exp_tri_wide$crt_76)
#table(IM_exp_tri_wide$crt_326)
table(IM_exp_tri_wide$dhps_436)
table(IM_exp_tri_wide$dhps_437)
table(IM_exp_tri_wide$dhps_540)
table(IM_exp_tri_wide$dhps_581)
table(IM_exp_tri_wide$dhps_613)
table(IM_exp_tri_wide$k13_675)
table(IM_exp_tri_wide$k13_622)
table(IM_exp_tri_wide$k13_580)
table(IM_exp_tri_wide$k13_574)
table(IM_exp_tri_wide$k13_568)
table(IM_exp_tri_wide$k13_561)
table(IM_exp_tri_wide$k13_553)
table(IM_exp_tri_wide$k13_543)
table(IM_exp_tri_wide$k13_539)
table(IM_exp_tri_wide$k13_538)
table(IM_exp_tri_wide$k13_537)
table(IM_exp_tri_wide$k13_527)
table(IM_exp_tri_wide$k13_515)
table(IM_exp_tri_wide$k13_493)
table(IM_exp_tri_wide$k13_481)
table(IM_exp_tri_wide$k13_476)
table(IM_exp_tri_wide$k13_469)
table(IM_exp_tri_wide$k13_458)
table(IM_exp_tri_wide$k13_449)
table(IM_exp_tri_wide$k13_446)
table(IM_exp_tri_wide$k13_441)
table(IM_exp_tri_wide$exo_415)
table(IM_exp_tri_wide$fd_193)
table(IM_exp_tri_wide$mdr2_484)
table(IM_exp_tri_wide$arps10_127)
table(IM_exp_tri_wide$dhps_431)

table(IM_exp_tri_wide$k13_662)
table(IM_exp_tri_wide$k13_661)
table(IM_exp_tri_wide$k13_659)
table(IM_exp_tri_wide$k13_656)
table(IM_exp_tri_wide$k13_641)
table(IM_exp_tri_wide$k13_631)
table(IM_exp_tri_wide$k13_578)
table(IM_exp_tri_wide$k13_546)
table(IM_exp_tri_wide$k13_544)
table(IM_exp_tri_wide$k13_520)
table(IM_exp_tri_wide$k13_494)
table(IM_exp_tri_wide$k13_485)
table(IM_exp_tri_wide$k13_483)
table(IM_exp_tri_wide$k13_480)
table(IM_exp_tri_wide$k13_454)
table(IM_exp_tri_wide$k13_607)
table(IM_exp_tri_wide$k13_605)
table(IM_exp_tri_wide$k13_599)
table(IM_exp_tri_wide$k13_598)
table(IM_exp_tri_wide$k13_597)
table(IM_exp_tri_wide$k13_596)

bi_meta <- merge(IM_exp_bi_wide, metadata, by.x="SampleName", by.y="barcode", all = TRUE)
tri_meta <- merge(IM_exp_tri_wide, metadata, by.x="SampleName", by.y="barcode", all = TRUE)


bi_meta_long <- pivot_longer(bi_meta, cols = arps10_127:mdr2_484) 
names(bi_meta_long) <- c("SampleName", "site", "Sample", "SNP", "Genotype")

tab_bi <- as.data.frame(xtabs(~SNP + Genotype + site, bi_meta_long))
tab_bi_wide <- spread(tab_bi, Genotype, Freq)
names(tab_bi_wide) <- c("SNP", "Site", "WT_N", "Mix_N", "Mut_N")
tab_bi_wide$N <- tab_bi_wide$WT_N + tab_bi_wide$Mix_N + tab_bi_wide$Mut_N
tab_bi_wide$WT_f <- round(tab_bi_wide$WT_N/tab_bi_wide$N * 100, 1)
tab_bi_wide$Mix_f <- round(tab_bi_wide$Mix_N/tab_bi_wide$N * 100, 1)
tab_bi_wide$Mut_f <- round(tab_bi_wide$Mut_N/tab_bi_wide$N * 100, 1)
tab_bi_wide <- subset(tab_bi_wide, select = c("SNP", "Site", "N", "WT_N", "Mix_N", "Mut_N", "WT_f", "Mix_f", "Mut_f"))
tab_bi_long <- pivot_longer(tab_bi_wide, cols = WT_f:Mut_f, names_to = "genotype", values_to = "prev")



tri_meta_long <- pivot_longer(tri_meta, cols = arps10_127:mdr2_484) 
names(tri_meta_long) <- c("SampleName", "site", "Sample", "SNP", "Genotype")

tab_tri <- as.data.frame(xtabs(~SNP + Genotype + site, tri_meta_long))
tab_tri_wide <- spread(tab_tri, Genotype, Freq)
names(tab_tri_wide) <- c("SNP", "Site", "WT_N", "ALT_N", "ALT2_N", "Mix_WT_ALT", "Mix_WT_ALT2")
tab_tri_wide$N <- tab_tri_wide$WT_N + tab_tri_wide$ALT_N + tab_tri_wide$ALT2_N + tab_tri_wide$Mix_WT_ALT + tab_tri_wide$Mix_WT_ALT2

tab_tri_wide$WT_f <- round(tab_tri_wide$WT_N/tab_tri_wide$N * 100, 1)
tab_tri_wide$ALT_f <- round(tab_tri_wide$ALT_N/tab_tri_wide$N * 100, 1)
tab_tri_wide$ALT2_f <- round(tab_tri_wide$ALT2_N/tab_tri_wide$N * 100, 1)
tab_tri_wide$Mix_WT_ALT_f <- round(tab_tri_wide$Mix_WT_ALT/tab_tri_wide$N * 100, 1)
tab_tri_wide$Mix_WT_ALT2_f <- round(tab_tri_wide$Mix_WT_ALT2/tab_tri_wide$N * 100, 1)


tab_tri_wide <- subset(tab_tri_wide, select = c("SNP", "Site", "N", "WT_N", "ALT_N", "ALT2_N", "Mix_WT_ALT", "Mix_WT_ALT2", "WT_f", "ALT_f", "ALT2_f", "Mix_WT_ALT_f", "Mix_WT_ALT2_f"))
tab_tri_long <- pivot_longer(tab_tri_wide, cols = WT_f:Mix_WT_ALT2_f, names_to = "genotype", values_to = "prev")


plot_tri <- 
  ggplot(tab_tri_long, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "ALT_f", "ALT2_f", "Mix_WT_ALT_f", "Mix_WT_ALT2_f", "Mix_ALT_ALT2_f",  "Mix_all3_f")))) +
  geom_col() +
  facet_wrap(facets = tab_tri_long$SNP) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4", "blue", "red"), labels = c("WT", "ALT", "ALT2", "Mix (WT+ALT)", "Mix (WT+ALT2", "Mix (ALT+ALT2)", "All 3 Variants")) +
  ylab("Prevalence (%)") +
  xlab("Year") +
  theme(axis.text.x = element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5))
plot_tri


k13 <- subset(tab_tri_long, str_detect(SNP, "k13_"))
all_k13 <-
ggplot(k13, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "ALT_f", "ALT2_f", "Mix_WT_ALT_f", "Mix_WT_ALT2_f", "Mix_ALT_ALT2_f",  "Mix_all3_f")))) +
  geom_col() +
  facet_wrap(facets = k13$SNP) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4", "blue", "red"), labels = c("WT", "ALT", "ALT2", "Mix (WT+ALT)", "Mix (WT+ALT2", "Mix (ALT+ALT2)", "All 3 Variants")) +
  ylab("Prevalence (%)") +
  xlab("Year") +
  theme(axis.text.x = element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5))

known_k13 <- subset(tab_tri_long, SNP=="k13_675" | SNP=="k13_469" | SNP=="k13_561")
known_k13 <-
  ggplot(known_k13, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "ALT_f", "ALT2_f", "Mix_WT_ALT_f", "Mix_WT_ALT2_f", "Mix_ALT_ALT2_f",  "Mix_all3_f")))) +
  geom_col() +
  facet_wrap(facets = known_k13$SNP) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4", "blue", "red"), labels = c("WT", "ALT", "ALT2", "Mix (WT+ALT)", "Mix (WT+ALT2", "Mix (ALT+ALT2)", "All 3 Variants")) +
  ylab("Prevalence (%)") +
  xlab("Year") +
  theme(axis.text.x = element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5))
known_k13

other_k13 <- subset(tab_tri_long, SNP=="k13_580" | SNP=="k13_441" | SNP=="k13_578")
other_k13 <-
  ggplot(other_k13, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "ALT_f", "ALT2_f", "Mix_WT_ALT_f", "Mix_WT_ALT2_f", "Mix_ALT_ALT2_f",  "Mix_all3_f")))) +
  geom_col() +
  facet_wrap(facets = other_k13$SNP) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4", "blue", "red"), labels = c("WT", "ALT", "ALT2", "Mix (WT+ALT)", "Mix (WT+ALT2", "Mix (ALT+ALT2)", "All 3 Variants")) +
  ylab("Prevalence (%)") +
  xlab("Year") +
  theme(axis.text.x = element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5))
other_k13

dhfr <- subset(tab_tri_long, str_detect(SNP, "dhfr_"))
dhfr <-
  ggplot(dhfr, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "ALT_f", "ALT2_f", "Mix_WT_ALT_f", "Mix_WT_ALT2_f", "Mix_ALT_ALT2_f",  "Mix_all3_f")))) +
  geom_col() +
  facet_wrap(facets = dhfr$SNP) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4", "blue", "red"), labels = c("WT", "ALT", "ALT2", "Mix (WT+ALT)", "Mix (WT+ALT2", "Mix (ALT+ALT2)", "All 3 Variants")) +
  ylab("Prevalence (%)") +
  xlab("Year") +
  theme(axis.text.x = element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5))
dhfr

dhps <- subset(tab_tri_long, str_detect(SNP, "dhps_"))
dhps <-
  ggplot(dhps, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "ALT_f", "ALT2_f", "Mix_WT_ALT_f", "Mix_WT_ALT2_f", "Mix_ALT_ALT2_f",  "Mix_all3_f")))) +
  geom_col() +
  facet_wrap(facets = dhps$SNP) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4", "blue", "red"), labels = c("WT", "ALT", "ALT2", "Mix (WT+ALT)", "Mix (WT+ALT2", "Mix (ALT+ALT2)", "All 3 Variants")) +
  ylab("Prevalence (%)") +
  xlab("Year") +
  theme(axis.text.x = element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5))
dhps


transporters <- subset(tab_tri_long, str_detect(SNP, "crt_") | str_detect(SNP, "mdr1_"))
transporters <-
  ggplot(transporters, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "ALT_f", "ALT2_f", "Mix_WT_ALT_f", "Mix_WT_ALT2_f", "Mix_ALT_ALT2_f",  "Mix_all3_f")))) +
  geom_col() +
  facet_wrap(facets = transporters$SNP) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4", "blue", "red"), labels = c("WT", "ALT", "ALT2", "Mix (WT+ALT)", "Mix (WT+ALT2", "Mix (ALT+ALT2)", "All 3 Variants")) +
  ylab("Prevalence (%)") +
  xlab("Year") +
  theme(axis.text.x = element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5))
transporters

k13_background <- subset(tab_tri_long, str_detect(SNP, "arps10_") | str_detect(SNP, "exo_") | str_detect(SNP, "fd_"))
k13_background <-
  ggplot(k13_background, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "ALT_f", "ALT2_f", "Mix_WT_ALT_f", "Mix_WT_ALT2_f", "Mix_ALT_ALT2_f",  "Mix_all3_f")))) +
  geom_col() +
  facet_wrap(facets = k13_background$SNP) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4", "blue", "red"), labels = c("WT", "ALT", "ALT2", "Mix (WT+ALT)", "Mix (WT+ALT2", "Mix (ALT+ALT2)", "All 3 Variants")) +
  ylab("Prevalence (%)") +
  xlab("Year") +
  theme(axis.text.x = element_text(angle=45, hjust = 1), plot.title = element_text(hjust = 0.5))
k13_background

ggsave(plot = known_k13, "figures/known_k13.tiff", height = 5, width = 5)
ggsave(plot = other_k13, "figures/other_k13.tiff", height = 5, width = 5)
ggsave(plot = dhfr, "figures/dhfr.tiff", height = 5, width = 5)
ggsave(plot = dhps, "figures/dhps.tiff", height = 5, width = 5)
ggsave(plot = transporters, "figures/transporters.tiff", height = 5, width = 5)
ggsave(plot = k13_background, "figures/k13_background.tiff", height = 5, width = 5)


