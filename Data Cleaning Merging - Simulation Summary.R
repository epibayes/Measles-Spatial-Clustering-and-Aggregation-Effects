# Data Cleaning and Merging for Simulation Data
# January 31, 2020
#############################

#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")

#############################################################################################
#############################################################################################
#make master data compilation of all different runs together 
library(tidyverse)
# read in first five simulations
#first, we want to only save 3 decimal places for previously saved files
#only do this once and save over old files, all new simulations run online will only contain 3 dec places for file size reasons
# sim_95_1 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95.csv")
# sim_95_1[] <- lapply(sim_95_1, function(x) {if (is.numeric(x)) round(x,3) else x})
# write.csv(sim_95_1, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_1.csv")
# 
# #now load in simulation 2
# sim_95_2 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_2.csv")
# sim_95_2[] <- lapply(sim_95_2, function(x) {if (is.numeric(x)) round(x,3) else x})
# write.csv(sim_95_2, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_2.csv")
# 
# #now load in simulation 3
# sim_95_3 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_3.csv")
# sim_95_3[] <- lapply(sim_95_3, function(x) {if (is.numeric(x)) round(x,3) else x})
# write.csv(sim_95_3, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_3.csv")
# 
# #now load in simulation 4
# sim_95_4 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_4.csv")
# sim_95_4[] <- lapply(sim_95_4, function(x) {if (is.numeric(x)) round(x,3) else x})
# write.csv(sim_95_4, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_4.csv")
# 
# #now load in simulation 5
# sim_95_5 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_5.csv")
# sim_95_5[] <- lapply(sim_95_5, function(x) {if (is.numeric(x)) round(x,3) else x})
# write.csv(sim_95_5, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_5.csv")

###########################################################################
#now just read in the 10 simulation files (3 dec places only) to merge into one large master file
sim_95_1 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_1.csv")
names(sim_95_1)
sim_95_1 <- plyr::rename(sim_95_1, c("X" = "X", "motif" = "motif", "run" = "run", "Prevalence" = "P1", "Recovered" = "R1", "Susceptibles" = "S1", "Moran_S" = "Moran_S_1", "Moran_Infected" = "Moran_Inf_1"))

sim_95_2 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_2.csv")
sim_95_2 <- plyr::rename(sim_95_2, c("X" = "X", "motif" = "motif_2", "run" = "run_2", "Prevalence" = "P2", "Recovered" = "R2", "Susceptibles" = "S2", "Moran_S" = "Moran_S_2", "Moran_Infected" = "Moran_Inf_2"))

#merge the first 2
sim_95_combo <- merge(sim_95_1, sim_95_2, by = "X")
drops <- c("motif_2","run_2")
sim_95_combo <- sim_95_combo[ , !(names(sim_95_combo) %in% drops)]

#merge in the third
sim_95_3 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_3.csv")
sim_95_3 <- plyr::rename(sim_95_3, c("X" = "X", "motif" = "motif_3", "run" = "run_3", "Prevalence" = "P3", "Recovered" = "R3", "Susceptibles" = "S3", "Moran_S" = "Moran_S_3", "Moran_Infected" = "Moran_Inf_3"))

sim_95_combo <- merge(sim_95_combo, sim_95_3, by = "X")
drops <- c("motif_3","run_3")
sim_95_combo <- sim_95_combo[ , !(names(sim_95_combo) %in% drops)]

#merge in the fourth
sim_95_4 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_4.csv")
sim_95_4 <- plyr::rename(sim_95_4, c("X" = "X", "motif" = "motif_4", "run" = "run_4", "Prevalence" = "P4", "Recovered" = "R4", "Susceptibles" = "S4", "Moran_S" = "Moran_S_4", "Moran_Infected" = "Moran_Inf_4"))

sim_95_combo <- merge(sim_95_combo, sim_95_4, by = "X")
drops <- c("motif_4","run_4")
sim_95_combo <- sim_95_combo[ , !(names(sim_95_combo) %in% drops)]

#merge in the 5th
sim_95_5 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_5.csv")
sim_95_5 <- plyr::rename(sim_95_5, c("X" = "X", "motif" = "motif_5", "run" = "run_5", "Prevalence" = "P5", "Recovered" = "R5", "Susceptibles" = "S5", "Moran_S" = "Moran_S_5", "Moran_Infected" = "Moran_Inf_5"))

sim_95_combo <- merge(sim_95_combo, sim_95_5, by = "X")
drops <- c("motif_5","run_5")
sim_95_combo <- sim_95_combo[ , !(names(sim_95_combo) %in% drops)]

#merge in the 6th
sim_95_6 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_6.csv")
sim_95_6 <- plyr::rename(sim_95_6, c("X" = "X", "motif" = "motif_6", "run" = "run_6", "Prevalence" = "P6", "Recovered" = "R6", "Susceptibles" = "S6", "Moran_S" = "Moran_S_6", "Moran_Infected" = "Moran_Inf_6"))

sim_95_combo <- merge(sim_95_combo, sim_95_6, by = "X")
drops <- c("motif_6","run_6")
sim_95_combo <- sim_95_combo[ , !(names(sim_95_combo) %in% drops)]

#merge in the 7th
sim_95_7 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_7.csv")
sim_95_7 <- plyr::rename(sim_95_7, c("X" = "X", "motif" = "motif_7", "run" = "run_7", "Prevalence" = "P7", "Recovered" = "R7", "Susceptibles" = "S7", "Moran_S" = "Moran_S_7", "Moran_Infected" = "Moran_Inf_7"))

sim_95_combo <- merge(sim_95_combo, sim_95_7, by = "X")
drops <- c("motif_7","run_7")
sim_95_combo <- sim_95_combo[ , !(names(sim_95_combo) %in% drops)]

#merge in the 8th
sim_95_8 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_8.csv")
sim_95_8 <- plyr::rename(sim_95_8, c("X" = "X", "motif" = "motif_8", "run" = "run_8", "Prevalence" = "P8", "Recovered" = "R8", "Susceptibles" = "S8", "Moran_S" = "Moran_S_8", "Moran_Infected" = "Moran_Inf_8"))

sim_95_combo <- merge(sim_95_combo, sim_95_8, by = "X")
drops <- c("motif_8","run_8")
sim_95_combo <- sim_95_combo[ , !(names(sim_95_combo) %in% drops)]

#merge in the 9th
sim_95_9 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_9.csv")
sim_95_9 <- plyr::rename(sim_95_9, c("X" = "X", "motif" = "motif_9", "run" = "run_9", "Prevalence" = "P9", "Recovered" = "R9", "Susceptibles" = "S9", "Moran_S" = "Moran_S_9", "Moran_Infected" = "Moran_Inf_9"))

sim_95_combo <- merge(sim_95_combo, sim_95_9, by = "X")
drops <- c("motif_9","run_9")
sim_95_combo <- sim_95_combo[ , !(names(sim_95_combo) %in% drops)]

#merge in the 10th
sim_95_10 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_10.csv")
sim_95_10 <- plyr::rename(sim_95_10, c("X" = "X", "motif" = "motif_10", "run" = "run_10", "Prevalence" = "P10", "Recovered" = "R10", "Susceptibles" = "S10", "Moran_S" = "Moran_S_10", "Moran_Infected" = "Moran_Inf_10"))

sim_95_combo <- merge(sim_95_combo, sim_95_10, by = "X")
drops <- c("motif_10","run_10")
sim_95_combo <- sim_95_combo[ , !(names(sim_95_combo) %in% drops)]

#write combo data out so we can just use that one moving forward and save satellite, smaller files on external hard drive
write.csv(sim_95_combo, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_combo.csv")

#############################################################################################
#now let's remove extraneous info in the first column and turn into numeric time after merging
library(stringr)
sim_95 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_95_combo.csv")

#drop extra counter column
drops <- c("X.1")
sim_95 <- sim_95[ , !(names(sim_95) %in% drops)]

#change X into time variable
sim_95$time <- as.numeric(substring(sim_95$X, regexpr("time ", sim_95$X) + 5))

#now let's reorder the columns and remove that first column
sim_95 <- sim_95[c(54,1:53)]
#now drop the cariable w motif, run, and time...
drops <- c("X")
sim_95 <- sim_95[ , !(names(sim_95) %in% drops)]

#now let's read in key for motif to Moran's I for 95% vaccination data
moran_key <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/moran_list_motif_key_95.csv")

#now let's remove extraneous info in the first column and turn into numeric motif identifier
moran_key$motif <- substring(moran_key$motif, regexpr("motif ", moran_key$motif) + 6)
moran_key$motif <- as.numeric(moran_key$motif)

#now merge with relevant information from the Moran key file
sim_95 <- merge(sim_95, moran_key, by.x = "motif")

#now write this to a file 
write.csv(sim_95, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_data_95_combined_with_key.csv")

###################################################################################################################
###################################################################################################################

### Now repeat cleaning for 94% dataset

###################################################################################################################
###################################################################################################################


library(tidyverse)

#now just read in the 10 simulation files (3 dec places only) to merge into one large master file
sim_94_1 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Run 1 94% Vaccination/simulation_summary_data_94_1.csv")
sim_94_1 <- plyr::rename(sim_94_1, c("X" = "X", "motif" = "motif", "run" = "run", "Prevalence" = "P1", "Recovered" = "R1", "Susceptibles" = "S1", "Moran_S" = "Moran_S_1", "Moran_Infected" = "Moran_Inf_1"))

sim_94_2 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Run 1 94% Vaccination/simulation_summary_data_94_2.csv")
sim_94_2 <- plyr::rename(sim_94_2, c("X" = "X", "motif" = "motif_2", "run" = "run_2", "Prevalence" = "P2", "Recovered" = "R2", "Susceptibles" = "S2", "Moran_S" = "Moran_S_2"))

#merge the first 2
sim_94_combo <- merge(sim_94_1, sim_94_2, by = "X")
drops <- c("motif_2","run_2", "Moran_Inf_1")
sim_94_combo <- sim_94_combo[ , !(names(sim_94_combo) %in% drops)]

#merge in the third
sim_94_3 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Run 1 94% Vaccination/simulation_summary_data_94_3.csv")
sim_94_3 <- plyr::rename(sim_94_3, c("X" = "X", "motif" = "motif_3", "run" = "run_3", "Prevalence" = "P3", "Recovered" = "R3", "Susceptibles" = "S3", "Moran_S" = "Moran_S_3"))

sim_94_combo <- merge(sim_94_combo, sim_94_3, by = "X")
drops <- c("motif_3","run_3")
sim_94_combo <- sim_94_combo[ , !(names(sim_94_combo) %in% drops)]

#merge in the fourth
sim_94_4 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Run 1 94% Vaccination/simulation_summary_data_94_4.csv")
sim_94_4 <- plyr::rename(sim_94_4, c("X" = "X", "motif" = "motif_4", "run" = "run_4", "Prevalence" = "P4", "Recovered" = "R4", "Susceptibles" = "S4", "Moran_S" = "Moran_S_4"))

sim_94_combo <- merge(sim_94_combo, sim_94_4, by = "X")
drops <- c("motif_4","run_4")
sim_94_combo <- sim_94_combo[ , !(names(sim_94_combo) %in% drops)]

#merge in the 5th
sim_94_5 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Run 1 94% Vaccination/simulation_summary_data_94_5.csv")
sim_94_5 <- plyr::rename(sim_94_5, c("X" = "X", "motif" = "motif_5", "run" = "run_5", "Prevalence" = "P5", "Recovered" = "R5", "Susceptibles" = "S5", "Moran_S" = "Moran_S_5"))

sim_94_combo <- merge(sim_94_combo, sim_94_5, by = "X")
drops <- c("motif_5","run_5")
sim_94_combo <- sim_94_combo[ , !(names(sim_94_combo) %in% drops)]

#merge in the 6th
sim_94_6 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Run 1 94% Vaccination/simulation_summary_data_94_6.csv")
sim_94_6 <- plyr::rename(sim_94_6, c("X" = "X", "motif" = "motif_6", "run" = "run_6", "Prevalence" = "P6", "Recovered" = "R6", "Susceptibles" = "S6", "Moran_S" = "Moran_S_6"))

sim_94_combo <- merge(sim_94_combo, sim_94_6, by = "X")
drops <- c("motif_6","run_6")
sim_94_combo <- sim_94_combo[ , !(names(sim_94_combo) %in% drops)]

#merge in the 7th
sim_94_7 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Run 1 94% Vaccination/simulation_summary_data_94_7.csv")
sim_94_7 <- plyr::rename(sim_94_7, c("X" = "X", "motif" = "motif_7", "run" = "run_7", "Prevalence" = "P7", "Recovered" = "R7", "Susceptibles" = "S7", "Moran_S" = "Moran_S_7"))

sim_94_combo <- merge(sim_94_combo, sim_94_7, by = "X")
drops <- c("motif_7","run_7")
sim_94_combo <- sim_94_combo[ , !(names(sim_94_combo) %in% drops)]

#merge in the 8th
sim_94_8 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Run 1 94% Vaccination/simulation_summary_data_94_8.csv")
sim_94_8 <- plyr::rename(sim_94_8, c("X" = "X", "motif" = "motif_8", "run" = "run_8", "Prevalence" = "P8", "Recovered" = "R8", "Susceptibles" = "S8", "Moran_S" = "Moran_S_8"))

sim_94_combo <- merge(sim_94_combo, sim_94_8, by = "X")
drops <- c("motif_8","run_8")
sim_94_combo <- sim_94_combo[ , !(names(sim_94_combo) %in% drops)]

#merge in the 9th
sim_94_9 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Run 1 94% Vaccination/simulation_summary_data_94_9.csv")
sim_94_9 <- plyr::rename(sim_94_9, c("X" = "X", "motif" = "motif_9", "run" = "run_9", "Prevalence" = "P9", "Recovered" = "R9", "Susceptibles" = "S9", "Moran_S" = "Moran_S_9"))

sim_94_combo <- merge(sim_94_combo, sim_94_9, by = "X")
drops <- c("motif_9","run_9")
sim_94_combo <- sim_94_combo[ , !(names(sim_94_combo) %in% drops)]

#merge in the 10th
sim_94_10 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Run 1 94% Vaccination/simulation_summary_data_94_10.csv")
sim_94_10 <- plyr::rename(sim_94_10, c("X" = "X", "motif" = "motif_10", "run" = "run_10", "Prevalence" = "P10", "Recovered" = "R10", "Susceptibles" = "S10", "Moran_S" = "Moran_S_10"))

sim_94_combo <- merge(sim_94_combo, sim_94_10, by = "X")
drops <- c("motif_10","run_10")
sim_94_combo <- sim_94_combo[ , !(names(sim_94_combo) %in% drops)]

#write combo data out so we can just use that one moving forward and save satellite, smaller files on external hard drive
write.csv(sim_94_combo, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/simulation_summary_data_94_combo.csv")

#############################################################################################
#now let's remove extraneous info in the first column and turn into numeric time after merging
library(stringr)
sim_94 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Full Datasets/simulation_summary_data_94_combo.csv")

#drop extra counter column
drops <- c("X.1")
sim_94 <- sim_94[ , !(names(sim_94) %in% drops)]

#change X into time variable
sim_94$time <- as.numeric(substring(sim_94$X, regexpr("time ", sim_94$X) + 5))

#now let's reorder the columns and remove that first column
sim_94 <- sim_94[c(44,1:43)]
#now drop the cariable w motif, run, and time...
drops <- c("X")
sim_94 <- sim_94[ , !(names(sim_94) %in% drops)]

#now let's read in key for motif summaries (30 runs) to Moran's I for 94% vaccination data
moran_key <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_summary_94.csv")
drops <- c("X")
moran_key <- moran_key[ , !(names(moran_key) %in% drops)]

#now merge with relevant information from the Moran key file
sim_94 <- merge(sim_94, moran_key, by= "motif")

#now write this to a file 
write.csv(sim_94, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Full Datasets/simulation_data_94_combined_with_key.csv")

######################################################################################################################
###################################    Now repeat for 99% dataset ###################################################
#now just read in the 10 simulation files (3 dec places only) to merge into one large master file
sim_99_1 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/simulation_summary_data_99_1.csv")
names(sim_99_1)
sim_99_1 <- plyr::rename(sim_99_1, c("X" = "X", "motif" = "motif", "run" = "run", "Prevalence" = "P1", "Recovered" = "R1", "Susceptibles" = "S1", "Moran_S" = "Moran_S_1"))

sim_99_2 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/simulation_summary_data_99_2.csv")
sim_99_2 <- plyr::rename(sim_99_2, c("X" = "X", "motif" = "motif_2", "run" = "run_2", "Prevalence" = "P2", "Recovered" = "R2", "Susceptibles" = "S2", "Moran_S" = "Moran_S_2"))

#merge the first 2
sim_99_combo <- merge(sim_99_1, sim_99_2, by = "X")
drops <- c("motif_2","run_2")
sim_99_combo <- sim_99_combo[ , !(names(sim_99_combo) %in% drops)]

#merge in the third
sim_99_3 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/simulation_summary_data_99_3.csv")
sim_99_3 <- plyr::rename(sim_99_3, c("X" = "X", "motif" = "motif_3", "run" = "run_3", "Prevalence" = "P3", "Recovered" = "R3", "Susceptibles" = "S3", "Moran_S" = "Moran_S_3"))

sim_99_combo <- merge(sim_99_combo, sim_99_3, by = "X")
drops <- c("motif_3","run_3")
sim_99_combo <- sim_99_combo[ , !(names(sim_99_combo) %in% drops)]

#merge in the fourth
sim_99_4 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/simulation_summary_data_99_4.csv")
sim_99_4 <- plyr::rename(sim_99_4, c("X" = "X", "motif" = "motif_4", "run" = "run_4", "Prevalence" = "P4", "Recovered" = "R4", "Susceptibles" = "S4", "Moran_S" = "Moran_S_4"))

sim_99_combo <- merge(sim_99_combo, sim_99_4, by = "X")
drops <- c("motif_4","run_4")
sim_99_combo <- sim_99_combo[ , !(names(sim_99_combo) %in% drops)]

#merge in the 5th
sim_99_5 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/simulation_summary_data_99_5.csv")
sim_99_5 <- plyr::rename(sim_99_5, c("X" = "X", "motif" = "motif_5", "run" = "run_5", "Prevalence" = "P5", "Recovered" = "R5", "Susceptibles" = "S5", "Moran_S" = "Moran_S_5"))

sim_99_combo <- merge(sim_99_combo, sim_99_5, by = "X")
drops <- c("motif_5","run_5")
sim_99_combo <- sim_99_combo[ , !(names(sim_99_combo) %in% drops)]

#merge in the 6th
sim_99_6 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/simulation_summary_data_99_6.csv")
sim_99_6 <- plyr::rename(sim_99_6, c("X" = "X", "motif" = "motif_6", "run" = "run_6", "Prevalence" = "P6", "Recovered" = "R6", "Susceptibles" = "S6", "Moran_S" = "Moran_S_6"))

sim_99_combo <- merge(sim_99_combo, sim_99_6, by = "X")
drops <- c("motif_6","run_6")
sim_99_combo <- sim_99_combo[ , !(names(sim_99_combo) %in% drops)]

#merge in the 7th
sim_99_7 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/simulation_summary_data_99_7.csv")
sim_99_7 <- plyr::rename(sim_99_7, c("X" = "X", "motif" = "motif_7", "run" = "run_7", "Prevalence" = "P7", "Recovered" = "R7", "Susceptibles" = "S7", "Moran_S" = "Moran_S_7"))

sim_99_combo <- merge(sim_99_combo, sim_99_7, by = "X")
drops <- c("motif_7","run_7")
sim_99_combo <- sim_99_combo[ , !(names(sim_99_combo) %in% drops)]

#merge in the 8th
sim_99_8 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/simulation_summary_data_99_8.csv")
sim_99_8 <- plyr::rename(sim_99_8, c("X" = "X", "motif" = "motif_8", "run" = "run_8", "Prevalence" = "P8", "Recovered" = "R8", "Susceptibles" = "S8", "Moran_S" = "Moran_S_8"))

sim_99_combo <- merge(sim_99_combo, sim_99_8, by = "X")
drops <- c("motif_8","run_8")
sim_99_combo <- sim_99_combo[ , !(names(sim_99_combo) %in% drops)]

#merge in the 9th
sim_99_9 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/simulation_summary_data_99_9.csv")
sim_99_9 <- plyr::rename(sim_99_9, c("X" = "X", "motif" = "motif_9", "run" = "run_9", "Prevalence" = "P9", "Recovered" = "R9", "Susceptibles" = "S9", "Moran_S" = "Moran_S_9"))

sim_99_combo <- merge(sim_99_combo, sim_99_9, by = "X")
drops <- c("motif_9","run_9")
sim_99_combo <- sim_99_combo[ , !(names(sim_99_combo) %in% drops)]

#merge in the 10th
sim_99_10 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/simulation_summary_data_99_10.csv")
sim_99_10 <- plyr::rename(sim_99_10, c("X" = "X", "motif" = "motif_10", "run" = "run_10", "Prevalence" = "P10", "Recovered" = "R10", "Susceptibles" = "S10", "Moran_S" = "Moran_S_10"))

sim_99_combo <- merge(sim_99_combo, sim_99_10, by = "X")
drops <- c("motif_10","run_10")
sim_99_combo <- sim_99_combo[ , !(names(sim_99_combo) %in% drops)]

#write combo data out so we can just use that one moving forward and save satellite, smaller files on external hard drive
write.csv(sim_99_combo, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Full Datasets/simulation_summary_data_99_combo.csv")




