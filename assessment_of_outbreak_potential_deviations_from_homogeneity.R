############## summary of how many outbreak/non outbreak cases at each vaccination level
############## outbreak defined as 5 + cases over 1 year (w sensitivity of 10+ and 20+ cases)
library(tidyverse)

############## run sensitivity analysis with outbreak defined as 10 or 20 cases over 1 year
###### 94% Vaccination ######
# read in final time threshold subset file
sim_94_f <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Full Datasets/simulation_data_94_with_key_threshold_subset_final_time.csv")
drops <- c("X", "sd_Moran", "var_Moran", "mean_Moran", "Moran_S_1", "Moran_S_2", "Moran_S_3", "Moran_S_4", "Moran_S_5", "Moran_S_6", "Moran_S_7", "Moran_S_8", "Moran_S_9", "Moran_S_10")
sim_94_f <- sim_94_f[ , !(names(sim_94_f) %in% drops)]

#load in isolation index file
iso_94 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_94.csv")
drops <- c("X", "level_1", "level_2", "level_3", "level_4")
iso_94 <- iso_94[ , !(names(iso_94) %in% drops)]

sim_94_f <- merge(sim_94_f, iso_94, by = "motif")

#create new column for cumulative incidence
sim_94_f$CI1 <- 15359- sim_94_f$S1
sim_94_f$CI2 <- 15359 - sim_94_f$S2
sim_94_f$CI3 <- 15359 - sim_94_f$S3
sim_94_f$CI4 <- 15359 - sim_94_f$S4
sim_94_f$CI5 <- 15359 - sim_94_f$S5
sim_94_f$CI6 <- 15359 - sim_94_f$S6
sim_94_f$CI7 <- 15359 - sim_94_f$S7
sim_94_f$CI8 <- 15359 - sim_94_f$S8
sim_94_f$CI9 <- 15359 - sim_94_f$S9
sim_94_f$CI10 <- 15359 - sim_94_f$S10

#let's reshape the data so we have all the S, P, CI in long form using reshape2()
library(reshape2)
drops <- c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10","S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10")
sim_94_f<- sim_94_f[ , !(names(sim_94_f) %in% drops)]

#get rid of isolation mean and variance for CI long form- will merge after converting
drops <- c("var_iso", "level_1", "level_2", "level_3", "level_4")
sim_94_f_long <- sim_94_f[ , !(names(sim_94_f) %in% drops)]

sim_94_f_long <- melt(sim_94_f_long, id.vars = c("motif", "run", "mean_iso"), variable.name = "CI", value.name = "CI_val")  
sim_94_f_long$outbreak_5_threshold <- ifelse(sim_94_f_long$CI_val > 5,1,0)
sim_94_f_long$outbreak_10_threshold <- ifelse(sim_94_f_long$CI_val > 10,1,0)
sim_94_f_long$outbreak_20_threshold <- ifelse(sim_94_f_long$CI_val > 20,1,0)

#now write out
write_csv(sim_94_f_long, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Outbreak_Probability_94_Pct_Vax_Coverage.csv")

###### 95% Vaccination ######
# read in final time threshold subset file
sim_95_f <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Full Datasets/simulation_data_95_with_key_threshold_subset_final_time.csv")
drops <- c("X", "level_1", "level_2", "level_3", "level_4",
           "sd_Moran", "var_Moran", "mean_Moran", "Moran_S_1", "Moran_S_2", "Moran_S_3", "Moran_S_4", "Moran_S_5", "Moran_S_6", "Moran_S_7", "Moran_S_8", "Moran_S_9", "Moran_S_10",
           "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10",
           "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", 
           "Moran_Inf_1", "Moran_Inf_2", "Moran_Inf_3", "Moran_Inf_4", "Moran_Inf_5", "Moran_Inf_6", "Moran_Inf_7", "Moran_Inf_8", "Moran_Inf_9", "Moran_Inf_10")
sim_95_f <- sim_95_f[ , !(names(sim_95_f) %in% drops)]

#load in isolation index file
iso_95 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_95.csv")
drops <- c("X", "var_iso", "level_1", "level_2", "level_3", "level_4")
iso_95 <- iso_95[ , !(names(iso_95) %in% drops)]

sim_95_f <- merge(sim_95_f, iso_95, by = "motif")

#let's reshape the data so we have all the S, P, CI in long form using reshape2()
sim_95_f_long <- melt(sim_95_f, id.vars = c("motif", "run", "mean_iso"), variable.name = "CI", value.name = "CI_val")  
sim_95_f_long$outbreak_5_threshold <- ifelse(sim_95_f_long$CI_val > 5,1,0)
sim_95_f_long$outbreak_10_threshold <- ifelse(sim_95_f_long$CI_val > 10,1,0)
sim_95_f_long$outbreak_20_threshold <- ifelse(sim_95_f_long$CI_val > 20,1,0)

#now write out
write_csv(sim_95_f_long, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Outbreak_Probability_95_Pct_Vax_Coverage.csv")

###### 99% Vaccination ######
sim_99_f <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Full Datasets/simulation_data_99_with_key_threshold_subset_final_time.csv")
    sim_99_f$CI1 <- 2560 - sim_99_f$S1
    sim_99_f$CI2 <- 2560 - sim_99_f$S2
    sim_99_f$CI3 <- 2560 - sim_99_f$S3
    sim_99_f$CI4 <- 2560 - sim_99_f$S4
    sim_99_f$CI5 <- 2560 - sim_99_f$S5
    sim_99_f$CI6 <- 2560 - sim_99_f$S6
    sim_99_f$CI7 <- 2560 - sim_99_f$S7
    sim_99_f$CI8 <- 2560 - sim_99_f$S8
    sim_99_f$CI9 <- 2560 - sim_99_f$S9
    sim_99_f$CI10 <- 2560 - sim_99_f$S10

drops <- c("X", "level_1", "level_2", "level_3", "level_4",
           "sd_Moran", "var_Moran", "mean_Moran", "Moran_S_1", "Moran_S_2", "Moran_S_3", "Moran_S_4", "Moran_S_5", "Moran_S_6", "Moran_S_7", "Moran_S_8", "Moran_S_9", "Moran_S_10",
           "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10",
           "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10") 

sim_99_f <- sim_99_f[ , !(names(sim_99_f) %in% drops)]

#load in isolation index file
iso_99 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_99.csv")
drops <- c("X", "var_iso", "level_1", "level_2", "level_3", "level_4")
iso_99 <- iso_99[ , !(names(iso_99) %in% drops)]

sim_99_f <- merge(sim_99_f, iso_99, by = "motif")

#let's reshape the data so we have all the S, P, CI in long form using reshape2()
sim_99_f_long <- melt(sim_99_f, id.vars = c("motif", "run", "mean_iso"), variable.name = "CI", value.name = "CI_val")  
sim_99_f_long$outbreak_5_threshold <- ifelse(sim_99_f_long$CI_val > 5,1,0)
sim_99_f_long$outbreak_10_threshold <- ifelse(sim_99_f_long$CI_val > 10,1,0)
sim_99_f_long$outbreak_20_threshold <- ifelse(sim_99_f_long$CI_val > 20,1,0)

#now write out
write_csv(sim_99_f_long, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Outbreak_Probability_99_Pct_Vax_Coverage.csv")

########### now evaluate outbreak probability based upon these different scenarios

# number of individual runs that yielded an outbreak defined by 3 thresholds for 94% vaccination
outbreak_potential_by_vaccination_level <- data.frame(vaccination_level = c(94, 94, 94, 95, 95, 95, 99, 99, 99), outbreak_threshold = c(5,10,20,5,10,20,5,10,20), number_of_runs_exceeding_outbreak_potential = rep(0,9), total_number_of_runs = c(11840, 11840, 11840, 13440, 13440, 13440, 24800, 24800, 24800), percent_runs_with_outbreak = rep(0,9))

# #calculate number of runs at/above/below outbreak potential
outbreak_potential_by_vaccination_level$number_of_runs_exceeding_outbreak_potential[1] <- sum(sim_94_f_long$outbreak_5_threshold)
outbreak_potential_by_vaccination_level$number_of_runs_exceeding_outbreak_potential[2] <- sum(sim_94_f_long$outbreak_10_threshold)
outbreak_potential_by_vaccination_level$number_of_runs_exceeding_outbreak_potential[3] <- sum(sim_94_f_long$outbreak_20_threshold)
outbreak_potential_by_vaccination_level$number_of_runs_exceeding_outbreak_potential[4] <- sum(sim_95_f_long$outbreak_5_threshold)
outbreak_potential_by_vaccination_level$number_of_runs_exceeding_outbreak_potential[5] <- sum(sim_95_f_long$outbreak_10_threshold)
outbreak_potential_by_vaccination_level$number_of_runs_exceeding_outbreak_potential[6] <- sum(sim_95_f_long$outbreak_20_threshold)
outbreak_potential_by_vaccination_level$number_of_runs_exceeding_outbreak_potential[7] <- sum(sim_99_f_long$outbreak_5_threshold)
outbreak_potential_by_vaccination_level$number_of_runs_exceeding_outbreak_potential[8] <- sum(sim_99_f_long$outbreak_10_threshold)
outbreak_potential_by_vaccination_level$number_of_runs_exceeding_outbreak_potential[9] <- sum(sim_99_f_long$outbreak_20_threshold)

outbreak_potential_by_vaccination_level$percent_runs_with_outbreak = outbreak_potential_by_vaccination_level$number_of_runs_exceeding_outbreak_potential/outbreak_potential_by_vaccination_level$total_number_of_runs

outbreak_potential <- as.data.frame(outbreak_potential_by_vaccination_level)
readr::write_csv(outbreak_potential, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/outbreak_potential_by_vaccination_level.csv")

################### plots #############
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

ggplot(sim_94_f_long, aes(x=mean_iso, y=outbreak_5_threshold)) + binomial_smooth(color = "blue") + geom_rug() + theme_minimal() +
  binomial_smooth(aes(x=mean_iso, y=outbreak_10_threshold), color = "red") + binomial_smooth(aes(x=mean_iso, y=outbreak_20_threshold), color = "black") +
  xlab("Isolation Index of Starting Motif") + ylab("Outbreak Probability") + ggtitle("94% Overall Vaccination: Outbreak Probability by Isolation Index")

ggplot(sim_95_f_long, aes(x=mean_iso, y=outbreak_5_threshold)) + binomial_smooth(color = "blue") + geom_rug() + theme_minimal() +
  binomial_smooth(aes(x=mean_iso, y=outbreak_10_threshold), color = "red") + binomial_smooth(aes(x=mean_iso, y=outbreak_20_threshold), color = "black") +
  xlab("Isolation Index of Starting Motif") + ylab("Outbreak Probability") + ggtitle("95 Overall Vaccination: Outbreak Probability by Isolation Index")

ggplot(sim_99_f_long, aes(x=mean_iso, y=outbreak_5_threshold)) + binomial_smooth(color = "blue") + geom_rug() + theme_minimal() +
  binomial_smooth(aes(x=mean_iso, y=outbreak_10_threshold), color = "red") + binomial_smooth(aes(x=mean_iso, y=outbreak_20_threshold), color = "black") +
  xlab("Isolation Index of Starting Motif") + ylab("Outbreak Probability") + ggtitle("99% Overall Vaccination: Outbreak Probability by Isolation Index")

#now make plots broken out by location of seed case
ggplot(sim_94_f_long, aes(x=mean_iso, y=outbreak_5_threshold)) + binomial_smooth() + facet_wrap(~run) + theme_minimal() +
  xlab("Isolation Index of Starting Motif") + ylab("Outbreak Probability (Using 5 Cases as Threshold)") + ggtitle("94% Overall Vaccination: Outbreak Probability by Isolation Index")

ggplot(sim_95_f_long, aes(x=mean_iso, y=outbreak_5_threshold)) + binomial_smooth() + facet_wrap(~run) + theme_minimal() +
  xlab("Isolation Index of Starting Motif") + ylab("Outbreak Probability (Using 5 Cases as Threshold)") + ggtitle("95% Overall Vaccination: Outbreak Probability by Isolation Index")

ggplot(sim_99_f_long, aes(x=mean_iso, y=outbreak_5_threshold)) + binomial_smooth() + facet_wrap(~run) + theme_minimal() +
  xlab("Isolation Index of Starting Motif") + ylab("Outbreak Probability (Using 5 Cases as Threshold)") + ggtitle("99% Overall Vaccination: Outbreak Probability by Isolation Index")
