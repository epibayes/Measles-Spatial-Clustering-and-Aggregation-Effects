###################### Sensitivity analysis of looking at impact of varying between transmission percent
###################### look at impact on total cases and also what the regression equations are...


####################### 99% #####################
#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")

sim_99_10_percent <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/simulation_summary_data_99_10_percent.csv")

#change X into time variable
sim_99_10_percent$time <- as.numeric(substring(sim_99_10_percent$X, regexpr("time ", sim_99_10_percent$X) + 5))

# #clear out X" column
drops <- c("X")
sim_99_10_percent <- sim_99_10_percent[ , !(names(sim_99_10_percent) %in% drops)]

# # #read in threshold file for key motifs
threshold_motifs_99 <- c(1,2,6,26,126)

# # #restrict dataset to only include the motifs that don't violate 1000 people per cell mark
sim_99_threshold <- sim_99_10_percent[!(sim_99_10_percent$motif %in% threshold_motifs_99),]

#merge with correct file for summary moran stats
#first - upload the moran key file
moran_list_99<- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_summary_99.csv")
drops <- c("X", "sd_Moran")
moran_list_99 <- moran_list_99[ , !(names(moran_list_99) %in% drops)]

#now merge 
sim_99_threshold_10 <- merge(sim_99_threshold, moran_list_99, by = "motif")

sim_99_f_10 <- sim_99_threshold_10[sim_99_threshold_10$time == 365,]

#load in isolation index file
iso_99 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_99.csv")
drops <- c("X", "level_1", "level_2", "level_3", "level_4")
iso_99 <- iso_99[ , !(names(iso_99) %in% drops)]

sim_99_f_10 <- merge(sim_99_f_10, iso_99, by = "motif")

#clear out two unnamed "X" columns and all "R" columns
drops <- c("X", "time")
sim_99_f_10 <- sim_99_f_10[ , !(names(sim_99_f_10) %in% drops)]

#write out file
write.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/final_time_dataset_99%_coverage_10%_between_percent.csv")


sim_99_f_10$CI <- 2560 - sim_99_f_10$Susceptibles

#total number of cases 
summary(sim_99_f_10$CI)


## now let's plot variance in Moran's I by cum incidence
ggplot(sim_99_f_10, aes(x = mean_iso, y = CI)) + geom_point(colour = "blue") + theme_classic() +
  labs(title = "Cumulative Incidence after 1 year across 620 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 99% (10% between cell transmission)", x = "Variance of Moran's I", y = "Cumulative Incidence", color = "Legend") 

#now we have restricted long-form dataset ot work with.
#first summarize total cases per quadrant of seed case
sim_99_Q1 <- sim_99_f_10[sim_99_f_10$run == 1,]
sim_99_Q2 <- sim_99_f_10[sim_99_f_10$run == 2,]
sim_99_Q3 <- sim_99_f_10[sim_99_f_10$run == 3,]
sim_99_Q4 <- sim_99_f_10[sim_99_f_10$run == 4,]
Q1 <- c(1, min(sim_99_Q1$CI), median(sim_99_Q1$CI), mean(sim_99_Q1$CI), max(sim_99_Q1$CI)) 
Q2 <- c(2, min(sim_99_Q2$CI), median(sim_99_Q2$CI), mean(sim_99_Q2$CI), max(sim_99_Q2$CI)) 
Q3 <- c(3, min(sim_99_Q3$CI), median(sim_99_Q3$CI), mean(sim_99_Q3$CI), max(sim_99_Q3$CI)) 
Q4 <- c(4, min(sim_99_Q4$CI), median(sim_99_Q4$CI), mean(sim_99_Q4$CI), max(sim_99_Q4$CI)) 

sum_cases <- as.data.frame(rbind(Q1, Q2, Q3, Q4))
names(sum_cases) <- c("run", "min", "median", "mean", "max")

########################################################################################################################
###############################################  25%       #############################################################

#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")

sim_99_25_percent <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/simulation_summary_data_99_25_percent.csv")

#change X into time variable
sim_99_25_percent$time <- as.numeric(substring(sim_99_25_percent$X, regexpr("time ", sim_99_25_percent$X) + 5))

# #clear out X" column
drops <- c("X")
sim_99_25_percent <- sim_99_25_percent[ , !(names(sim_99_25_percent) %in% drops)]

# # #read in threshold file for key motifs
threshold_motifs_99 <- c(1,2,6,26,126)

# # #restrict dataset to only include the motifs that don't violate 1000 people per cell mark
sim_99_threshold <- sim_99_25_percent[!(sim_99_25_percent$motif %in% threshold_motifs_99),]

#merge with correct file for summary moran stats
#first - upload the moran key file
moran_list_99<- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_summary_99.csv")
drops <- c("X", "sd_Moran")
moran_list_99 <- moran_list_99[ , !(names(moran_list_99) %in% drops)]

#now merge 
sim_99_threshold_25 <- merge(sim_99_threshold, moran_list_99, by = "motif")

sim_99_f_25 <- sim_99_threshold_25[sim_99_threshold_25$time == 365,]

#load in isolation index file
iso_99 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_99.csv")
drops <- c("X", "level_1", "level_2", "level_3", "level_4")
iso_99 <- iso_99[ , !(names(iso_99) %in% drops)]

sim_99_f_25 <- merge(sim_99_f_25, iso_99, by = "motif")

#clear out two unnamed "X" columns and all "R" columns
drops <- c("X", "time")
sim_99_f_25 <- sim_99_f_25[ , !(names(sim_99_f_25) %in% drops)]

#write out file
write.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/final_time_dataset_99%_coverage_25%_between_percent.csv")


sim_99_f_25$CI <- 2560 - sim_99_f_25$Susceptibles

#total number of cases 
summary(sim_99_f_25$CI)


## now let's plot variance in Moran's I by cum incidence
ggplot(sim_99_f_25, aes(x = mean_iso, y = CI)) + geom_point(colour = "blue") + theme_classic() +
  labs(title = "Cumulative Incidence after 1 year across 620 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 99% (10% between cell transmission)", x = "Variance of Moran's I", y = "Cumulative Incidence", color = "Legend") 

#now we have restricted long-form dataset ot work with.
#first summarize total cases per quadrant of seed case
sim_99_Q1 <- sim_99_f_25[sim_99_f_25$run == 1,]
sim_99_Q2 <- sim_99_f_25[sim_99_f_25$run == 2,]
sim_99_Q3 <- sim_99_f_25[sim_99_f_25$run == 3,]
sim_99_Q4 <- sim_99_f_25[sim_99_f_25$run == 4,]
Q1 <- c(1, min(sim_99_Q1$CI), median(sim_99_Q1$CI), mean(sim_99_Q1$CI), max(sim_99_Q1$CI)) 
Q2 <- c(2, min(sim_99_Q2$CI), median(sim_99_Q2$CI), mean(sim_99_Q2$CI), max(sim_99_Q2$CI)) 
Q3 <- c(3, min(sim_99_Q3$CI), median(sim_99_Q3$CI), mean(sim_99_Q3$CI), max(sim_99_Q3$CI)) 
Q4 <- c(4, min(sim_99_Q4$CI), median(sim_99_Q4$CI), mean(sim_99_Q4$CI), max(sim_99_Q4$CI)) 

sum_cases <- as.data.frame(rbind(Q1, Q2, Q3, Q4))
names(sum_cases) <- c("run", "min", "median", "mean", "max")

########################################################################################################################
###############################################  75%       #############################################################

#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")

sim_99_75_percent <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/simulation_summary_data_99_75_percent.csv")

#change X into time variable
sim_99_75_percent$time <- as.numeric(substring(sim_99_75_percent$X, regexpr("time ", sim_99_75_percent$X) + 5))

# #clear out X" column
drops <- c("X")
sim_99_75_percent <- sim_99_75_percent[ , !(names(sim_99_75_percent) %in% drops)]

# # #read in threshold file for key motifs
threshold_motifs_99 <- c(1,2,6,26,126)

# # #restrict dataset to only include the motifs that don't violate 1000 people per cell mark
sim_99_threshold <- sim_99_75_percent[!(sim_99_75_percent$motif %in% threshold_motifs_99),]

#merge with correct file for summary moran stats
#first - upload the moran key file
moran_list_99<- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_summary_99.csv")
drops <- c("X", "sd_Moran")
moran_list_99 <- moran_list_99[ , !(names(moran_list_99) %in% drops)]

#now merge 
sim_99_threshold_75 <- merge(sim_99_threshold, moran_list_99, by = "motif")

sim_99_f_75 <- sim_99_threshold_75[sim_99_threshold_75$time == 365,]

#load in isolation index file
iso_99 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_99.csv")
drops <- c("X", "level_1", "level_2", "level_3", "level_4")
iso_99 <- iso_99[ , !(names(iso_99) %in% drops)]

sim_99_f_75 <- merge(sim_99_f_75, iso_99, by = "motif")

#clear out two unnamed "X" columns and all "R" columns
drops <- c("X", "time")
sim_99_f_75 <- sim_99_f_75[ , !(names(sim_99_f_75) %in% drops)]

#write out file
write.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/final_time_dataset_99%_coverage_75%_between_percent.csv")


sim_99_f_75$CI <- 2560 - sim_99_f_75$Susceptibles

#total number of cases 
summary(sim_99_f_75$CI)


## now let's plot variance in Moran's I by cum incidence
ggplot(sim_99_f_75, aes(x = mean_iso, y = CI)) + geom_point(colour = "blue") + theme_classic() +
  labs(title = "Cumulative Incidence after 1 year across 620 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 99% (10% between cell transmission)", x = "Variance of Moran's I", y = "Cumulative Incidence", color = "Legend") 

#now we have restricted long-form dataset ot work with.
#first summarize total cases per quadrant of seed case
sim_99_Q1 <- sim_99_f_75[sim_99_f_75$run == 1,]
sim_99_Q2 <- sim_99_f_75[sim_99_f_75$run == 2,]
sim_99_Q3 <- sim_99_f_75[sim_99_f_75$run == 3,]
sim_99_Q4 <- sim_99_f_75[sim_99_f_75$run == 4,]
Q1 <- c(1, min(sim_99_Q1$CI), median(sim_99_Q1$CI), mean(sim_99_Q1$CI), max(sim_99_Q1$CI)) 
Q2 <- c(2, min(sim_99_Q2$CI), median(sim_99_Q2$CI), mean(sim_99_Q2$CI), max(sim_99_Q2$CI)) 
Q3 <- c(3, min(sim_99_Q3$CI), median(sim_99_Q3$CI), mean(sim_99_Q3$CI), max(sim_99_Q3$CI)) 
Q4 <- c(4, min(sim_99_Q4$CI), median(sim_99_Q4$CI), mean(sim_99_Q4$CI), max(sim_99_Q4$CI)) 

sum_cases <- as.data.frame(rbind(Q1, Q2, Q3, Q4))
names(sum_cases) <- c("run", "min", "median", "mean", "max")

quad <- c(1,2,3,4)
mean_10 <- c(469.2, 141, 148.8, 71.9)
mean_25 <- c(419.6, 117.9, 123.9, 70.7)
mean_50 <- c(295.9, 56.3, 53.9, 18.8)
mean_75 <- c(122.6, 0.014, 0.014, 0.086)

means <- as.data.frame(cbind(quad, mean_10, mean_25, mean_50, mean_75))
names(means) <- c("quad", "ten", "twentyfive", "fifty", "seventyfive")

library(ggplot2)
ggplot(means, aes(x = as.factor(quad))) +
  geom_point(aes(y = ten), col = "green", size = 3.5) +
  geom_point(aes(y = twentyfive), col = "yellow", size = 3.5) +
  geom_point(aes(y = fifty), col = "orange", size = 3.5) +
  geom_point(aes(y = seventyfive), col = "red", size = 3.5) + theme_dark() +
labs(title = "Effect of varying between cell transmission percent on case burden: Overall Vaccination Coverage at 99%", x = "Quadrant of Seed Case", y = "Cumulative Incidence", color = "Legend") +
ggplot2::scale_color_manual(values = c("ten percent" = "green", "twenty five percent" =  "yellow", "fifty percent" = "orange","seventyfive percent" =  "red"), labels = c(10, 25, 50, 75))

##################################################################################################################################
##################################################################################################################################
######################################### Sensitivity analysis of looking at impact of varying between transmission percent ######
######################################### look at impact on total cases and also what the regression equations are...  ###########
######################################### look now at dataset with 95% overall vaccination   #####################################
##################################################################################################################################
##################################################################################################################################


#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")

sim_95_10_percent <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/95% Vaccination/simulation_summary_data_95_10_percent.csv")

#change X into time variable
sim_95_10_percent$time <- as.numeric(substring(sim_95_10_percent$X, regexpr("time ", sim_95_10_percent$X) + 5))

# #clear out X" column
drops <- c("X")
sim_95_10_percent <- sim_95_10_percent[ , !(names(sim_95_10_percent) %in% drops)]

# these have already been restricted to threshold

#merge with correct file for summary moran stats
#first - upload the moran key file
moran_list_95<- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_summary_95.csv")
drops <- c("X", "sd_Moran")
moran_list_95 <- moran_list_95[ , !(names(moran_list_95) %in% drops)]

#now merge 
sim_95_10_percent <- merge(sim_95_10_percent, moran_list_95, by = "motif")

sim_95_f_10 <- sim_95_10_percent[sim_95_10_percent$time == 365,]

#load in isolation index file
iso_95 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_95.csv")
drops <- c("X", "level_1", "level_2", "level_3", "level_4")
iso_95 <- iso_95[ , !(names(iso_95) %in% drops)]

sim_95_f_10 <- merge(sim_95_f_10, iso_95, by = "motif")

#clear out two unnamed "X" columns and all "R" columns
drops <- c("X", "time")
sim_95_f_10 <- sim_95_f_10[ , !(names(sim_95_f_10) %in% drops)]

#write out file
write.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/95% Vaccination/final_time_dataset_95%_coverage_10%_between_percent.csv")


sim_95_f_10$CI <- 12799 - sim_95_f_10$Susceptibles

#total number of cases 
summary(sim_95_f_10$CI)


## now let's plot variance in Moran's I by cum incidence
ggplot(sim_95_f_10, aes(x = mean_iso, y = CI)) + geom_point(colour = "blue") + theme_classic() +
  labs(title = "Cumulative Incidence after 1 year across 620 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 95% (10% between cell transmission)", x = "Variance of Moran's I", y = "Cumulative Incidence", color = "Legend") 

#now we have restricted long-form dataset ot work with.
#first summarize total cases per quadrant of seed case
sim_95_Q1 <- sim_95_f_10[sim_95_f_10$run == 1,]
sim_95_Q2 <- sim_95_f_10[sim_95_f_10$run == 2,]
sim_95_Q3 <- sim_95_f_10[sim_95_f_10$run == 3,]
sim_95_Q4 <- sim_95_f_10[sim_95_f_10$run == 4,]
Q1 <- c(1, min(sim_95_Q1$CI), median(sim_95_Q1$CI), mean(sim_95_Q1$CI), max(sim_95_Q1$CI)) 
Q2 <- c(2, min(sim_95_Q2$CI), median(sim_95_Q2$CI), mean(sim_95_Q2$CI), max(sim_95_Q2$CI)) 
Q3 <- c(3, min(sim_95_Q3$CI), median(sim_95_Q3$CI), mean(sim_95_Q3$CI), max(sim_95_Q3$CI)) 
Q4 <- c(4, min(sim_95_Q4$CI), median(sim_95_Q4$CI), mean(sim_95_Q4$CI), max(sim_95_Q4$CI)) 

sum_cases <- as.data.frame(rbind(Q1, Q2, Q3, Q4))
names(sum_cases) <- c("run", "min", "median", "mean", "max")

########################################################################################################################
###############################################  25%       #############################################################

#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")

sim_95_25_percent <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/95% Vaccination/simulation_summary_data_95_25_percent.csv")

#change X into time variable
sim_95_25_percent$time <- as.numeric(substring(sim_95_25_percent$X, regexpr("time ", sim_95_25_percent$X) + 5))

# #clear out X" column
drops <- c("X")
sim_95_25_percent <- sim_95_25_percent[ , !(names(sim_95_25_percent) %in% drops)]


#merge with correct file for summary moran stats
#first - upload the moran key file
moran_list_95<- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_summary_95.csv")
drops <- c("X", "sd_Moran")
moran_list_95 <- moran_list_95[ , !(names(moran_list_95) %in% drops)]

#now merge 
sim_95_25_percent <- merge(sim_95_25_percent, moran_list_95, by = "motif")

sim_95_f_25 <- sim_95_25_percent[sim_95_25_percent$time == 365,]

#load in isolation index file
iso_95 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_95.csv")
drops <- c("X", "level_1", "level_2", "level_3", "level_4")
iso_95 <- iso_95[ , !(names(iso_95) %in% drops)]

sim_95_f_25 <- merge(sim_95_f_25, iso_95, by = "motif")

#clear out two unnamed "X" columns and all "R" columns
drops <- c("X", "time")
sim_95_f_25 <- sim_95_f_25[ , !(names(sim_95_f_25) %in% drops)]

#write out file
write.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/95% Vaccination/final_time_dataset_95%_coverage_25%_between_percent.csv")


sim_95_f_25$CI <- 12799 - sim_95_f_25$Susceptibles

#total number of cases 
summary(sim_95_f_25$CI)


## now let's plot variance in Moran's I by cum incidence
ggplot(sim_95_f_25, aes(x = mean_iso, y = CI)) + geom_point(colour = "blue") + theme_classic() +
  labs(title = "Cumulative Incidence after 1 year across 620 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 95% (10% between cell transmission)", x = "Variance of Moran's I", y = "Cumulative Incidence", color = "Legend") 

#now we have restricted long-form dataset ot work with.
#first summarize total cases per quadrant of seed case
sim_95_Q1 <- sim_95_f_25[sim_95_f_25$run == 1,]
sim_95_Q2 <- sim_95_f_25[sim_95_f_25$run == 2,]
sim_95_Q3 <- sim_95_f_25[sim_95_f_25$run == 3,]
sim_95_Q4 <- sim_95_f_25[sim_95_f_25$run == 4,]
Q1 <- c(1, min(sim_95_Q1$CI), median(sim_95_Q1$CI), mean(sim_95_Q1$CI), max(sim_95_Q1$CI)) 
Q2 <- c(2, min(sim_95_Q2$CI), median(sim_95_Q2$CI), mean(sim_95_Q2$CI), max(sim_95_Q2$CI)) 
Q3 <- c(3, min(sim_95_Q3$CI), median(sim_95_Q3$CI), mean(sim_95_Q3$CI), max(sim_95_Q3$CI)) 
Q4 <- c(4, min(sim_95_Q4$CI), median(sim_95_Q4$CI), mean(sim_95_Q4$CI), max(sim_95_Q4$CI)) 

sum_cases <- as.data.frame(rbind(Q1, Q2, Q3, Q4))
names(sum_cases) <- c("run", "min", "median", "mean", "max")

########################################################################################################################
###############################################  75%       #############################################################

#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")

sim_95_75_percent <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/95% Vaccination/simulation_summary_data_95_75_percent.csv")

#change X into time variable
sim_95_75_percent$time <- as.numeric(substring(sim_95_75_percent$X, regexpr("time ", sim_95_75_percent$X) + 5))

# #clear out X" column
drops <- c("X")
sim_95_75_percent <- sim_95_75_percent[ , !(names(sim_95_75_percent) %in% drops)]

#merge with correct file for summary moran stats
#first - upload the moran key file
moran_list_95<- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_summary_95.csv")
drops <- c("X", "sd_Moran")
moran_list_95 <- moran_list_95[ , !(names(moran_list_95) %in% drops)]

#now merge 
sim_95_75_percent <- merge(sim_95_75_percent, moran_list_95, by = "motif")

sim_95_f_75 <- sim_95_75_percent[sim_95_75_percent$time == 365,]

#load in isolation index file
iso_95 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_95.csv")
drops <- c("X", "level_1", "level_2", "level_3", "level_4")
iso_95 <- iso_95[ , !(names(iso_95) %in% drops)]

sim_95_f_75 <- merge(sim_95_f_75, iso_95, by = "motif")

#clear out two unnamed "X" columns and all "R" columns
drops <- c("X", "time")
sim_95_f_75 <- sim_95_f_75[ , !(names(sim_95_f_75) %in% drops)]

#write out file
write.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/95% Vaccination/final_time_dataset_95%_coverage_75%_between_percent.csv")


sim_95_f_75$CI <- 12799 - sim_95_f_75$Susceptibles

#total number of cases 
summary(sim_95_f_75$CI)


## now let's plot variance in Moran's I by cum incidence
ggplot(sim_95_f_75, aes(x = mean_iso, y = CI)) + geom_point(colour = "blue") + theme_classic() +
  labs(title = "Cumulative Incidence after 1 year across 620 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 95% (10% between cell transmission)", x = "Variance of Moran's I", y = "Cumulative Incidence", color = "Legend") 

#now we have restricted long-form dataset ot work with.
#first summarize total cases per quadrant of seed case
sim_95_Q1 <- sim_95_f_75[sim_95_f_75$run == 1,]
sim_95_Q2 <- sim_95_f_75[sim_95_f_75$run == 2,]
sim_95_Q3 <- sim_95_f_75[sim_95_f_75$run == 3,]
sim_95_Q4 <- sim_95_f_75[sim_95_f_75$run == 4,]
Q1 <- c(1, min(sim_95_Q1$CI), median(sim_95_Q1$CI), mean(sim_95_Q1$CI), max(sim_95_Q1$CI)) 
Q2 <- c(2, min(sim_95_Q2$CI), median(sim_95_Q2$CI), mean(sim_95_Q2$CI), max(sim_95_Q2$CI)) 
Q3 <- c(3, min(sim_95_Q3$CI), median(sim_95_Q3$CI), mean(sim_95_Q3$CI), max(sim_95_Q3$CI)) 
Q4 <- c(4, min(sim_95_Q4$CI), median(sim_95_Q4$CI), mean(sim_95_Q4$CI), max(sim_95_Q4$CI)) 

sum_cases <- as.data.frame(rbind(Q1, Q2, Q3, Q4))
names(sum_cases) <- c("run", "min", "median", "mean", "max")

quad <- c(1,2,3,4)
mean_10 <- c(6213, 5393, 5380, 5003)
mean_25 <- c(5876, 5100, 5109, 4742)
mean_50 <- c(4748, 3779, 3778, 3248)
mean_75 <- c(2652, 1172, 1192, 649)

means <- as.data.frame(cbind(quad, mean_10, mean_25, mean_50, mean_75))
names(means) <- c("quad", "ten", "twentyfive", "fifty", "seventyfive")

library(ggplot2)
ggplot(means, aes(x = as.factor(quad))) +
  geom_point(aes(y = ten), col = "green", size = 3.5) + geom_smooth(method = "loess", aes(x = quad, y = ten), col = "green") +
  geom_point(aes(y = twentyfive), col = "yellow", size = 3.5) +geom_smooth(method = "loess", aes(x = quad, y = twentyfive), col = "yellow") +
  geom_point(aes(y = fifty), col = "orange", size = 3.5) + geom_smooth(method = "loess", aes(x = quad, y = fifty), col = "orange") +
  geom_point(aes(y = seventyfive), col = "red", size = 3.5) + geom_smooth(method = "loess", aes(x = quad, y = seventyfive), col = "red") + theme_dark() +
  labs(title = "Effect of varying between cell transmission percent on case burden: Overall Vaccination Coverage at 95%", x = "Quadrant of Seed Case", y = "Cumulative Incidence", color = "Legend") +
  ggplot2::scale_color_manual(values = c("ten percent" = "green", "twenty five percent" =  "yellow", "fifty percent" = "orange","seventyfive percent" =  "red"), labels = c(10, 25, 50, 75))

##################################################################################################################################
##################################################################################################################################
######################################### Sensitivity analysis of looking at impact of varying between transmission percent ######
######################################### look at impact on total cases and also what the regression equations are...  ###########
######################################### look now at dataset with 94% overall vaccination   #####################################
##################################################################################################################################
##################################################################################################################################


#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")

sim_94_10_percent <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Run 1 94% Vaccination/simulation_summary_data_94_10_percent.csv")

#change X into time variable
sim_94_10_percent$time <- as.numeric(substring(sim_94_10_percent$X, regexpr("time ", sim_94_10_percent$X) + 5))

# #clear out X" column
drops <- c("X")
sim_94_10_percent <- sim_94_10_percent[ , !(names(sim_94_10_percent) %in% drops)]

# these have already been restricted to threshold

#merge with correct file for summary moran stats
#first - upload the moran key file
moran_list_94<- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_summary_94.csv")
drops <- c("X", "sd_Moran")
moran_list_94 <- moran_list_94[ , !(names(moran_list_94) %in% drops)]

#now merge 
sim_94_10_percent <- merge(sim_94_10_percent, moran_list_94, by = "motif")

sim_94_f_10 <- sim_94_10_percent[sim_94_10_percent$time == 365,]

#load in isolation index file
iso_94 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_94.csv")
drops <- c("X", "level_1", "level_2", "level_3", "level_4")
iso_94 <- iso_94[ , !(names(iso_94) %in% drops)]

sim_94_f_10 <- merge(sim_94_f_10, iso_94, by = "motif")

#clear out two unnamed "X" columns and all "R" columns
drops <- c("X", "time")
sim_94_f_10 <- sim_94_f_10[ , !(names(sim_94_f_10) %in% drops)]

#write out file
write.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/94% Vaccination/final_time_dataset_94%_coverage_10%_between_percent.csv")


sim_94_f_10$CI <- 15359 - sim_94_f_10$Susceptibles

#total number of cases 
summary(sim_94_f_10$CI)


## now let's plot variance in Moran's I by cum incidence
ggplot(sim_94_f_10, aes(x = mean_iso, y = CI)) + geom_point(colour = "blue") + theme_classic() +
  labs(title = "Cumulative Incidence after 1 year across 620 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 94% (10% between cell transmission)", x = "Variance of Moran's I", y = "Cumulative Incidence", color = "Legend") 

#now we have restricted long-form dataset ot work with.
#first summarize total cases per quadrant of seed case
sim_94_Q1 <- sim_94_f_10[sim_94_f_10$run == 1,]
sim_94_Q2 <- sim_94_f_10[sim_94_f_10$run == 2,]
sim_94_Q3 <- sim_94_f_10[sim_94_f_10$run == 3,]
sim_94_Q4 <- sim_94_f_10[sim_94_f_10$run == 4,]
Q1 <- c(1, min(sim_94_Q1$CI), median(sim_94_Q1$CI), mean(sim_94_Q1$CI), max(sim_94_Q1$CI)) 
Q2 <- c(2, min(sim_94_Q2$CI), median(sim_94_Q2$CI), mean(sim_94_Q2$CI), max(sim_94_Q2$CI)) 
Q3 <- c(3, min(sim_94_Q3$CI), median(sim_94_Q3$CI), mean(sim_94_Q3$CI), max(sim_94_Q3$CI)) 
Q4 <- c(4, min(sim_94_Q4$CI), median(sim_94_Q4$CI), mean(sim_94_Q4$CI), max(sim_94_Q4$CI)) 

sum_cases <- as.data.frame(rbind(Q1, Q2, Q3, Q4))
names(sum_cases) <- c("run", "min", "median", "mean", "max")

########################################################################################################################
###############################################  25%       #############################################################

#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")

sim_94_25_percent <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Run 1 94% Vaccination/simulation_summary_data_94_25_percent.csv")

#change X into time variable
sim_94_25_percent$time <- as.numeric(substring(sim_94_25_percent$X, regexpr("time ", sim_94_25_percent$X) + 5))

# #clear out X" column
drops <- c("X")
sim_94_25_percent <- sim_94_25_percent[ , !(names(sim_94_25_percent) %in% drops)]


#merge with correct file for summary moran stats
#first - upload the moran key file
moran_list_94<- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_summary_94.csv")
drops <- c("X", "sd_Moran")
moran_list_94 <- moran_list_94[ , !(names(moran_list_94) %in% drops)]

#now merge 
sim_94_25_percent <- merge(sim_94_25_percent, moran_list_94, by = "motif")

sim_94_f_25 <- sim_94_25_percent[sim_94_25_percent$time == 365,]

#load in isolation index file
iso_94 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_94.csv")
drops <- c("X", "level_1", "level_2", "level_3", "level_4")
iso_94 <- iso_94[ , !(names(iso_94) %in% drops)]

sim_94_f_25 <- merge(sim_94_f_25, iso_94, by = "motif")

#clear out two unnamed "X" columns and all "R" columns
drops <- c("X", "time")
sim_94_f_25 <- sim_94_f_25[ , !(names(sim_94_f_25) %in% drops)]

#write out file
#write.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/94% Vaccination/final_time_dataset_94%_coverage_25%_between_percent.csv")


sim_94_f_25$CI <- 15359 - sim_94_f_25$Susceptibles

#total number of cases 
summary(sim_94_f_25$CI)


## now let's plot variance in Moran's I by cum incidence
ggplot(sim_94_f_25, aes(x = mean_iso, y = CI)) + geom_point(colour = "blue") + theme_classic() +
  labs(title = "Cumulative Incidence after 1 year across 620 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 94% (10% between cell transmission)", x = "Variance of Moran's I", y = "Cumulative Incidence", color = "Legend") 

#now we have restricted long-form dataset ot work with.
#first summarize total cases per quadrant of seed case
sim_94_Q1 <- sim_94_f_25[sim_94_f_25$run == 1,]
sim_94_Q2 <- sim_94_f_25[sim_94_f_25$run == 2,]
sim_94_Q3 <- sim_94_f_25[sim_94_f_25$run == 3,]
sim_94_Q4 <- sim_94_f_25[sim_94_f_25$run == 4,]
Q1 <- c(1, min(sim_94_Q1$CI), median(sim_94_Q1$CI), mean(sim_94_Q1$CI), max(sim_94_Q1$CI)) 
Q2 <- c(2, min(sim_94_Q2$CI), median(sim_94_Q2$CI), mean(sim_94_Q2$CI), max(sim_94_Q2$CI)) 
Q3 <- c(3, min(sim_94_Q3$CI), median(sim_94_Q3$CI), mean(sim_94_Q3$CI), max(sim_94_Q3$CI)) 
Q4 <- c(4, min(sim_94_Q4$CI), median(sim_94_Q4$CI), mean(sim_94_Q4$CI), max(sim_94_Q4$CI)) 

sum_cases <- as.data.frame(rbind(Q1, Q2, Q3, Q4))
names(sum_cases) <- c("run", "min", "median", "mean", "max")

########################################################################################################################
###############################################  75%       #############################################################

#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")

sim_94_75_percent <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Run 1 94% Vaccination/simulation_summary_data_94_75_percent.csv")

#change X into time variable
sim_94_75_percent$time <- as.numeric(substring(sim_94_75_percent$X, regexpr("time ", sim_94_75_percent$X) + 5))

# #clear out X" column
drops <- c("X")
sim_94_75_percent <- sim_94_75_percent[ , !(names(sim_94_75_percent) %in% drops)]

#merge with correct file for summary moran stats
#first - upload the moran key file
moran_list_94<- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_summary_94.csv")
drops <- c("X", "sd_Moran")
moran_list_94 <- moran_list_94[ , !(names(moran_list_94) %in% drops)]

#now merge 
sim_94_75_percent <- merge(sim_94_75_percent, moran_list_94, by = "motif")

sim_94_f_75 <- sim_94_75_percent[sim_94_75_percent$time == 365,]

#load in isolation index file
iso_94 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_94.csv")
drops <- c("X", "level_1", "level_2", "level_3", "level_4")
iso_94 <- iso_94[ , !(names(iso_94) %in% drops)]

sim_94_f_75 <- merge(sim_94_f_75, iso_94, by = "motif")

#clear out two unnamed "X" columns and all "R" columns
drops <- c("X", "time")
sim_94_f_75 <- sim_94_f_75[ , !(names(sim_94_f_75) %in% drops)]

#write out file
#write.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/94% Vaccination/final_time_dataset_94%_coverage_75%_between_percent.csv")


sim_94_f_75$CI <- 15359 - sim_94_f_75$Susceptibles

#total number of cases 
summary(sim_94_f_75$CI)


## now let's plot variance in Moran's I by cum incidence
ggplot(sim_94_f_75, aes(x = mean_iso, y = CI)) + geom_point(colour = "blue") + theme_classic() +
  labs(title = "Cumulative Incidence after 1 year across 620 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 94% (10% between cell transmission)", x = "Variance of Moran's I", y = "Cumulative Incidence", color = "Legend") 

#now we have restricted long-form dataset ot work with.
#first summarize total cases per quadrant of seed case
sim_94_Q1 <- sim_94_f_75[sim_94_f_75$run == 1,]
sim_94_Q2 <- sim_94_f_75[sim_94_f_75$run == 2,]
sim_94_Q3 <- sim_94_f_75[sim_94_f_75$run == 3,]
sim_94_Q4 <- sim_94_f_75[sim_94_f_75$run == 4,]
Q1 <- c(1, min(sim_94_Q1$CI), median(sim_94_Q1$CI), mean(sim_94_Q1$CI), max(sim_94_Q1$CI)) 
Q2 <- c(2, min(sim_94_Q2$CI), median(sim_94_Q2$CI), mean(sim_94_Q2$CI), max(sim_94_Q2$CI)) 
Q3 <- c(3, min(sim_94_Q3$CI), median(sim_94_Q3$CI), mean(sim_94_Q3$CI), max(sim_94_Q3$CI)) 
Q4 <- c(4, min(sim_94_Q4$CI), median(sim_94_Q4$CI), mean(sim_94_Q4$CI), max(sim_94_Q4$CI)) 

sum_cases <- as.data.frame(rbind(Q1, Q2, Q3, Q4))
names(sum_cases) <- c("run", "min", "median", "mean", "max")

quad <- c(1,2,3,4)
mean_10 <- c(8042, 7238, 7241, 6793)
mean_25 <- c(7594, 6855, 6855, 6459)
mean_50 <- c(6250, 5256, 5466, 4748)
mean_75 <- c(3610, 1931, 1967, 1168)

means <- as.data.frame(cbind(quad, mean_10, mean_25, mean_50, mean_75))
names(means) <- c("quad", "ten", "twentyfive", "fifty", "seventyfive")

library(ggplot2)
ggplot(means, aes(x = as.factor(quad))) +
  geom_point(aes(y = ten), col = "green", size = 5) + geom_hline(yintercept = 7239, col = "green", size = 2) +
  geom_point(aes(y = twentyfive), col = "yellow", size = 5) +geom_hline(yintercept = 6941, col = "yellow", size = 2) +
  geom_point(aes(y = fifty), col = "orange", size = 5) + geom_hline(yintercept= 5382, col = "orange", size = 2) +
  geom_point(aes(y = seventyfive), col = "red", size = 5) + geom_hline(yintercept = 2169, col = "red", size = 2) + theme_dark() +
  labs(title = "Effect of varying between cell transmission percent on case burden: Overall Vaccination Coverage at 94%", x = "Quadrant of Seed Case", y = "Cumulative Incidence", color = "Legend") +
  ggplot2::scale_color_manual(values = c("ten percent" = "green", "twenty five percent" =  "yellow", "fifty percent" = "orange","seventyfive percent" =  "red"), labels = c(10, 25, 50, 75))

