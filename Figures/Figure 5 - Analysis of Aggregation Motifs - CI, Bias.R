# Analysis of Aggregation Files for Final Figure
#Just read in final time threshold subset file
agg_95_f <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/aggregation_end_time_data_95%_vax.csv")

# drop x
drops <- c("X")
agg_95_f <- agg_95_f[ , !(names(agg_95_f) %in% drops)]

#load in isolation index file
iso_95 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_95%_vax_aggregation_summary.csv")

agg_95_f <- merge(agg_95_f, iso_95, by = "motif")

#make holder dataset
aggregation_case_counts_block <- data.frame(vaccination = c("94", "95", "98", "99"), min = c(0,0,0,0), Q1 = c(0,0,0,0), median = c(0,0,0,0), mean= c(0,0,0,0), Q3= c(0,0,0,0), max = c(0,0,0,0))
aggregation_case_counts_tract <- data.frame(vaccination = c("94", "95", "98", "99"), min = c(0,0,0,0), Q1 = c(0,0,0,0), median = c(0,0,0,0), mean= c(0,0,0,0), Q3= c(0,0,0,0), max = c(0,0,0,0))
aggregation_case_counts_neighborhood <- data.frame(vaccination = c("94", "95", "98", "99"), min = c(0,0,0,0), Q1 = c(0,0,0,0), median = c(0,0,0,0), mean= c(0,0,0,0), Q3= c(0,0,0,0), max = c(0,0,0,0))
aggregation_case_counts_quadrant <- data.frame(vaccination = c("94", "95", "98", "99"), min = c(0,0,0,0), Q1 = c(0,0,0,0), median = c(0,0,0,0), mean= c(0,0,0,0), Q3= c(0,0,0,0), max = c(0,0,0,0))

# create new columns for cumulative incidence
agg_95_f$CI_l1 <- 12799- agg_95_f$S1
agg_95_f$CI_l2 <- 12799- agg_95_f$S2
agg_95_f$CI_l3 <- 12799- agg_95_f$S3
agg_95_f$CI_l4 <- 12799- agg_95_f$S4

#fill in table - for block level
aggregation_case_counts_block$min[2] <- summary(agg_95_f$CI_l1)[1]
aggregation_case_counts_block$Q1[2] <- summary(agg_95_f$CI_l1)[2]
aggregation_case_counts_block$median[2] <- summary(agg_95_f$CI_l1)[3]
aggregation_case_counts_block$mean[2] <- summary(agg_95_f$CI_l1)[4]
aggregation_case_counts_block$Q3[2] <- summary(agg_95_f$CI_l1)[5]
aggregation_case_counts_block$max[2] <- summary(agg_95_f$CI_l1)[6]

#fill in table - for tract level
aggregation_case_counts_tract$min[2] <- summary(agg_95_f$CI_l2)[1]
aggregation_case_counts_tract$Q1[2] <- summary(agg_95_f$CI_l2)[2]
aggregation_case_counts_tract$median[2] <- summary(agg_95_f$CI_l2)[3]
aggregation_case_counts_tract$mean[2] <- summary(agg_95_f$CI_l2)[4]
aggregation_case_counts_tract$Q3[2] <- summary(agg_95_f$CI_l2)[5]
aggregation_case_counts_tract$max[2] <- summary(agg_95_f$CI_l2)[6]

#fill in table - for neighborhood level
aggregation_case_counts_neighborhood$min[2] <- summary(agg_95_f$CI_l3)[1]
aggregation_case_counts_neighborhood$Q1[2] <- summary(agg_95_f$CI_l3)[2]
aggregation_case_counts_neighborhood$median[2] <- summary(agg_95_f$CI_l3)[3]
aggregation_case_counts_neighborhood$mean[2] <- summary(agg_95_f$CI_l3)[4]
aggregation_case_counts_neighborhood$Q3[2] <- summary(agg_95_f$CI_l3)[5]
aggregation_case_counts_neighborhood$max[2] <- summary(agg_95_f$CI_l3)[6]

#fill in table - for quadrant level
aggregation_case_counts_quadrant$min[2] <- summary(agg_95_f$CI_l4)[1]
aggregation_case_counts_quadrant$Q1[2] <- summary(agg_95_f$CI_l4)[2]
aggregation_case_counts_quadrant$median[2] <- summary(agg_95_f$CI_l4)[3]
aggregation_case_counts_quadrant$mean[2] <- summary(agg_95_f$CI_l4)[4]
aggregation_case_counts_quadrant$Q3[2] <- summary(agg_95_f$CI_l4)[5]
aggregation_case_counts_quadrant$max[2] <- summary(agg_95_f$CI_l4)[6]

#create colums for difference
agg_95_f$CI_l2_diff <- -(agg_95_f$CI_l1-agg_95_f$CI_l2)
agg_95_f$CI_l3_diff <- -(agg_95_f$CI_l1-agg_95_f$CI_l3)
agg_95_f$CI_l4_diff <- -(agg_95_f$CI_l1-agg_95_f$CI_l4)


agg_95_f$CI_l2_diff_pct <- (agg_95_f$CI_l1-agg_95_f$CI_l2)/agg_95_f$CI_l1
agg_95_f$CI_l3_diff_pct <- (agg_95_f$CI_l1-agg_95_f$CI_l3)/agg_95_f$CI_l1
agg_95_f$CI_l4_diff_pct <- (agg_95_f$CI_l1-agg_95_f$CI_l4)/agg_95_f$CI_l1

#create column for percentage of cases identified
agg_95_f$CI_l2_prop <- agg_95_f$CI_l2/agg_95_f$CI_l1
agg_95_f$CI_l3_prop <- agg_95_f$CI_l3/agg_95_f$CI_l1
agg_95_f$CI_l4_prop <- agg_95_f$CI_l4/agg_95_f$CI_l1

######################
#now turn level into factor variables
agg_95_f$l_1 <- factor(agg_95_f$level_1, order = "TRUE", levels = c("0.25", "0.4", "0.58", "0.7", "0.85"))
agg_95_f$l_2 <- factor(agg_95_f$level_2, order = "TRUE", levels = c("0.25", "0.4", "0.58", "0.7", "0.85"))
agg_95_f$l_3 <- factor(agg_95_f$level_3, order = "TRUE", levels = c("0.25", "0.4", "0.58", "0.7", "0.85"))
agg_95_f$l_4 <- factor(agg_95_f$level_4, order = "TRUE", levels = c("0.25", "0.4", "0.58", "0.7", "0.85"))

library(ggplot2)
library(ggthemes)

## now let's plot Isolation by cum incidence
ggplot(agg_95_f, aes(x = mean_iso_1, y = CI_l1)) + geom_point(colour = "blue") + theme_classic() +
  geom_point(aes(x = mean_iso_1+.1, y = CI_l2), colour = "green") + 
  geom_point(aes(x = mean_iso_1+.2, y = CI_l3), colour = "red") +
  geom_point(aes(x = mean_iso_1+.3, y = CI_l4), colour = "purple") +
  labs(title = "Cumulative Incidence across 336 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 95%: Aggregation Effects", x = "Isolation Index of Starting Motif (true data/level 1)", y = "Cumulative Incidence", color = "Legend") 

## now let's plot Isolation by cum incidence
ggplot(agg_95_f, aes(x = mean_iso_1, y = CI_l1)) + geom_smooth(method = "loess", colour = "black", se = FALSE) + theme_clean() +
  geom_smooth(method = "loess", aes(x = mean_iso_1, y = CI_l2), colour = "red", se = FALSE) + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, y = CI_l3), colour = "blueviolet", se = FALSE) +
  geom_smooth(method = "loess", aes(x = mean_iso_1, y = CI_l4), colour = "deepskyblue", se = FALSE) +
  geom_point(aes(x = mean_iso_1, y = CI_l1), colour = "black", alpha = 1/4) +
  geom_point(aes(x = mean_iso_1, y = CI_l2), colour = "red", alpha = 1/4) + 
  geom_point(aes(x = mean_iso_1, y = CI_l3), colour = "blueviolet", alpha = 1/4) +
  geom_point(aes(x = mean_iso_1, y = CI_l4), colour = "deepskyblue", alpha = 1/4) +
  labs(title = "Cumulative Incidence across 336 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 95%: Aggregation Effects", x = "Isolation Index of Starting Motif (true data/level 1)", y = "Cumulative Incidence", color = "Legend") 


## now let's plot by bias
ggplot(agg_95_f, aes(x = mean_iso_1, y = CI_l2_diff_pct)) + geom_smooth(method = "loess", colour = "red") + theme_classic() +
  geom_abline(aes(slope = 0, intercept = 0), colour = "grey") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l3_diff_pct), colour = "blueviolet") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l4_diff_pct), colour = "deepskyblue") +
  geom_point(aes(x = mean_iso_1, y = CI_l2_diff_pct), colour = "red", alpha = 1/10) + 
  geom_point(aes(x = mean_iso_1, y = CI_l3_diff_pct), colour = "blueviolet", alpha = 1/10) +
  geom_point(aes(x = mean_iso_1, y = CI_l4_diff_pct), colour = "deepskyblue", alpha = 1/10) +
  labs(title = "Cumulative Incidence across 336 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 95%: Aggregation Effects", x = "Isolation Index of Starting Motif (true data/level 1)", y = "Percent Difference in Estimated Cumulative Incidence", color = "Legend") 


########################### FIGURE 5
## now let's plot proportion of cases detected
ggplot(agg_95_f, aes(x = mean_iso_1, y = CI_l2_prop)) + geom_smooth(method = "loess", colour = "red", se = FALSE) + theme_classic() +
  geom_abline(aes(slope = 0, intercept = 1), colour = "grey") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l3_prop),se = FALSE, colour = "blueviolet") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l4_prop),se = FALSE, colour = "deepskyblue") +
  geom_point(aes(x = mean_iso_1, y = CI_l2_prop), colour = "red", alpha = 1/6, size = 2.5) + 
  geom_point(aes(x = mean_iso_1, y = CI_l3_prop), colour = "blueviolet", alpha = 1/6, size = 2.5) +
  geom_point(aes(x = mean_iso_1, y = CI_l4_prop), colour = "deepskyblue", alpha = 1/6, size = 2.5) +
  labs(title = "Proportion of Cases Predicted across 336 Motifs: Overall Vaccination Coverage at 95%: Aggregation Effects", x = "Isolation Index of Starting Motif (true data/level 1)", y = "Prportion of Estimated Cumulative Incidence Predicted", color = "Legend") 

## now let's plot by bias
ggplot(agg_95_f, aes(x = mean_iso_1, y = CI_l2_diff)) + geom_smooth(method = "loess", colour = "red", se = FALSE) + theme_clean() +
  geom_abline(aes(slope = 0, intercept = 0), colour = "grey") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l3_diff), se = FALSE, colour = "blueviolet") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l4_diff), se = FALSE, colour = "deepskyblue") +
  geom_point(aes(x = mean_iso_1, y = CI_l2_diff), colour = "red", alpha = 1/6, size = 2.5) + 
  geom_point(aes(x = mean_iso_1, y = CI_l3_diff), colour = "blueviolet", alpha = 1/6, size = 2.5) +
  geom_point(aes(x = mean_iso_1, y = CI_l4_diff), colour = "deepskyblue", alpha = 1/6, size =2.5) +
  labs(title = "Cumulative Incidence across 336 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 95%: Aggregation Effects", x = "Isolation Index of Starting Motif (true data/level 1)", y = "Difference in Estimated Cumulative Incidence", color = "Legend") 

# Figure S11
ggplot(agg_95_f, aes(x = mean_iso_1, y = mean_iso_1)) + geom_smooth(method = "loess", se = FALSE, colour = "black") + theme_classic() +
  geom_smooth(method = "loess", aes(x = mean_iso_1, mean_iso_2), se = FALSE, colour = "red") +
  geom_smooth(method = "loess", aes(x = mean_iso_1, mean_iso_3), se = FALSE, colour = "blueviolet") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, mean_iso_4), se = FALSE, colour = "deepskyblue") +
  #geom_point(aes(x = mean_iso_1, y = mean_iso_2), colour = "red", alpha = 1/10) + 
  #geom_point(aes(x = mean_iso_1, y = mean_iso_3), colour = "blueviolet", alpha = 1/10) +
  #geom_point(aes(x = mean_iso_1, y = mean_iso_4), colour = "deepskyblue", alpha = 1/10) +
  labs(title = "Cumulative Incidence across 336 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 95%: Aggregation Effects", x = "Isolation Index of Starting Motif (true data/level 1)", y = "Isolation Index of Aggregated Motifs", color = "Legend") 


#################################### ANALYSIS OF 94% ##################################################
#Just read in final time threshold subset file
agg_94_f <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/aggregation_end_time_data_94%_vax.csv")

# drop x
drops <- c("X")
agg_94_f <- agg_94_f[ , !(names(agg_94_f) %in% drops)]

#load in isolation index file
iso_94 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_94%_vax_aggregation_summary.csv")

agg_94_f <- merge(agg_94_f, iso_94, by = "motif")

# create new columns for cumulative incidence
agg_94_f$CI_l1 <- 15359- agg_94_f$S1
agg_94_f$CI_l2 <- 15359- agg_94_f$S2
agg_94_f$CI_l3 <- 15359- agg_94_f$S3
agg_94_f$CI_l4 <- 15359- agg_94_f$S4


#fill in table - for block level
aggregation_case_counts_block$min[1] <- summary(agg_94_f$CI_l1)[1]
aggregation_case_counts_block$Q1[1] <- summary(agg_94_f$CI_l1)[2]
aggregation_case_counts_block$median[1] <- summary(agg_94_f$CI_l1)[3]
aggregation_case_counts_block$mean[1] <- summary(agg_94_f$CI_l1)[4]
aggregation_case_counts_block$Q3[1] <- summary(agg_94_f$CI_l1)[5]
aggregation_case_counts_block$max[1] <- summary(agg_94_f$CI_l1)[6]

#fill in table - for tract level
aggregation_case_counts_tract$min[1] <- summary(agg_94_f$CI_l2)[1]
aggregation_case_counts_tract$Q1[1] <- summary(agg_94_f$CI_l2)[2]
aggregation_case_counts_tract$median[1] <- summary(agg_94_f$CI_l2)[3]
aggregation_case_counts_tract$mean[1] <- summary(agg_94_f$CI_l2)[4]
aggregation_case_counts_tract$Q3[1] <- summary(agg_94_f$CI_l2)[5]
aggregation_case_counts_tract$max[1] <- summary(agg_94_f$CI_l2)[6]

#fill in table - for neighborhood level
aggregation_case_counts_neighborhood$min[1] <- summary(agg_94_f$CI_l3)[1]
aggregation_case_counts_neighborhood$Q1[1] <- summary(agg_94_f$CI_l3)[2]
aggregation_case_counts_neighborhood$median[1] <- summary(agg_94_f$CI_l3)[3]
aggregation_case_counts_neighborhood$mean[1] <- summary(agg_94_f$CI_l3)[4]
aggregation_case_counts_neighborhood$Q3[1] <- summary(agg_94_f$CI_l3)[5]
aggregation_case_counts_neighborhood$max[1] <- summary(agg_94_f$CI_l3)[6]

#fill in table - for quadrant level
aggregation_case_counts_quadrant$min[1] <- summary(agg_94_f$CI_l4)[1]
aggregation_case_counts_quadrant$Q1[1] <- summary(agg_94_f$CI_l4)[2]
aggregation_case_counts_quadrant$median[1] <- summary(agg_94_f$CI_l4)[3]
aggregation_case_counts_quadrant$mean[1] <- summary(agg_94_f$CI_l4)[4]
aggregation_case_counts_quadrant$Q3[1] <- summary(agg_94_f$CI_l4)[5]
aggregation_case_counts_quadrant$max[1] <- summary(agg_94_f$CI_l4)[6]

#create colums for difference
agg_94_f$CI_l2_diff <- -(agg_94_f$CI_l1-agg_94_f$CI_l2)
agg_94_f$CI_l3_diff <- -(agg_94_f$CI_l1-agg_94_f$CI_l3)
agg_94_f$CI_l4_diff <- -(agg_94_f$CI_l1-agg_94_f$CI_l4)

agg_94_f$CI_l2_diff_pct <- (agg_94_f$CI_l1-agg_94_f$CI_l2)/agg_94_f$CI_l1
agg_94_f$CI_l3_diff_pct <- (agg_94_f$CI_l1-agg_94_f$CI_l3)/agg_94_f$CI_l1
agg_94_f$CI_l4_diff_pct <- (agg_94_f$CI_l1-agg_94_f$CI_l4)/agg_94_f$CI_l1

#create column for percentage of cases identified
agg_94_f$CI_l2_prop <- agg_94_f$CI_l2/agg_94_f$CI_l1
agg_94_f$CI_l3_prop <- agg_94_f$CI_l3/agg_94_f$CI_l1
agg_94_f$CI_l4_prop <- agg_94_f$CI_l4/agg_94_f$CI_l1

######################
#now turn level into factor variables
agg_94_f$l_1 <- factor(agg_94_f$level_1, order = "TRUE", levels = c("0.25", "0.4", "0.58", "0.7", "0.85"))
agg_94_f$l_2 <- factor(agg_94_f$level_2, order = "TRUE", levels = c("0.25", "0.4", "0.58", "0.7", "0.85"))
agg_94_f$l_3 <- factor(agg_94_f$level_3, order = "TRUE", levels = c("0.25", "0.4", "0.58", "0.7", "0.85"))
agg_94_f$l_4 <- factor(agg_94_f$level_4, order = "TRUE", levels = c("0.25", "0.4", "0.58", "0.7", "0.85"))

library(ggplot2)
library(ggthemes)

####################### FIG 5 Sensitivity
## now let's plot proportion of cases detected
ggplot(agg_94_f, aes(x = mean_iso_1, y = CI_l2_prop)) + geom_smooth(method = "loess", colour = "red", se = FALSE) + theme_classic() +
  geom_abline(aes(slope = 0, intercept = 1), colour = "grey") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l3_prop),se = FALSE, colour = "blueviolet") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l4_prop),se = FALSE, colour = "deepskyblue") +
  geom_point(aes(x = mean_iso_1, y = CI_l2_prop), colour = "red", alpha = 1/6, size = 2.5) + 
  geom_point(aes(x = mean_iso_1, y = CI_l3_prop), colour = "blueviolet", alpha = 1/6, size = 2.5) +
  geom_point(aes(x = mean_iso_1, y = CI_l4_prop), colour = "deepskyblue", alpha = 1/6, size = 2.5) +
  labs(title = "Proportion of Cases Predicted across 296 Motifs: Overall Vaccination Coverage at 94%: Aggregation Effects", x = "Isolation Index of Starting Motif (true data/level 1)", y = "Prportion of Estimated Cumulative Incidence Predicted", color = "Legend") 

## now let's plot by bias
ggplot(agg_94_f, aes(x = mean_iso_1, y = CI_l2_diff)) + geom_smooth(method = "loess", colour = "red", se = FALSE) + theme_clean() +
  geom_abline(aes(slope = 0, intercept = 0), colour = "grey") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l3_diff), se = FALSE, colour = "blueviolet") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l4_diff), se = FALSE, colour = "deepskyblue") +
  geom_point(aes(x = mean_iso_1, y = CI_l2_diff), colour = "red", alpha = 1/6, size = 2.5) + 
  geom_point(aes(x = mean_iso_1, y = CI_l3_diff), colour = "blueviolet", alpha = 1/6, size = 2.5) +
  geom_point(aes(x = mean_iso_1, y = CI_l4_diff), colour = "deepskyblue", alpha = 1/6, size =2.5) +
  labs(title = "Cumulative Incidence across 296 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 94%: Aggregation Effects", x = "Isolation Index of Starting Motif (true data/level 1)", y = "Difference in Estimated Cumulative Incidence", color = "Legend") 

#################################### ANALYSIS OF 98% ##################################################

#Just read in final time threshold subset file
agg_98_f <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/aggregation_end_time_data_98%_vax.csv")

# drop x
drops <- c("X")
agg_98_f <- agg_98_f[ , !(names(agg_98_f) %in% drops)]

#load in isolation index file
iso_98 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_98%_vax_aggregation_summary.csv")

agg_98_f <- merge(agg_98_f, iso_98, by = "motif")

# create new columns for cumulative incidence
agg_98_f$CI_l1 <- 5120- agg_98_f$S1
agg_98_f$CI_l2 <- 5120- agg_98_f$S2
agg_98_f$CI_l3 <- 5120- agg_98_f$S3
agg_98_f$CI_l4 <- 5120- agg_98_f$S4

#fill in table - for block level
aggregation_case_counts_block$min[3] <- summary(agg_98_f$CI_l1)[1]
aggregation_case_counts_block$Q1[3] <- summary(agg_98_f$CI_l1)[2]
aggregation_case_counts_block$median[3] <- summary(agg_98_f$CI_l1)[3]
aggregation_case_counts_block$mean[3] <- summary(agg_98_f$CI_l1)[4]
aggregation_case_counts_block$Q3[3] <- summary(agg_98_f$CI_l1)[5]
aggregation_case_counts_block$max[3] <- summary(agg_98_f$CI_l1)[6]

#fill in table - for tract level
aggregation_case_counts_tract$min[3] <- summary(agg_98_f$CI_l2)[1]
aggregation_case_counts_tract$Q1[3] <- summary(agg_98_f$CI_l2)[2]
aggregation_case_counts_tract$median[3] <- summary(agg_98_f$CI_l2)[3]
aggregation_case_counts_tract$mean[3] <- summary(agg_98_f$CI_l2)[4]
aggregation_case_counts_tract$Q3[3] <- summary(agg_98_f$CI_l2)[5]
aggregation_case_counts_tract$max[3] <- summary(agg_98_f$CI_l2)[6]

#fill in table - for neighborhood level
aggregation_case_counts_neighborhood$min[3] <- summary(agg_98_f$CI_l3)[1]
aggregation_case_counts_neighborhood$Q1[3] <- summary(agg_98_f$CI_l3)[2]
aggregation_case_counts_neighborhood$median[3] <- summary(agg_98_f$CI_l3)[3]
aggregation_case_counts_neighborhood$mean[3] <- summary(agg_98_f$CI_l3)[4]
aggregation_case_counts_neighborhood$Q3[3] <- summary(agg_98_f$CI_l3)[5]
aggregation_case_counts_neighborhood$max[3] <- summary(agg_98_f$CI_l3)[6]

#fill in table - for quadrant level
aggregation_case_counts_quadrant$min[3] <- summary(agg_98_f$CI_l4)[1]
aggregation_case_counts_quadrant$Q1[3] <- summary(agg_98_f$CI_l4)[2]
aggregation_case_counts_quadrant$median[3] <- summary(agg_98_f$CI_l4)[3]
aggregation_case_counts_quadrant$mean[3] <- summary(agg_98_f$CI_l4)[4]
aggregation_case_counts_quadrant$Q3[3] <- summary(agg_98_f$CI_l4)[5]
aggregation_case_counts_quadrant$max[3] <- summary(agg_98_f$CI_l4)[6]


#create colums for difference
agg_98_f$CI_l2_diff <- -(agg_98_f$CI_l1-agg_98_f$CI_l2)
agg_98_f$CI_l3_diff <- -(agg_98_f$CI_l1-agg_98_f$CI_l3)
agg_98_f$CI_l4_diff <- -(agg_98_f$CI_l1-agg_98_f$CI_l4)


agg_98_f$CI_l2_diff_pct <- (agg_98_f$CI_l1-agg_98_f$CI_l2)/agg_98_f$CI_l1
agg_98_f$CI_l3_diff_pct <- (agg_98_f$CI_l1-agg_98_f$CI_l3)/agg_98_f$CI_l1
agg_98_f$CI_l4_diff_pct <- (agg_98_f$CI_l1-agg_98_f$CI_l4)/agg_98_f$CI_l1

#create column for percentage of cases identified
agg_98_f$CI_l2_prop <- agg_98_f$CI_l2/agg_98_f$CI_l1
agg_98_f$CI_l3_prop <- agg_98_f$CI_l3/agg_98_f$CI_l1
agg_98_f$CI_l4_prop <- agg_98_f$CI_l4/agg_98_f$CI_l1

####################### FIG 5 Sensitivity
## now let's plot proportion of cases detected
ggplot(agg_98_f, aes(x = mean_iso_1, y = CI_l2_prop)) + geom_smooth(method = "loess", colour = "red", se = FALSE) + theme_classic() +
  geom_abline(aes(slope = 0, intercept = 1), colour = "grey") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l3_prop),se = FALSE, colour = "blueviolet") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l4_prop),se = FALSE, colour = "deepskyblue") +
  geom_point(aes(x = mean_iso_1, y = CI_l2_prop), colour = "red", alpha = 1/6, size = 2.5) + 
  geom_point(aes(x = mean_iso_1, y = CI_l3_prop), colour = "blueviolet", alpha = 1/6, size = 2.5) +
  geom_point(aes(x = mean_iso_1, y = CI_l4_prop), colour = "deepskyblue", alpha = 1/6, size = 2.5) +
  ylim(c(0,2)) +
  labs(title = "Proportion of Cases Predicted across 543 Motifs: Overall Vaccination Coverage at 98%: Aggregation Effects", x = "Isolation Index of Starting Motif (true data/level 1)", y = "Prportion of Estimated Cumulative Incidence Predicted", color = "Legend") 

## now let's plot by bias
ggplot(agg_98_f, aes(x = mean_iso_1, y = CI_l2_diff)) + geom_smooth(method = "loess", colour = "red", se = FALSE) + theme_clean() +
  geom_abline(aes(slope = 0, intercept = 0), colour = "grey") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l3_diff), se = FALSE, colour = "blueviolet") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l4_diff), se = FALSE, colour = "deepskyblue") +
  geom_point(aes(x = mean_iso_1, y = CI_l2_diff), colour = "red", alpha = 1/6, size = 2.5) + 
  geom_point(aes(x = mean_iso_1, y = CI_l3_diff), colour = "blueviolet", alpha = 1/6, size = 2.5) +
  geom_point(aes(x = mean_iso_1, y = CI_l4_diff), colour = "deepskyblue", alpha = 1/6, size =2.5) +
  labs(title = "Cumulative Incidence across 543 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 98%: Aggregation Effects", x = "Isolation Index of Starting Motif (true data/level 1)", y = "Difference in Estimated Cumulative Incidence", color = "Legend") 


#################################### ANALYSIS OF 99% ##################################################

#Just read in final time threshold subset file
agg_99_f <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/aggregation_end_time_data_99%_vax.csv")

# drop x
drops <- c("X")
agg_99_f <- agg_99_f[ , !(names(agg_99_f) %in% drops)]

#load in isolation index file
iso_99 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_99%_vax_aggregation_summary.csv")

agg_99_f <- merge(agg_99_f, iso_99, by = "motif")

# create new columns for cumulative incidence
agg_99_f$CI_l1 <- 2560- agg_99_f$S1
agg_99_f$CI_l2 <- 2560- agg_99_f$S2
agg_99_f$CI_l3 <- 2560- agg_99_f$S3
agg_99_f$CI_l4 <- 2560- agg_99_f$S4

#fill in table - for block level
aggregation_case_counts_block$min[4] <- summary(agg_99_f$CI_l1)[1]
aggregation_case_counts_block$Q1[4] <- summary(agg_99_f$CI_l1)[2]
aggregation_case_counts_block$median[4] <- summary(agg_99_f$CI_l1)[3]
aggregation_case_counts_block$mean[4] <- summary(agg_99_f$CI_l1)[4]
aggregation_case_counts_block$Q3[4] <- summary(agg_99_f$CI_l1)[5]
aggregation_case_counts_block$max[4] <- summary(agg_99_f$CI_l1)[6]

#fill in table - for tract level
aggregation_case_counts_tract$min[4] <- summary(agg_99_f$CI_l2)[1]
aggregation_case_counts_tract$Q1[4] <- summary(agg_99_f$CI_l2)[2]
aggregation_case_counts_tract$median[4] <- summary(agg_99_f$CI_l2)[3]
aggregation_case_counts_tract$mean[4] <- summary(agg_99_f$CI_l2)[4]
aggregation_case_counts_tract$Q3[4] <- summary(agg_99_f$CI_l2)[5]
aggregation_case_counts_tract$max[4] <- summary(agg_99_f$CI_l2)[6]

#fill in table - for neighborhood level
aggregation_case_counts_neighborhood$min[4] <- summary(agg_99_f$CI_l3)[1]
aggregation_case_counts_neighborhood$Q1[4] <- summary(agg_99_f$CI_l3)[2]
aggregation_case_counts_neighborhood$median[4] <- summary(agg_99_f$CI_l3)[3]
aggregation_case_counts_neighborhood$mean[4] <- summary(agg_99_f$CI_l3)[4]
aggregation_case_counts_neighborhood$Q3[4] <- summary(agg_99_f$CI_l3)[5]
aggregation_case_counts_neighborhood$max[4] <- summary(agg_99_f$CI_l3)[6]

#fill in table - for quadrant level
aggregation_case_counts_quadrant$min[4] <- summary(agg_99_f$CI_l4)[1]
aggregation_case_counts_quadrant$Q1[4] <- summary(agg_99_f$CI_l4)[2]
aggregation_case_counts_quadrant$median[4] <- summary(agg_99_f$CI_l4)[3]
aggregation_case_counts_quadrant$mean[4] <- summary(agg_99_f$CI_l4)[4]
aggregation_case_counts_quadrant$Q3[4] <- summary(agg_99_f$CI_l4)[5]
aggregation_case_counts_quadrant$max[4] <- summary(agg_99_f$CI_l4)[6]

#create colums for difference
agg_99_f$CI_l2_diff <- -(agg_99_f$CI_l1-agg_99_f$CI_l2)
agg_99_f$CI_l3_diff <- -(agg_99_f$CI_l1-agg_99_f$CI_l3)
agg_99_f$CI_l4_diff <- -(agg_99_f$CI_l1-agg_99_f$CI_l4)

agg_99_f$CI_l2_diff_pct <- (agg_99_f$CI_l1-agg_99_f$CI_l2)/agg_99_f$CI_l1
agg_99_f$CI_l3_diff_pct <- (agg_99_f$CI_l1-agg_99_f$CI_l3)/agg_99_f$CI_l1
agg_99_f$CI_l4_diff_pct <- (agg_99_f$CI_l1-agg_99_f$CI_l4)/agg_99_f$CI_l1


#create column for percentage of cases identified
agg_99_f$CI_l2_prop <- agg_99_f$CI_l2/agg_99_f$CI_l1
agg_99_f$CI_l3_prop <- agg_99_f$CI_l3/agg_99_f$CI_l1
agg_99_f$CI_l4_prop <- agg_99_f$CI_l4/agg_99_f$CI_l1


####################### FIG 5 Sensitivity
## now let's plot proportion of cases detected
ggplot(agg_99_f, aes(x = mean_iso_1, y = CI_l2_prop)) + geom_smooth(method = "loess", colour = "red", se = FALSE) + theme_classic() +
  geom_abline(aes(slope = 0, intercept = 1), colour = "grey") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l3_prop),se = FALSE, colour = "blueviolet") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l4_prop),se = FALSE, colour = "deepskyblue") +
  geom_point(aes(x = mean_iso_1, y = CI_l2_prop), colour = "red", alpha = 1/6, size = 2.5) + 
  geom_point(aes(x = mean_iso_1, y = CI_l3_prop), colour = "blueviolet", alpha = 1/6, size = 2.5) +
  geom_point(aes(x = mean_iso_1, y = CI_l4_prop), colour = "deepskyblue", alpha = 1/6, size = 2.5) +
  labs(title = "Proportion of Cases Predicted across 620 Motifs: Overall Vaccination Coverage at 99%: Aggregation Effects", x = "Isolation Index of Starting Motif (true data/level 1)", y = "Prportion of Estimated Cumulative Incidence Predicted", color = "Legend") 

## now let's plot by bias
ggplot(agg_99_f, aes(x = mean_iso_1, y = CI_l2_diff)) + geom_smooth(metzhod = "loess", colour = "red", se = FALSE) + theme_clean() +
  geom_abline(aes(slope = 0, intercept = 0), colour = "grey") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l3_diff), se = FALSE, colour = "blueviolet") + 
  geom_smooth(method = "loess", aes(x = mean_iso_1, CI_l4_diff), se = FALSE, colour = "deepskyblue") #+
  geom_point(aes(x = mean_iso_1, y = CI_l2_diff), colour = "red", alpha = 1/6, size = 2.5) + 
  geom_point(aes(x = mean_iso_1, y = CI_l3_diff), colour = "blueviolet", alpha = 1/6, size = 2.5) +
  geom_point(aes(x = mean_iso_1, y = CI_l4_diff), colour = "deepskyblue", alpha = 1/6, size =2.5) +
  labs(title = "Cumulative Incidence across 620 Motifs that do not exceed cell-level population of 1000: Overall Vaccination Coverage at 99%: Aggregation Effects", x = "Isolation Index of Starting Motif (true data/level 1)", y = "Difference in Estimated Cumulative Incidence", color = "Legend") 




########### finalize + manipulate case output table
aggregation_case_counts_block$level <- c("block", "block", "block", "block")
aggregation_case_counts_tract$level <- c("tract", "tract", "tract", "tract")
aggregation_case_counts_neighborhood$level <- c("neighborhood", "neighborhood", "neighborhood", "neighborhood")
aggregation_case_counts_quadrant$level <- c("quadrant", "quadrant", "quadrant", "quadrant")

aggregation_case_counts <- rbind(aggregation_case_counts_block, aggregation_case_counts_tract, aggregation_case_counts_neighborhood, aggregation_case_counts_quadrant)
write.csv(aggregation_case_counts, "/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /aggregation_case_counts.csv")
