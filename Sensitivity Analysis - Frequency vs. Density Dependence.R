###################### Sensitivity analysis of looking at impact of varying frequency vs. density dependent transmission
###################### look at impact on total cases and also what the regression equations are...


####################### 94% #####################
#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")

sim_94_density <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Run 1 94% Vaccination/simulation_summary_data_94_density_dependent_foi.csv")

#change X into time variable
sim_94_density$time <- as.numeric(substring(sim_94_density$X, regexpr("time ", sim_94_density$X) + 5))

# #clear out X" column
drops <- c("X")
sim_94_density <- sim_94_density[ , !(names(sim_94_density) %in% drops)]

sim_94_density_f <- sim_94_density[sim_94_density$time == 365,]

sim_94_density_f$CI <- 15359 - sim_94_density_f$Susceptibles

#total number of cases 
summary(sim_94_density_f$CI)

#now we have restricted long-form dataset ot work with.
#first summarize total cases per quadrant of seed case
Q1 <- sim_94_density_f[sim_94_density_f$run == 1,]
Q2 <- sim_94_density_f[sim_94_density_f$run == 2,]
Q3 <- sim_94_density_f[sim_94_density_f$run == 3,]
Q4 <- sim_94_density_f[sim_94_density_f$run == 4,]
Q1 <- c(1, median(Q1$CI), mean(Q1$CI)) 
Q2 <- c(2, median(Q2$CI), mean(Q2$CI)) 
Q3 <- c(3, median(Q3$CI), mean(Q3$CI)) 
Q4 <- c(4, median(Q4$CI), mean(Q4$CI)) 


sum_cases <- as.data.frame(rbind(Q1, Q2, Q3, Q4))
names(sum_cases) <- c("run", "median", "mean")

sum_cases$pct <- sum_cases$mean/15359*100


####################### 95% #####################
#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")

sim_95_density <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/95% Vaccination/simulation_summary_data_95_density_dependent_foi.csv")

#change X into time variable
sim_95_density$time <- as.numeric(substring(sim_95_density$X, regexpr("time ", sim_95_density$X) + 5))

# #clear out X" column
drops <- c("X")
sim_95_density <- sim_95_density[ , !(names(sim_95_density) %in% drops)]

sim_95_density_f <- sim_95_density[sim_95_density$time == 365,]

sim_95_density_f$CI <- 12799 - sim_95_density_f$Susceptibles

#total number of cases 
summary(sim_95_density_f$CI)

#now we have restricted long-form dataset ot work with.
#first summarize total cases per quadrant of seed case
Q1 <- sim_95_density_f[sim_95_density_f$run == 1,]
Q2 <- sim_95_density_f[sim_95_density_f$run == 2,]
Q3 <- sim_95_density_f[sim_95_density_f$run == 3,]
Q4 <- sim_95_density_f[sim_95_density_f$run == 4,]
Q1 <- c(1, median(Q1$CI), mean(Q1$CI)) 
Q2 <- c(2, median(Q2$CI), mean(Q2$CI)) 
Q3 <- c(3, median(Q3$CI), mean(Q3$CI)) 
Q4 <- c(4, median(Q4$CI), mean(Q4$CI)) 


sum_cases <- as.data.frame(rbind(Q1, Q2, Q3, Q4))
names(sum_cases) <- c("run", "median", "mean")

sum_cases$pct <- sum_cases$mean/12799*100

####################### 98% #####################
#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")

sim_98_density <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/98% Vaccination/simulation_summary_data_98_density_dependent_foi.csv")

#change X into time variable
sim_98_density$time <- as.numeric(substring(sim_98_density$X, regexpr("time ", sim_98_density$X) + 5))

# #clear out X" column
drops <- c("X")
sim_98_density <- sim_98_density[ , !(names(sim_98_density) %in% drops)]

sim_98_density_f <- sim_98_density[sim_98_density$time == 365,]

sim_98_density_f$CI <- 5120 - sim_98_density_f$Susceptibles

#total number of cases 
summary(sim_98_density_f$CI)

#now we have restricted long-form dataset ot work with.
#first summarize total cases per quadrant of seed case
Q1 <- sim_98_density_f[sim_98_density_f$run == 1,]
Q2 <- sim_98_density_f[sim_98_density_f$run == 2,]
Q3 <- sim_98_density_f[sim_98_density_f$run == 3,]
Q4 <- sim_98_density_f[sim_98_density_f$run == 4,]
Q1 <- c(1, median(Q1$CI), mean(Q1$CI)) 
Q2 <- c(2, median(Q2$CI), mean(Q2$CI)) 
Q3 <- c(3, median(Q3$CI), mean(Q3$CI)) 
Q4 <- c(4, median(Q4$CI), mean(Q4$CI)) 


sum_cases <- as.data.frame(rbind(Q1, Q2, Q3, Q4))
names(sum_cases) <- c("run", "median", "mean")

sum_cases$pct <- sum_cases$mean/5120*100


####################### 99% #####################
#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")

sim_99_density <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/99% Vaccination/simulation_summary_data_99_density_dependent_foi.csv")

#change X into time variable
sim_99_density$time <- as.numeric(substring(sim_99_density$X, regexpr("time ", sim_99_density$X) + 5))

# #clear out X" column
drops <- c("X")
sim_99_density <- sim_99_density[ , !(names(sim_99_density) %in% drops)]

sim_99_density_f <- sim_99_density[sim_99_density$time == 365,]

sim_99_density_f$CI <- 2560 - sim_99_density_f$Susceptibles

#total number of cases 
summary(sim_99_density_f$CI)

#now we have restricted long-form dataset ot work with.
#first summarize total cases per quadrant of seed case
Q1 <- sim_99_density_f[sim_99_density_f$run == 1,]
Q2 <- sim_99_density_f[sim_99_density_f$run == 2,]
Q3 <- sim_99_density_f[sim_99_density_f$run == 3,]
Q4 <- sim_99_density_f[sim_99_density_f$run == 4,]
Q1 <- c(1, median(Q1$CI), mean(Q1$CI)) 
Q2 <- c(2, median(Q2$CI), mean(Q2$CI)) 
Q3 <- c(3, median(Q3$CI), mean(Q3$CI)) 
Q4 <- c(4, median(Q4$CI), mean(Q4$CI)) 


sum_cases <- as.data.frame(rbind(Q1, Q2, Q3, Q4))
names(sum_cases) <- c("run", "median", "mean")

sum_cases$pct <- sum_cases$mean/2560*100
