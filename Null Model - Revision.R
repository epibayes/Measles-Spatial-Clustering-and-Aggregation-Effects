######## Sensitivity Analysis: Null Model

############# July 13, 2020
############# Response to reviewer comments

#load necessary packages
library(ggplot2)
library(maps)
library(mapdata)
library(maptools)
library(gridExtra)
library(here)
library(tidyverse)
library(spdep)
require(dplyr)
library(RANN)
library(rgdal)
library(RColorBrewer)
library(knitr)
library(magrittr)

require(readr)
require(dplyr)

#pull in neighbor defining function
source("/Users/ninamasters/measles-spatial-model/Neighbor_Functions.R")

setwd("/Users/ninamasters/measles-spatial-model")

#### #make simplified environment ####
#define function to make nxn grid of spatial points data frame and spatial polygon data frame
make_me_a_grid <- function(n){
  x <- seq(1,n, by = 1)
  y <- seq(1,n, by = 1)
  num_cells = length(x)*length(y)
  xy <<- expand.grid(x=x, y=y)
  
  grid.pts<<-SpatialPointsDataFrame(coords= xy, data=xy)
  #make points a gridded object
  gridded(grid.pts) <- TRUE
  #plot(grid.pts)
  
  grid <- as(grid.pts, "SpatialPolygons") #encode as spatial polygons
  #str(grid)
  
  gridspdf <<- SpatialPolygonsDataFrame(grid, data=data.frame(ID=row.names(grid), row.names=row.names(grid)))
  return(grid.pts)
}

#make grid of 16 x 16 to get 256-square cell grid
make_me_a_grid(16)

gridspdf$ID <- sapply(slot(gridspdf, "polygons"), function(x) slot(x, "ID"))

# set up adjacencies / boundaries using both queen and rook style boundaries
queen_boundaries_grid <- poly2nb(gridspdf, queen = T)

####################################################################################
# Make Randomly generated 'null model' for 94% (15,360 non-vaccinators), 95% (12,800 non-vaccinators), 98% (5120 non-vaccinators) and 99% (2560 non-vaccinators)
# swap in the relevant code to generate each and then run SIR model on it for null model results

##### 94% ######
#create distribution of non-vaccinators in each quadrant at the individual cell level (n = 256) 
vax_motif <- data.frame(cells = c(1:256), probability = rep(0.06/256, 256)) %>%
  mutate(nonvax_percell = rmultinom(1, 15360, probability)) %>%
  mutate(percent_nonvax = nonvax_percell/1000*100, vax_percell = 1000-nonvax_percell, infected_percell = 0) %>%
  select(cells, percent_nonvax, nonvax_percell, vax_percell, infected_percell) %>%
  mutate(ID = paste("g", cells, sep=""))

##### 95% ######
vax_motif <- data.frame(cells = c(1:256), probability = rep(0.05/256, 256)) %>%
  mutate(nonvax_percell = rmultinom(1, 12800, probability)) %>%
  mutate(percent_nonvax = nonvax_percell/1000*100, vax_percell = 1000-nonvax_percell, infected_percell = 0) %>%
  select(cells, percent_nonvax, nonvax_percell, vax_percell, infected_percell) %>%
  mutate(ID = paste("g", cells, sep=""))

##### 98% ######
vax_motif <- data.frame(cells = c(1:256), probability = rep(0.02/256, 256)) %>%
  mutate(nonvax_percell = rmultinom(1, 5120, probability)) %>%
  mutate(percent_nonvax = nonvax_percell/1000*100, vax_percell = 1000-nonvax_percell, infected_percell = 0) %>%
  select(cells, percent_nonvax, nonvax_percell, vax_percell, infected_percell) %>%
  mutate(ID = paste("g", cells, sep=""))

##### 99% ######
vax_motif <- data.frame(cells = c(1:256), probability = rep(0.01/256, 256)) %>%
  mutate(nonvax_percell = rmultinom(1, 2560, probability)) %>%
  mutate(percent_nonvax = nonvax_percell/1000*100, vax_percell = 1000-nonvax_percell, infected_percell = 0) %>%
  select(cells, percent_nonvax, nonvax_percell, vax_percell, infected_percell) %>%
  mutate(ID = paste("g", cells, sep=""))

#rename for vax motif
names(vax_motif) <- c("cell", "percent_nonvax", "So", "Ro", "Io", "ID")

#check that we have the correct total
sum(vax_motif$So)
sum(vax_motif$Io)  
sum(vax_motif$Ro)

# now merge with gridspdf, the spatial polygon data frame containing our grid information
vax_motif_initial <<- merge(gridspdf,vax_motif, by="ID")


##################################################### 
################  run model
#####################################################

summary_list = list()

for (j in 1:4){
vax_motif_grid <- vax_motif_initial
summary_list[[paste("quadrant", j)]] <- list()
  
#now randomly drop in one case into the four different test locations: one in each quadrant
#quadrant_infection vector has quad 1, 2, 3, 4
quadrant_infection <- c(196, 205, 52, 61)
names(quadrant_infection) <- c("quad1", "quad2", "quad3", "quad4")

#add in one infected case in that cell
vax_motif_grid@data$Io[quadrant_infection[j]] <- 1

#now remove one susceptible in that cell (i.e. the susceptible person got infected)
vax_motif_grid@data$So[quadrant_infection[j]] <- vax_motif_grid@data$So[quadrant_infection[j]]-1


# Define Model Parameters and initialize states

R0 <- 16 #generate within-cell Ro (16 for measles)
between_trans_percent <- 0.5  #how much transmission is happening between cells
gamma <- 1/14 #2 week recovery, in days
within_beta <- R0*(1-between_trans_percent)*gamma #internal force of infection (w/in cell)
between_beta <- R0*between_trans_percent*gamma #external force of infection (b/t cells)

#generate neighbors list and neighbors matrix (both with spatial weights) from the queen neighbors generated in clustering function
queen_matrix <- nb2mat(queen_boundaries_grid, zero.policy = TRUE, style = "W")
queen_list <- nb2listw(queen_boundaries_grid)

#format this into classic adjacency list with "from' and "to" ids
queen_adj <- purrr::map_dfr(queen_list$neighbours, function(x) data.frame(to = as.numeric(x)), .id = "from_id")
queen_adj <- mutate(queen_adj, "from_id" = as.numeric(queen_adj$from_id), "to_id" = as.numeric(queen_adj$to))
queen_adj <- select(queen_adj, from_id, to_id)

## Get mean number of neighbors for weight (i.e. those on edges, corners will have fewer neighbors)
num_neighbors <- queen_adj %>%
  group_by(from_id) %>%
  summarize(Num_neighbs = n()) 

#order 
num_neighbors[order(num_neighbors$from_id, decreasing = FALSE),] 

#join queen_adj with the number of neighbors
queen_adj <- dplyr::left_join(queen_adj, num_neighbors, by = "from_id")

#### set time for the simulation to run, 1 year, in days ####
t <- 365

total_I <- rep(0, t+1)
total_I[1] <- sum(vax_motif_grid@data$Io)

total_S <- rep(sum(vax_motif_grid@data$So), t+1) # begin with total susceptible population
total_R <- rep(sum(vax_motif_grid@data$Ro), t+1) # begin with total susceptible population

#define rho which depends upon the number of neighbors so we are scaling transmission by how many neighbors
rho <- 1 / (queen_adj$Num_neighbs)

out_A <- matrix(0, nrow = 256, ncol = 256)

#for loop to generate neighbor weights matrix
for (k in 1:nrow(queen_adj)) {
  out_A[queen_adj$from_id[k], queen_adj$to_id[k]] <- rho[k]
}

for (i in 1:t) {
  
  ## Sample potential imports
  outside_exposure <- (out_A %*% as.numeric(vax_motif_grid@data$Io)) 
  
  ## Get household outside exposure
  cell_foi <- (between_beta * outside_exposure / (1000*num_neighbors$Num_neighbs)) + (within_beta * as.numeric(vax_motif_grid@data$Io)/1000) 
  #print(cell_foi)
  
  ## New infections
  
  vax_motif_grid@data$So <- vax_motif_grid@data$So - cell_foi*vax_motif_grid@data$So
  vax_motif_grid@data$Io <- vax_motif_grid@data$Io - gamma*vax_motif_grid@data$Io + cell_foi*vax_motif_grid@data$So 
  vax_motif_grid@data$Ro <- vax_motif_grid@data$Ro + gamma*vax_motif_grid@data$Io
  
  total_S[i+1] <- sum(vax_motif_grid@data$So)  
  total_I[i+1] <- sum(vax_motif_grid@data$Io)
  total_R[i+1] <- sum(vax_motif_grid@data$Ro)
  
  # calculate the Moran's I
  w <- poly2nb(vax_motif_grid, queen = T, row.names = vax_motif_grid$ID)
  class(w)
  
  #use style = "W" for row-standardized values of Moran's I
  wm <- nb2mat(w, style='W', zero.policy = TRUE)
  ww <- nb2listw(w, style = 'W', zero.policy =TRUE)
  
  m_S <- moran(vax_motif_grid$So, ww, length(ww$neighbours), S0=Szero(ww))
  m_I <- moran(vax_motif_grid$Io, ww, length(ww$neighbours), S0=Szero(ww))
  
  #save output as nested list with p motifs, j runs, i time points
  summary_list[[paste("quadrant",j)]][[paste("time",i)]]  <-  list(run = j, Prevalence = total_I[i], Recovered = total_R[i], Susceptibles = total_S[i], Moran_S = m_S$I, Moran_Infected = m_I$I)
  
}
}

library(rlist)
sim_94 <- unlist(summary_list, recursive = FALSE, use.names = TRUE)
sim_94 <- as.data.frame(do.call(Map, c(f = rbind, sim_94)))
write.csv(sim_94, "/Users/ninamasters/measles-spatial-model/simulation_summary_data_null_model_94_percent_overall.csv")

sim_95 <- unlist(summary_list, recursive = FALSE, use.names = TRUE)
sim_95 <- as.data.frame(do.call(Map, c(f = rbind, sim_95)))
write.csv(sim_95, "/Users/ninamasters/measles-spatial-model/simulation_summary_data_null_model_95_percent_overall.csv")

sim_98 <- unlist(summary_list, recursive = FALSE, use.names = TRUE)
sim_98 <- as.data.frame(do.call(Map, c(f = rbind, sim_98)))
write.csv(sim_98, "/Users/ninamasters/measles-spatial-model/simulation_summary_data_null_model_98_percent_overall.csv")

sim_99 <- unlist(summary_list, recursive = FALSE, use.names = TRUE)
sim_99 <- as.data.frame(do.call(Map, c(f = rbind, sim_99)))
write.csv(sim_99, "/Users/ninamasters/measles-spatial-model/simulation_summary_data_null_model_99_percent_overall.csv")

######################################################################
##################### Load and Clean/Analyze Data
sim_94 <- read_csv("/Users/ninamasters/measles-spatial-model/simulation_summary_data_null_model_94_percent_overall.csv") 
#change X1 into time variable
sim_94 <- mutate(sim_94, time = as.numeric(substring(sim_94$X1, regexpr("time ", sim_94$X1) + 5))) %>%
select(run, Prevalence, Recovered, Susceptibles, time)

sim_95 <- read_csv("/Users/ninamasters/measles-spatial-model/simulation_summary_data_null_model_95_percent_overall.csv") 
sim_95 <- mutate(sim_95, time = as.numeric(substring(sim_95$X1, regexpr("time ", sim_95$X1) + 5))) %>%
  select(run, Prevalence, Recovered, Susceptibles, time)

sim_98 <- read_csv("/Users/ninamasters/measles-spatial-model/simulation_summary_data_null_model_98_percent_overall.csv") 
sim_98 <- mutate(sim_98, time = as.numeric(substring(sim_98$X1, regexpr("time ", sim_98$X1) + 5))) %>%
  select(run, Prevalence, Recovered, Susceptibles, time)

sim_99 <- read_csv("/Users/ninamasters/measles-spatial-model/simulation_summary_data_null_model_99_percent_overall.csv")
sim_99 <- mutate(sim_99, time = as.numeric(substring(sim_99$X1, regexpr("time ", sim_99$X1) + 5))) %>%
  select(run, Prevalence, Recovered, Susceptibles, time)

#check that time = 1 all start over / no carryover thru loops
sim_94_1 <- sim_94[sim_94$time == 1,]
sim_95_1 <- sim_95[sim_95$time == 1,]
sim_98_1 <- sim_98[sim_98$time == 1,]
sim_99_1 <- sim_99[sim_99$time == 1,]

sim_94_f <- sim_94[sim_94$time == 365,]
sim_95_f <- sim_95[sim_95$time == 365,]
sim_98_f <- sim_98[sim_98$time == 365,]
sim_99_f <- sim_99[sim_99$time == 365,]

#tabulate total final case counts into a dataframe
null_model_counts <- data.frame(vax_level = rep(c(94,95,98,99), each = 4), run = rep(c(1,2,3,4),4), S_end = rep(0,16), S_o = rep(c(15359,12799, 5119, 2559), each = 4))
null_model_counts$S_end[1] <- sim_94_f$Susceptibles[1]
null_model_counts$S_end[2] <- sim_94_f$Susceptibles[2]
null_model_counts$S_end[3] <- sim_94_f$Susceptibles[3]
null_model_counts$S_end[4] <- sim_94_f$Susceptibles[4]

null_model_counts$S_end[5] <- sim_95_f$Susceptibles[1]
null_model_counts$S_end[6] <- sim_95_f$Susceptibles[2]
null_model_counts$S_end[7] <- sim_95_f$Susceptibles[3]
null_model_counts$S_end[8] <- sim_95_f$Susceptibles[4]

null_model_counts$S_end[9] <- sim_98_f$Susceptibles[1]
null_model_counts$S_end[10] <- sim_98_f$Susceptibles[2]
null_model_counts$S_end[11] <- sim_98_f$Susceptibles[3]
null_model_counts$S_end[12] <- sim_98_f$Susceptibles[4]

null_model_counts$S_end[13] <- sim_99_f$Susceptibles[1]
null_model_counts$S_end[14] <- sim_99_f$Susceptibles[2]
null_model_counts$S_end[15] <- sim_99_f$Susceptibles[3]
null_model_counts$S_end[16] <- sim_99_f$Susceptibles[4]

null_model_counts$CI <- null_model_counts$S_o - null_model_counts$S_end

write_csv(null_model_counts, "/Users/ninamasters/measles-spatial-model/null_model_final_output_results.csv")
