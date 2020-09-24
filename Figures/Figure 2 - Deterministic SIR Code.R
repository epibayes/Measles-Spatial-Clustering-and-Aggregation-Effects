# Write R Code for deterministic SIR model in spatial environment 
# Aim 1
# Last Updated Jan 7, 2020 
# Nina Masters
# NB: THis code is useful for running particular loops of code and generating plots when there is a specific motif desired

##################################################################################
#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
#pull in necessary functions to make initial grids
source("/Users/ninamasters/measles-spatial-model/Clustering_Functions_95.R")
#pull in neighbor defining function
source("/Users/ninamasters/measles-spatial-model/Neighbor_Functions.R")

##################################################################################

# load in initial conditions and generate a spatial polygon data frame
make_me_a_motif(x4,y,z4,q4)

#initialize states: S, I, R
#vaccinated people in each cell (i.e. initially in "R" state)
vax_motif_grid@data$Ro
#set to 0 for test case
#vax_motif_grid@data$Ro <- rep(0,256)

#non-vaccinated people in each cell (i.e. initially in "S" state)
vax_motif_grid@data$So
#set to 1000 for test case
#vax_motif_grid@data$So <- rep(1000, 256)

#initially have no infecteds
#vax_motif_grid@data$Io 
#set sequentially for test case
vax_motif_grid@data$Io

#now randomly drop in one case
# val_quad1 <- 196
# val_quad2 <- 205
# val_quad3 <- 52
# val_quad4 <- 61

val_quad <- 196
#alternatively, could do: sample(1:256, 1) # to randomly select cell in which we drop one case

#add in one infected case in that cell
vax_motif_grid@data$Io[val_quad] <- 1

#now remove one susceptible in that cell (i.e. the susceptible person got infected)
#don't do this for 98% and 99%, can create negatives - in these cases just add an infected in from outside
vax_motif_grid@data$So[val_quad] <- vax_motif_grid@data$So[val_quad]-1

#check that total population is 256,000
total_pop = sum(vax_motif_grid@data$Ro, vax_motif_grid@data$So, vax_motif_grid@data$Io) #256,000

###################################################################################################

# Define Model Parameters 
R0 <- 16 #generate within-cell Ro (16 for measles)
between_trans_percent <- 0.5  #how much transmission is happening between cells
gamma <- 1/14 #2 week recovery, in days
within_beta <- R0*(1-between_trans_percent)*gamma #internal force of infection (w/in cell)
between_beta <- R0*between_trans_percent*gamma #external force of infection (b/t cells)

#beta <- R0*(gamma) #in weeks, beta ~ 8

#generate neighbors list and neighbors matrix (both with spatial weights) from the queen neighbors generated in clustering function
queen_matrix <- nb2mat(queen_boundaries_grid, zero.policy = TRUE, style = "W")
queen_list <- nb2listw(queen_boundaries_grid)

#format this into classic adjacency list with "from' and "to" ids
queen_adj <- purrr::map_dfr(queen_list$neighbours, function(x) data.frame(to = as.numeric(x)), .id = "from_id")
queen_adj <- mutate(queen_adj, "from_id" = as.numeric(queen_adj$from_id), "to_id" = as.numeric(queen_adj$to))
queen_adj <- dplyr::select(queen_adj, from_id, to_id)

## Get mean number of neighbors for weight (i.e. those on edges, corners will have fewer neighbors)
num_neighbors <- queen_adj %>%
  group_by(from_id) %>%
  summarize(Num_neighbs = n()) 

#order 
num_neighbors[order(num_neighbors$from_id, decreasing = FALSE),] 

#join queen_adj with the number of neighbors
queen_adj <- left_join(queen_adj, num_neighbors, by = "from_id")

#### set time for the simulation to run ####
# time step in days
t <- 365

incidence <- rep(0, t+1)
total_I <- rep(0, t+1)
total_I[1] <- sum(vax_motif_grid@data$Io)

total_S <- rep(sum(vax_motif_grid@data$So), t+1) # begin with total susceptible population
total_R <- rep(sum(vax_motif_grid@data$Ro), t+1) # begin with total susceptible population

output_list <-list()

#define rho which depends upon the number of neighbors so we are scaling transmission by how many neighbors
rho <- 1 / (queen_adj$Num_neighbs)

out_A <- matrix(0, nrow = 256, ncol = 256)

# start for loop
for (i in 1:nrow(queen_adj)) {
  out_A[queen_adj$from_id[i], queen_adj$to_id[i]] <- rho[i]
}

#Check for NA in matrix
#View(out_A)
#isNA <- is.na(out_A)


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
  
  tmp <- list(ID = vax_motif_grid@data$ID, S = vax_motif_grid@data$So, I= vax_motif_grid@data$Io, R = vax_motif_grid@data$Ro)
  output_list[[i]] <- tmp
 
  
  #number of incident cases
  #cumulative_incidence <- 15359 - total_S #94%
  cumulative_incidence <- 12799 - total_S #95%
  #cumulative_incidence <- 5120 - total_S #98%
  #cumulative_incidence <- 2560 - total_S #99%
  
  incidence[i+1] <- total_S[i] - total_S[i+1]
  
  #plotvar <- vax_motif_grid@data$Io
  
  #define color scheme for plot
  #plotclr <- brewer.pal(9,"PuBu")
  
  
  # #second color code for the number of cases / number of susceptibles in each cell
  # colcode <- ifelse((plotvar <=   0), plotclr[1],
  #                   ifelse((plotvar >   0 & plotvar <=  0.5), plotclr[2],
  #                          ifelse((plotvar >   0.5 & plotvar <=   1), plotclr[3],
  #                                 ifelse((plotvar >   1 & plotvar <=   3), plotclr[4],
  #                                        ifelse((plotvar >   3 & plotvar <=  5), plotclr[5],
  #                                               ifelse((plotvar >  5 & plotvar <=  10), plotclr[6],
  #                                                      ifelse((plotvar >  10 & plotvar <=  15), plotclr[7],
  #                                                             ifelse((plotvar > 15 & plotvar <=20), plotclr[8],
  #                                                                    plotclr[9]))))))))
  # 
  # 
  # # #define plot margins
  # par(mar=c(2, 3, 3, 2))
  # # #for x, y4, z4, q4
  # png(filename = paste("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Clustering Motifs/Simulations/Figure 2 Day", i,".png"), width = 900, height = 650, units = "px", pointsize = 12, bg = "white")
  # plot(vax_motif_grid, col = colcode)
  # legend("left", legend = c("0", "0-0.5", "0.5-1", "1-3", "3-5", "5-10", "10-15", "15-20", ">20"),  col = brewer.pal(9,"PuBu")[1:9], lty = 1, lwd = 8, cex = 1)
  # title(paste("Motif with 85% Vaccination at Quadrant Level, 94% Overall Vaccination: Number of Infections, Day", i))
  # # #add legend for test case with homogeneous vaccination
  # # 
  # dev.off()


}

#get final # of cases
cumulative_incidence[366]

#plot map of those who got sick at the end vs. those who were susceptible at the beginning (Attack Rate Map)

# start by freezing initial conditions at t1
t1 <- output_list[[1]]
t1 <- as.data.frame(do.call(Map, c(f = rbind, t1)))
#transpose
library(data.table)
t1 <- transpose(t1)
t1 <- plyr::rename(t1, c("V1" = "ID", "V2" = "S1", "V3" = "I1", "V4" = "R1"))

# repeat by freezing final conditions at t365
t365 <- output_list[[365]]
t365 <- as.data.frame(do.call(Map, c(f = rbind, t365)))
#transpose
t365 <- transpose(t365)
t365 <- plyr::rename(t365, c("V1" = "ID", "V2" = "S365", "V3" = "I365", "V4" = "R365"))

start_end_map <- merge(t1, t365, by = "ID")
start_end_map$tot_cases <- as.numeric(start_end_map$S1) - as.numeric(start_end_map$S365)
start_end_map$inf_pct <- (as.numeric(start_end_map$S1) - as.numeric(start_end_map$S365))/as.numeric(start_end_map$S1)*100
start_end_map$surv_pct <- as.numeric(start_end_map$S365)/as.numeric(start_end_map$S1)*100
# 

#merge in this image of start and end to vax_motif_grid
vax_motif_grid_end <- merge(vax_motif_grid, start_end_map, by = "ID")


# now we plot it
plotvar <- vax_motif_grid_end@data$tot_cases
summary(vax_motif_grid_end@data$inf_pct)

#define color scheme for plot
plotclr <- rev(brewer.pal(9,"RdBu"))

#second color code for the number of cases / number of susceptibles in each cell
colcode <- ifelse((plotvar <=   0), plotclr[1],
                  ifelse((plotvar >   0 & plotvar <=   1), plotclr[2],
                         ifelse((plotvar >   1 & plotvar <=   5), plotclr[3],
                                ifelse((plotvar >  5 & plotvar <=   10), plotclr[4],
                                       ifelse((plotvar >   10 & plotvar <=  20), plotclr[5],
                                              ifelse((plotvar >  20 & plotvar <=  50), plotclr[6],
                                                     ifelse((plotvar >  50 & plotvar <=  100), plotclr[7],
                                                            ifelse((plotvar > 100 & plotvar <= 150), plotclr[8],
                                                                   plotclr[9]))))))))


# # PLOT SURVIVOR MAP
par(mar=c(5, 8, 3, 2))
plot(vax_motif_grid_end, col = colcode)
legend(3,16, legend = c("0", "0-1", "1-5", "5-10", "10-20", "20-50", "50-100", "100-150", ">150"),  col = rev(brewer.pal(9,"RdBu"))[1:9], lty = 1, lwd = 12, cex = 0.6)
title("Year-End Infection Rate Per Cell (95% Overall Vaccination Coverage)")


#plot final cumulative incidence curve
t2 <- 1:(t+1)
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(t2, incidence, type ="l", col = "blue", xlim = c(0,400), ylim= c(0,100), bty = "n", xlab = "time (days)", ylab = "# hosts", main = "Incident infections over simulation time: drop in quadrant 1 (94% Total Vax)")
