# Automated Aggregation code

# Aim 1
# Last Updated Feb 14, 2020 
# Nina Masters

##################################################################################
source("~/packages_clustering.R")
#pull in necessary functions to make initial grids
source("~/Clustering_Functions_99.R")
#pull in neighbor defining function
source("~/Neighbor_Functions.R")
#pull in new aggregation function 
source("~/aggregation_function.R")

##################################################################################

# outer for loop: interate through all the different motif combinations to generate them all 
#then for each motif, we want to iterate through the four test sites and save output. 

#do this 20 times and save the output
#everything inside the p loop and turn into a function... will make life easier

library(foreach)
library(doParallel)

#### Setup parallel backend ####
cores = detectCores()

mycluster = makeCluster(cores[1]-2) 
# leaving 2 cores out so I can run other things and not overwhelm the 
# computer, you could leave just one out and use more cores for this

registerDoParallel(mycluster)

# p's to create 625 motifs

names <- seq(1:625)
#read in motifs that violate condition of 1,000 per cell
violate <- read.csv("~/motif_violate_94_threshold.csv") # for 94%
#violate <- read.csv("~/motif_violate_95_threshold.csv") # for 95%
#violate <- read.csv("~/motif_violate_98_threshold.csv") # for 98%
#violate <- c(1,2,6,26,126) # for 99%

names <- names[!(names %in% violate$motif)] #for 94%, 95%, 98%
#names <- names[!(names %in% violate)] # for 99%

Sim_Data = foreach(p = unique(names)) %dopar% {

#now simulate the model..... (would swap out different aggregate So variables)
  
  summary_list = list()
  
  summary_list[[paste("motif",p)]] <- list()
  
  source("~/packages_clustering.R")
  #pull in necessary functions to make initial grids
  source("~/Clustering_Functions_95.R")
  #pull in neighbor defining function
  source("~/Neighbor_Functions.R")
  #pull in new aggregation function 
  source("~/aggregation_function_95.R")
  
  # j's to  iterate through dropping test case in the four different locations
  
    for (j in 1:4){
    
        summary_list[[paste("motif",p)]][[paste("quadrant", j)]] <- list()
        
        
        #now randomly drop in one case into the four different test locations: one in each quadrant
        #quadrant_infection vector has quad 1, 2, 3, 4
        quadrant_infection <- c(196, 205, 52, 61)
        names(quadrant_infection) <- c("quad1", "quad2", "quad3", "quad4")
        
        name <- paste0("motif",names[p],": ", combo_list$v1[[p]],combo_list$v2[[p]],combo_list$v3[[p]],combo_list$v4[[p]])
        temp <- motif_df(combo_list$v1[[p]],combo_list$v2[[p]],combo_list$v3[[p]],combo_list$v4[[p]])
        assign(name, temp)
                
        
        #add in one infected case in that cell (add for all four levels of aggregation)
        vax_motif_grid@data$Io[val_quad] <- 1
        vax_motif_grid@data$level_2_aggregate_Io[val_quad] <- 1
        vax_motif_grid@data$level_3_aggregate_Io[val_quad] <- 1
        vax_motif_grid@data$level_4_aggregate_Io[val_quad] <- 1
        
        #now remove one susceptible in that cell (i.e. the susceptible person got infected)
        #don't do this for 98% and 99%
        #can run in tandem for all levels of aggregation
        vax_motif_grid@data$So[val_quad] <- vax_motif_grid@data$So[val_quad]-1
        vax_motif_grid@data$level_2_aggregate_So[val_quad] <- vax_motif_grid@data$level_2_aggregate_So[val_quad] -1
        vax_motif_grid@data$level_3_aggregate_So[val_quad] <- vax_motif_grid@data$level_3_aggregate_So[val_quad] -1
        vax_motif_grid@data$level_4_aggregate_So[val_quad] <- vax_motif_grid@data$level_4_aggregate_So[val_quad] -1
        
        # Define Model Parameters 
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
        total_I_2 <- rep(0, t+1)
        total_I_3 <- rep(0, t+1)
        total_I_4 <- rep(0, t+1)
        
        total_I[1] <- sum(vax_motif_grid@data$Io)
        total_I_2[1] <- sum(vax_motif_grid@data$level_2_aggregate_Io)
        total_I_3[1] <- sum(vax_motif_grid@data$level_3_aggregate_Io)
        total_I_4[1] <- sum(vax_motif_grid@data$level_4_aggregate_Io)
        
        total_S <- rep(sum(vax_motif_grid@data$So), t+1) # begin with total susceptible population
        total_S_2 <- rep(sum(vax_motif_grid@data$level_2_aggregate_So), t+1) # begin with total susceptible population
        total_S_3 <- rep(sum(vax_motif_grid@data$level_3_aggregate_So), t+1) # begin with total susceptible population
        total_S_4 <- rep(sum(vax_motif_grid@data$level_4_aggregate_So), t+1) # begin with total susceptible population
        
        total_R <- rep(sum(vax_motif_grid@data$Ro), t+1) # begin with total susceptible population
        total_R_2 <- rep(sum(vax_motif_grid@data$level_2_aggregate_Ro), t+1) # begin with total susceptible population
        total_R_3 <- rep(sum(vax_motif_grid@data$level_3_aggregate_Ro), t+1) # begin with total susceptible population
        total_R_4 <- rep(sum(vax_motif_grid@data$level_4_aggregate_Ro), t+1) # begin with total susceptible population
        
        #define rho which depends upon the number of neighbors so we are scaling transmission by how many neighbors
        rho <- 1 / (queen_adj$Num_neighbs)

        out_A <- matrix(0, nrow = 256, ncol = 256)

                  # start for loop
                  for (i in 1:nrow(queen_adj)) {
                    out_A[queen_adj$from_id[i], queen_adj$to_id[i]] <- rho[i]
                  }

            for (i in 1:t) {
              
              ## Sample potential imports for all four scales
              outside_exposure   <- (out_A %*% as.numeric(vax_motif_grid@data$Io)) 
              outside_exposure_2 <- (out_A %*% as.numeric(vax_motif_grid@data$level_2_aggregate_Io)) 
              outside_exposure_3 <- (out_A %*% as.numeric(vax_motif_grid@data$level_3_aggregate_Io)) 
              outside_exposure_4 <- (out_A %*% as.numeric(vax_motif_grid@data$level_4_aggregate_Io)) 
              
              
              ## Get household outside exposure
              cell_foi   <- (between_beta * outside_exposure /   (1000*num_neighbors$Num_neighbs)) + (within_beta * as.numeric(vax_motif_grid@data$Io)/1000) 
              cell_foi_2 <- (between_beta * outside_exposure_2 / (1000*num_neighbors$Num_neighbs)) + (within_beta * as.numeric(vax_motif_grid@data$level_2_aggregate_Io)/1000) 
              cell_foi_3 <- (between_beta * outside_exposure_3 / (1000*num_neighbors$Num_neighbs)) + (within_beta * as.numeric(vax_motif_grid@data$level_3_aggregate_Io)/1000) 
              cell_foi_4 <- (between_beta * outside_exposure_4 / (1000*num_neighbors$Num_neighbs)) + (within_beta * as.numeric(vax_motif_grid@data$level_4_aggregate_Io)/1000) 
              
              
              ## New infections: level 1
              vax_motif_grid@data$So <- vax_motif_grid@data$So - cell_foi*vax_motif_grid@data$So
              vax_motif_grid@data$Io <- vax_motif_grid@data$Io - gamma*vax_motif_grid@data$Io + cell_foi*vax_motif_grid@data$So 
              vax_motif_grid@data$Ro <- vax_motif_grid@data$Ro + gamma*vax_motif_grid@data$Io
              
              # new infections: level 2
              vax_motif_grid@data$level_2_aggregate_So <- vax_motif_grid@data$level_2_aggregate_So - cell_foi_2*vax_motif_grid@data$level_2_aggregate_So
              vax_motif_grid@data$level_2_aggregate_Io <- vax_motif_grid@data$level_2_aggregate_Io - gamma*vax_motif_grid@data$level_2_aggregate_Io + cell_foi_2*vax_motif_grid@data$level_2_aggregate_So 
              vax_motif_grid@data$level_2_aggregate_Ro <- vax_motif_grid@data$level_2_aggregate_Ro + gamma*vax_motif_grid@data$level_2_aggregate_Io
              
              # new infections: level 3
              vax_motif_grid@data$level_3_aggregate_So <- vax_motif_grid@data$level_3_aggregate_So - cell_foi_3*vax_motif_grid@data$level_3_aggregate_So
              vax_motif_grid@data$level_3_aggregate_Io <- vax_motif_grid@data$level_3_aggregate_Io - gamma*vax_motif_grid@data$level_3_aggregate_Io + cell_foi_3*vax_motif_grid@data$level_3_aggregate_So 
              vax_motif_grid@data$level_3_aggregate_Ro <- vax_motif_grid@data$level_3_aggregate_Ro + gamma*vax_motif_grid@data$level_3_aggregate_Io
              
              # new infections: level 4
              vax_motif_grid@data$level_4_aggregate_So <- vax_motif_grid@data$level_4_aggregate_So - cell_foi_4*vax_motif_grid@data$level_4_aggregate_So
              vax_motif_grid@data$level_4_aggregate_Io <- vax_motif_grid@data$level_4_aggregate_Io - gamma*vax_motif_grid@data$level_4_aggregate_Io + cell_foi_4*vax_motif_grid@data$level_4_aggregate_So 
              vax_motif_grid@data$level_4_aggregate_Ro <- vax_motif_grid@data$level_4_aggregate_Ro + gamma*vax_motif_grid@data$level_4_aggregate_Io
              
              #now get totals for level 1
              total_S[i+1] <- sum(vax_motif_grid@data$So)  
              total_I[i+1] <- sum(vax_motif_grid@data$Io)
              total_R[i+1] <- sum(vax_motif_grid@data$Ro)
              
              #level 2
              total_S_2[i+1] <- sum(vax_motif_grid@data$level_2_aggregate_So)  
              total_I_2[i+1] <- sum(vax_motif_grid@data$level_2_aggregate_Io)
              total_R_2[i+1] <- sum(vax_motif_grid@data$level_2_aggregate_Ro)
              
              #level 3
              total_S_3[i+1] <- sum(vax_motif_grid@data$level_3_aggregate_So)  
              total_I_3[i+1] <- sum(vax_motif_grid@data$level_3_aggregate_Io)
              total_R_3[i+1] <- sum(vax_motif_grid@data$level_3_aggregate_Ro)
              
              #level 4
              total_S_4[i+1] <- sum(vax_motif_grid@data$level_4_aggregate_So)  
              total_I_4[i+1] <- sum(vax_motif_grid@data$level_4_aggregate_Io)
              total_R_4[i+1] <- sum(vax_motif_grid@data$level_4_aggregate_Ro)
              
              
              tmp <- list(ID = vax_motif_grid@data$ID, S = vax_motif_grid@data$So, I= vax_motif_grid@data$Io, R = vax_motif_grid@data$Ro, 
                          S2 = vax_motif_grid@data$level_2_aggregate_So, I2 = vax_motif_grid@data$level_2_aggregate_Io, R2 = vax_motif_grid@data$level_2_aggregate_Ro,
                          S3 = vax_motif_grid@data$level_3_aggregate_So, I3 = vax_motif_grid@data$level_3_aggregate_Io, R3 = vax_motif_grid@data$level_3_aggregate_Ro,
                          S4 = vax_motif_grid@data$level_4_aggregate_So, I4 = vax_motif_grid@data$level_4_aggregate_Io, R4 = vax_motif_grid@data$level_4_aggregate_Ro)
              
              summary_list[[paste("motif",p)]][[paste("quadrant",j)]][[paste("time",i)]]  <-  tmp
  
            }
    }

        final_output = summary_list
        
}        
