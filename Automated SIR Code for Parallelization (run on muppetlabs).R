# Write R Code for deterministic SIR model in spatial environment : test to generate automated loop that can iterate through different motifs, and then four different test drop sites for each and output necessary data
# Aim 1
# Last Updated January 22, 2020
# Nina Masters

##################################################################################
#set working directory
setwd("/Users/ninamasters/measles-spatial-model")
#pull in necessary packages
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
#pull in necessary functions to make initial grids
source("/Users/ninamasters/measles-spatial-model/Clustering_Functions.R")
#pull in neighbor defining function
source("/Users/ninamasters/measles-spatial-model/Neighbor_Functions.R")

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

names <- seq(1:2)
Simulation_Data = foreach(p = 1:length(names)) %dopar% {
  
  output_list = list()
  summary_list = list()
  output_list[[paste("motif",p)]] <- list()
  summary_list[[paste("motif",p)]] <- list()
  
  source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
  #pull in necessary functions to make initial grids
  source("/Users/ninamasters/measles-spatial-model/Clustering_Functions.R")
  #pull in neighbor defining function
  source("/Users/ninamasters/measles-spatial-model/Neighbor_Functions.R")
  
  # j's to  iterate through dropping test case in the four different locations

  
  for (j in 1:4){
  output_list[[paste("motif",p)]][[paste("quadrant", j)]] <- list()
  summary_list[[paste("motif",p)]][[paste("quadrant", j)]] <- list()
  
      
      #now randomly drop in one case into the four different test locations: one in each quadrant
      #quadrant_infection vector has quad 1, 2, 3, 4
      quadrant_infection <- c(196, 205, 52, 61)
      names(quadrant_infection) <- c("quad1", "quad2", "quad3", "quad4")
      
    
      name <- paste0("motif",names[p],": ", combo_list$v1[[p]],combo_list$v2[[p]],combo_list$v3[[p]],combo_list$v4[[p]])
      temp <- motif_df(combo_list$v1[[p]],combo_list$v2[[p]],combo_list$v3[[p]],combo_list$v4[[p]])
      assign(name, temp)
      
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
              
              #desired output
              tmp = list(motif = p, run = j, ID = vax_motif_grid@data$ID, S = vax_motif_grid@data$So, Prevalence = vax_motif_grid@data$Io, R = vax_motif_grid@data$Ro)
              output_list[[paste("motif",p)]][[paste("quadrant",j)]][[paste("time",i)]] <- tmp
              
              ###################################################
              #need to make summary output not indexed by cell
              summary_list[[paste("motif",p)]][[paste("quadrant",j)]][[paste("time",i)]]  <-  list(motif = p, run = j, Prevalence = total_I[i], Recovered = total_R[i], Susceptibles = total_S[i], Moran_S = m_S$I, Moran_Infected = m_I$I)
              
            }
      
  }
 
  write.csv(output_list, paste("individual_output_motif",p, ".csv"))
  final_output = summary_list
  
  
}

#individual output
library(rlist)

# sim_try <- unlist(Simulation_Data, recursive = FALSE)[ c(TRUE,FALSE) ]
# sim_try <- unlist(sim_try, recursive = FALSE, use.names = TRUE)
# output_df <- as.data.frame(sim_try) %>%
# select(c("individual.motif.1.quadrant.1.time.1.ID", contains(".S"), contains(".Prevalence"), contains(".R")), -contains(".run"))

#now do for summary output
sim_try <- unlist(Simulation_Data, recursive = FALSE, use.names = TRUE)
sim_try <- unlist(sim_try, recursive = FALSE, use.names = TRUE)
sim_try <- unlist(sim_try, recursive = FALSE, use.names = TRUE)

sim_try <- as.data.frame(do.call(Map, c(f = rbind, sim_try)))
write.csv(sim_try, "simulation_summary_data.csv")



