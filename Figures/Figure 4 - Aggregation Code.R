# Automated Aggregation code

# Aim 1
# Last Updated Feb 14, 2020 
# Nina Masters

##################################################################################
setwd("/Users/ninamasters/measles-spatial-model")

source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
#pull in necessary functions to make initial grids
source("/Users/ninamasters/measles-spatial-model/Clustering_Functions_95.R")
#pull in neighbor defining function
source("/Users/ninamasters/measles-spatial-model/Neighbor_Functions.R")
#pull in new aggregation function 
source("/Users/ninamasters/measles-spatial-model/aggregation_function.R")

##################################################################################

#aggregate_my_motif(x,y,z4,q4)
#aggregate_my_motif(x, y3, z, q4)
aggregate_my_motif(x,y1,z3,q4)

############################################ LEVEL 1 ####################################
plotvar <- vax_motif_grid$So/10 #So/1000 = proportion nonvax, So/1000*100 = So/10 = percent nonvax per cell

#define color scheme for plot
plotclr <- brewer.pal(9,"PuRd")

colcode <- ifelse((plotvar <=   0), plotclr[1],
                  ifelse((plotvar >   0 & plotvar <=   1), plotclr[2],
                         ifelse((plotvar >   1 & plotvar <=   2), plotclr[3],
                                ifelse((plotvar >   2 & plotvar <=   5), plotclr[4],
                                       ifelse((plotvar >   5 & plotvar <=  10), plotclr[5],
                                              ifelse((plotvar >  10 & plotvar <=  15), plotclr[6],
                                                     ifelse((plotvar >  15 & plotvar <=  20), plotclr[7],
                                                            ifelse((plotvar > 20 & plotvar <= 25), plotclr[8],
                                                                   plotclr[9]))))))))

#define plot margins
par(mar=c(3, 8, 5, 1))

# # calculate the Moran's I
w <- poly2nb(vax_motif_grid, queen = T, row.names = vax_motif_grid$ID)
class(w)

#use style = "W" for row-standardized values of Moran's I
wm <- nb2mat(w, style='W', zero.policy = TRUE)
ww <- nb2listw(w, style = 'W', zero.policy =TRUE)

m <- moran(vax_motif_grid$So/10, ww, length(ww$neighbours), S0=Szero(ww))

#define plot margins
par(mar=c(3, 8, 5, 1))

#for x, y4, z4, q4
plot(vax_motif_grid, col = colcode) 
#legend(-7,14, legend = c("0", "0-1%", "1-2%", "2-5%", "5-10%", "10-15%", "15-20%", "20-25%", ">25%"),  col = brewer.pal(9,"PuRd")[1:9], lty = 1, lwd = 10, cex = 0.5)
title(paste("Moran =", round(m$I,3),", Iso =", round(isolation_list$isolation_index_l1,2)),"Level 1 Aggregation (Truth)", line = 0.6)


################################################ NOW RUN THE MODEL FOR CI
  summary_list = list()
  
    for (j in 1:4){
    
      summary_list[[paste("quadrant",j)]] <- list()
      
        #now randomly drop in one case into the four different test locations: one in each quadrant
        #quadrant_infection vector has quad 1, 2, 3, 4
        quadrant_infection <- c(196, 205, 52, 61)
        names(quadrant_infection) <- c("quad1", "quad2", "quad3", "quad4")
        
        #add in one infected case in that cell (add for all four levels of aggregation)
        vax_motif_grid@data$Io[quadrant_infection[j]] <- 1
        vax_motif_grid@data$level_2_aggregate_Io[quadrant_infection[j]] <- 1
        vax_motif_grid@data$level_3_aggregate_Io[quadrant_infection[j]] <- 1
        vax_motif_grid@data$level_4_aggregate_Io[quadrant_infection[j]] <- 1
        
        #now remove one susceptible in that cell (i.e. the susceptible person got infected)
        #don't do this for 98% and 99%
        #can run in tandem for all levels of aggregation
        vax_motif_grid@data$So[quadrant_infection[j]] <- vax_motif_grid@data$So[quadrant_infection[j]]-1
        vax_motif_grid@data$level_2_aggregate_So[quadrant_infection[j]] <- vax_motif_grid@data$level_2_aggregate_So[quadrant_infection[j]] -1
        vax_motif_grid@data$level_3_aggregate_So[quadrant_infection[j]] <- vax_motif_grid@data$level_3_aggregate_So[quadrant_infection[j]] -1
        vax_motif_grid@data$level_4_aggregate_So[quadrant_infection[j]] <- vax_motif_grid@data$level_4_aggregate_So[quadrant_infection[j]] -1
        
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
              
              
          
            }
        
        #just save last image at the end of the simulation
        summary_list[[paste("quadrant",j)]][[paste("time",i)]]  <-  list(S1 = total_S[i], S2 = total_S_2[i], S3 = total_S_3[i], S4 = total_S_4[i])
        
    }


#now convert nested list output to dataframe to export as csv
#do this for summary data
  
sim_try <- unlist(summary_list, recursive = FALSE, use.names = TRUE)
sim_try <- as.data.frame(do.call(Map, c(f = rbind, sim_try)))

sim_try$CI1 <- 12799 - sim_try$S1
sim_try$CI2 <- 12799 - sim_try$S2
sim_try$CI3 <- 12799 - sim_try$S3
sim_try$CI4 <- 12799 - sim_try$S4

mean(sim_try$CI1/12799)
mean(sim_try$CI2/12799)
mean(sim_try$CI3/12799)
mean(sim_try$CI4/12799)

mean(sim_try$CI1)
mean(sim_try$CI2)
mean(sim_try$CI3)
mean(sim_try$CI4)
