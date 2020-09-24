# Make motif function that includes aggregation at four subsequent scales
# NOTE  - THIS IS CURRENTLY SET FOR MOTIFS of 98% VAX
# February 14, 2020

setwd("/Users/ninamasters/measles-spatial-model")

library(ggplot2)
library(maps)
library(mapdata)
library(maptools)
library(gridExtra)
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


# set up adjacencies / boundaries using queen boundaries
queen_boundaries_grid <- poly2nb(gridspdf, queen = T)
#rook_boundaries_grid <- poly2nb(gridspdf, queen = F)

############################### Pull in unique codes on aggregation scales ###############################
############################### here we have 1,2,3,4 for quadrants, etc. #################################
level_unique <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Clustering Motifs/level_unique.csv")
motifsu <- merge(gridspdf,level_unique, by="ID")

#calculate isolation at the non-aggregate scale
isolation <- function(vax_motif_grid){
  a <- vax_motif_grid$So[(vax_motif_grid$So+vax_motif_grid$Ro) > 0]
  n <- (vax_motif_grid$So+vax_motif_grid$Ro)[(vax_motif_grid$So+vax_motif_grid$Ro) > 0]
  total_a <- sum(a)
  d <<- sum((a/total_a)*(a/n))
  
  return(d)
}

#calculate isolation at agg level 2
isolation_2 <- function(vax_motif_grid){
  a <- vax_motif_grid$level_2_aggregate_So[(vax_motif_grid$level_2_aggregate_So+vax_motif_grid$level_2_aggregate_Ro) > 0]
  n <- (vax_motif_grid$level_2_aggregate_So+vax_motif_grid$level_2_aggregate_Ro)[(vax_motif_grid$level_2_aggregate_So+vax_motif_grid$level_2_aggregate_Ro) > 0]
  total_a <- sum(a)
  d <<- sum((a/total_a)*(a/n))
  
  return(d)
}

#calculate isolation at agg level 3
isolation_3 <- function(vax_motif_grid){
  a <- vax_motif_grid$level_3_aggregate_So[(vax_motif_grid$level_3_aggregate_So+vax_motif_grid$level_3_aggregate_Ro) > 0]
  n <- (vax_motif_grid$level_3_aggregate_So+vax_motif_grid$level_3_aggregate_Ro)[(vax_motif_grid$level_3_aggregate_So+vax_motif_grid$level_3_aggregate_Ro) > 0]
  total_a <- sum(a)
  d <<- sum((a/total_a)*(a/n))
  
  return(d)
}

#calculate isolation at agg level 4
isolation_4 <- function(vax_motif_grid){
  a <- vax_motif_grid$level_4_aggregate_So[(vax_motif_grid$level_4_aggregate_So+vax_motif_grid$level_4_aggregate_Ro) > 0]
  n <- (vax_motif_grid$level_4_aggregate_So+vax_motif_grid$level_4_aggregate_Ro)[(vax_motif_grid$level_4_aggregate_So+vax_motif_grid$level_4_aggregate_Ro) > 0]
  total_a <- sum(a)
  d <<- sum((a/total_a)*(a/n))
  
  return(d)
}

aggregate_my_motif <- function(x,y,z,q){
  
  # load in initial conditions and generate a spatial polygon data frame: this is the 'true' data:
  #create distribution of the non-vaccinators in each quadrant
  quadrant <- c(1,2,3,4)
  probability_quadrant <- as.data.frame(cbind(quadrant,x))
  
  #create distribution of the non-vaccinators in each quadrant at the neighborhood level
  neighborhoods <- c(1:16)
  neighb_quad <- c(1,1,2,2,1,1,2,2,3,3,4,4,3,3,4,4)
  y_neighb <- rep(y, 2)
  
  probability_neighb <- as.data.frame(cbind(neighborhoods, neighb_quad))
  probability_neighb <- as.data.frame(cbind(probability_neighb, y_neighb))
  names(probability_neighb) <- c("neighb", "quadrant", "y")
  
  prob_2 <- dplyr::left_join(probability_quadrant, probability_neighb, by="quadrant")
  
  #create distribution of non-vaccinators in each quadrant at the block level (n = 64)
  blocks <- c(1:64)
  block_neighb <- c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,13,13,14,14,15,15,16,16)
  y_block <- rep(z, 4)
  probability_block <- as.data.frame(cbind(blocks, block_neighb))
  probability_block <- as.data.frame(cbind(probability_block, y_block))
  names(probability_block) <- c("block", "neighb", "z")
  
  prob_3 <- dplyr::left_join(prob_2, probability_block, by="neighb")
  
  #create distribution of non-vaccinators in each quadrant at the individual cell level (n = 256)
  cells <- c(1:256)
  cell_block <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,21,21,22,22,23,23,24,24,17,17,18,18,19,19,20,20,21,21,22,22,23,23,24,24,25,25,26,26,27,27,28,28,29,29,30,30,31,31,32,32,25,25,26,26,27,27,28,28,29,29,30,30,31,31,32,32,33,33,34,34,35,35,36,36,37,37,38,38,39,39,40,40,33,33,34,34,35,35,36,36,37,37,38,38,39,39,40,40,41,41,42,42,43,43,44,44,45,45,46,46,47,47,48,48,41,41,42,42,43,43,44,44,45,45,46,46,47,47,48,48,49,49,50,50,51,51,52,52,53,53,54,54,55,55,56,56,49,49,50,50,51,51,52,52,53,53,54,54,55,55,56,56,57,57,58,58,59,59,60,60,61,61,62,62,63,63,64,64,57,57,58,58,59,59,60,60,61,61,62,62,63,63,64,64)
  y_cell <- rep(q, 8)
  
  probability_cell <- as.data.frame(cbind(cells, cell_block))
  probability_cell <- as.data.frame(cbind(probability_cell, y_cell))
  names(probability_cell) <- c("cell", "block", "q")
  
  prob_4 <- dplyr::left_join(prob_3, probability_cell, by="block")
  
  #now calculate cumulative probability
  prob_4 <- dplyr::mutate(prob_4,cumulative_prob = (x*y*z*q))
  
  #sort data by cell
  prob_4 <- prob_4[order(prob_4$cell),]
  
  #### multinomial draw according to distributions in prob_4
  #add the number of nonvaccinators per cell to prob_4 data frame
  prob_4 <- dplyr::mutate(prob_4, nonvax_percell =rmultinom(1, 5120, prob_4$cumulative_prob))
  
  #create percent of non-vaccinators (/1000, multiplied by 100 to get percent)
  prob_4 <- dplyr::mutate(prob_4, percent_nonvax =nonvax_percell/1000*100, vax_percell = 1000-nonvax_percell, infected_percell = 0) #final prob_4 dataset
  
  vax_motif <- dplyr::select(prob_4, cell, percent_nonvax, nonvax_percell, vax_percell, infected_percell)
  vax_motif <- dplyr::mutate(vax_motif, ID = paste("g", cell, sep=""))
  names(vax_motif) <- c("cell", "percent_nonvax", "So", "Ro", "Io", "ID")
  
  # now merge with gridspdf, the spatial polygon data frame containing our grid information
  vax_motif_grid <<- merge(gridspdf,vax_motif, by="ID")
  
  # #clear out additional "cell" column from level_unique file
  drops <- c("cell")
  level_unique <- level_unique[ , !(names(level_unique) %in% drops)]
  
  vax_motif_grid <<- merge(vax_motif_grid,level_unique, by="ID")
  
  #now we want to aggregate this up a level...
  #first let's take 'true' base-level data and aggregate / average up to block group level
  agg_block <- aggregate(list(vax_motif_grid@data$So), by = list(vax_motif_grid@data$level2u), FUN = "sum", drop = TRUE) 
  colnames(agg_block) <- c("level2u", "So_agg_block")
  
  vax_motif_grid <<- merge(vax_motif_grid, agg_block, by = "level2u")
  
  #now let's take 'true' base-level data and aggregate / average up to neighborhood level
  agg_neighb <- aggregate(list(vax_motif_grid@data$So), by = list(vax_motif_grid@data$level3u), FUN = "sum", drop = TRUE) 
  colnames(agg_neighb) <- c("level3u", "So_agg_neighb")
  
  vax_motif_grid <<- merge(vax_motif_grid, agg_neighb, by = "level3u")
  
  #now let's take 'true' base-level data and aggregate / average up to quadrant level
  agg_quad <- aggregate(list(vax_motif_grid@data$So), by = list(vax_motif_grid@data$level4u), FUN = "sum", drop = TRUE) 
  colnames(agg_quad) <- c("level4u", "So_agg_quad")
  
  vax_motif_grid <<- merge(vax_motif_grid, agg_quad, by = "level4u")
  
  #now divide total aggregated So_agg_block by 4
  vax_motif_grid@data$level_2_aggregate_So <<- vax_motif_grid@data$So_agg_block/4
  
  #now divide total aggregated So_agg_neighb by 16
  vax_motif_grid@data$level_3_aggregate_So <<- vax_motif_grid@data$So_agg_neighb/16
  
  #now divide total aggregated So_agg_quad by 64
  vax_motif_grid@data$level_4_aggregate_So <<- vax_motif_grid@data$So_agg_quad/64
  
  #now create similar level 2 agg Ro
  vax_motif_grid@data$level_2_aggregate_Ro <<- 1000-vax_motif_grid@data$level_2_aggregate_So
  
  #now create similar level 3 agg Ro
  vax_motif_grid@data$level_3_aggregate_Ro <<- 1000-vax_motif_grid@data$level_3_aggregate_So
  
  #now create similar level 4 agg Ro
  vax_motif_grid@data$level_4_aggregate_Ro <<- 1000-vax_motif_grid@data$level_4_aggregate_So
  
  #now create similar level 2 agg Io
  vax_motif_grid@data$level_2_aggregate_Io <<- vax_motif_grid@data$Io
  
  #now create similar level 3 agg Io
  vax_motif_grid@data$level_3_aggregate_Io <<- vax_motif_grid@data$Io
  
  #now create similar level 4 agg Io
  vax_motif_grid@data$level_4_aggregate_Io <<- vax_motif_grid@data$Io
  
  isolation_list <<- list(level_1 = x[1], level_2 = y[1], level_3 = z[1], level_4 = q[1], isolation_index_l1 = isolation(vax_motif_grid), isolation_index_l2 = isolation_2(vax_motif_grid), isolation_index_l3 = isolation_3(vax_motif_grid), isolation_index_l4 = isolation_4(vax_motif_grid))
  
  ###################### now add the Moran's I stuff
  # calculate the Moran's I
  w <- poly2nb(vax_motif_grid, queen = T, row.names = vax_motif_grid$ID)
  class(w)
  
  #use style = "W" for row-standardized values of Moran's I
  wm <- nb2mat(w, style='W', zero.policy = TRUE)
  ww <- nb2listw(w, style = 'W', zero.policy =TRUE)
  
  m1 <- moran(vax_motif_grid$So/10, ww, length(ww$neighbours), S0=Szero(ww))
  m2 <- moran(vax_motif_grid$level_2_aggregate_So/10, ww, length(ww$neighbours), S0=Szero(ww))
  m3 <- moran(vax_motif_grid$level_3_aggregate_So/10, ww, length(ww$neighbours), S0=Szero(ww))
  m4 <- moran(vax_motif_grid$level_4_aggregate_So/10, ww, length(ww$neighbours), S0=Szero(ww))
  
  # save moran list as well
  moran_list <<- list(level_1 = x[1], level_2 = y[1], level_3 = z[1], level_4 = q[1], Moran_l1 = m1$I, Moran_l2 = m2$I, Moran_l3 = m3$I, Moran_l4 = m4$I)

}



#can create different "x" motifs here for outermost layer of clustering
x <- c(0.85, 0.05, 0.05, 0.05)
x1 <- c(0.7, 0.1, 0.1, 0.1)
x2 <- c(0.58, 0.14, 0.14, 0.14)
x3 <- c(0.4, 0.2, 0.2, 0.2)
x4 <- c(0.25, 0.25, 0.25, 0.25)

#can create different "y" motifs for neighborhood level
y <- c(0.85, 0.05, 0.85, 0.05, 0.05, 0.05, 0.05, 0.05)
y1 <- c(0.7, 0.1, 0.7, 0.1, 0.1, 0.1, 0.1, 0.1)
y2 <- c(0.58, 0.14, 0.58, 0.14, 0.14, 0.14, 0.14, 0.14)
y3 <- c(0.4, 0.2, 0.4, 0.2, 0.2, 0.2, 0.2, 0.2)
y4 <- c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25) #homogeneous quadrants at neighborhood level

#can create different "y" motifs for block level
z <- c(0.85, 0.05, 0.85, 0.05, 0.85, 0.05, 0.85, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
z1 <- c(0.7, 0.1, 0.7, 0.1, 0.7, 0.1, 0.7, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
z2 <- c(0.58, 0.14, 0.58, 0.14, 0.58, 0.14, 0.58, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14)
z3 <- c(0.4, 0.2, 0.4, 0.2, 0.4, 0.2, 0.4, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
z4 <- c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25) #homogeneous quadrants at neighborhood level

#can create different "y" motifs for individual level
q <- c(0.85, 0.05, 0.85, 0.05, 0.85, 0.05, 0.85, 0.05, 0.85, 0.05, 0.85, 0.05, 0.85, 0.05, 0.85, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
q1 <- c(0.7, 0.1, 0.7, 0.1, 0.7, 0.1, 0.7, 0.1, 0.7, 0.1, 0.7, 0.1, 0.7, 0.1, 0.7, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
q2 <- c(0.58, 0.14, 0.58, 0.14, 0.58, 0.14, 0.58, 0.14, 0.58, 0.14, 0.58, 0.14, 0.58, 0.14, 0.58, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14)
q3 <- c(0.4, 0.2, 0.4, 0.2, 0.4, 0.2, 0.4, 0.2, 0.4, 0.2, 0.4, 0.2, 0.4, 0.2, 0.4, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
q4 <- c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25) #homogeneous quadrants at neighborhood level

#combination of grids:
x_grid <- list(x, x1, x2, x3, x4)
names(x_grid) <- c("x", "x1", "x2", "x3", "x4")
y_grid <- list(y, y1, y2, y3, y4)
names(y_grid) <- c("y", "y1", "y2", "y3", "y4")
z_grid <- list(z, z1, z2, z3, z4)
names(z_grid) <- c("z", "z1", "z2", "z3", "z4")
q_grid <- list(q, q1, q2, q3, q4)
names(q_grid) <- c("q", "q1", "q2", "q3", "q4")

combo_list <- as.list(expand.grid(x_grid, y_grid, z_grid, q_grid))
names(combo_list) <- c("v1", "v2", "v3", "v4")
