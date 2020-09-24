## Saving motifs (images) for each vaccination level

setwd("/Users/ninamasters/measles-spatial-model/")
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
#pull in necessary functions to make initial grids
source("/Users/ninamasters/measles-spatial-model/Clustering_Functions_94.R")

#for 94% vaccination coverage

  #read in violating motifs (i.e. motifs that exceed 1,000 people per cell)
  violate <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/motif_violate_94_threshold.csv")
  names <- seq(1:625)
  names <- names[!(names %in% violate$motif)]
  
  # for for the 625 motifs, run loop to spit out the image
  for(i in names){
    
    make_me_a_motif(combo_list$v1[[i]],combo_list$v2[[i]],combo_list$v3[[i]],combo_list$v4[[i]])
  }
    
###### Clear environment 
# for 95% vaccination coverage
  
  source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
  #pull in necessary functions to make initial grids
  source("/Users/ninamasters/measles-spatial-model/Clustering_Functions_95.R")
  
  #read in violating motifs (i.e. motifs that exceed 1,000 people per cell)
  violate <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/motif_violate_95_threshold.csv")
  names <- seq(1:625)
  names <- names[!(names %in% violate$motif)]
  
  # for for the 625 motifs, run loop to spit out the image
  for(i in names){
    
    make_me_a_motif_plot(combo_list$v1[[i]],combo_list$v2[[i]],combo_list$v3[[i]],combo_list$v4[[i]])
  }
  
  ###### Clear environment 
  # for 99% vaccination coverage
  
  source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
  #pull in necessary functions to make initial grids
  source("/Users/ninamasters/measles-spatial-model/Clustering_Functions_99.R")
  
  names <- seq(1:625)
  
  #violating motifs for 99%
  violate <- c(1,2,6,26,126)
  
  names <- names[!(names %in% violate)]
  
  # for for the 625 motifs, run loop to spit out the image
  for(i in names){
    
    make_me_a_motif(combo_list$v1[[i]],combo_list$v2[[i]],combo_list$v3[[i]],combo_list$v4[[i]])
  }
  
  