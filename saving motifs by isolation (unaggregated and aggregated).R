
setwd("/Users/ninamasters/measles-spatial-model/")
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
#pull in necessary functions to make initial grids
source("/Users/ninamasters/measles-spatial-model/Clustering_Functions_94.R")
source("/Users/ninamasters/measles-spatial-model/aggregation_function_94.R")

########### Saving Motifs in the for loop: save the isolation of each motif!

#############################################################################################################
########################################### 94% coverage ####################################################
#############################################################################################################

isolation_list_motifs = list()

#run this outer loop 30 times to get replicates of Moran's I
for (p in 1:30){
  isolation_list_motifs[[paste("run", p)]] <- list()
  
  #read in violating motifs (i.e. motifs that exceed 1,000 people per cell)
  violate <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/motif_violate_94_threshold.csv")
  names <- seq(1:625)
  names <- names[!(names %in% violate$motif)]
  
  # for for the 625 motifs, run loop to spit out moran's I and clustering at each level
  for(i in unique(names)){
    isolation_list_motifs[[paste("run", p)]][[paste("motif",i)]] <- list()
    
    make_me_an_isolation(combo_list$v1[[i]],combo_list$v2[[i]],combo_list$v3[[i]],combo_list$v4[[i]])
    
    isolation_list_motifs[[paste("run", p)]][[paste("motif", i)]] <- isolation_list
    
  }
}

#now turn list into data frame for 94% coverage
isolation_list <- unlist(isolation_list_motifs, recursive = FALSE, use.names = TRUE)
isolation_list_94 <- as.data.frame(do.call(Map, c(f = rbind, isolation_list)))

#write to csv
write.csv(isolation_list_94, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_motifs_30_94.csv")

########### Saving Motifs in the for loop: save the isolation of each motif!

#############################################################################################################
########################################### 95% coverage ####################################################
#############################################################################################################

setwd("/Users/ninamasters/measles-spatial-model/")
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
#pull in necessary functions to make initial grids
source("/Users/ninamasters/measles-spatial-model/Clustering_Functions_95.R")
source("/Users/ninamasters/measles-spatial-model/aggregation_function_95.R")


isolation_list_motifs = list()

#run this outer loop 30 times to get replicates of Moran's I
for (p in 1:30){
  isolation_list_motifs[[paste("run", p)]] <- list()
 
  #read in violating motifs (i.e. motifs that exceed 1,000 people per cell)
  violate <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/motif_violate_95_threshold.csv")
  names <- seq(1:625)
  names <- names[!(names %in% violate$motif)]
  
  # for for the 625 motifs, run loop to spit out moran's I and clustering at each level
  for(i in unique(names)){
    isolation_list_motifs[[paste("run", p)]][[paste("motif",i)]] <- list()
    
    make_me_an_isolation(combo_list$v1[[i]],combo_list$v2[[i]],combo_list$v3[[i]],combo_list$v4[[i]])
    
    isolation_list_motifs[[paste("run", p)]][[paste("motif", i)]] <- isolation_list
    
  }
}

#now turn list into data frame for 95% coverage
isolation_list <- unlist(isolation_list_motifs, recursive = FALSE, use.names = TRUE)
isolation_list_95 <- as.data.frame(do.call(Map, c(f = rbind, isolation_list)))

#write to csv
write.csv(isolation_list_95, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_motifs_30_95.csv")


#############################################################################################################
########################################### 99% coverage ####################################################
#############################################################################################################

setwd("/Users/ninamasters/measles-spatial-model/")
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
#pull in necessary functions to make initial grids
source("/Users/ninamasters/measles-spatial-model/Clustering_Functions_99.R")
source("/Users/ninamasters/measles-spatial-model/aggregation_function_99.R")


isolation_list_motifs = list()

#run this outer loop 30 times to get replicates of Moran's I
for (p in 1:30){
  isolation_list_motifs[[paste("run", p)]] <- list()
  
  names <- seq(1:625)

  #violating motifs for 99%
  violate <- c(1,2,6,26,126)
  
  names <- names[!(names %in% violate)]
  
  # for for the 625 motifs, run loop to spit out moran's I and clustering at each level
  for(i in unique(names)){
    isolation_list_motifs[[paste("run", p)]][[paste("motif",i)]] <- list()
    
    make_me_an_isolation(combo_list$v1[[i]],combo_list$v2[[i]],combo_list$v3[[i]],combo_list$v4[[i]])
    
    isolation_list_motifs[[paste("run", p)]][[paste("motif", i)]] <- isolation_list
    
  }
}

#now turn list into data frame for 99% coverage
isolation_list <- unlist(isolation_list_motifs, recursive = FALSE, use.names = TRUE)
isolation_list_99 <- as.data.frame(do.call(Map, c(f = rbind, isolation_list)))

#write to csv
write.csv(isolation_list_99, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_motifs_30_99.csv")


############################################################################################################
############################################################################################################
####################@@@@@@@## Now generate summary of isolation indices ####################################
############################################################################################################
############################################################################################################

#############################################################################################################
########################################### 94% coverage ####################################################
#############################################################################################################

iso_94 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_motifs_30_94.csv")

library(stringr)

#clean up run and motif columns (remove text)
iso_94$run <- as.numeric(substring(iso_94$run, regexpr("run ", iso_94$run) + 4))
iso_94$motif <- as.numeric(substring(iso_94$motif, regexpr("motif ", iso_94$motif) + 6))

motif_iso_summary <- iso_94 %>%
  group_by(motif) %>%
  summarize(mean_iso = mean(isolation_index), var_iso = var(isolation_index))

#now drop unnecessary columns from iso_94 file because we want the level information in there
drops <- c("run", "isolation_index")
iso_94 <- iso_94[ , !(names(iso_94) %in% drops)]
#now select only distinct rows here
iso_94<-iso_94[!duplicated(iso_94$motif), ]

iso_94_summary <- merge(iso_94, motif_iso_summary, by = "motif")

#write to csv file w summary info (mean moran, sd, variance)
write.csv(iso_94_summary, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_94.csv")


#############################################################################################################
########################################### 95% coverage ####################################################
#############################################################################################################

iso_95 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_motifs_30_95.csv")

library(stringr)

#clean up run and motif columns (remove text)
iso_95$run <- as.numeric(substring(iso_95$run, regexpr("run ", iso_95$run) + 4))
iso_95$motif <- as.numeric(substring(iso_95$motif, regexpr("motif ", iso_95$motif) + 6))

motif_iso_summary <- iso_95 %>%
  group_by(motif) %>%
  summarize(mean_iso = mean(isolation_index), var_iso = var(isolation_index))

#now drop unnecessary columns from iso_95 file because we want the level information in there
drops <- c("run", "isolation_index")
iso_95 <- iso_95[ , !(names(iso_95) %in% drops)]
#now select only distinct rows here
iso_95<-iso_95[!duplicated(iso_95$motif), ]

iso_95_summary <- merge(iso_95, motif_iso_summary, by = "motif")

#write to csv file w summary info (mean moran, sd, variance)
write.csv(iso_95_summary, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_95.csv")


#############################################################################################################
########################################### 99% coverage ####################################################
#############################################################################################################

iso_99 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_motifs_30_99.csv")

library(stringr)

#clean up run and motif columns (remove text)
iso_99$run <- as.numeric(substring(iso_99$run, regexpr("run ", iso_99$run) + 4))
iso_99$motif <- as.numeric(substring(iso_99$motif, regexpr("motif ", iso_99$motif) + 6))

motif_iso_summary <- iso_99 %>%
  group_by(motif) %>%
  summarize(mean_iso = mean(isolation_index), var_iso = var(isolation_index))

#now drop unnecessary columns from iso_99 file because we want the level information in there
drops <- c("run", "isolation_index")
iso_99 <- iso_99[ , !(names(iso_99) %in% drops)]
#now select only distinct rows here
iso_99<-iso_99[!duplicated(iso_99$motif), ]

iso_99_summary <- merge(iso_99, motif_iso_summary, by = "motif")

#write to csv file w summary info (mean moran, sd, variance)
write.csv(iso_99_summary, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Isolation Keys/isolation_list_summary_99.csv")

#############################################################################################################
#############################################################################################################

#############################################################################################################
####################################### RUN FOR AGGREGATION #################################################
########################################### 95% coverage ####################################################
#############################################################################################################

#############################################################################################################
#############################################################################################################

setwd("/Users/ninamasters/measles-spatial-model/")
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
#pull in necessary functions to make initial grids
source("/Users/ninamasters/measles-spatial-model/Clustering_Functions_95.R")
source("/Users/ninamasters/measles-spatial-model/aggregation_function_95.R")


isolation_list_motifs = list()


#run this outer loop 30 times to get replicates of Moran's I
for (p in 1:30){
  isolation_list_motifs[[paste("run", p)]] <- list()
  
  #read in violating motifs (i.e. motifs that exceed 1,000 people per cell)
  violate <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/motif_violate_95_threshold.csv")
  names <- seq(1:625)
  names <- names[!(names %in% violate$motif)]
  
  # for for the 625 motifs, run loop to spit out moran's I and clustering at each level
  for(i in unique(names)){
    isolation_list_motifs[[paste("run", p)]][[paste("motif",i)]] <- list()
    
    aggregate_my_motif(combo_list$v1[[i]],combo_list$v2[[i]],combo_list$v3[[i]],combo_list$v4[[i]])
    
    isolation_list_motifs[[paste("run", p)]][[paste("motif", i)]] <- isolation_list
    
  }
}

#now turn list into data frame for 95% coverage with aggregation at all levels
isolation_list <- unlist(isolation_list_motifs, recursive = FALSE, use.names = TRUE)
isolation_list_95_aggregation <- as.data.frame(do.call(Map, c(f = rbind, isolation_list)))

#write to csv
write.csv(isolation_list_95_aggregation, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_95%_vax_aggregation_levels.csv")


iso_95_agg <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_95%_vax_aggregation_levels.csv")

library(stringr)

#clean up run and motif columns (remove text)
iso_95_agg$run <- as.numeric(substring(iso_95_agg$run, regexpr("run ", iso_95_agg$run) + 4))
iso_95_agg$motif <- as.numeric(substring(iso_95_agg$motif, regexpr("motif ", iso_95_agg$motif) + 6))

motif_iso_summary <- iso_95_agg %>%
  group_by(motif) %>%
  summarize(mean_iso_1 = mean(isolation_index_l1), mean_iso_2 = mean(isolation_index_l2), mean_iso_3 = mean(isolation_index_l3), mean_iso_4 = mean(isolation_index_l4))

#now drop unnecessary columns from iso_95 file because we want the level information in there
drops <- c("run", "isolation_index_l1", "isolation_index_l2", "isolation_index_l3", "isolation_index_l4")
iso_95_agg <- iso_95_agg[ , !(names(iso_95_agg) %in% drops)]
#now select only distinct rows here
iso_95_agg<-iso_95_agg[!duplicated(iso_95_agg$motif), ]

iso_95_agg_summary <- merge(iso_95_agg, motif_iso_summary, by = "motif")

#write to csv file w summary info (mean moran, sd, variance)
write.csv(iso_95_agg_summary, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_95%_vax_aggregation_summary.csv")


#############################################################################################################
##########################################    94% Vax    ###################################################
#############################################################################################################

setwd("/Users/ninamasters/measles-spatial-model/")
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
#pull in necessary functions to make initial grids
source("/Users/ninamasters/measles-spatial-model/Clustering_Functions_94.R")
source("/Users/ninamasters/measles-spatial-model/aggregation_function_94.R")


isolation_list_motifs = list()


#run this outer loop 30 times to get replicates of Moran's I
for (p in 1:30){
  isolation_list_motifs[[paste("run", p)]] <- list()
  
  #read in violating motifs (i.e. motifs that exceed 1,000 people per cell)
  violate <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/motif_violate_94_threshold.csv")
  names <- seq(1:625)
  names <- names[!(names %in% violate$motif)]
  
  # for for the 625 motifs, run loop to spit out moran's I and clustering at each level
  for(i in unique(names)){
    isolation_list_motifs[[paste("run", p)]][[paste("motif",i)]] <- list()
    
    aggregate_my_motif(combo_list$v1[[i]],combo_list$v2[[i]],combo_list$v3[[i]],combo_list$v4[[i]])
    
    isolation_list_motifs[[paste("run", p)]][[paste("motif", i)]] <- isolation_list
    
  }
}

#now turn list into data frame for 95% coverage with aggregation at all levels
isolation_list <- unlist(isolation_list_motifs, recursive = FALSE, use.names = TRUE)
isolation_list_94_aggregation <- as.data.frame(do.call(Map, c(f = rbind, isolation_list)))

#write to csv
write.csv(isolation_list_94_aggregation, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_94%_vax_aggregation_levels.csv")


iso_94_agg <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_94%_vax_aggregation_levels.csv")

library(stringr)

#clean up run and motif columns (remove text)
iso_94_agg$run <- as.numeric(substring(iso_94_agg$run, regexpr("run ", iso_94_agg$run) + 4))
iso_94_agg$motif <- as.numeric(substring(iso_94_agg$motif, regexpr("motif ", iso_94_agg$motif) + 6))

motif_iso_summary <- iso_94_agg %>%
  group_by(motif) %>%
  summarize(mean_iso_1 = mean(isolation_index_l1), mean_iso_2 = mean(isolation_index_l2), mean_iso_3 = mean(isolation_index_l3), mean_iso_4 = mean(isolation_index_l4))

#now drop unnecessary columns from iso_94 file because we want the level information in there
drops <- c("run", "isolation_index_l1", "isolation_index_l2", "isolation_index_l3", "isolation_index_l4")
iso_94_agg <- iso_94_agg[ , !(names(iso_94_agg) %in% drops)]
#now select only distinct rows here
iso_94_agg<-iso_94_agg[!duplicated(iso_94_agg$motif), ]

iso_94_agg_summary <- merge(iso_94_agg, motif_iso_summary, by = "motif")

#write to csv file w summary info (mean moran, sd, variance)
write.csv(iso_94_agg_summary, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_94%_vax_aggregation_summary.csv")


#############################################################################################################
##########################################    98% Vax    ###################################################
#############################################################################################################

setwd("/Users/ninamasters/measles-spatial-model/")
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
#pull in necessary functions to make initial grids
source("/Users/ninamasters/measles-spatial-model/Clustering_Functions_98.R")
source("/Users/ninamasters/measles-spatial-model/aggregation_function_98.R")


isolation_list_motifs = list()


#run this outer loop 30 times to get replicates of Moran's I
for (p in 1:30){
  isolation_list_motifs[[paste("run", p)]] <- list()
  
  #read in violating motifs (i.e. motifs that exceed 1,000 people per cell)
  violate <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/motif_violate_98_threshold.csv")
  names <- seq(1:625)
  names <- names[!(names %in% violate$motif)]
  
  # for for the 625 motifs, run loop to spit out moran's I and clustering at each level
  for(i in unique(names)){
    isolation_list_motifs[[paste("run", p)]][[paste("motif",i)]] <- list()
    
    aggregate_my_motif(combo_list$v1[[i]],combo_list$v2[[i]],combo_list$v3[[i]],combo_list$v4[[i]])
    
    isolation_list_motifs[[paste("run", p)]][[paste("motif", i)]] <- isolation_list
    
  }
}

#now turn list into data frame for 95% coverage with aggregation at all levels
isolation_list <- unlist(isolation_list_motifs, recursive = FALSE, use.names = TRUE)
isolation_list_98_aggregation <- as.data.frame(do.call(Map, c(f = rbind, isolation_list)))

#write to csv
write.csv(isolation_list_98_aggregation, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_98%_vax_aggregation_levels.csv")


iso_98_agg <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_98%_vax_aggregation_levels.csv")

library(stringr)

#clean up run and motif columns (remove text)
iso_98_agg$run <- as.numeric(substring(iso_98_agg$run, regexpr("run ", iso_98_agg$run) + 4))
iso_98_agg$motif <- as.numeric(substring(iso_98_agg$motif, regexpr("motif ", iso_98_agg$motif) + 6))

motif_iso_summary <- iso_98_agg %>%
  group_by(motif) %>%
  summarize(mean_iso_1 = mean(isolation_index_l1), mean_iso_2 = mean(isolation_index_l2), mean_iso_3 = mean(isolation_index_l3), mean_iso_4 = mean(isolation_index_l4))

#now drop unnecessary columns from iso_98 file because we want the level information in there
drops <- c("run", "isolation_index_l1", "isolation_index_l2", "isolation_index_l3", "isolation_index_l4")
iso_98_agg <- iso_98_agg[ , !(names(iso_98_agg) %in% drops)]
#now select only distinct rows here
iso_98_agg<-iso_98_agg[!duplicated(iso_98_agg$motif), ]

iso_98_agg_summary <- merge(iso_98_agg, motif_iso_summary, by = "motif")

#write to csv file w summary info (mean moran, sd, variance)
write.csv(iso_98_agg_summary, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_98%_vax_aggregation_summary.csv")


#############################################################################################################
##########################################    99% Vax    ###################################################
#############################################################################################################

setwd("/Users/ninamasters/measles-spatial-model/")
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
#pull in necessary functions to make initial grids
source("/Users/ninamasters/measles-spatial-model/Clustering_Functions_99.R")
source("/Users/ninamasters/measles-spatial-model/aggregation_function_99.R")


isolation_list_motifs = list()


#run this outer loop 30 times to get replicates of Moran's I
for (p in 1:30){
  isolation_list_motifs[[paste("run", p)]] <- list()
  
  #read in violating motifs (i.e. motifs that exceed 1,000 people per cell)
  violate <- c(1,2,6,26,126)
  names <- seq(1:625)
  names <- names[!(names %in% violate)]
  
  # for for the 625 motifs, run loop to spit out moran's I and clustering at each level
  for(i in unique(names)){
    isolation_list_motifs[[paste("run", p)]][[paste("motif",i)]] <- list()
    
    aggregate_my_motif(combo_list$v1[[i]],combo_list$v2[[i]],combo_list$v3[[i]],combo_list$v4[[i]])
    
    isolation_list_motifs[[paste("run", p)]][[paste("motif", i)]] <- isolation_list
    
  }
}

#now turn list into data frame for 95% coverage with aggregation at all levels
isolation_list <- unlist(isolation_list_motifs, recursive = FALSE, use.names = TRUE)
isolation_list_99_aggregation <- as.data.frame(do.call(Map, c(f = rbind, isolation_list)))

#write to csv
write.csv(isolation_list_99_aggregation, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_99%_vax_aggregation_levels.csv")


iso_99_agg <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_99%_vax_aggregation_levels.csv")

library(stringr)

#clean up run and motif columns (remove text)
iso_99_agg$run <- as.numeric(substring(iso_99_agg$run, regexpr("run ", iso_99_agg$run) + 4))
iso_99_agg$motif <- as.numeric(substring(iso_99_agg$motif, regexpr("motif ", iso_99_agg$motif) + 6))

motif_iso_summary <- iso_99_agg %>%
  group_by(motif) %>%
  summarize(mean_iso_1 = mean(isolation_index_l1), mean_iso_2 = mean(isolation_index_l2), mean_iso_3 = mean(isolation_index_l3), mean_iso_4 = mean(isolation_index_l4))

#now drop unnecessary columns from iso_99 file because we want the level information in there
drops <- c("run", "isolation_index_l1", "isolation_index_l2", "isolation_index_l3", "isolation_index_l4")
iso_99_agg <- iso_99_agg[ , !(names(iso_99_agg) %in% drops)]
#now select only distinct rows here
iso_99_agg<-iso_99_agg[!duplicated(iso_99_agg$motif), ]

iso_99_agg_summary <- merge(iso_99_agg, motif_iso_summary, by = "motif")

#write to csv file w summary info (mean moran, sd, variance)
write.csv(iso_99_agg_summary, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/isolation_keys_99%_vax_aggregation_summary.csv")
