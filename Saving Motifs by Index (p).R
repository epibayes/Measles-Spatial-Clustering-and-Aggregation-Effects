########### Saving Motifs in the for loop

setwd("/Users/ninamasters/measles-spatial-model/")
source("/Users/ninamasters/measles-spatial-model/packages_clustering.R")
#pull in necessary functions to make initial grids
source("/Users/ninamasters/measles-spatial-model/Clustering_Functions_95.R")


# Need to save the 625 motifs by their index in the combo_list file for reference (as to their Moran's I, etc.)

moran_list_motifs = list()

#run this outer loop 30 times to get replicates of Moran's I
for (p in 1:30){
    moran_list_motifs[[paste("run", p)]] <- list()
    
    # for for the 625 motifs, run loop to spit out moran's I and clustering at each level
    for(i in 1:625){
      moran_list_motifs[[paste("run", p)]][[paste("motif",i)]] <- list()
        
      make_me_a_moran(combo_list$v1[[i]],combo_list$v2[[i]],combo_list$v3[[i]],combo_list$v4[[i]])
    
      moran_list_motifs[[paste("run", p)]][[paste("motif", i)]] <- moran_list
  
    }
}

#now turn list into data frame for 95% coverage
moran_list <- unlist(moran_list_motifs, recursive = FALSE, use.names = TRUE)
moran_list_95 <- as.data.frame(do.call(Map, c(f = rbind, moran_list)))

#write to csv
write.csv(moran_list_95, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/moran_list_motifs_30_95.csv")



#####################################################################################################################
# clean up / get spread measures on the different motifs

motif_key_95 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_30_95.csv")
library(stringr)

#clean up run and motif columns (remove text)
motif_key_95$run <- as.numeric(substring(motif_key_95$run, regexpr("run ", motif_key_95$run) + 4))
motif_key_95$motif <- as.numeric(substring(motif_key_95$motif, regexpr("motif ", motif_key_95$motif) + 6))

motif_moran_summary <- motif_key_95 %>%
    group_by(motif) %>%
    summarize(mean_Moran = mean(Moran), sd_Moran = sd(Moran), var_Moran = var(Moran))

#now drop unnecessary columns from motif_key_94 file because we want the level information in there
drops <- c("run", "Moran")
motif_key_95 <- motif_key_95[ , !(names(motif_key_95) %in% drops)]
#now select only distinct rows here
motif_key_95<-motif_key_95[!duplicated(motif_key_95$motif), ]

motif_moran_summary <- merge(motif_moran_summary, motif_key_95, by = "motif")

#write to csv file w summary info (mean moran, sd, variance)
write.csv(motif_moran_summary, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_summary_95.csv")

#####################################################################################################################
# 94%

motif_key_94 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_30_94.csv")
library(stringr)

#clean up run and motif columns (remove text)
motif_key_94$run <- as.numeric(substring(motif_key_94$run, regexpr("run ", motif_key_94$run) + 4))
motif_key_94$motif <- as.numeric(substring(motif_key_94$motif, regexpr("motif ", motif_key_94$motif) + 6))

motif_moran_summary <- motif_key_94 %>%
  group_by(motif) %>%
  summarize(mean_Moran = mean(Moran), sd_Moran = sd(Moran), var_Moran = var(Moran))

#now drop unnecessary columns from motif_key_94 file because we want the level information in there
drops <- c("run", "Moran")
motif_key_94 <- motif_key_94[ , !(names(motif_key_94) %in% drops)]
#now select only distinct rows here
motif_key_94<-motif_key_94[!duplicated(motif_key_94$motif), ]

motif_moran_summary <- merge(motif_moran_summary, motif_key_94, by = "motif")

#write to csv file w summary info (mean moran, sd, variance)
write.csv(motif_moran_summary, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_summary_94.csv")

##################################################################################################################
motif_key_99 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_30_99.csv")

library(stringr)

#clean up run and motif columns (remove text)
motif_key_99$run <- as.numeric(substring(motif_key_99$run, regexpr("run ", motif_key_99$run) + 4))
motif_key_99$motif <- as.numeric(substring(motif_key_99$motif, regexpr("motif ", motif_key_99$motif) + 6))

motif_moran_summary <- motif_key_99 %>%
  group_by(motif) %>%
  summarize(mean_Moran = mean(Moran), sd_Moran = sd(Moran), var_Moran = var(Moran))

#now drop unnecessary columns from motif_key_94 file because we want the level information in there
drops <- c("run", "Moran")
motif_key_99 <- motif_key_99[ , !(names(motif_key_99) %in% drops)]
#now select only distinct rows here
motif_key_99<-motif_key_99[!duplicated(motif_key_99$motif), ]

motif_moran_summary <- merge(motif_moran_summary, motif_key_99, by = "motif")

#write to csv file w summary info (mean moran, sd, variance)
write.csv(motif_moran_summary, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_summary_99.csv")


##################################################################################################################
motif_key_98 <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_30_98.csv")

library(stringr)

#clean up run and motif columns (remove text)
motif_key_98$run <- as.numeric(substring(motif_key_98$run, regexpr("run ", motif_key_98$run) + 4))
motif_key_98$motif <- as.numeric(substring(motif_key_98$motif, regexpr("motif ", motif_key_98$motif) + 6))

motif_moran_summary <- motif_key_98 %>%
  group_by(motif) %>%
  summarize(mean_Moran = mean(Moran), sd_Moran = sd(Moran), var_Moran = var(Moran))

#now drop unnecessary columns from motif_key_94 file because we want the level information in there
drops <- c("run", "Moran")
motif_key_98 <- motif_key_98[ , !(names(motif_key_98) %in% drops)]
#now select only distinct rows here
motif_key_98<-motif_key_98[!duplicated(motif_key_98$motif), ]

motif_moran_summary <- merge(motif_moran_summary, motif_key_98, by = "motif")

#write to csv file w summary info (mean moran, sd, variance)
write.csv(motif_moran_summary, "/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Moran Keys/moran_list_motifs_summary_98.csv")


