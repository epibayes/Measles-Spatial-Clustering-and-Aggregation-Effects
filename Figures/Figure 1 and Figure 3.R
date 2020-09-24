# Figure 1 Aim 1 Paper
# Need Washtenaw county MI info

library(ggplot2)
library(maptools)
library(maps)
library(rgdal)
library(sp)
library(raster)
library(rgeos)
library(sp)
library(spdep)
library(viridis)
library(hrbrthemes)

require("rgdal") # requires sp, will use proj.4 if installed
require("maptools")
require("ggplot2")
require("plyr")
require("dplyr")
#set working directory
setwd("/Users/ninamasters/measles-spatial-model")

################################################################################################################
######################################### TOP ROW OF FIGURE 1 ##################################################
################################################################################################################

#read in data files
mdhhs_kindergarten <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /Final Cleaned Data/MDHHS_KG_FINAL_DEDUP.csv") 
mdhhs_kindergarten <- mutate(mdhhs_kindergarten, complete_percent = COMPLETE_QY/STUDENT_QY) 
mdhhs_kindergarten <- mutate(mdhhs_kindergarten, waiv_percent = WAIV_TOT_QY / STUDENT_QY) 
mdhhs_kindergarten <- mutate(mdhhs_kindergarten, med_waiv_percent = NUM_MED_WAIV / WAIV_TOT_QY) 
mdhhs_kindergarten <- mutate(mdhhs_kindergarten, phil_waiv_percent = NUM_PHIL_WAIV / WAIV_TOT_QY) 
mdhhs_kindergarten <- mutate(mdhhs_kindergarten, relig_waiv_percent = NUM_RELIG_WAIV / WAIV_TOT_QY) 
mdhhs_kindergarten <- mutate(mdhhs_kindergarten, year = RPT_YR - 1)

#remove few schools w NA in name - couldn't find record
mdhhs_kindergarten <- mdhhs_kindergarten[-grep("N/A", mdhhs_kindergarten$NAME_FOR_ADDRESS_MATCH), ]

kinder_summary <- select(mdhhs_kindergarten, c("NAME_FOR_ADDRESS_MATCH", "MATCH_COUNTY", "year", "STUDENT_QY", "COMPLETE_QY", "WAIV_TOT_QY", "NUM_MED_WAIV", "NUM_RELIG_WAIV", "NUM_PHIL_WAIV"))

#summarize kindergarten overall annual vaccination data to the county level
kinder_county_summary <- aggregate(list(kinder_summary$STUDENT_QY, kinder_summary$COMPLETE_QY,kinder_summary$WAIV_TOT_QY,kinder_summary$NUM_PHIL_WAIV,kinder_summary$NUM_RELIG_WAIV, kinder_summary$NUM_MED_WAIV), by = list(kinder_summary$year, kinder_summary$MATCH_COUNTY), FUN = sum) 
names(kinder_county_summary) = c("YEAR", "COUNTY", "STUDENT_QY", "COMPLETE_QY", "WAIV_TOT_QY", "NUM_PHIL_WAIV", "NUM_RELIG_WAIV", "NUM_MED_WAIV")
kinder_county_summary <- mutate(kinder_county_summary, complete_percent = COMPLETE_QY/STUDENT_QY)
kinder_county_summary <- mutate(kinder_county_summary, waiv_percent = WAIV_TOT_QY / STUDENT_QY) 
kinder_county_summary <- mutate(kinder_county_summary, phil_county_percent = NUM_PHIL_WAIV/STUDENT_QY) 
kinder_county_summary <- mutate(kinder_county_summary, relig_county_percent = NUM_RELIG_WAIV/STUDENT_QY) 
kinder_county_summary <- mutate(kinder_county_summary, med_county_percent = NUM_MED_WAIV/STUDENT_QY) 

### add county-level vaccination data from 2008 - 2018

#now subset to just be washtenaw county
washtenaw_vax <- kinder_county_summary[kinder_county_summary$COUNTY == "washtenaw",]
washtenaw_vax <- subset(washtenaw_vax, YEAR == 2018,)
wayne_vax <- kinder_county_summary[kinder_county_summary$COUNTY == "wayne",]
wayne_vax <- subset(wayne_vax, YEAR == 2018,)
oakland_vax <- kinder_county_summary[kinder_county_summary$COUNTY == "oakland",]
oakland_vax <- subset(oakland_vax, YEAR == 2018,)

#now get all data for oakland county

mdhhs_oakland <- mdhhs_kindergarten[mdhhs_kindergarten$MATCH_COUNTY == "oakland",]
mdhhs_oakland_2018 <- mdhhs_oakland[mdhhs_oakland$RPT_YR == '2019',]
summary(mdhhs_oakland_2018$waiv_percent)

#subset to > 50% waivers
mdhhs_oak_2018_high <- mdhhs_oakland_2018[mdhhs_oakland_2018$waiv_percent >= .5,]

#get school district statistics
oakland_sd_summary <- aggregate(list(mdhhs_oakland$STUDENT_QY, mdhhs_oakland$COMPLETE_QY,mdhhs_oakland$WAIV_TOT_QY,mdhhs_oakland$NUM_PHIL_WAIV,mdhhs_oakland$NUM_RELIG_WAIV, mdhhs_oakland$NUM_MED_WAIV), by = list(mdhhs_oakland$year, mdhhs_oakland$MDHHS_DISTRICT), FUN = sum) 
names(oakland_sd_summary) = c("YEAR", "SCHOOL_DISTRICT", "STUDENT_QY", "COMPLETE_QY", "WAIV_TOT_QY", "NUM_PHIL_WAIV", "NUM_RELIG_WAIV", "NUM_MED_WAIV")
oakland_sd_summary_18_19 <- oakland_sd_summary[oakland_sd_summary$YEAR == 2018,]
oakland_sd_summary_18_19 <- mutate(oakland_sd_summary_18_19, waiv_percent = WAIV_TOT_QY / STUDENT_QY) 

summary(oakland_sd_summary_18_19$waiv_percent)

#get tract-level statistics
oakland_tract_summary <- aggregate(list(mdhhs_oakland$STUDENT_QY, mdhhs_oakland$COMPLETE_QY,mdhhs_oakland$WAIV_TOT_QY,mdhhs_oakland$NUM_PHIL_WAIV,mdhhs_oakland$NUM_RELIG_WAIV, mdhhs_oakland$NUM_MED_WAIV), by = list(mdhhs_oakland$year, mdhhs_oakland$MDHHS_DISTRICT), FUN = sum) 
names(oakland_tract_summary) = c("YEAR", "TRACT", "STUDENT_QY", "COMPLETE_QY", "WAIV_TOT_QY", "NUM_PHIL_WAIV", "NUM_RELIG_WAIV", "NUM_MED_WAIV")
oakland_tract_summary_18_19 <- oakland_tract_summary[oakland_tract_summary$YEAR == 2018,]
oakland_tract_summary_18_19 <- mutate(oakland_tract_summary_18_19, waiv_percent = WAIV_TOT_QY / STUDENT_QY) 

summary(mdhhs_oakland$waiv_percent)

#bring in county file for boundary, subset to specific counties
county_boundary <-rgdal::readOGR("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /MI-Boundary-Data/County/Counties_v17a.shp", "Counties_v17a")
washtenaw <- county_boundary[county_boundary$NAME  == "Washtenaw",]
wayne <- county_boundary[county_boundary$NAME  == "Wayne",]
oakland <- county_boundary[county_boundary$NAME  == "Oakland",]


#for 2018
plotvar <- washtenaw_vax$waiv_percent*100
plotvar <- wayne_vax$waiv_percent*100
plotvar <- oakland_vax$waiv_percent*100

#define color scheme for plot
plotclr <- viridis(9)

colcode <- ifelse((plotvar <=   0), plotclr[1],
                  ifelse((plotvar >   0 & plotvar <=   2), plotclr[2],
                         ifelse((plotvar >   2 & plotvar <=   5), plotclr[3],
                                ifelse((plotvar >   5 & plotvar <=   10), plotclr[4],
                                       ifelse((plotvar >   10 & plotvar <=  15), plotclr[5],
                                              ifelse((plotvar >  15 & plotvar <=  20), plotclr[6],
                                                     ifelse((plotvar >  20 & plotvar <=  40), plotclr[7],
                                                            ifelse((plotvar > 40 & plotvar <= 80), plotclr[8],
                                                                   plotclr[9]))))))))

par(mar=c(0,0,0,0))
plot(washtenaw, col =colcode)

plot(washtenaw, col = colcode)
plot(wayne, col = colcode)
plot(oakland, col = colcode)  



#############################################################################################################
## now we want to overlay on this subset our subdivisions of vaccination coverage: CONGRESSIONAL DISTRICTS
#boundary data
congress_boundary <-rgdal::readOGR("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /MI-Boundary-Data/House_Districts/Michigan_State_House_Districts_v17a.shp", "Michigan_State_House_Districts_v17a")

#vax data
#read in data files
mdhhs_kindergarten <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /Final Cleaned Data/MDHHS_KG_FINAL_DEDUP.csv") 

#merge into data with spatial boundaries
mi_schools_congress <-rgdal::readOGR("/Users/ninamasters/Desktop/Spatial_Join_County/Spatial_Join_Congress_Districts/Spatial_Join_Congress_Districts.shp", "Spatial_Join_Congress_Districts")

mi_congress <- mi_schools_congress@data
mi_congress <- select(mi_congress, c("NAME_FOR_A", "NAME_2", "PARTY"))
names(mi_congress) <- c("NAME_FOR_ADDRESS_MATCH", "CONGRESS_DISTRICT", "PARTY")

# #now merge 
mi_schools_congress <- sp::merge(mdhhs_kindergarten, mi_congress, by = "NAME_FOR_ADDRESS_MATCH", all = FALSE)
mi_schools_congress$year <- mi_schools_congress$RPT_YR - 1

#now summarize kindergarten overall annual vaccination data at the congressional district level
vax_congress_summary <- aggregate(list(mi_schools_congress$STUDENT_QY, mi_schools_congress$COMPLETE_QY,mi_schools_congress$WAIV_TOT_QY,mi_schools_congress$NUM_PHIL_WAIV,mi_schools_congress$NUM_RELIG_WAIV, mi_schools_congress$NUM_MED_WAIV), by = list(mi_schools_congress$year, mi_schools_congress$CONGRESS_DISTRICT), FUN = sum) 
names(vax_congress_summary) = c("YEAR", "CONGRESS_DISTRICT", "STUDENT_QY", "COMPLETE_QY", "WAIV_TOT_QY", "NUM_PHIL_WAIV", "NUM_RELIG_WAIV", "NUM_MED_WAIV")
vax_congress_summary <- mutate(vax_congress_summary, complete_percent = COMPLETE_QY/STUDENT_QY)
vax_congress_summary <- mutate(vax_congress_summary, waiv_percent = WAIV_TOT_QY / STUDENT_QY) 
vax_congress_summary <- mutate(vax_congress_summary, phil_county_percent = NUM_PHIL_WAIV/STUDENT_QY) 
vax_congress_summary <- mutate(vax_congress_summary, relig_county_percent = NUM_RELIG_WAIV/STUDENT_QY) 
vax_congress_summary <- mutate(vax_congress_summary, med_county_percent = NUM_MED_WAIV/STUDENT_QY) 

#restrict to 2018
vax_congress_summ_18 <- dplyr::filter(vax_congress_summary, YEAR == 2018)

# #now merge back in with congressional boundaries
vax_congress <- sp::merge(congress_boundary, vax_congress_summ_18, by.x = "NAME", by.y = "CONGRESS_DISTRICT")


#for 2018
plotvar <- vax_congress$waiv_percent*100

#define color scheme for plot
plotclr <- viridis(9)

colcode <- ifelse((plotvar <=   0), plotclr[1],
                  ifelse((plotvar >   0 & plotvar <=   2), plotclr[2],
                         ifelse((plotvar >   2 & plotvar <=   5), plotclr[3],
                                ifelse((plotvar >   5 & plotvar <=   10), plotclr[4],
                                       ifelse((plotvar >   10 & plotvar <=  15), plotclr[5],
                                              ifelse((plotvar >  15 & plotvar <=  20), plotclr[6],
                                                     ifelse((plotvar >  20 & plotvar <=  40), plotclr[7],
                                                            ifelse((plotvar > 40 & plotvar <= 80), plotclr[8],
                                                                   plotclr[9]))))))))

par(mar=c(0,0,0,0))

#plot overall boundaries and fill ins for whole state of MI
plot(congress_boundary, col = "white")
plot(vax_congress, col = colcode, add = TRUE) 


#bring in county file for boundary, subset to specific counties
county_boundary <-rgdal::readOGR("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /MI-Boundary-Data/County/Counties_v17a.shp", "Counties_v17a")
washtenaw <- county_boundary[county_boundary$NAME  == "Washtenaw",]
wayne <- county_boundary[county_boundary$NAME  == "Wayne",]
oakland <- county_boundary[county_boundary$NAME  == "Oakland",]

#crop to county boundary
### washtenaw
overall_washtenaw_congress <- crop(x= congress_boundary, y = washtenaw)

plot(overall_washtenaw_congress, col = "white")
plot(vax_congress, col = colcode, add = TRUE)
plot(washtenaw, add = TRUE)

plot(washtenaw, col = "black")

### wayne
overall_wayne_congress <- crop(x= congress_boundary, y = wayne)

plot(overall_wayne_congress, col = "white")
plot(vax_congress, col = colcode, add = TRUE)
plot(wayne, add = TRUE)

plot(wayne, col = "black")
#raster image
library(raster)

wayne_outline <- raster::mask(wayne)

### oakland
overall_oakland_congress <- crop(x= congress_boundary, y = oakland)

plot(overall_oakland_congress, col = "white", bg = "black")
plot(vax_congress, col = colcode, add = TRUE)

plot(oakland, col = "black")

#############################################################################################################
## now we want to overlay on this subset our subdivisions of vaccination coverage: CITIES/TOWNSHIPS
## CITIES AND TOWNSHIPS

#boundary data
city_boundary <-rgdal::readOGR("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /MI-Boundary-Data/Cities/Minor_Civil_Divisions_Cities__Townships_v17a.shp", "Minor_Civil_Divisions_Cities__Townships_v17a")
summary(city_boundary$NAME)

#vax data
vax_data <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /Final Linked Vaccine and Census Data/MI full data table with spatial join and vaccination data.csv") 

drops <- c("City", "Score", "DisplayX", "DisplayY", "REVISED_SCHOOL_NAME", "LINK_BG", "LINK_TRACT", "SCHOOL_DISTRICT", "DCODE", "SCHOOL_TYPE")
vax_data <- vax_data[, !(names(vax_data) %in% drops)]

#summarize kindergarten overall annual vaccination data to the county level for 2018
city_vax <- aggregate(list(vax_data$STUDENT_QY_2019, vax_data$COMPLETE_QY_2019,vax_data$WAIV_TOT_QY_2019,vax_data$NUM_PHIL_WAIV_2019,vax_data$NUM_RELIG_WAIV_2019, vax_data$NUM_MED_WAIV_2019), by = list(vax_data$City_Caps), FUN = sum) 
names(city_vax) = c("CITY", "STUDENT_QY", "COMPLETE_QY", "WAIV_TOT_QY", "NUM_PHIL_WAIV", "NUM_RELIG_WAIV", "NUM_MED_WAIV")
city_vax <- mutate(city_vax, complete_percent = COMPLETE_QY/STUDENT_QY)
city_vax <- mutate(city_vax, waiv_percent = WAIV_TOT_QY / STUDENT_QY) 
city_vax <- mutate(city_vax, phil_county_percent = NUM_PHIL_WAIV/STUDENT_QY) 
city_vax <- mutate(city_vax, relig_county_percent = NUM_RELIG_WAIV/STUDENT_QY) 
city_vax <- mutate(city_vax, med_county_percent = NUM_MED_WAIV/STUDENT_QY) 

#now merge 
city_vax <- sp::merge(city_boundary, city_vax, by.x = "NAME", by.y = "CITY", all = FALSE)

#for 2018
plotvar <- city_vax$waiv_percent*100

#define color scheme for plot
plotclr <- viridis(9)

colcode <- ifelse((plotvar <=   0), plotclr[1],
                  ifelse((plotvar >   0 & plotvar <=   2), plotclr[2],
                         ifelse((plotvar >   2 & plotvar <=   5), plotclr[3],
                                ifelse((plotvar >   5 & plotvar <=   10), plotclr[4],
                                       ifelse((plotvar >   10 & plotvar <=  15), plotclr[5],
                                              ifelse((plotvar >  15 & plotvar <=  20), plotclr[6],
                                                     ifelse((plotvar >  20 & plotvar <=  40), plotclr[7],
                                                            ifelse((plotvar > 40 & plotvar <= 80), plotclr[8],
                                                                   plotclr[9]))))))))

par(mar=c(0,0,0,0))

#plot overall boundaries and fill ins for whole state of MI
plot(city_boundary, col = "white")
plot(city_vax, col = colcode, add = TRUE) 


#bring in county file for boundary, subset to specific counties
county_boundary <-rgdal::readOGR("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /MI-Boundary-Data/County/Counties_v17a.shp", "Counties_v17a")
washtenaw <- county_boundary[county_boundary$NAME  == "Washtenaw",]
wayne <- county_boundary[county_boundary$NAME  == "Wayne",]
oakland <- county_boundary[county_boundary$NAME  == "Oakland",]

#crop to county boundary
### washtenaw
overall_washtenaw_city <- crop(x= city_boundary, y = washtenaw)
city_washtenaw <- raster::crop(x = city_vax, y = washtenaw)

plot(overall_washtenaw_city, col = "white")
plot(city_washtenaw, col = colcode, add = TRUE)


### wayne
overall_wayne_city <- crop(x= city_boundary, y = wayne)
city_wayne <- raster::crop(x = city_vax, y = wayne)

plot(overall_wayne_city, col = "white")
plot(city_wayne, col = colcode, add = TRUE)


### oakland
overall_oakland_city <- crop(x= city_boundary, y = oakland)
city_oakland <- crop(city_vax, oakland)

plot(overall_oakland_city, col = "white")
plot(city_oakland, col = colcode, add = TRUE)

#############################################################################################################
## now we want to overlay on this subset our subdivisions of vaccination coverage
## school district level waiver summary info

#vax data
sd_waiv <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /Final Linked Vaccine and Census Data/Kindergarten School District Level Waiver Summary.csv")
drops <- c("X")
sd_waiv <- sd_waiv[ , !(names(sd_waiv) %in% drops)]

#spatial data
sd_boundary <-rgdal::readOGR("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /MI-Boundary-Data/School District/School_Districts_v17a.shp", "School_Districts_v17a")

#now merge 
sd_vax <- sp::merge(sd_boundary, sd_waiv, by = "DCODE", all = FALSE)


#for 2018
plotvar <- sd_vax$waiv_percent_2018

#define color scheme for plot
plotclr <- viridis(9)

colcode <- ifelse((plotvar <=   0), plotclr[1],
                  ifelse((plotvar >   0 & plotvar <=   2), plotclr[2],
                         ifelse((plotvar >   2 & plotvar <=   5), plotclr[3],
                                ifelse((plotvar >   5 & plotvar <=   10), plotclr[4],
                                       ifelse((plotvar >   10 & plotvar <=  15), plotclr[5],
                                              ifelse((plotvar >  15 & plotvar <=  20), plotclr[6],
                                                     ifelse((plotvar >  20 & plotvar <=  40), plotclr[7],
                                                            ifelse((plotvar > 40 & plotvar <= 80), plotclr[8],
                                                                   plotclr[9]))))))))

par(mar=c(0,0,0,0))

#plot overall boundaries and fill ins for whole state of MI
plot(sd_boundary, col = "white")
plot(sd_vax, col = colcode, add = TRUE) 


#bring in county file for boundary, subset to specific counties
county_boundary <-rgdal::readOGR("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /MI-Boundary-Data/County/Counties_v17a.shp", "Counties_v17a")
washtenaw <- county_boundary[county_boundary$NAME  == "Washtenaw",]
wayne <- county_boundary[county_boundary$NAME  == "Wayne",]
oakland <- county_boundary[county_boundary$NAME  == "Oakland",]

#crop to county boundary
### washtenaw
overall_washtenaw_sd <- crop(x= sd_boundary, y = washtenaw)
sd_washtenaw <- raster::crop(x = sd_vax, y = washtenaw)

plot(overall_washtenaw_sd, col = "white")
plot(sd_washtenaw, col = colcode, add = TRUE)


### wayne
overall_wayne_sd <- crop(x= sd_boundary, y = wayne)
sd_wayne <- raster::crop(x = sd_vax, y = wayne)

plot(overall_wayne_sd, col = "white")
plot(sd_wayne, col = colcode, add = TRUE)


### oakland
overall_oakland_sd <- crop(x= sd_boundary, y = oakland)
sd_oakland <- crop(sd_vax, oakland)

plot(overall_oakland_sd, col = "white")
plot(sd_oakland, col = colcode, add = TRUE)


###################################################################################
### now bring in census tracts
#vax data
tract_waiv <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /Final Linked Vaccine and Census Data/Kindergarten Tract Level Waiver Summary.csv")
drops <- c("X")
tract_waiv <- tract_waiv[ , !(names(tract_waiv) %in% drops)]

#spatial data
tract_boundary <-rgdal::readOGR("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /MI-Boundary-Data/Tract/2010_Census_Tracts_v17a.shp", "2010_Census_Tracts_v17a")

#now merge 
tract_vax <- sp::merge(tract_boundary, tract_waiv, by.x = "LINK", by.y = "LINK_TRACT", all = FALSE)

#for 2018
plotvar <- tract_vax$waiv_percent_2018

#define color scheme for plot
plotclr <- viridis(9)

colcode <- ifelse((plotvar <=   0), plotclr[1],
                  ifelse((plotvar >   0 & plotvar <=   2), plotclr[2],
                         ifelse((plotvar >   2 & plotvar <=   5), plotclr[3],
                                ifelse((plotvar >   5 & plotvar <=   10), plotclr[4],
                                       ifelse((plotvar >   10 & plotvar <=  15), plotclr[5],
                                              ifelse((plotvar >  15 & plotvar <=  20), plotclr[6],
                                                     ifelse((plotvar >  20 & plotvar <=  40), plotclr[7],
                                                            ifelse((plotvar > 40 & plotvar <= 80), plotclr[8],
                                                                   plotclr[9]))))))))

par(mar=c(0,0,0,0))

plot(tract_boundary, col = "white")
plot(tract_vax, col = colcode, add = TRUE)

#bring in county file for boundary
county_boundary <-rgdal::readOGR("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /MI-Boundary-Data/County/Counties_v17a.shp", "Counties_v17a")
washtenaw <- county_boundary[county_boundary$NAME  == "Washtenaw",]
wayne <- county_boundary[county_boundary$NAME  == "Wayne",]
oakland <- county_boundary[county_boundary$NAME  == "Oakland",]

common.crs <- CRS(proj4string(washtenaw))
tract_vax_reprojected <- spTransform(tract_vax, common.crs)

CRS(proj4string(washtenaw))
CRS(proj4string(tract_vax_reprojected))

require(rgdal)
require(rgeos)
require(maps)
#crop to county boundary

washtenaw@proj4string
tract_vax@proj4string

plot(tract_vax)
dim(washtenaw)
tract_washtenaw <- crop(x = tract_vax, y = washtenaw)

tract_wayne <- crop(x = tract_vax, y = wayne)
par(mar=c(0,0,0,0))

### washtenaw
overall_washtenaw_tract <- crop(x= tract_boundary, y = washtenaw)
tract_washtenaw <- crop(tract_vax, washtenaw)

plot(overall_washtenaw_tract, col = "white")
plot(tract_washtenaw, col = colcode, add = TRUE)
legend("right", legend = c("<=0%", "0-2%", "2-5%", "5-10%", "10-15%", "15-20%", "20-40%", "40-80%", ">80%"), col = viridis(9), lty = 1, lwd =30, cex = 1)


### wayne
overall_wayne_tract <- crop(x= tract_boundary, y = wayne)
tract_wayne<- crop(tract_vax, wayne)

plot(overall_wayne_tract, col = "white")
plot(tract_wayne, col = colcode, add = TRUE)
legend("right", legend = c("<=0%", "0-2%", "2-5%", "5-10%", "10-15%", "15-20%", "20-40%", "40-80%", ">80%"), col = viridis(9), lty = 1, lwd =30, cex = 1)

### oakland
overall_oakland_tract <- crop(x= tract_boundary, y = oakland)
tract_oakland <- crop(tract_vax, oakland)

plot(overall_oakland_tract, col = "white")
plot(tract_oakland, col = colcode, add = TRUE)

plot(tract_washtenaw, col = colcode)
text(coordinates(sd_washtenaw$ID),  sd_washtenaw$NAME, pos = 4)


plot(sd_washtenaw, col = colcode)


######################################################################
### now bring in census block groups
#vax data
bg_waiv <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /Final Linked Vaccine and Census Data/Kindergarten Block-Group Level Waiver Summary.csv")
drops <- c("X")
bg_waiv <- bg_waiv[ , !(names(bg_waiv) %in% drops)]

#spatial data
bg_boundary <-rgdal::readOGR("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /MI-Boundary-Data/Block-Group/2010_Block_Groups_v17a.shp", "2010_Block_Groups_v17a")

#now merge 
bg_vax <- sp::merge(bg_boundary, bg_waiv, by.x = "LINK",  by.y = "LINK_BG", all = FALSE)

#for 2018
plotvar <- bg_vax$waiv_percent_2018

#define color scheme for plot
plotclr <- viridis(9)

colcode <- ifelse((plotvar <=   0), plotclr[1],
                  ifelse((plotvar >   0 & plotvar <=   2), plotclr[2],
                         ifelse((plotvar >   2 & plotvar <=   5), plotclr[3],
                                ifelse((plotvar >   5 & plotvar <=   10), plotclr[4],
                                       ifelse((plotvar >   10 & plotvar <=  15), plotclr[5],
                                              ifelse((plotvar >  15 & plotvar <=  20), plotclr[6],
                                                     ifelse((plotvar >  20 & plotvar <=  40), plotclr[7],
                                                            ifelse((plotvar > 40 & plotvar <= 80), plotclr[8],
                                                                   plotclr[9]))))))))

par(mar=c(0,0,0,0))

plot(bg_boundary, col = "white")
plot(bg_vax, col = colcode, add = TRUE) 


#bring in county file for boundary
county_boundary <-rgdal::readOGR("/Users/ninamasters/Desktop/Dissertation/Aim 3 - MDHHS /MI-Boundary-Data/County/Counties_v17a.shp", "Counties_v17a")
washtenaw <- county_boundary[county_boundary$NAME  == "Washtenaw",]
wayne <- county_boundary[county_boundary$NAME  == "Wayne",]
oakland <- county_boundary[county_boundary$NAME  == "Oakland",]

require(rgdal)
require(rgeos)

#crop to county boundary
#first make sure all in same projection
# 
common.crs <- CRS(proj4string(county_boundary))
bg_vax_reprojected <- spTransform(bg_vax, common.crs)

bg_washtenaw_boundary <- crop(x = bg_boundary, y = washtenaw)
bg_washtenaw <- crop(x = bg_vax_reprojected, y = washtenaw)


bg_oakland_boundary <- crop(x = bg_boundary, y = oakland)
bg_oakland <- crop(x = bg_vax_reprojected, y = oakland)

plot(bg_oakland_boundary, col = "white")
plot(bg_oakland, col = colcode, add = TRUE)
legend("right", legend = c("<=0%", "0-2%", "2-5%", "5-10%", "10-15%", "15-20%", "20-40%", "40-80%", ">80%"), col = viridis(9), lty = 1, lwd =30, cex = 1)



par(mar=c(0,0,0,0))
plot(bg_washtenaw_boundary, col = "white")
plot(bg_washtenaw, col = colcode, add = TRUE)
legend("bottom", legend = c("Waiver Percentage, 2018"), col = colcode, lty = 1, lwd = 15, cex = 1)

overall_wayne_bg <- crop(x= bg_boundary, y = wayne)
bg_wayne<- crop(bg_vax, wayne)

plot(overall_wayne_bg, col = "white")
plot(bg_wayne, col = colcode, add = TRUE)
legend("right", legend = c("<=0%", "0-2%", "2-5%", "5-10%", "10-15%", "15-20%", "20-40%", "40-80%", ">80%"), col = viridis(9), lty = 1, lwd =30, cex = 1)

################################################################################################################
######################################### MID ROW OF FIGURE 1 ##################################################
################################################################################################################

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

par(mfrow = c(1,1))
par(mar=c(1,1,1,1))
#plot grid
plot(gridspdf, xlim = c(0,17), ylim = c(0,17), main = "16 x 16 Grid of cells")
gridspdf$ID <- sapply(slot(gridspdf, "polygons"), function(x) slot(x, "ID"))
#overlay the cell ID onto the plot
text(coordinates(gridspdf),labels=sapply(slot(gridspdf,"polygons"),function(i) slot(i,"ID")), cex=0.5)


# set up adjacencies / boundaries using both queen and rook style boundaries
queen_boundaries_grid <- poly2nb(gridspdf, queen = T)
rook_boundaries_grid <- poly2nb(gridspdf, queen = F)
par(mfrow = c(1,1))
par(mar=c(1,1,1,1))
#plot both the rook and queen contiguity to ensure everything worked
plot(queen_boundaries_grid, coords = xy, col = "red", xlim = c(0,17), ylim = c(0,17), main = "Queen Contiguity Boundaries", xlab = "x coordinates", ylab = "y coordinates")
plot(rook_boundaries_grid, coords = xy, col = "blue", add = T)

###########################################################################################################
############ Merge in levels of different aggregation to overlay on the grid ##############################
###########################################################################################################
require(sp) 
level <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Clustering Motifs/level.csv")
motifs <- merge(gridspdf,level, by="ID")

plotvar <- motifs$level1
plotvar2 <- motifs$level2
plotvar3 <- motifs$level3
plotvar4 <- motifs$level4


library("viridis")

c1 <- viridis(20)
col <- c1[5]

plotclr <- c("white", col)

colcode <- ifelse((plotvar == 0), plotclr[1], plotclr[2])
colcode2 <- ifelse((plotvar2 == 0), plotclr[1], plotclr[2])
colcode3 <- ifelse((plotvar3 == 0), plotclr[1], plotclr[2])
colcode4 <- ifelse((plotvar4 == 0), plotclr[1], plotclr[2])

par(mar=c(1,1,1,1))
#plot level 1
plot(motifs, col = colcode, main = "Level 1 Grid Resolution")
lines(motifs, col = "white", lwd = .25, lty = 1)
#plot level 2
plot(motifs, col = colcode2, main = "Level 2 Grid Resolution")
lines(motifs, col = "white", lwd = .25, lty = 1)
#plot level 3
plot(motifs, col = colcode3, main = "Level 3 Grid Resolution")
lines(motifs, col = "white", lwd = .25, lty = 1)
#plot level 4
plot(motifs, col = colcode4, main = "Level 4 Grid Resolution")
lines(motifs, col = "white", lwd = .25, lty = 1)

######################################################################################################################
############################### Now third component of the figure is aggregation stuff ###############################
######################################################################################################################

# FIRST READ IN FINAL TIME AGGREGATION SUBSET FOR 95%
agg_95_f <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/aggregation_end_time_data_95%_vax.csv")

# drop x
drops <- c("X")
agg_95_f <- agg_95_f[ , !(names(agg_95_f) %in% drops)]

# create new columns for cumulative incidence
agg_95_f$CI_l1 <- 12799- agg_95_f$S1
agg_95_f$CI_l2 <- 12799- agg_95_f$S2
agg_95_f$CI_l3 <- 12799- agg_95_f$S3
agg_95_f$CI_l4 <- 12799- agg_95_f$S4

summary(agg_95_f$CI_l1)
sd(agg_95_f$CI_l1)
summary(agg_95_f$CI_l2)
sd(agg_95_f$CI_l2)
summary(agg_95_f$CI_l3)
sd(agg_95_f$CI_l3)
summary(agg_95_f$CI_l4)
sd(agg_95_f$CI_l4)

#create colums for difference
agg_95_f$CI_l2_diff <- -(agg_95_f$CI_l1-agg_95_f$CI_l2)
agg_95_f$CI_l3_diff <- -(agg_95_f$CI_l1-agg_95_f$CI_l3)
agg_95_f$CI_l4_diff <- -(agg_95_f$CI_l1-agg_95_f$CI_l4)

agg_95_f$CI_l2_diff_pct <- (agg_95_f$CI_l1-agg_95_f$CI_l2)/agg_95_f$CI_l1
agg_95_f$CI_l3_diff_pct <- (agg_95_f$CI_l1-agg_95_f$CI_l3)/agg_95_f$CI_l1
agg_95_f$CI_l4_diff_pct <- (agg_95_f$CI_l1-agg_95_f$CI_l4)/agg_95_f$CI_l1

summary(agg_95_f$CI_l2_diff_pct)
summary(agg_95_f$CI_l3_diff_pct)
summary(agg_95_f$CI_l4_diff_pct)

#create column for percentage of cases identified
agg_95_f$CI_l2_prop <- agg_95_f$CI_l2/agg_95_f$CI_l1
agg_95_f$CI_l3_prop <- agg_95_f$CI_l3/agg_95_f$CI_l1
agg_95_f$CI_l4_prop <- agg_95_f$CI_l4/agg_95_f$CI_l1

######################################## 94% ##################################################

#Just read in final time threshold subset file
agg_94_f <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/aggregation_end_time_data_94%_vax.csv")

# drop x
drops <- c("X")
agg_94_f <- agg_94_f[ , !(names(agg_94_f) %in% drops)]

# create new columns for cumulative incidence
agg_94_f$CI_l1 <- 15359- agg_94_f$S1
agg_94_f$CI_l2 <- 15359- agg_94_f$S2
agg_94_f$CI_l3 <- 15359- agg_94_f$S3
agg_94_f$CI_l4 <- 15359- agg_94_f$S4

summary(agg_94_f$CI_l1)
sd(agg_94_f$CI_l1)
summary(agg_94_f$CI_l2)
sd(agg_94_f$CI_l2)
summary(agg_94_f$CI_l3)
sd(agg_94_f$CI_l3)
summary(agg_94_f$CI_l4)
sd(agg_94_f$CI_l4)

#create colums for difference
agg_94_f$CI_l2_diff <- -(agg_94_f$CI_l1-agg_94_f$CI_l2)
agg_94_f$CI_l3_diff <- -(agg_94_f$CI_l1-agg_94_f$CI_l3)
agg_94_f$CI_l4_diff <- -(agg_94_f$CI_l1-agg_94_f$CI_l4)

agg_94_f$CI_l2_diff_pct <- (agg_94_f$CI_l1-agg_94_f$CI_l2)/agg_94_f$CI_l1
agg_94_f$CI_l3_diff_pct <- (agg_94_f$CI_l1-agg_94_f$CI_l3)/agg_94_f$CI_l1
agg_94_f$CI_l4_diff_pct <- (agg_94_f$CI_l1-agg_94_f$CI_l4)/agg_94_f$CI_l1

summary(agg_94_f$CI_l2_diff_pct)
summary(agg_94_f$CI_l3_diff_pct)
summary(agg_94_f$CI_l4_diff_pct)

#create column for percentage of cases identified
agg_94_f$CI_l2_prop <- agg_94_f$CI_l2/agg_94_f$CI_l1
agg_94_f$CI_l3_prop <- agg_94_f$CI_l3/agg_94_f$CI_l1
agg_94_f$CI_l4_prop <- agg_94_f$CI_l4/agg_94_f$CI_l1

########################################### 98% ##################################################

#Just read in final time threshold subset file
agg_98_f <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/aggregation_end_time_data_98%_vax.csv")

# drop x
drops <- c("X")
agg_98_f <- agg_98_f[ , !(names(agg_98_f) %in% drops)]

# create new columns for cumulative incidence
agg_98_f$CI_l1 <- 5120- agg_98_f$S1
agg_98_f$CI_l2 <- 5120- agg_98_f$S2
agg_98_f$CI_l3 <- 5120- agg_98_f$S3
agg_98_f$CI_l4 <- 5120- agg_98_f$S4

summary(agg_98_f$CI_l1)
sd(agg_98_f$CI_l1)
summary(agg_98_f$CI_l2)
sd(agg_98_f$CI_l2)
summary(agg_98_f$CI_l3)
sd(agg_98_f$CI_l3)
summary(agg_98_f$CI_l4)
sd(agg_98_f$CI_l4)

#create colums for difference
agg_98_f$CI_l2_diff <- -(agg_98_f$CI_l1-agg_98_f$CI_l2)
agg_98_f$CI_l3_diff <- -(agg_98_f$CI_l1-agg_98_f$CI_l3)
agg_98_f$CI_l4_diff <- -(agg_98_f$CI_l1-agg_98_f$CI_l4)


agg_98_f$CI_l2_diff_pct <- (agg_98_f$CI_l1-agg_98_f$CI_l2)/agg_98_f$CI_l1
agg_98_f$CI_l3_diff_pct <- (agg_98_f$CI_l1-agg_98_f$CI_l3)/agg_98_f$CI_l1
agg_98_f$CI_l4_diff_pct <- (agg_98_f$CI_l1-agg_98_f$CI_l4)/agg_98_f$CI_l1


summary(agg_98_f$CI_l2_diff_pct)
summary(agg_98_f$CI_l3_diff_pct)
summary(agg_98_f$CI_l4_diff_pct)

#create column for percentage of cases identified
agg_98_f$CI_l2_prop <- agg_98_f$CI_l2/agg_98_f$CI_l1
agg_98_f$CI_l3_prop <- agg_98_f$CI_l3/agg_98_f$CI_l1
agg_98_f$CI_l4_prop <- agg_98_f$CI_l4/agg_98_f$CI_l1

#################################### ANALYSIS OF 98% ##################################################

#Just read in final time threshold subset file
agg_99_f <- read.csv("/Users/ninamasters/Desktop/Dissertation/Aim 1 - Spatial Model/Simulation Output/Aggregation/aggregation_end_time_data_99%_vax.csv")

# drop x
drops <- c("X")
agg_99_f <- agg_99_f[ , !(names(agg_99_f) %in% drops)]

# create new columns for cumulative incidence
agg_99_f$CI_l1 <- 2560- agg_99_f$S1
agg_99_f$CI_l2 <- 2560- agg_99_f$S2
agg_99_f$CI_l3 <- 2560- agg_99_f$S3
agg_99_f$CI_l4 <- 2560- agg_99_f$S4

summary(agg_99_f$CI_l1)
sd(agg_99_f$CI_l1)
summary(agg_99_f$CI_l2)
sd(agg_99_f$CI_l2)
summary(agg_99_f$CI_l3)
sd(agg_99_f$CI_l3)
summary(agg_99_f$CI_l4)
sd(agg_99_f$CI_l4)

#create colums for difference
agg_99_f$CI_l2_diff <- -(agg_99_f$CI_l1-agg_99_f$CI_l2)
agg_99_f$CI_l3_diff <- -(agg_99_f$CI_l1-agg_99_f$CI_l3)
agg_99_f$CI_l4_diff <- -(agg_99_f$CI_l1-agg_99_f$CI_l4)


agg_99_f$CI_l2_diff_pct <- (agg_99_f$CI_l1-agg_99_f$CI_l2)/agg_99_f$CI_l1
agg_99_f$CI_l3_diff_pct <- (agg_99_f$CI_l1-agg_99_f$CI_l3)/agg_99_f$CI_l1
agg_99_f$CI_l4_diff_pct <- (agg_99_f$CI_l1-agg_99_f$CI_l4)/agg_99_f$CI_l1


summary(agg_99_f$CI_l2_diff_pct)
summary(agg_99_f$CI_l3_diff_pct)
summary(agg_99_f$CI_l4_diff_pct)

#create column for percentage of cases identified
agg_99_f$CI_l2_prop <- agg_99_f$CI_l2/agg_99_f$CI_l1
agg_99_f$CI_l3_prop <- agg_99_f$CI_l3/agg_99_f$CI_l1
agg_99_f$CI_l4_prop <- agg_99_f$CI_l4/agg_99_f$CI_l1

######################################################################################################################
######################## Now create a master file that has summary data for all vax levels ###########################
######################################################################################################################
require(dplyr)
agg_94_summary <- dplyr::select(agg_94_f,c(CI_l1, CI_l2, CI_l3, CI_l4))
agg_94_summary$pct <- rep(94,1184)

agg_95_summary <- dplyr::select(agg_95_f,c(CI_l1, CI_l2, CI_l3, CI_l4))
agg_95_summary$pct <- rep(95,1344)

agg_98_summary <- dplyr::select(agg_98_f,c(CI_l1, CI_l2, CI_l3, CI_l4))
agg_98_summary$pct <- rep(98,2172)

agg_99_summary <- dplyr::select(agg_99_f,c(CI_l1, CI_l2, CI_l3, CI_l4))
agg_99_summary$pct <- rep(99,2480)

agg_overall_summary <- rbind(agg_94_summary, agg_95_summary, agg_98_summary, agg_99_summary)

library(reshape2)

agg_melt <- melt(agg_overall_summary, id.vars = c("pct"), variable.name = "aggregation_level", value.name = "CI")
agg_melt$pct <- as.factor(agg_melt$pct)

agg_melt$key <- paste(agg_melt$pct, agg_melt$aggregation_level)

#take means of the CI to use as fill
agg_means_94 <- as.data.frame(0)
agg_means_94$pct <- 94
agg_means_94$CI_l1 <- mean(agg_94_summary$CI_l1)
agg_means_94$CI_l2 <- mean(agg_94_summary$CI_l2)
agg_means_94$CI_l3 <- mean(agg_94_summary$CI_l3)
agg_means_94$CI_l4 <- mean(agg_94_summary$CI_l4)

agg_means_95 <- as.data.frame(0)
agg_means_95$pct <- 95
agg_means_95$CI_l1 <- mean(agg_95_summary$CI_l1)
agg_means_95$CI_l2 <- mean(agg_95_summary$CI_l2)
agg_means_95$CI_l3 <- mean(agg_95_summary$CI_l3)
agg_means_95$CI_l4 <- mean(agg_95_summary$CI_l4)

agg_means_98 <- as.data.frame(0)
agg_means_98$pct <- 98
agg_means_98$CI_l1 <- mean(agg_98_summary$CI_l1)
agg_means_98$CI_l2 <- mean(agg_98_summary$CI_l2)
agg_means_98$CI_l3 <- mean(agg_98_summary$CI_l3)
agg_means_98$CI_l4 <- mean(agg_98_summary$CI_l4)

agg_means_99 <- as.data.frame(0)
agg_means_99$pct <- 99
agg_means_99$CI_l1 <- mean(agg_99_summary$CI_l1)
agg_means_99$CI_l2 <- mean(agg_99_summary$CI_l2)
agg_means_99$CI_l3 <- mean(agg_99_summary$CI_l3)
agg_means_99$CI_l4 <- mean(agg_99_summary$CI_l4)

agg_means <- rbind(agg_means_94, agg_means_95, agg_means_98, agg_means_99) %>%
  dplyr::select(c(pct, CI_l1, CI_l2, CI_l3, CI_l4))

agg_means <- melt(agg_means, id.vars = c("pct"), variable.name = "aggregation_level", value.name = "mean_CI")
agg_means$key <- paste(agg_means$pct, agg_means$aggregation_level)
agg_means <- dplyr::select(agg_means, c(key, mean_CI))

# now merge in aggregate means to melted data
agg_melt <- merge(agg_melt, agg_means, by = "key")

library(ggplot2)
library(viridis)
library(hrbrthemes)
library(dplyr)
library(ggthemes)

#bubbleplot
ggplot(data = agg_melt, aes(x = aggregation_level, y = pct, size = CI, fill = CI)) + geom_point(alpha = 0.5, shape = 21, color = "black") +
  scale_size(range = c(0.1, 50), name = "Outbreak Size") +
  scale_fill_viridis() +
  theme_minimal() 


#bubbleplot
ggplot(data = agg_melt) + geom_point(aes(x = aggregation_level, y = pct, size = CI), alpha = 0.15, shape = 21, fill = "NA", col = "darkgrey") +
  geom_point(aes(x = aggregation_level, y = pct, size = mean_CI, fill = mean_CI), alpha = 0.5, shape = 21, color = "black") +
  scale_size(range = c(0.1, 45), name = "Outbreak Size") +
  scale_fill_viridis(option ="B") +
  theme_clean() 



# try boxplot
ggplot(data = agg_melt) + geom_boxplot(aes(y = CI, fill = mean)) + facet_wrap(pct~aggregation_level) +
  theme_void() 

#subset to 94
ggplot(data = agg_melt[agg_melt$pct == "94",], aes(x = aggregation_level, y = CI, fill = mean)) + geom_boxplot() 

View(agg_melt[agg_melt$pct == "94",])
