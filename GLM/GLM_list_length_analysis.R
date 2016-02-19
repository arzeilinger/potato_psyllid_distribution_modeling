#### GLMM model with list length analysis and spatial random effect
#### Inspired by Isaac et al. 2014 and using Gio's code

#### Preliminaries
my_packages<-c('lme4', 'data.table', 'tidyr', 'lattice', 'dplyr')
lapply(my_packages, require, character.only=T)

# Load California polygon, object named CAsp
#load("CA_SpatialPolygons.Rdata")
# Load Occupancy functions
rm(list = ls())
## Set working directory
setwd("C:/Users/Adam/Documents/UC Berkeley post doc/BIGCB/Pest Project/Potato psyllid")
source("C:/Users/Adam/Documents/GitHub/potato_psyllid_distribution_modeling/museum_specimen_analysis_functions.R")

# Load species lists data set with climate data 
AllLists <- readRDS("Potato psyllid data/All_Hemip_Lists_Climate_2016-02-18.rds")
# Transform to data frame with pp detection
detectData <- detectDataFunc(AllLists) # All lists

# Exploration for dynamic occupancy model
# Split data by collecting event (year+cell)
# repeatVisits <- split(AllData,AllData$eventID)
# numVisits <- as.numeric(unlist(lapply(1:length(repeatVisits), function(x) length(unique(repeatVisits[[x]]$ymo)))))

# Exploring the data set
hist(detectData[detectData$list_length >= 3,]$list_length, breaks = seq(3,40,1))
hist(detectData$collectors_length) # very few lists with multiple collectors
table(detectData$detection, detectData$list_length) # Still have a problem of declining detections with greater list length

# GLMM with all species lists
glmMod <- glmer(detection ~ standardize(log(list_length)) + standardize(year) + standardize(month) + I(standardize(month)^2) + 
                  standardize(aet) + standardize(cwd) + standardize(tmn) + standardize(tmx) + (1|cellID), 
                family = "binomial", data = detectData)              

# GLMM with "long" species lists (length > 2)
longListsData <- detectData[detectData$list_length >= 3,]
glmModLL <- glmer(detection ~ standardize(log(list_length)) + standardize(year) + standardize(month) + I(standardize(month)^2) + 
                    standardize(aet) + standardize(cwd) + standardize(tmn) + standardize(tmx) + (1|cellID), 
                  family = "binomial", data = longListsData)              
