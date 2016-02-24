#### GLMM model with list length analysis and spatial random effect
#### Inspired by Isaac et al. 2014 and using Gio's code

rm(list = ls())
#### Preliminaries
my_packages<-c('lme4', 'data.table', 'tidyr', 'lattice', 'dplyr')
lapply(my_packages, require, character.only=T)

## Set working directory and load Occupancy functions
setwd("C:/Users/Adam/Documents/UC Berkeley post doc/BIGCB/Pest Project/Potato psyllid")
source("C:/Users/Adam/Documents/GitHub/potato_psyllid_distribution_modeling/museum_specimen_analysis_functions.R")

# Load species lists data set with climate data 
AllLists <- readRDS("Potato psyllid data/All_Hemip_Lists_Climate_2016-02-18.rds")

# Collectors of potato psyllids, from RawRecords data set in making_species_lists
ppCollectors <- readRDS("Potato psyllid data/potato_psyllid_collectors.rds")

# Keep only long lists (ll > 2)
longLists <- AllLists %>% rbindlist() %>% as.data.frame() %>% make_lists(., min.list.length = 3)
longData <- longLists %>% rbindlist() %>% as.data.frame()

# lists that contain collectors of potato psyllids
ppcCollections <- longData[longData$Collector %in% ppCollectors, "collectionID"]
ppcData <- longData[longData$collectionID %in% ppcCollections,]
ppcLists <- ppcData %>% make_lists(., min.list.length = 3)

# Transform to data frame with pp detection
longListsData <- detectDataFunc(ppcLists) # Long lists
ppLists.i <- which(longListsData$detection == 1)


# Exploration for dynamic occupancy model
# Split data by collecting event (year+cell)
# repeatVisits <- split(AllData,AllData$eventID)
# numVisits <- as.numeric(unlist(lapply(1:length(repeatVisits), function(x) length(unique(repeatVisits[[x]]$ymo)))))

# Exploring the data set
hist(detectData$collectors_length) # very few lists with multiple collectors
table(detectData$detection, detectData$list_length) # Still have a problem of declining detections with greater list length
# Distribution of list lengths that contain potato psyllid occurrences
ppData <- detectData[detectData$detection == 1,]
hist(ppData$list_length)
# Most potato psyllid occurrences are in very short lists, with a declining number of long lists containing psyllids


# Exploring long lists (ll >= 3)
longListsData <- detectData[detectData$list_length >= 3,]
hist(longListsData$list_length, breaks = seq(3,40,1))
length(unique(longListsData$cellID)) # 858 different cells -- too many for the site random effect
cellTable <- longListsData %>% group_by(cellID) %>% summarise(nvisits = length(cellID))

#################################################################################
#### GLMMs

# GLMM with all species lists
glmMod <- glmer(detection ~ standardize(log(list_length)) + standardize(year) + standardize(month) + I(standardize(month)^2) + 
                  standardize(aet) + standardize(cwd) + standardize(tmn) + standardize(tmx) + (1|cellID), 
                family = "binomial", data = detectData)              
summary(glmMod)

# GLMM with "long" species lists (length > 2)
# Full model
glmModLL1 <- glmer(detection ~ standardize(log(list_length)) + standardize(year) + standardize(month) + I(standardize(month)^2) + 
                    standardize(aet) + standardize(cwd) + standardize(tmn) + standardize(tmx) + (1|cellID), 
                  family = "binomial", data = longListsData,
                  control = glmerControl(optimizer = "bobyqa"))
# Model without month
glmModLL2 <- glmer(detection ~ standardize(log(list_length)) + standardize(year) +  
                    standardize(aet) + standardize(cwd) + standardize(tmn) + standardize(tmx) + (1|cellID), 
                  family = "binomial", data = longListsData,
                  control = glmerControl(optimizer = "bobyqa"))
# Model with only tmn and tmx
glmModLL3 <- glmer(detection ~ standardize(log(list_length)) + standardize(year) + 
                    standardize(tmn) + standardize(tmx) + (1|cellID), 
                  family = "binomial", data = longListsData,
                  control = glmerControl(optimizer = "bobyqa"))
# Model with only climate variables
glmModLL4 <- glmer(detection ~ standardize(log(list_length)) +  
                    standardize(aet) + standardize(cwd) + standardize(tmn) + standardize(tmx) + (1|cellID), 
                  family = "binomial", data = longListsData,
                  control = glmerControl(optimizer = "bobyqa"))
# Model with only year and LL
glmModLL5 <- glmer(detection ~ standardize(log(list_length)) + standardize(year) +
                    (1|cellID), 
                  family = "binomial", data = longListsData,
                  control = glmerControl(optimizer = "bobyqa"))
AIC(glmModLL1,glmModLL2,glmModLL3,glmModLL4,glmModLL5)
summary(glmModLL4)
