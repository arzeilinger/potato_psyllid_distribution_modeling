#### GLMM model with list length analysis and spatial random effect
#### Inspired by Isaac et al. 2014 and using Gio's code

rm(list = ls())
#### Preliminaries
my_packages<-c('lme4', 'lmerTest', 'data.table', 'tidyr', 'lattice', 'dplyr', 'bbmle', 'optimx')
lapply(my_packages, require, character.only=T)

## Set working directory and load Occupancy functions
setwd("C:/Users/Adam/Documents/UC Berkeley post doc/BIGCB/Pest Project/Potato psyllid")
source("C:/Users/Adam/Documents/GitHub/potato_psyllid_distribution_modeling/museum_specimen_analysis_functions.R")

# Load species lists data set with climate data 
AllLists <- readRDS("Potato psyllid data/All_Hemip_Lists_Climate_15km_Cells_2016-02-24.rds")

# Collectors of potato psyllids, from RawRecords data set in making_species_lists
ppCollectors <- readRDS("Potato psyllid data/potato_psyllid_collectors.rds")

# Keep only long lists (ll > 2)
longLists <- AllLists %>% rbindlist() %>% as.data.frame() %>% make_lists(., min.list.length = 3)
longListsDF <- longLists %>% rbindlist() %>% as.data.frame()

# lists that contain collectors of potato psyllids
ppcCollections <- longListsDF[longListsDF$Collector %in% ppCollectors, "collectionID"]
ppcData <- longListsDF[longListsDF$collectionID %in% ppcCollections,]
ppcLists <- ppcData %>% make_lists(., min.list.length = 3)

# Transform to data frame with pp detection
detectData <- detectDataFunc(ppcLists) # Long lists


# Exploring data
table(detectData$detection)
with(detectData, hist(list_length, breaks = seq(min(list_length),max(list_length),1)))
ppData <- detectData[detectData$detection == 1,]
hist(ppData$list_length)

length(unique(detectData$cellID)) # 858 different cells -- too many for the site random effect
detectData %>% group_by(cellID) %>% summarise(nvisits = length(cellID))

#################################################################################
#### GLMMs
require(optimx)
# GLMM with "long" species lists (length > 2)
# Full model
glmModLL1 <- glmer(detection ~ standardize(log(list_length)) + standardize(year) + standardize(month) + I(standardize(month)^2) + 
                    standardize(aet) + standardize(cwd) + standardize(tmn) + standardize(tmx) + (1|cellID), 
                  family = "binomial", data = detectData,
                  control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
# Model without month
glmModLL2 <- glmer(detection ~ standardize(log(list_length)) + standardize(year) +  
                    standardize(aet) + standardize(cwd) + standardize(tmn) + standardize(tmx) + (1|cellID), 
                  family = "binomial", data = detectData,
                  control = glmerControl(optimizer = "Nelder_Mead"))
# Model with only tmn and tmx
glmModLL3 <- glmer(detection ~ standardize(log(list_length)) + standardize(year) + 
                    standardize(tmn) + standardize(tmx) + (1|cellID), 
                  family = "binomial", data = detectData,
                  control = glmerControl(optimizer = "bobyqa"))
# Model with only climate variables
glmModLL4 <- glmer(detection ~ standardize(log(list_length)) +  
                    standardize(aet) + standardize(cwd) + standardize(tmn) + standardize(tmx) + (1|cellID), 
                  family = "binomial", data = detectData,
                  control = glmerControl(optimizer = "bobyqa"))
# Model with only year and LL
glmModLL5 <- glmer(detection ~ standardize(log(list_length)) + standardize(year) +
                    (1|cellID), 
                  family = "binomial", data = detectData,
                  control = glmerControl(optimizer = "bobyqa"))
# Model with year*month interaction
glmModLL6 <- glmer(detection ~ standardize(log(list_length)) + standardize(year)*standardize(month) +
                     standardize(year)*I(standardize(month)^2) +
                    (1|cellID), 
                  family = "binomial", data = detectData,
                  control = glmerControl(optimizer = "bobyqa"))
AICtab(glmModLL1,glmModLL2,glmModLL3,glmModLL4,glmModLL5,glmModLL6, base = TRUE)
summary(glmModLL1)
