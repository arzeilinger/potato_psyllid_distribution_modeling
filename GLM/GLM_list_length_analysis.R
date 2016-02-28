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

# Keep only long lists (ll > 3)
longLists <- AllLists %>% rbindlist() %>% as.data.frame() %>% make_lists(., min.list.length = 3)
longListsDF <- longLists %>% rbindlist() %>% as.data.frame()

# select only lists that contain collectors of potato psyllids
ppcCollections <- longListsDF[longListsDF$Collector %in% ppCollectors, "collectionID"]
ppcData <- longListsDF[longListsDF$collectionID %in% ppcCollections,]
ppcLists <- ppcData %>% make_lists(., min.list.length = 3)

# Transform to data frame with pp detection
detectData <- detectDataFunc(ppcLists) 
detectData$lnlist_length <- log(detectData$list_length)

# standardize numeric covariates, include as new variables in data frame
covars <- c("year", "month", "lnlist_length", "aet", "cwd", "tmn", "tmx")
covars.i <- as.numeric(sapply(covars, function(x) which(names(detectData) == x), simplify = TRUE))
for(i in covars.i){
  var.i <- names(detectData)[i]
  stdname.i <- paste("std", var.i, sep = "")
  stdvar.i <- standardize(detectData[,var.i])
  detectData[,stdname.i] <- stdvar.i
}


# Exploring data
table(detectData$detection)
with(detectData, hist(list_length, breaks = seq(min(list_length),max(list_length),1)))
ppData <- detectData[detectData$detection == 1,]
hist(ppData$list_length)

length(unique(detectData$cellID)) # 858 different cells -- too many for the site random effect
detectData %>% group_by(cellID) %>% summarise(nvisits = length(cellID))


#### Subsample nondetections
nabsence <- ppData %>% nrow()*10 # Set number of absences to 10 x number of presences (SDM rule of thumb)
# Subsample nondetections and recombine with detections
subabsenceData <- detectData[detectData$detection == 0,] %>% sample_n(., size = nabsence, replace = FALSE) %>% rbind(.,ppData)


#################################################################################
#### GLMMs
# GLMM with "long" species lists (length > 3)
# Using previously standardized covariates
# Full model
glmModFull <- glmer(detection ~ stdlnlist_length + stdyear + stdmonth + I(stdmonth^2) + 
                    stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                  family = "binomial", data = subabsenceData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
summary(glmModFull)
# Model with season instead of month (season is a factor)
glmModSeason <- glmer(detection ~ stdlnlist_length + season*stdyear + 
                    stdaet + stdcwd + season*stdtmn + season*stdtmx + (1|cellID), 
                  family = "binomial", data = subabsenceData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
# Model without month
glmModLL2 <- glmer(detection ~ stdlnlist_length + stdyear +  
                    stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                  family = "binomial", data = subabsenceData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
# Model with only tmn and tmx
glmModLL3 <- glmer(detection ~ stdlnlist_length + stdyear + 
                    stdtmn + stdtmx + (1|cellID), 
                  family = "binomial", data = subabsenceData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
# Model with only climate variables
glmModLL4 <- glmer(detection ~ stdlnlist_length +  
                    stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                  family = "binomial", data = subabsenceData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
# Model with only year and LL
glmModLL5 <- glmer(detection ~ stdlnlist_length + stdyear +
                    (1|cellID), 
                  family = "binomial", data = subabsenceData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
# Model with year*month interaction
glmModLL6 <- glmer(detection ~ stdlnlist_length + stdyear*stdmonth +
                     standardize(year)*I(standardize(month)^2) +
                    (1|cellID), 
                  family = "binomial", data = subabsenceData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
ms <- AICtab(glmModFull, glmModSeason, glmModLL2, glmModLL3, glmModLL4, glmModLL5, glmModLL6, base = TRUE)
# glmModLL2 (full model without month terms) seems to be best
summary(glmModLL2)



###########################################################################################################
#### Run subsampling of nondetections and model selection many times
# Subsample nondetections
nabsence <- ppData %>% nrow()*10 # Set number of absences to 10 x number of presences (SDM rule of thumb)

sampleAndSelect <- function(x){
  # Subsample nondetections and recombine with detections
  subabsence <- detectData[detectData$detection == 0,] %>% sample_n(., size = x, replace = FALSE) %>% rbind(.,ppData)
  #### Run subabsence data set through models
  # Full model
  glmModFull <- glmer(detection ~ stdlnlist_length + stdyear + stdmonth + I(stdmonth^2) + 
                        stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                      family = "binomial", data = subabsence,
                      #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                      control = glmerControl(optimizer = "bobyqa"))
  # Model without month
  glmMod2 <- glmer(detection ~ stdlnlist_length + stdyear +  
                       stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                     family = "binomial", data = subabsence,
                     #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                     control = glmerControl(optimizer = "bobyqa"))
  # Model with only climate variables
  glmMod3 <- glmer(detection ~ stdlnlist_length +  
                       stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                     family = "binomial", data = subabsence,
                     #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                     control = glmerControl(optimizer = "bobyqa"))
  # Model with only year and LL
  glmMod4 <- glmer(detection ~ stdlnlist_length + stdyear +
                       (1|cellID), 
                     family = "binomial", data = subabsence,
                     #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                     control = glmerControl(optimizer = "bobyqa"))  
  selectModel <- AICtab(glmModFull, glmMod2, glmMod3, glmMod4, base = TRUE)
  bestModel <- attr(selectModel, "row.names")[1]
  return(bestModel)
}

times <- 50
runs <- rep(nabsence, times) 
sampleSelectRuns <- sapply(runs, sampleAndSelect, simplify = TRUE)


#### Get coefficients from the best model, run multiple times
modelCoefFunc <- function(x){
  # Subsample nondetections and recombine with detections
  subabsence <- detectData[detectData$detection == 0,] %>% sample_n(., size = x, replace = FALSE) %>% rbind(.,ppData)
  # Model without month
  glmMod2 <- glmer(detection ~ stdlnlist_length + stdyear +  
                       stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                     family = "binomial", data = subabsence,
                     #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                     control = glmerControl(optimizer = "bobyqa"))
  return(summary(glmMod2))
}
extractModelRuns <- lapply(runs, modelCoefFunc)

#### Model predictions
# Using the full model
predictData <- dplyr::filter(detectData, !is.na(aet) & !is.na(cwd) & !is.na(tmn) & !is.na(tmx))
predictData$predOcc <- predict(glmModLL2, type = "response", re.form = NA)

plot(x = predictData$year, y = predictData$predOcc)
plot(x = predictData$lnlist_length, y = predictData$predOcc)
plot(x = predictData$tmn, y = predictData$predOcc)
plot(x = predictData$aet, y = predictData$predOcc)
