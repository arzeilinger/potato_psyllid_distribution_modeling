#### GLMM model with list length analysis and spatial random effect
#### Inspired by Isaac et al. 2014 and using Gio's code

rm(list = ls())
#### Preliminaries
my_packages<-c('lme4', 'lmerTest', 'data.table', 'tidyr', 'lattice', 'dplyr', 'bbmle', 'optimx')
lapply(my_packages, require, character.only=T)

## Load functions
source("R_functions/museum_specimen_analysis_functions.R")

# Load species lists data set with climate data 
AllLists <- readRDS("Potato psyllid data/All_Hemip_Lists_Climate_15km_Cells_2016-03-25.rds")

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
detectData <- dplyr::filter(detectData, !is.na(aet) & !is.na(cwd) & !is.na(tmn) & !is.na(tmx))
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


#################################################################################
#### GLMMs
## GLMM with "long" species lists (length > 3)
## Using previously standardized covariates and full data set

# Full model
glmModFull <- glmer(detection ~ stdlnlist_length + stdyear + stdmonth + I(stdmonth^2) + 
                    stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                  family = "binomial", data = detectData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
summary(glmModFull)
# Model with season instead of month (season is a factor)
glmModSeason <- glmer(detection ~ stdlnlist_length + season*stdyear + 
                    stdaet + stdcwd + season*stdtmn + season*stdtmx + (1|cellID), 
                  family = "binomial", data = detectData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
# Model without month
glmModLL2 <- glmer(detection ~ stdlnlist_length + stdyear +  
                    stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                  family = "binomial", data = detectData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
# Model with only tmn and tmx
glmModLL3 <- glmer(detection ~ stdlnlist_length + stdyear + 
                    stdtmn + stdtmx + (1|cellID), 
                  family = "binomial", data = detectData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
# Model with only climate variables
glmModLL4 <- glmer(detection ~ stdlnlist_length +  
                    stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                  family = "binomial", data = detectData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
# Model with only year and LL
glmModLL5 <- glmer(detection ~ stdlnlist_length + stdyear +
                    (1|cellID), 
                  family = "binomial", data = detectData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
# Model with year*month interaction
glmModLL6 <- glmer(detection ~ stdlnlist_length + stdyear*stdmonth +
                     standardize(year)*I(standardize(month)^2) +
                    (1|cellID), 
                  family = "binomial", data = detectData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
AICtab(glmModFull, glmModSeason, glmModLL2, glmModLL3, glmModLL4, glmModLL5, glmModLL6, base = TRUE)
# glmModLL2 (full model without month terms) seems to be best
summary(glmModLL2)


###########################################################################################################
#### Subsampling nondetections and running model selection many times
# Remove NAs from covariates
detectData <- dplyr::filter(detectData, !is.na(aet) & !is.na(cwd) & !is.na(tmn) & !is.na(tmx))

# Number of nondetections
nabsence <- ppData %>% nrow()*10 # Set number of absences to 10 x number of presences (Gio's SDM rule of thumb)

sampleSelect <- function(x){
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

times <- 10 # Number of times to subsample nondetections and fit to model
runs <- rep(nabsence, times) 
runSampleSelect <- sapply(runs, sampleSelect, simplify = TRUE)
table(runSampleSelect)
# glmMod2 is most often the model selected, but not always


#### Get coefficients from the best model (glmMod2), run multiple times
modelCoefFunc <- function(x){
  # Subsample nondetections and recombine with detections
  subabsence <- detectData[detectData$detection == 0,] %>% sample_n(., size = x, replace = FALSE) %>% rbind(.,ppData)
  # Model without month
  glmMod2 <- glmer(detection ~ stdlnlist_length + stdyear +  
                       stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                     family = "binomial", data = subabsence,
                     #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                     control = glmerControl(optimizer = "bobyqa"))
  return(glmMod2)
}

extractModelRuns <- lapply(runs, modelCoefFunc)
modelSummaries <- lapply(extractModelRuns, summary)


#### Model predictions
# Using one subsample iteration of best model
bestModel <- extractModelRuns[[1]]
detectData$predOcc <- predict(bestModel, type = "response", re.form = NA)

plot(x = detectData$year, y = detectData$predOcc)
plot(x = detectData$lnlist_length, y = detectData$predOcc)
plot(x = detectData$tmn, y = detectData$predOcc)
plot(x = detectData$aet, y = detectData$predOcc)
