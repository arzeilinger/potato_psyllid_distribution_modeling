#### GLMM model with list length analysis and spatial random effect
#### Inspired by Isaac et al. 2014 and using Gio's code

rm(list = ls())
#### Preliminaries
my_packages<-c('lme4', 'lmerTest', 'data.table', 'tidyr', 'lattice', 'dplyr', 'bbmle', 'optimx')
lapply(my_packages, require, character.only=T)

## Load functions
source("R_functions/museum_specimen_analysis_functions.R")

## Load detection data set, set up from making_jags_glmm_data.R script
detectData <- readRDS("output/potato_psyllid_detection_dataset.rds")

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
glmModFull <- glmer(detection ~ stdlnlist_length + stdyear + stdlnlist_length*stdyear + 
                      stdmonth + I(stdmonth^2) + 
                    stdaet + stdtmn + stdtmx + (1|cellID), 
                  family = "binomial", data = detectData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
summary(glmModFull)
# Model without month
glmModLL2 <- glmer(detection ~ stdlnlist_length + stdyear + stdlnlist_length*stdyear + 
                    stdaet + stdtmn + stdtmx + (1|cellID), 
                  family = "binomial", data = detectData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
# Model with only tmn and tmx
glmModLL3 <- glmer(detection ~ stdlnlist_length + stdyear + stdlnlist_length*stdyear + 
                    stdtmn + stdtmx + (1|cellID), 
                  family = "binomial", data = detectData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
# Model with only climate variables
glmModLL4 <- glmer(detection ~ stdlnlist_length +
                    stdaet + stdtmn + stdtmx + (1|cellID), 
                  family = "binomial", data = detectData,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
# Model with only year and LL
glmModLL5 <- glmer(detection ~ stdlnlist_length + stdyear + stdlnlist_length*stdyear + 
                    (1|cellID), 
                  family = "binomial", data = detectData,
                  control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  #control = glmerControl(optimizer = "bobyqa"))

AICtab(glmModFull, glmModLL2, glmModLL3, glmModLL4, glmModLL5, base = TRUE)
# glmModLL5 has lowest AIC but doesn't converge well
# glmModFull performs nearly as well as converges successfully
summary(glmModFull)


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


########################################################################################################################
#### Model predictions
# Using one subsample iteration of best model
bestModel <- glmModFull
detectData$predOcc <- predict(bestModel, type = "response", re.form = NA)

plot(x = detectData$year, y = detectData$predOcc)
plot(x = detectData$lnlist_length, y = detectData$predOcc)
plot(x = detectData$tmn, y = detectData$predOcc)
plot(x = detectData$aet, y = detectData$predOcc)
plot(x = detectData$month, y = detectData$predOcc)
