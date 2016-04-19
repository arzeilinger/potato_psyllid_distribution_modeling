#### Single-visit conditional likelihood model
#### From Solymos et al. 2012

rm(list = ls())
#### Preliminaries
my_packages<-c('lme4', 'lmerTest', 'data.table', 'tidyr', 'lattice', 'dplyr', 'bbmle', 
               'detect', 'optimx')
lapply(my_packages, require, character.only=T)

## load functions
source("R_functions/museum_specimen_analysis_functions.R")

# Load species lists data set with climate data 
AllLists <- readRDS("abundance_modeling/All_Hemip_Counts_Climate_15km_Cells_2016-02-26.rds")

# Collectors of potato psyllids, from RawRecords data set in making_species_lists
ppCollectors <- readRDS("potato_psyllid_collectors.rds")
AllListsDF <- AllLists %>% rbindlist() %>% as.data.frame()
# lists that contain collectors of potato psyllids
ppcCollections <- AllListsDF[AllListsDF$Collector %in% ppCollectors, "collectionID"]
ppcData <- AllListsDF[AllListsDF$collectionID %in% ppcCollections,]
ppcLists <- ppcData %>% make_lists(., min.list.length = 1)

# Transform to data frame with pp detection
countData <- countDataFunc(ppcLists, "Bactericera.cockerelli")
# Select only long lists (list_length >= 3)
countData <- countData[countData$list_length > 2,]
countData$lnlist_length <- log(countData$list_length)
countData$lntotal_n <- log(countData$total_n)

# standardize numeric covariates, include as new variables in data frame
covars <- c("year", "month", "lnlist_length", "lntotal_n", "diversity", "aet", "cwd", "tmn", "tmx")
covars.i <- as.numeric(sapply(covars, function(x) which(names(countData) == x), simplify = TRUE))
for(i in covars.i){
  var.i <- names(countData)[i]
  stdname.i <- paste("std", var.i, sep = "")
  stdvar.i <- standardize(countData[,var.i])
  countData[,stdname.i] <- stdvar.i
}

# Just potato psyllid data points
ppData <- countData[countData$count >= 1,]

# Exploring data
table(countData$count)
with(countData, hist(list_length, breaks = seq(min(list_length),max(list_length),1)))
hist(ppData$list_length)
hist(ppData$count)
length(unique(countData$cellID)) # 858 different cells -- too many for the site random effect
cellTable <- countData %>% group_by(cellID) %>% summarise(nvisits = length(cellID))


# Making data set for HM seminar GLMM group
countData <- countData[,c("year", "cellID", "month", "list_length", "aet", "cwd", "tmn", "tmx",
                          "stdlnlist_length", "stdaet", "stdcwd", "stdtmn", "stdtmx", "count")]
metaData <- c("count = count of potato psyllid museum specimens collected in a given collecting event.",
              "year = year collected. cellID = spatial cell. month = month collected.",
              "list_length = total number of insect species collected.",
              "Climate variables: aet = actual evapotransporation; cwd = climate water deficit;", 
              "tmn = minimum annual temperature; tmx = maximum annual temperature.",
              "Variables with prefix 'std' are standardized based on mean and SD.",
              "Data compiled by Adam Zeilinger 2016-04-02.")

countData2 <- list(metaData, countData)
saveRDS(countData2, file = "potato_psyllid_count_data_for_GLMM.rds")

###########################################################################################################
#### Subsampling nondetections and running model selection many times
# Remove NAs from covariates
countData <- dplyr::filter(countData, !is.na(aet) & !is.na(cwd) & !is.na(tmn) & !is.na(tmx))

# Number of nondetections
nabsence <- ppData %>% nrow()*10 # Set number of absences to 10 x number of presences (Gio's SDM rule of thumb)
# Subsample nondetections
subabsence <- countData[countData$count == 0,] %>% sample_n(., size = nabsence, replace = FALSE) %>% rbind(.,ppData)

#################################################################################
#### GLMM
# GLMM with "long" species lists (length > 2)
# Full model
glmModFull <- glmer(count ~ stdlnlist_length + stdyear + stdmonth + I(stdmonth^2) + 
                    stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                  family = "poisson", data = subabsence,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
summary(glmModFull)
# List_length Model without month
glmModLL2 <- glmer(count ~ stdlnlist_length + stdyear +  
                    stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                  family = "poisson", data = subabsence,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
summary(glmModLL2)

# Diversity model
glmDiversity <- glmer(count ~ stddiversity + stdyear +  
                    stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                  family = "poisson", data = subabsence,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))

# Total N model
glmN <- glmer(count ~ stdlntotal_n + stdyear + #stdmonth + I(stdmonth^2) +
                    stdaet + stdcwd + stdtmn + stdtmx + (1|cellID), 
                  family = "poisson", data = subabsence,
                  #control = glmerControl(optCtrl = list(method = "spg", ftol = 1e-20, maxit = 100000)))
                  control = glmerControl(optimizer = "bobyqa"))
AICtab(glmModLL2, glmDiversity, glmN, base = TRUE)
summary(glmN)

#######################################################################################
#### Solymos' conditional likelihood model
#######################################################################################

clLLMod <- svabu(count ~ stdyear + stdaet + stdcwd + stdtmn + stdtmx | 
                 stdlnlist_length + stdmonth + I(stdmonth^2), 
               data = countData)
clDiversityMod <- svabu(count ~ stdyear + stdaet + stdcwd + stdtmn + stdtmx | 
                 stddiversity + stdmonth + I(stdmonth^2), 
               data = countData)
clnMod <- svabu(count ~ stdyear + stdaet + stdcwd + stdtmn + stdtmx | 
                 stdtotal_n + stdmonth + I(stdmonth^2), 
               data = countData)
AICtab(clLLMod, clDiversityMod, clnMod, base = TRUE)
summary(clnMod)
summary(clDiversityMod)
