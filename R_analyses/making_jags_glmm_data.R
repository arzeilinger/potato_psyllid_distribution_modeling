#### Bayesian binomial GLMM for list length analysis of potato psyllid occupancy
#### Inspired by Isaac et al. 2014
#### Specifying spatial random effect as in Van Strien et al. 2015 (I think) 

rm(list = ls())
#### Preliminaries
my_packages<-c('data.table', 'tidyr', 'lattice', 'dplyr')
lapply(my_packages, require, character.only=T)

## load functions
source("R_functions/museum_specimen_analysis_functions.R")

# Load species lists data set with climate data 
AllLists <- readRDS("output/All_Hemip_Lists_Climate_15km_Cells_2016-04-14.rds")

# Collectors of potato psyllids, from RawRecords data set in making_species_lists
ppCollectors <- readRDS("output/potato_psyllid_collectors.rds")

# Keep only long lists (ll >= 3)
longLists <- AllLists %>% rbindlist() %>% as.data.frame() %>% make_lists(., min.list.length = 4)
longListsDF <- longLists %>% rbindlist() %>% as.data.frame()

# select only lists that contain collectors of potato psyllids
ppcCollections <- longListsDF[longListsDF$Collector %in% ppCollectors, "collectionID"]
ppcData <- longListsDF[longListsDF$collectionID %in% ppcCollections,]
ppcLists <- ppcData %>% make_lists(., min.list.length = 4)

# Transform to data frame with pp detection
detectData <- detectDataFunc(ppcLists) 
detectData <- dplyr::filter(detectData, !is.na(aet) & !is.na(cwd) & !is.na(tmn) & !is.na(tmx))
detectData$lnlist_length <- log(detectData$list_length)
str(detectData)

#### Exploring intercorrelations among climate variables
climateVars <- detectData[,c("aet", "cwd", "tmn", "tmx")]
pairs(climateVars)
# CWD is highly correlated with AET, probably should drop CWD


# standardize numeric covariates, include as new variables in data frame
covars <- c("year", "month", "lnlist_length", "aet", "cwd", "tmn", "tmx")
covars.i <- as.numeric(sapply(covars, function(x) which(names(detectData) == x), simplify = TRUE))
for(i in covars.i){
  var.i <- names(detectData)[i]
  stdname.i <- paste("std", var.i, sep = "")
  stdvar.i <- standardize(detectData[,var.i])
  detectData[,stdname.i] <- stdvar.i
}

# Save detectData
saveRDS(detectData, file = "output/potato_psyllid_detection_dataset.rds")
write.csv(detectData, file = "output/potato_psyllid_detection_dataset.csv", row.names = FALSE)

# Just potato psyllid occurrences
ppData <- detectData[detectData$detection == 1,]

#### Make JAGS data set, with lists as rows (i) and sites as columns (j)
detectionMatrix <- makeEcoDataMatrix("detection")
stdcovars <- paste("std", covars, sep = "")
jagsData <- lapply(stdcovars, function(x) makeEcoDataMatrix(x, fill = 0))
names(jagsData) <- stdcovars

jagsGLMMdata <- list(detectionMatrix = detectionMatrix,
                     year = jagsData$stdyear,
                     month = jagsData$stdmonth,
                     list_length = jagsData$stdlnlist_length,
                     aet = jagsData$stdaet,
                     cwd = jagsData$stdcwd,
                     tmn = jagsData$stdtmn,
                     tmx = jagsData$stdtmx,
                     nlist = nrow(detectionMatrix),
                     nsite = ncol(detectionMatrix))
saveRDS(jagsGLMMdata, file = "output/Data_JAGS_GLMM.rds")

# Generate test data, a subset of full data set, to make sure model runs
# testData <- detectData[1:100,]
# detectionMatrix <- makeEcoDataMatrix("detection", testData)
# jagsData <- lapply(stdcovars, function(x) makeEcoDataMatrix(x, testData, fill = 0))
# names(jagsData) <- stdcovars
# jagsTestdata <- list(detectionMatrix = detectionMatrix,
#                      year = jagsData$stdyear,
#                      list_length = jagsData$stdlnlist_length,
#                      aet = jagsData$stdaet,
#                      cwd = jagsData$stdcwd,
#                      tmn = jagsData$stdtmn,
#                      tmx = jagsData$stdtmx,
#                      nlist = nrow(detectionMatrix),
#                      nsite = ncol(detectionMatrix))
# saveRDS(jagsTestdata, file = "C:/Users/Adam/Documents/GitHub/potato_psyllid_distribution_modeling/GLM/Data_JAGS_GLMM_test.rds")
