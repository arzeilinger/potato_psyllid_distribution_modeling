#### Static occupancy modeling of potato psyllid occurrences
#### Using 2 composite seasons as "seasons" and months as repeat surveys 
#### Following Van Strien et al. 2015 

#### Preliminaries
my_packages<-c('snow', 'rjags', 'R2jags', 'dclone', 
               'data.table', 'tidyr', 'lattice')
lapply(my_packages, require, character.only=T)

# Load California polygon, object named CAsp
#load("CA_SpatialPolygons.Rdata")
# Load Occupancy functions
source("C:/Users/Adam/Documents/UC Berkeley post doc/BIGCB/Pest Project/Potato psyllid/Potato psyllid analyses/Seasonal climate occupancy model/Occupancy_functions_and_polygons.R")


#### Load and structure species lists data set with climate data 
AllData <- readRDS("C:/Users/Adam/Documents/UC Berkeley post doc/BIGCB/Pest Project/Potato psyllid/Potato psyllid data/SA_spp_lists_climate_data_2015-10-27.rds")
AllData$MonthCollected <- as.numeric(levels(factor(AllData$MonthCollected)))[factor(AllData$MonthCollected)]
AllLists <- split(AllData,AllData$collectionID)
PALists <- make_lists(AllData, 4) # "Presence-absence" lists

# Exploration for dynamic occupancy model
# Split data by collecting event (year+cell)
# repeatVisits <- split(AllData,AllData$eventID)
# numVisits <- as.numeric(unlist(lapply(1:length(repeatVisits), function(x) length(unique(repeatVisits[[x]]$ymo)))))

#####################################################################################################
#### Make data set with year, cellID, season, list length, and pp occurrence

detectDataFunc <- function(spLists){
  detectData <- data.frame(year = as.numeric(unlist(lapply(strsplit(names(spLists), ", "), function(x) x[1]))), 
                           cellID = unlist(lapply(strsplit(names(spLists), ", "), function(x) x[2])),
                           month = unlist(lapply(strsplit(names(spLists), ", "), function(x) x[3])),
                           season = unlist(lapply(1:length(spLists), function(x) unique(spLists[[x]]$Season))),
                           eventID = unlist(lapply(1:length(spLists), function(x) unique(spLists[[x]][,"eventID"]))),
                           list_length = as.numeric(unlist(lapply(spLists, nrow))),
                           cwd = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$cwd, na.rm = TRUE))),
                           aet = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$aet, na.rm = TRUE))),
                           tmn = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$tmn, na.rm = TRUE))),
                           tmx = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$tmx, na.rm = TRUE))),
                           detection = unlist(lapply(spLists, function(x) as.numeric(any(x$Species == "Bactericera.cockerelli")))))
  detectData$month <- as.numeric(levels(detectData$month))[detectData$month]
  detectData$cellID <- as.numeric(levels(detectData$cellID))[detectData$cellID]
  return(detectData)
}
#detectPAData <- detectDataFunc(PALists) # Lists with length > 3 
detectData <- detectDataFunc(AllLists) # All lists

## Make new composite season variable
detectData$newSeason <- ifelse(detectData$season == "winter" | detectData$season == "spring", 
                               "winterspring", "summerfall")
detectData$season.year <- paste(detectData$newSeason, detectData$year, sep="-")


#### Making new month numbers
## Note: When I use intra-year period as "season" in occ model, 
## need to repeat month numbers within each season.
## This makes reshaping occurrence data possible
## Months within "season" are visits
old.month <- 1:12
new.month <- rep(c(2,3,4,5,6,1), 2)
detectData$newMonth <- -99
for(i in old.month){
  set.i <- which(detectData$month == i)
  detectData$newMonth[set.i] <- new.month[i]
}


#############################################################################################
#### Make data sets for JAGS model
#############################################################################################
#### Reshape occurrence data for occupancy model
# Requires assumption of closure of sites
occData <- spread(detectData[,c("season.year", "cellID", "newMonth", "detection")], 
                  key = newMonth, value = detection)
# Get dynamic occ model indices-- year(t), site(i), month(j)
season.years <- unique(occData$season.year)
season.years <- season.years[order(season.years)]
nseason.year <- length(season.years)
nrep <- length(unique(detectData$newMonth))
sites <- unique(occData$cellID)
sites <- sites[order(sites)]
nsite <- length(sites)

#### 3-dimensional arrays for detection submodel
# Generate 3-dimensional array with site(1), month(2), and year(3) in the dimensions
# Following structure of Royle and Dorazio 2008 Chap 9
pparray <- array(NA, dim = c(nsite, nrep, nseason.year))
for(i in 1:nseason.year){
  syr.i <- unique(detectData$season.year)[i]
  occData.i <- occData[occData$season.year == syr.i,]
  for(j in unique(occData.i$cellID)){
    siteIndex.j <- which(sites == j)
    pparray[siteIndex.j,,i] <- as.numeric(occData.i[occData.i$cellID == j,3:ncol(occData)])
  }
}


#### Construct array for list length data
llData <- spread(detectData[,c("season.year", "cellID", "newMonth", "list_length")], 
                 key = newMonth, value = list_length)
# Generate 3-dimensional array with site(1), month(2), and year(3) in the dimensions
# Following structure of Royle and Dorazio 2008 Chap 9
llarray <- array(0, dim = c(nsite, nrep, nseason.year))
for(i in 1:nseason.year){
  yr.i <- unique(detectData$season.year)[i]
  llData.i <- llData[llData$season.year == yr.i,]
  for(j in unique(llData.i$cellID)){
    siteIndex.j <- which(sites == j)
    llarray[siteIndex.j,,i] <- as.numeric(llData.i[llData.i$cellID == j,3:ncol(occData)])
  }
}
llarray[is.na(llarray)] <- 0

#### Month array can be filled in with all month numbers, even though no occurrence data is associated with it
marray <- array(NA, dim = c(nsite, nrep, nseason.year))
months1 <- matrix(seq(1,6,by=1), nrow = nsite, ncol = nrep, byrow = TRUE) 
for(i in grep("winterspring", season.years)){
  season.year.i <- grep("winterspring", season.years)[i]
  marray[,,i] <- months1
}
months2 <- matrix(seq(7,12,by=1), nrow = nsite, ncol = nrep, byrow = TRUE) 
for(i in grep("summerfall", season.years)){
  season.year.i <- grep("summerfall", season.years)[i]
  marray[,,i] <- months2
}

#### 2-dimensional dataframes for ecological submodel
#### Construct data.frames for each climate variable
tmnData <- makeDataEcoModel("tmn", detectData)
aetData <- makeDataEcoModel("aet", detectData)
cwdData <- makeDataEcoModel("cwd", detectData)
tmxData <- makeDataEcoModel("tmx", detectData)

# Generate data for derived occupancy regression
# From Van Strien et al. 2015
facsyear <- as.numeric(factor(season.years))
sumX <- sum(facsyear)
sumX2 <- sum(facsyear*facsyear)

### Save data sets for dynamic occupancy model
# Standardize covariates
llarray <- standardize(llarray)
marray <- standardize(marray)
tmnData <- standardize(tmnData)
aetData <- standardize(aetData)
cwdData <- standardize(cwdData)
tmxData <- standardize(tmxData)
occDataList <- list(y = pparray, 
                       month = marray,
                       all_lists = llarray,
                       tmn = tmnData,
                       aet = aetData,
                       cwd = cwdData,
                       tmx = tmxData,
                       nsite = dim(pparray)[1],
                       nrep = dim(pparray)[2],
                       nseason.year = dim(pparray)[3],
                       sumX = sumX,
                       sumX2 = sumX2)

saveRDS(occDataList, file = "seasonal_climate_occupancy_model_data.rds")
saveRDS(occDataList, file = "Cluster runs/seasonal_climate_occupancy_model_data.rds")
saveRDS(occDataList, file = "C:/Users/Adam/Documents/GitHub/potato_psyllid_distribution_modeling/seasonal_climate_occupancy_model_data.rds")

