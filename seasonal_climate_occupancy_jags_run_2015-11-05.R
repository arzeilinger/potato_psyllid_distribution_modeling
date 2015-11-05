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
#### Construct data.frame with annual minimum temperature (tmn)
tmnData <- aggregate(detectData$tmn, by = list(detectData$season.year, detectData$cellID), mean, na.rm = TRUE)
names(tmnData) <- c("season.year","cellID","tmn")
tmnData <- as.matrix(spread(tmnData, key = season.year, value = tmn)[,-1])
tmnData[is.nan(tmnData)] <- NA

#### Construct data.frame with mean seasonal aet
aetData <- aggregate(detectData$aet, by = list(detectData$season.year, detectData$cellID), mean, na.rm = TRUE)
names(aetData) <- c("season.year","cellID","aet")
aetData <- as.matrix(spread(aetData, key = season.year, value = aet)[,-1])
aetData[is.nan(aetData)] <- NA

#### Construct data.frame with mean seasonal cwd
cwdData <- aggregate(detectData$cwd, by = list(detectData$season.year, detectData$cellID), mean, na.rm = TRUE)
names(cwdData) <- c("season.year","cellID","cwd")
cwdData <- as.matrix(spread(cwdData, key = season.year, value = cwd)[,-1])
cwdData[is.nan(cwdData)] <- NA

#### Construct data.frame with mean seasonal tmx
cwdData <- aggregate(detectData$cwd, by = list(detectData$season.year, detectData$cellID), mean, na.rm = TRUE)
names(cwdData) <- c("season.year","cellID","cwd")
cwdData <- as.matrix(spread(cwdData, key = season.year, value = cwd)[,-1])
cwdData[is.nan(cwdData)] <- NA


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
occDataList <- list(y = pparray, 
                       month = marray,
                       all_lists = llarray,
                       tmn = tmnData,
                       aet = aetData,
                       cwd = cwdData,
                       nsite = dim(pparray)[1],
                       nrep = dim(pparray)[2],
                       nseason.year = dim(pparray)[3],
                       sumX = sumX,
                       sumX2 = sumX2)

saveRDS(occDataList, file = "seasonal_climate_occupancy_model_data.rds")
saveRDS(occDataList, file = "Cluster runs/seasonal_climate_occupancy_model_data.rds")


##################################################################################
#### Dynaimc occupancy model with imperfect detection
#### Based on Royle and Dorazio 2008, Panel 9.2 and Van Strien et al. 2013

# setwd("C:/Users/Adam/Documents/UC Berkeley post doc/BIGCB/Pest Project/Potato psyllid/Potato psyllid analyses")

sink("seasonal_climate_occ_model.jags")
cat("
model {

## Priors
# Ecological submodel priors
for(t in 1:nseason.year){
  a[t]~dnorm(0,0.01)
}
b1~dnorm(0,0.01)
b2~dnorm(0,0.01)
b3~dnorm(0,0.01)
b4~dnorm(0,0.01)
# Random effect for site
for(i in 1:nsite){
  eta[i]~dnorm(0,tau)
}
tau <- 1/(sigma*sigma)
sigma~dunif(0,5)
# Observation sub-model priors
for(t in 1:nseason.year){
  alpha.p[t]~dnorm(0,0.01)
}
beta1.p~dnorm(0,0.01)
beta2.p~dnorm(0,0.01)
delta1.p~dnorm(0,0.01)
delta2.p~dnorm(0,0.01)
# Prior for missing values of tmn, aet, and cwd
for(t in 1:nseason.year){
  for(i in 1:nsite){
    tmn[i,t]~dnorm(0,0.01)
    aet[i,t]~dnorm(0,0.01)
    cwd[i,t]~dnorm(0,0.01)
  }
}

# Ecological process sub-model
for(i in 1:nsite){
  #z[i,1]~dbern(psi)
  for(t in 1:nseason.year){
     logit(muZ[i,t]) <- a[t] + b1*tmn[i,t] + b2*pow(tmn[i,t],2) + 
                        b3*aet[i,t] + b4*cwd[i,t] + eta[i]
     z[i,t]~dbern(muZ[i,t])
 }
}

# Observation sub-model
# Adapted from Kery and Schaub 2012
# Long list lengths separated from short and relate to p as Michaelis-Menten function
# Following Van Strien et al. 2013, 2015
for (t in 1:nseason.year){
  for(i in 1:nsite){
    for(j in 1:nrep){
      y[i,j,t] ~ dbern(mu.p[i,j,t])
      mu.p[i,j,t] <- z[i,t]*p[i,j,t]
      p[i,j,t] <- 1 / (1 + exp(-lp.lim[i,j,t]))
      lp.lim[i,j,t] <- min(999, max(-999, lp[i,j,t]))
      lp[i,j,t] <- alpha.p[t] + beta1.p*month[i,j,t] + beta1.p*pow(month[i,j,t], 2) + (delta1.p*all_lists[i,j,t])/(all_lists[i,j,t] + delta2.p)
    } #j
  } #i
} #t

## Derived parameters
# Finite sample occupancy and mean detection
for(t in 1:nseason.year){
  psi.fs[t]<-sum(z[1:nsite,t])/nsite  
  mean.p[t] <- exp(alpha.p[t]) / (1 + exp(alpha.p[t]))    # Sort of average detection
}
# Overall trend in occupancy psi
sumY <- sum(psi.fs[1:nseason.year])
for(t in 1:nseason.year){
  sumxy[t] <- psi.fs[t]*t
}
sumXY <- sum(sumxy[1:nseason.year])
regres.psi <- (sumXY - ((sumX*sumY)/nseason.year))/(sumX2 - ((sumX*sumX)/nseason.year))

}
",fill=TRUE)
sink()

##################################################################################

## Constructing data list
occDataList <- readRDS("C:/Users/Adam/Documents/UC Berkeley post doc/BIGCB/Pest Project/Potato psyllid/Potato psyllid analyses/Dynamic occupancy model/seasonal_climate_occupancy_model_data.rds")
y <- occDataList$y
nseason.year <- occDataList$nseason.year
nrep <- occDataList$nrep
nsite <- occDataList$nsite
## Initial values
Zst <- sapply(1:dim(y)[3], function(x) rowSums(y[,,x], na.rm = TRUE))
Zst[Zst > 1] <- 1
apst <- rep(runif(1, -3, 3), nseason.year)
azst <- rep(runif(1, -3, 3), nseason.year)
b1zst <- runif(1, -3, 3)
b2zst <- runif(1, -3, 3)
b3zst <- runif(1, -3, 3)
b4zst <- runif(1, -3, 3)
b5zst <- runif(1, -3, 3)
b1pst <- runif(1, -3, 3)
b2pst <- runif(1, -3, 3)
d1pst <- runif(1, -3, 3)
d2pst <- runif(1, -3, 3)
d3pst <- runif(1, -3, 3)
inits <- function() list (z=Zst,alpha.p=apst,a=azst)#,
#       beta1.p=b1pst, beta2.p=b2pst, 
#       b1=b1zst, b2=b2zst, b3=b3zst, b4=b4zst, b5=b5zst,
#       delta1.p=d1pst, delta2.p=d2pst)


## Parameters to monitor
parameters <- c("mean.p", "psi.fs", "regres.psi", 
		"delta1.p", "delta2.p",
		"b1", "b2", "b3", "b4")

# MCMC parameters
ni=40100; nb=1000; nt=10; nc=3
# for jags.parfit(), burn-in iterations = n.adapt + n.update
n.adapt <- 500; n.update <- 500

# Make a SOCK cluster using snow
cl <- makeCluster(3, type = "SOCK")
date()
# Call to jags.parfit
occOutput <- jags.parfit(cl, data = occDataList,
                            params = parameters,
                            model = "seasonal_climate_occ_model.jags",
                            inits = inits,
                            n.adapt = n.adapt, n.update = n.update,
                            n.iter = ni, thin = nt, n.chains = nc)
date()
stopCluster(cl) # Close the cluster
# Compute statistics and save output
saveRDS(occOutput, file = "seasonal_climate_occ_jags_out_full.rds")
occdctab <- dctable(dynoccOutput)
occResults <- data.frame(rbindlist(occdctab))
occResults$param <- names(occdctab)
saveRDS(occResults, file = "seasonal_climate_occ_jags_out_params.rds")
print(occResults)
print("SUCCESS!")
