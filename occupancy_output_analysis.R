#### Analysis of seasonal climate occupancy model results
#### 2015-11-09

setwd("C:/Users/Adam/Documents/GitHub/potato_psyllid_distribution_modeling")

# Estimated parameters
occResults <- readRDS("seasonal_climate_occ_jags_out_params.rds")
dim(occResults)
head(occResults)
occResults$params <- row.names(occResults)
sum(is.na(occResults$r.hat)) # Multivariate Rhat calculated

# Just means, 95% CI, and Rhat
occSummary <- occResults[,c("mean","X2.5.","X97.5.","r.hat","params")]

# Loading original data
occDataList <- readRDS("seasonal_climate_occupancy_model_data.rds")

# Looking at Rhat
hist(occSummary$r.hat)
unconvParams <- occResults[occResults$r.hat < 0.9 | occResults$r.hat > 1.1, c("r.hat", "params")]

# looking at occupancy results over time
psiOut <- occSummary[grep("psi.fs", occSummary$params),]
psiOut$season.year <- occDataList$season.years
psiOut$year <- as.numeric(unlist(lapply(strsplit(psiOut$season.year, "-"), function(x) x[2])))

## Just winter-spring
ann1occ <- psiOut[grep("winterspring", psiOut$season.year),]
plot(ann1occ$year, ann1occ$mean)
# decline in occupancy in winter-spring over years

## Just summer-fall
ann2occ <- psiOut[grep("summerfall", psiOut$season.year),]
plot(ann2occ$year, ann2occ$mean)
# increase in occupancy in summer-fall over years


#### Detection probability
detectOut <- occSummary[grep("mean.p", occSummary$params),]
detectOut$season.year <- occDataList$season.years
detectOut$year <- as.numeric(unlist(lapply(strsplit(detectOut$season.year, "-"), function(x) x[2])))

## Just winter-spring
ann1occ <- detectOut[grep("winterspring", detectOut$season.year),]
plot(ann1occ$year, ann1occ$mean)
# decline in detection overall in winter-spring over years, but seems like two trends: one high p and one low p

## Just summer-fall
ann2occ <- detectOut[grep("summerfall", detectOut$season.year),]
plot(ann2occ$year, ann2occ$mean)
# increase in detection overall in summer-fall over years


## List length vs detection using logistic model
llv <- seq(0,15,by=1)
a <- occSummary[occSummary$params == "delta1.p","mean"]
b <- occSummary[occSummary$params == "delta2.p","mean"]
predp <- function(a, b, llv = llv){
  return(exp(a + b*llv)/(1 + exp(a + b*llv)))
}
est.pd <- predp(a = a, b = b, llv = llv)
plot(llv, est.pd)


## Exploring parameter values and finding the best (theoretical) parameters
ppx <- predp(a = -6, b = 1, llv)
plot(llv, ppx)
