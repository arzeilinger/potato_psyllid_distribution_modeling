#### Bayesian list length analysis of occupancy of pest species other than potato psyllids
#### Inspired by Isaac et al. 2014
#### Specifying spatial random effect as in Van Strien et al. 2015 (I think) 

rm(list = ls())
#### Preliminaries
my_packages<-c('data.table', 'tidyr', 'lattice', 'dplyr', 'ggplot2')
lapply(my_packages, require, character.only=T)

## load functions
source("R_functions/museum_specimen_analysis_functions.R")

# Load species lists data set with climate data 
longLists <- readRDS("output/Hemip_Long_Lists_Climate_15km_Cells_2016-06-14.rds")

# Make into dataframe
longListsDF <- longLists %>% rbindlist() %>% as.data.frame()
longListsDF$Species[longListsDF$Species == "Myzus.(Nectarosiphon).persicae"] <- "Myzus.persicae"

speciesFreq <- as.data.frame(table(longListsDF$Species))
speciesFreq <- speciesFreq[order(speciesFreq$Freq),]
tail(speciesFreq)
# Lygus hesperus and Myzus persicae are pests best represented in data set

longLists <- make_lists(longListsDF, min.list.length = 3)

######################################################################
#### Analysis of Lygus hesperus occupancy
######################################################################

# Collectors of potato psyllids, from RawRecords data set in making_species_lists
lygusCollectors <- readRDS("output/lygus_hesperus_collectors.rds")
# select only lists that contain collectors of potato psyllids
lygusLists <- onlyCollectors(lygusCollectors)

# Transform to data frame with pp detection
detectData <- detectDataFunc(lygusLists, focalSpecies = "Lygus.hesperus") 
detectData <- dplyr::filter(detectData, !is.na(aet) & !is.na(cwd) & !is.na(tmn) & !is.na(tmx))
detectData$lnlist_length <- log(detectData$list_length)
str(detectData)
table(detectData$detection)

# standardize numeric covariates, include as new variables in data frame
covars <- c("year", "month", "lnlist_length", "aet", "cwd", "tmn", "tmx")
covars.i <- as.numeric(sapply(covars, function(x) which(names(detectData) == x), simplify = TRUE))
for(i in covars.i){
  var.i <- names(detectData)[i]
  stdname.i <- paste("std", var.i, sep = "")
  stdvar.i <- standardize(detectData[,var.i])
  detectData[,stdname.i] <- stdvar.i
}

# Additional covariates for quadratic effects and interactions
detectData$stdmonth2 <- detectData$stdmonth^2
detectData$stdllyr <- detectData$stdlnlist_length*detectData$stdyear
str(detectData)

# Save detectData
saveRDS(detectData, file = "output/lygus_detection_dataset.rds")
write.csv(detectData, file = "output/lygus_detection_dataset.csv", row.names = FALSE)


#### Make list of "flat" data vectors for NIMBLE model
# make indices
N <- nrow(detectData)
siteID <- detectData$cellID %>% factor(., levels = unique(.)) %>% as.numeric()
nsite <- detectData$cellID %>% unique() %>% length()

nimbleData <- with(detectData, 
                   list(N = N,
                        nsite = nsite,
                        aet = stdaet,
                        tmn = stdtmn,
                        tmx = stdtmx,
                        year = stdyear,
                        month = stdmonth,
                        month2 = stdmonth2,
                        list_length = stdlnlist_length,
                        year_list_length = stdllyr,
                        y = detection,
                        siteID = siteID))
saveRDS(nimbleData, file = "output/lygus_data_nimble_occupancy.rds")


############################################################################################
#### Figures for lists

#### Histogram of list length
# Just potato psyllid occurrences
# lhData <- detectData[detectData$detection == 1,]
# 
# list_length_histogram <- ggplot(detectData,aes(x=list_length)) + 
#   geom_histogram(fill = "darkgrey", alpha = 1, binwidth = 1) +
#   geom_histogram(data=subset(detectData,detection == 1),fill = "black", alpha = 1, binwidth = 1) +
#   xlab("List length") + ylab("Frequency") + 
#   theme_bw() + 
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         #panel.border = element_blank(),
#         panel.background = element_blank()) 
# list_length_histogram

#ggsave(filename = "results/figures/list_length_histogram.tiff", plot = list_length_histogram)


#############################################################################################
#### Occupancy analysis of lygus hesperus
############################################################################################
my.packages <- c("nimble", "coda", "lattice", "akima")
lapply(my.packages, require, character.only = TRUE)

#### Load data set for model fitting
inputData <- readRDS('output/lygus_data_nimble_occupancy.rds')
source('R_functions/nimble_definitions.R')

# #### Increasing memory limit for R
# memory.limit()
# memory.limit(size = 7000) # Size in Mb
# memory.limit()


#####################################################
#### Define model in BUGS/NIMBLE language

code <- nimbleCode({
    mu_alpha ~ dnorm(0, 0.001)
    sigma_alpha ~ dunif(0, 1000)
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu_alpha, sd = sigma_alpha)  ## site random effect
    }
    for(i in 1:9) {
        beta[i] ~ dnorm(0, 0.001)
    }
    # for(i in 1:4) {
    #     betaseason[i] ~ dnorm(0, 0.001)    ## new fixed effects for each season
    #     #betaseasonyear[i] ~ dnorm(0, 0.001)
    # }
    for(i in 1:N) {
        logit(p_occ[i]) <- alpha[siteID[i]] + beta[4]*aet[i] + beta[5]*tmn[i] + beta[6]*tmx[i] + beta[7]*year[i] + beta[8]*month[i] + beta[9]*month2[i]
        logit(p_obs[i]) <- beta[1] + beta[2]*list_length[i] + beta[3]*year_list_length[i]
        y[i] ~ dOccupancy(p_occ[i], p_obs[i])
    }
})

constants <- with(inputData,
                  list(N=N, nsite=nsite, 
                       aet=aet, tmn=tmn, tmx=tmx, 
                       year=year, 
                       month=month,
                       month2=month2,
                       list_length=list_length, 
                       year_list_length=year_list_length, 
                       siteID=siteID))

data <- with(inputData, list(y=y))

inits <- list(mu_alpha=0, sigma_alpha=1, alpha=rep(0,inputData$nsite), beta=rep(0,9))#, betaseason=rep(0,4))

modelInfo_month <- list(code=code, constants=constants, data=data, inits=inits, name='month_model')


#### Set up model and samplers
Rmodel <- nimbleModel(modelInfo_month$code,
                      modelInfo_month$constants,
                      modelInfo_month$data,
                      modelInfo_month$inits)

Cmodel <- compileNimble(Rmodel)

spec <- configureMCMC(Rmodel)

#### Best configuration of samplers for random effect occupancy model
spec$removeSamplers('beta[1:9]')
spec$addSampler('beta[1:3]', 'RW_block') # detection sub-model sampler
spec$addSampler('beta[4:9]', 'RW_block') # occupancy sub-model sampler
spec$removeSamplers('sigma_alpha')
spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha')) # random effect sampler
spec$getSamplers() # Check samplers
spec$addMonitors('p_occ') # add a monitor to get p_occ in output

#### Compile MCMC in R and C++
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


#### Run MCMC with 150,000 iterations and 50,000 burn-in
niter <- 150000
burnin <- 50000

samplesList <- lapply(1:3, mcmcClusterFunction)

save(samplesList, file = 'output/lygus_MCMC_list.RData')

#### Loading saved MCMC run, sames as list, "samplesList"
load(file = 'output/lygus_MCMC_list.RData')


#######################################################################################
#### Analyzing results
###############################################################################################
#### Posterior Inferences
#### Mean and 95% Credible intervals
results <- as.data.frame(cbind(apply(samplesList[[1]], 2, mean),
                               apply(samplesList[[1]], 2, function(x) quantile(x, 0.025)),
                               apply(samplesList[[1]], 2, function(x) quantile(x, 0.975))))
names(results) <- c("mean", "cil", "ciu")
results$params <- row.names(results)
resultsPars <- results[-grep("p_occ", results$params), c("mean", "cil", "ciu", "params")]
resultsPars$covar <- c("det_intercept", "list_length", "year_list_length", "aet", "tmn", "tmx", "year", "month", "month2", NA, NA)
lygusResults <- resultsPars

# Make results table for ms
lygusTable <- lygusResults[, c("mean", "cil", "ciu")] %>% round(., digits = 2) %>%
  cbind(., lygusResults$covar)
lygusTable$summary <- with(lygusTable, paste(mean, " [", cil, ", ", ciu, "]", sep = ""))
write.csv(lygusTable, file = "results/lygus_occupancy_results_for_ms.csv", row.names = TRUE)
lygusTable



######################################################################
#### Analysis of Myzus persicae occupancy
######################################################################

# Collectors of potato psyllids, from RawRecords data set in making_species_lists
myzusCollectors <- readRDS("output/myzus_persicae_collectors.rds")
# select only lists that contain collectors of potato psyllids
myzusLists <- onlyCollectors(myzusCollectors)

# Transform to data frame with pp detection
detectData <- detectDataFunc(myzusLists, focalSpecies = "Myzus.persicae") 
detectData <- dplyr::filter(detectData, !is.na(aet) & !is.na(cwd) & !is.na(tmn) & !is.na(tmx))
detectData$lnlist_length <- log(detectData$list_length)
str(detectData)
table(detectData$detection)

# standardize numeric covariates, include as new variables in data frame
covars <- c("year", "month", "lnlist_length", "aet", "cwd", "tmn", "tmx")
covars.i <- as.numeric(sapply(covars, function(x) which(names(detectData) == x), simplify = TRUE))
for(i in covars.i){
  var.i <- names(detectData)[i]
  stdname.i <- paste("std", var.i, sep = "")
  stdvar.i <- standardize(detectData[,var.i])
  detectData[,stdname.i] <- stdvar.i
}

# Additional covariates for quadratic effects and interactions
detectData$stdmonth2 <- detectData$stdmonth^2
detectData$stdllyr <- detectData$stdlnlist_length*detectData$stdyear
str(detectData)

# Save detectData
saveRDS(detectData, file = "output/myzus_detection_dataset.rds")
write.csv(detectData, file = "output/myzus_detection_dataset.csv", row.names = FALSE)


#### Make list of "flat" data vectors for NIMBLE model
# make indices
N <- nrow(detectData)
siteID <- detectData$cellID %>% factor(., levels = unique(.)) %>% as.numeric()
nsite <- detectData$cellID %>% unique() %>% length()

nimbleData <- with(detectData, 
                   list(N = N,
                        nsite = nsite,
                        aet = stdaet,
                        tmn = stdtmn,
                        tmx = stdtmx,
                        year = stdyear,
                        month = stdmonth,
                        month2 = stdmonth2,
                        list_length = stdlnlist_length,
                        year_list_length = stdllyr,
                        y = detection,
                        siteID = siteID))
saveRDS(nimbleData, file = "output/myzus_data_nimble_occupancy.rds")


#############################################################################################
#### Occupancy analysis of Myzus persicae
############################################################################################
my.packages <- c("nimble", "coda", "lattice", "akima")
lapply(my.packages, require, character.only = TRUE)

#### Load data set for model fitting
inputData <- readRDS('output/myzus_data_nimble_occupancy.rds')
source('R_functions/nimble_definitions.R')

# #### Increasing memory limit for R
# memory.limit()
# memory.limit(size = 7000) # Size in Mb
# memory.limit()


#####################################################
#### Define model in BUGS/NIMBLE language

code <- nimbleCode({
    mu_alpha ~ dnorm(0, 0.001)
    sigma_alpha ~ dunif(0, 1000)
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu_alpha, sd = sigma_alpha)  ## site random effect
    }
    for(i in 1:9) {
        beta[i] ~ dnorm(0, 0.001)
    }
    # for(i in 1:4) {
    #     betaseason[i] ~ dnorm(0, 0.001)    ## new fixed effects for each season
    #     #betaseasonyear[i] ~ dnorm(0, 0.001)
    # }
    for(i in 1:N) {
        logit(p_occ[i]) <- alpha[siteID[i]] + beta[4]*aet[i] + beta[5]*tmn[i] + beta[6]*tmx[i] + beta[7]*year[i] + beta[8]*month[i] + beta[9]*month2[i]
        logit(p_obs[i]) <- beta[1] + beta[2]*list_length[i] + beta[3]*year_list_length[i]
        y[i] ~ dOccupancy(p_occ[i], p_obs[i])
    }
})

constants <- with(inputData,
                  list(N=N, nsite=nsite, 
                       aet=aet, tmn=tmn, tmx=tmx, 
                       year=year, 
                       month=month,
                       month2=month2,
                       list_length=list_length, 
                       year_list_length=year_list_length, 
                       siteID=siteID))

data <- with(inputData, list(y=y))

inits <- list(mu_alpha=0, sigma_alpha=1, alpha=rep(0,inputData$nsite), beta=rep(0,9))#, betaseason=rep(0,4))

modelInfo_month <- list(code=code, constants=constants, data=data, inits=inits, name='month_model')


#### Set up model and samplers
Rmodel <- nimbleModel(modelInfo_month$code,
                      modelInfo_month$constants,
                      modelInfo_month$data,
                      modelInfo_month$inits)

Cmodel <- compileNimble(Rmodel)

spec <- configureMCMC(Rmodel)

#### Best configuration of samplers for random effect occupancy model
spec$removeSamplers('beta[1:9]')
spec$addSampler('beta[1:3]', 'RW_block') # detection sub-model sampler
spec$addSampler('beta[4:9]', 'RW_block') # occupancy sub-model sampler
spec$removeSamplers('sigma_alpha')
spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha')) # random effect sampler
spec$getSamplers() # Check samplers
spec$addMonitors('p_occ') # add a monitor to get p_occ in output

#### Compile MCMC in R and C++
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


#### Run MCMC with 150,000 iterations and 50,000 burn-in
niter <- 150000
burnin <- 50000

samplesList <- lapply(1:3, mcmcClusterFunction)

save(samplesList, file = 'output/myzus_MCMC_list.RData')

#### Loading saved MCMC run, sames as list, "samplesList"
load(file = 'output/myzus_MCMC_list.RData')

#######################################################################################
#### Analyzing results
###############################################################################################
#### Posterior Inferences
#### Mean and 95% Credible intervals
results <- as.data.frame(cbind(apply(samplesList[[1]], 2, mean),
                               apply(samplesList[[1]], 2, function(x) quantile(x, 0.025)),
                               apply(samplesList[[1]], 2, function(x) quantile(x, 0.975))))
names(results) <- c("mean", "cil", "ciu")
results$params <- row.names(results)
resultsPars <- results[-grep("p_occ", results$params), c("mean", "cil", "ciu", "params")]
resultsPars$covar <- c("det_intercept", "list_length", "year_list_length", "aet", "tmn", "tmx", "year", "month", "month2", NA, NA)
myzusResults <- resultsPars

# Make results table for ms
myzusTable <- myzusResults[, c("mean", "cil", "ciu")] %>% round(., digits = 2) %>%
  cbind(., myzusResults$covar)
myzusTable$summary <- with(myzusTable, paste(mean, " [", cil, ", ", ciu, "]", sep = ""))

write.csv(myzusTable, file = "results/myzus_occupancy_results_for_ms.csv", row.names = TRUE)

lygusTable
myzusTable
