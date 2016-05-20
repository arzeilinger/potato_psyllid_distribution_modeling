#### NIMBLE occupancy model for ZIB GLMM model of potato psyllid occupancy
#### Originally developed by Daniel Turek

my.packages <- c("nimble", "coda", "lattice", "snow", "snowfall")
lapply(my.packages, require, character.only = TRUE)

#### Load data set for model fitting
inputData <- readRDS('output/data_nimble_zib.rds')
source('R_functions/nimble_definitions.R')

#### Increasing memory limit for R
memory.limit()
memory.limit(size = 7000) # Size in Mb
memory.limit()


#####################################################
#### Define model in BUGS/NIMBLE language

code <- nimbleCode({
    mu_alpha ~ dnorm(0, 0.001)
    sigma_alpha ~ dunif(0, 1000)
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu_alpha, sd = sigma_alpha)  ## site random effect
    }
    for(i in 1:7) {
        beta[i] ~ dnorm(0, 0.001)
    }
    for(i in 1:4) {
        betaseason[i] ~ dnorm(0, 0.001)    ## new fixed effects for each season
        #betaseasonyear[i] ~ dnorm(0, 0.001)
    }
    for(i in 1:N) {
        logit(p_occ[i]) <- alpha[siteID[i]] + beta[4]*aet[i] + beta[5]*tmn[i] + beta[6]*tmx[i] + beta[7]*year[i] + betaseason[season[i]] #+ betaseasonyear[season[i]]*year[i]
        logit(p_obs[i]) <- beta[1] + beta[2]*list_length[i] + beta[3]*year_list_length[i]
        y[i] ~ dOccupancy(p_occ[i], p_obs[i])
    }
})

constants <- with(inputData,
                  list(N=N, nsite=nsite, 
                       aet=aet, tmn=tmn, tmx=tmx, 
                       year=year, 
                       season=season, 
                       list_length=list_length, 
                       year_list_length=year_list_length, 
                       siteID=siteID))

data <- with(inputData, list(y=y))

inits <- list(mu_alpha=0, sigma_alpha=1, alpha=rep(0,inputData$nsite), beta=rep(0,7), betaseason=rep(0,4))

modelInfo_season <- list(code=code, constants=constants, data=data, inits=inits, name='season_model')


#### Set up model and samplers
Rmodel <- nimbleModel(modelInfo_season$code,
                      modelInfo_season$constants,
                      modelInfo_season$data,
                      modelInfo_season$inits)

Cmodel <- compileNimble(Rmodel)

spec <- configureMCMC(Rmodel)

#### Best configuration of samplers for random effect occupancy model
spec$removeSamplers('beta[1:7]')
spec$addSampler('beta[1:3]', 'RW_block') # detection sub-model sampler
spec$addSampler('beta[4:7]', 'RW_block') # occupancy sub-model sampler
spec$removeSamplers('sigma_alpha')
spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha')) # random effect sampler
spec$getSamplers() # Check samplers
spec$addMonitors('p_occ') # add a monitor to get p_occ in output

#### Compile MCMC in R and C++
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


#### Run MCMC with 500,000 iterations and 100,000 burn-in
niter <- 500000
burnin <- 100000

# Function for running MCMC multiple times
mcmcClusterFunction <- function(x){
  set.seed(x)
  Cmcmc$run(niter)
  samples <- as.matrix(Cmcmc$mvSamples)[(burnin+1):niter,]
  return(samples)
}

samplesList <- lapply(1:2, mcmcClusterFunction)



##############################################################################
#### Attempt at cluster
sfInit(parallel = TRUE, cpus = 2)
date()

sfLibrary(nimble)
sfExport("Cmcmc", "Rmodel", "niter", "burnin", "mcmcClusterFunction")
sfExport("spec")
sfExport("Cmodel")
sfExport("Rmcmc")
sfLapply(1:2, mcmcClusterFunction)

date()
sfStop() # Close the cluster
#################################################################################


mcmc1 <- coda::as.mcmc(samplesList[[1]])
mcmc2 <- coda::as.mcmc(samplesList[[2]])
mcmcs <- coda::mcmc.list(mcmc1, mcmc2)

save(samplesList, mcmcs, file = 'output/MCMC_season.RData')

#### Loading saved MCMC run
load(file = 'output/MCMC_season.RData')


#######################################################################
#### Assessing convergence and summarizing/plotting results

## Rhat
coda::gelman.diag(mcmcs, autoburnin = FALSE)
## Effective sample size
apply(samples1, 2, coda::effectiveSize)

## Posterior Density Plots
plot(mcmcs[[1]], ask = FALSE)
for(i in 1:3)
    plot(coda::as.mcmc(samples1[, (3*i-2):(3*i)]), ask = FALSE)


#### Posterior Inferences
#### Mean and 95% Credible intervals
results <- as.data.frame(cbind(apply(samplesList[[1]], 2, mean),
                               apply(samplesList[[1]], 2, function(x) quantile(x, 0.025)),
                               apply(samplesList[[1]], 2, function(x) quantile(x, 0.975))))
results$params <- row.names(results)
results[1:15,] # Coefficient results

#### Plotting P(occupancy) against covariates
pocc <- results[grep("p_occ", results$params),]
names(pocc) <- c("mean", "cil", "ciu")
# Load detection data set
detectData <- readRDS("../potato_psyllid_distribution_modeling/output/potato_psyllid_detection_dataset.rds")
detectData$pocc <- pocc$mean

plot(x = detectData$year, y = detectData$pocc)
plot(x = detectData$list_length, y = detectData$pocc)

xyplot(pocc ~ year|factor(season), data = detectData)














#### Wrap Up

##I think this last version of the model (Part II) is our best yet.  The handling of the month variable truly needed to be changed.

##Also, now we see significant effects from a number of the variables, with some nice temporal trends as you had hoped.

##Let me know when you've had a chance to look over this.

##Cheers!  - Daniel

