#### NIMBLE occupancy model for ZIB GLMM model of potato psyllid occupancy
#### Originally developed by Daniel Turek

my.packages <- c("nimble", "coda", "lattice", "akima")
lapply(my.packages, require, character.only = TRUE)

#### Load data set for model fitting
inputData <- readRDS('output/data_nimble_zib.rds')
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

save(samplesList, file = 'output/MCMC_month_list.RData')

#### Loading saved MCMC run, sames as list, "samplesList"
load(file = 'output/MCMC_month_list.RData')


#######################################################################
#### Assessing convergence and summarizing/plotting results
#### Assessing convergence of only covariates and mu.alpha. Memory requirements too great to assess for all p_occ[i]

# Make mcmc.list with only covariates and mu.alpha
mcmcs <- mcmc.list(lapply(1:length(samplesList), function(x) as.mcmc(samplesList[[x]][,1:12])))

## Rhat
coda::gelman.diag(mcmcs, autoburnin = FALSE)
## Effective sample size
effectiveSize(mcmcs)

## Posterior Density Plots
pdf("results/figures/trace_and_posterior_density_plots.pdf")
  plot(mcmcs[[1]], ask = FALSE)
dev.off()


#### Posterior Inferences
#### Mean and 95% Credible intervals
results <- as.data.frame(cbind(apply(samplesList[[1]], 2, mean),
                               apply(samplesList[[1]], 2, function(x) quantile(x, 0.025)),
                               apply(samplesList[[1]], 2, function(x) quantile(x, 0.975))))
names(results) <- c("mean", "cil", "ciu")
results$params <- row.names(results)
results[1:15,] # Coefficient results


##############################################################################################################
#### Plots

#### Plotting P(occupancy) against covariates
pocc <- results[grep("p_occ", results$params),]
# Load detection data set
detectData <- readRDS("output/potato_psyllid_detection_dataset.rds")
detectData$pocc <- pocc$mean

# Year vs P(occupancy)
tiff("results/figures/occupancy_year_vs_pocc.tif")
  plot(x = detectData$year, y = detectData$pocc)
  lines(smooth.spline(detectData$year, detectData$pocc, nknots = 4, tol = 1e-6), lwd = 2)
dev.off()
  
# List length vs P(occupancy)
tiff("results/figures/occupancy_list_length_vs_pocc.tif")
  plot(x = detectData$list_length, y = detectData$pocc)
  lines(smooth.spline(detectData$list_length, detectData$pocc, nknots = 4, tol = 1e-6), lwd = 2)
dev.off()

# List length vs month
tiff("results/figures/occupancy_month_vs_pocc.tif")
  plot(x = detectData$month, y = detectData$pocc)
  lines(smooth.spline(detectData$month, detectData$pocc, nknots = 4, tol = 1e-6), lwd = 2)
dev.off()


# trivariate plots with month and year
zz <- with(detectData, interp(x = year, y = month, z = pocc, duplicate = 'median'))
pdf("results/figures/year-month-occupancy_contourplot_nimble_occupancy.pdf")
  filled.contour(zz, col = topo.colors(32), xlab = "Year", ylab = "Month")
dev.off()



######################################################################################################
#### Extra code

# set.seed(1)
# Cmcmc$run(niter)
# samples1 <- as.matrix(Cmcmc$mvSamples)[(burnin+1):niter,]
# 
# set.seed(2)
# Cmcmc$run(niter)
# samples2 <- as.matrix(Cmcmc$mvSamples)[(burnin+1):niter,]
# 
# set.seed(3)
# Cmcmc$run(niter)
# samples3 <- as.matrix(Cmcmc$mvSamples)[(burnin+1):niter,]
# 
# samplesList <- list(samples1, samples2, samples3)
# 
# save(samples1, samples2, samples3, file = 'output/MCMC_season.RData')


# ##############################################################################
# #### Attempt at cluster
# sfInit(parallel = TRUE, cpus = 2)
# date()
# 
# sfLibrary(nimble)
# sfExport("Cmcmc", "Rmodel", "niter", "burnin", "mcmcClusterFunction")
# sfExport("spec")
# sfExport("Cmodel")
# sfExport("Rmcmc")
# sfLapply(1:2, mcmcClusterFunction)
# 
# date()
# sfStop() # Close the cluster
# #################################################################################
