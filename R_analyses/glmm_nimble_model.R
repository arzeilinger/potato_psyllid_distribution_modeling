#### NIMBLE GLMM model of potato psyllid occupancy
#### Adapted from code from Daniel Turek

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

glmmCode <- nimbleCode({
    mu_alpha ~ dnorm(0, 0.001)
    sigma_alpha ~ dunif(0, 1000)
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu_alpha, sd = sigma_alpha)  ## site random effect
    }
    for(i in 1:8) {
        beta[i] ~ dnorm(0, 0.001)
    }
    # for(i in 1:4) {
    #     betaseason[i] ~ dnorm(0, 0.001)    ## new fixed effects for each season
    # }
    for(i in 1:N) {
        logit(p_occ[i]) <- alpha[siteID[i]] + beta[1]*list_length[i] + beta[2]*year_list_length[i] + 
          beta[3]*aet[i] + beta[4]*tmn[i] + beta[5]*tmx[i] + beta[6]*year[i] + beta[7]*month[i] + beta[8]*month2[i]
        y[i] ~ dbin(size = 1, prob = p_occ[i])
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

inits <- list(mu_alpha=0, sigma_alpha=1, alpha=rep(0,inputData$nsite), beta=rep(0,8))

modelInfo_glmm <- list(code=glmmCode, constants=constants, data=data, inits=inits, name='glmm_month_model')


#### Set up model and samplers
Rmodel <- nimbleModel(modelInfo_glmm$code,
                      modelInfo_glmm$constants,
                      modelInfo_glmm$data,
                      modelInfo_glmm$inits)

Cmodel <- compileNimble(Rmodel)

spec <- configureMCMC(Rmodel)

#### Best configuration of samplers for random effect occupancy model
spec$removeSamplers('beta[1:8]')
spec$addSampler('beta[1:8]', 'RW_block') # linear coefficients
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

save(samplesList, file = 'output/MCMC_glmm_list2.RData')

#### Loading saved MCMC run, sames as list, "samplesList"
load(file = 'output/MCMC_glmm_list.RData')



#######################################################################
#### Assessing convergence and summarizing/plotting results
#### Assessing convergence of only covariates and mu.alpha. Memory requirements too great to assess for all p_occ[i]

# Make mcmc.list with only covariates and mu.alpha
mcmcs <- mcmc.list(lapply(1:length(samplesList), function(x) as.mcmc(samplesList[[x]][,1:20])))

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


#########################################################################################################
#### Plots

#### Plotting P(occupancy) against covariates
pocc <- results[grep("p_occ", results$params),]
# Load detection data set
detectData <- readRDS("output/potato_psyllid_detection_dataset.rds")
detectData$pocc <- pocc$mean

# Year vs P(occupancy)
tiff("results/figures/glmm_year_vs_pocc.tif")
  plot(x = detectData$year, y = detectData$pocc)
  lines(smooth.spline(detectData$year, detectData$pocc, nknots = 4, tol = 1e-6), lwd = 2)
dev.off()
  
# List length vs P(occupancy)
tiff("results/figures/glmm_list_length_vs_pocc.tif")
  plot(x = detectData$list_length, y = detectData$pocc)
  lines(smooth.spline(detectData$list_length, detectData$pocc, nknots = 4, tol = 1e-6), lwd = 2)
dev.off()

# List length vs month
tiff("results/figures/glmm_month_vs_pocc.tif")
  plot(x = detectData$month, y = detectData$pocc)
  lines(smooth.spline(detectData$month, detectData$pocc, nknots = 4, tol = 1e-6), lwd = 2)
dev.off()


# trivariate plots with month and year
zz <- with(detectData, interp(x = year, y = month, z = pocc, duplicate = 'median'))
pdf("results/figures/year-month-occupancy_contourplot_nimble_glmm.pdf")
  filled.contour(zz, col = topo.colors(32), xlab = "Year", ylab = "Month")
dev.off()

