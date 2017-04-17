#### NIMBLE occupancy model for ZIB GLMM model of potato psyllid occupancy
#### Originally developed by Daniel Turek

# Clear workspace
rm(list = ls())
# Load packages
my.packages <- c("nimble", "coda", "lattice", "akima")
lapply(my.packages, require, character.only = TRUE)

#### Load data set for model fitting
inputData <- readRDS('output/bactericera_data_nimble_occupancy.rds')
source('R_functions/nimble_definitions.R')


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
    # Priors for season factor; not used in current model
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
spec$addMonitors(c('p_occ')) # add a monitor to get p_occ in output
#spec$addMonitors(c('p_obs')) # add a monitor to get p_obs in output

#### Compile MCMC in R and C++
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


#### Run MCMC with 150,000 iterations and 50,000 burn-in
niter <- 150000
burnin <- 50000

ti <- Sys.time()
samplesList <- lapply(3, mcmcClusterFunction)
tf <- Sys.time()

# The time it took to run MCMC
tf-ti

save(samplesList, file = 'output/MCMC_list_climate_pocc.RData')

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
print(results[1:15,]) # Coefficient results

#### For more analysis, and generating figures, see R script "analyzing_nimble_models.R"
