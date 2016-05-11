#### Bayesian GLMM analysis

#### Preliminaries
#rm(list = ls())
my_packages<-c('data.table', 'snow', 'dclone', 'rjags', 'R2jags')
lapply(my_packages, require, character.only=T)

#setwd("C:/Users/Adam/Documents/GitHub/potato_psyllid_distribution_modeling/GLM")

jagsGLMMdata <- readRDS("Data_JAGS_GLMM.rds")

################################################################################
#### Binomial GLMM Model specification
sink("climate_glmm_model.jags")
cat("
    model {    
    ## Priors
    
    # For the site random effect
    for(j in 1:nsite) { 
      alpha[j] ~ dnorm(mu.alpha, tau.alpha) 
    }
    mu.alpha ~ dnorm(0, 0.001)
    tau.alpha <- 1 / (sigma.alpha * sigma.alpha)
    sigma.alpha ~ dunif(0, 5)
    
    # Grand mean
    #mu ~ dnorm(0, 0.01)

    # For the fixed effect coefficients
    beta1 ~ dnorm(0, 0.001)
    beta2 ~ dnorm(0, 0.001)
    beta3 ~ dnorm(0, 0.001)
    beta4 ~ dnorm(0, 0.001)
    beta5 ~ dnorm(0, 0.001)
    beta6 ~ dnorm(0, 0.001)
    beta7 ~ dnorm(0, 0.001)
    beta8 ~ dnorm(0, 0.001)
    beta9 ~ dnorm(0, 0.001)
    beta10 ~ dnorm(0, 0.001)

    # Likelihood
    for (i in 1:nlist){ # i = events (year-months)
      for(j in 1:nsite) { # j = sites
        detectionMatrix[i,j] ~ dbern(p[i,j])  # Distribution for random part
        logit(p[i,j]) <- beta1*list_length[i,j] + beta2*list_length[i,j]*year[i,j] + # List length effects
                          beta3*year[i,j] + beta4*month[i,j] + beta5*pow(month[i,j], 2) + # Month quadratic effect
                          beta6*month[i,j]*year[i,j] + beta7*pow(month[i,j], 2)*year[i,j] + # month-year interactions
                          beta8*aet[i,j] + beta9*tmn[i,j] + beta10*tmx[i,j] + # Climate effects
			  alpha[j] # Site random effects
      } #j
    } #i
    }",fill = TRUE)
sink()

#########################################################################################

#### Specifications of JAGS run
# Initial values
# Specify initial values for mu.alpha, sigma.alpha, and beta1
inits <- function() list(mu.alpha = runif(1, -3, 3),
                         sigma.alpha = runif(1, 0, 5))#,
                         #beta1 = runif(1, -3, 3),
                         #beta2 = runif(1, -3, 3),
                         #beta3 = runif(1, -3, 3),
                         #beta4 = runif(1, -3, 3),
                         #beta5 = runif(1, -3, 3),
                         #beta6 = runif(1, -3, 3),
                         #beta7 = runif(1, -3, 3))
# Monitored parameters
params <- c('beta1', 'beta2', 'beta3', 'beta4', 'beta5', 'beta6', 'beta7', 'beta8', 'beta9', 'beta10', 'alpha')
# MCMC specifications
ni=101000; nt=10; nc=3
# for jags.parfit(), burn-in iterations = n.adapt + n.update
n.adapt <- 500; n.update <- 500
nb <- n.adapt + n.update


#### Non-parallel jags()
# glmmOutput <- jags(data = jagsTestdata, 
#                 inits = inits, 
#                 parameters.to.save = params, 
#                 model.file = "climate_glmm_model.jags", 
#                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
#                 working.directory = getwd())    


#### Parallel JAGS
# Make a SOCK cluster using snow
cl <- makeCluster(3, type = "SOCK")
date()
# Call to jags.parfit
glmmOutput <- jags.parfit(cl, data = jagsGLMMdata,
                            params = params,
                            model = "climate_glmm_model.jags",
                            inits = inits,
                            n.adapt = n.adapt, n.update = n.update,
                            n.iter = ni, thin = nt, n.chains = nc)
date()
stopCluster(cl) # Close the cluster
#### Compute statistics and save output
saveRDS(glmmOutput, file = "climate_glmm_jags_out_full.rds")
glmmdctab <- dctable(glmmOutput)
glmmResults <- data.frame(rbindlist(glmmdctab))
glmmResults$params <- names(glmmdctab)
saveRDS(glmmResults, file = "climate_glmm_jags_out_params.rds")
print(glmmResults[grep("beta", glmmResults$params),])
print("SUCCESS!")
