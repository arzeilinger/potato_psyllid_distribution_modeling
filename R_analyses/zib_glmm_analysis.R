#### Zero-inflated (Bayesian) binomial GLMM analysis

#### Preliminaries
#rm(list = ls())
my_packages<-c('data.table', 'snow', 'dclone', 'rjags', 'R2jags')
lapply(my_packages, require, character.only=T)

#setwd("C:/Users/Adam/Documents/GitHub/potato_psyllid_distribution_modeling/GLM")

jagsGLMMdata <- readRDS("Data_JAGS_GLMM.rds")

################################################################################
#### Zero-Inflated Binomial GLMM Model specification
sink("climate_zib_glmm_model.jags")
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
    muq ~ dnorm(0, 0.001)
    
    # For the fixed effect coefficients
    betaq ~ dnorm(0, 0.001)
    betap1 ~ dnorm(0, 0.001)
    betap2 ~ dnorm(0, 0.001)
    betap3 ~ dnorm(0, 0.001)
    betap4 ~ dnorm(0, 0.001)
    betap5 ~ dnorm(0, 0.001)
    betap6 ~ dnorm(0, 0.001)
    betap7 ~ dnorm(0, 0.001)
    betap8 ~ dnorm(0, 0.001)

    # Likelihood
    for (i in 1:nlist){ # i = events (year-months)
      for(j in 1:nsite) { # j = sites
        detectionMatrix[i,j] ~ dbern(q[i,j])  # Distribution for random part; observed presences relating to detection probability
        Y[i,j] ~ dbern(p[i,j]) # Occupancy probability
        logit(q[i,j]) <- Y[i,j] + muq + betaq*list_length[i,j] # Logistic regression for detection
        logit(p[i,j]) <- betap1*year[i,j] + betap2*pow(year[i,j],2) + # Year quadratic effects
                          betap2*month[i,j] + betap3*pow(month[i,j],2) + # month quadratic effects
                          betap4*aet[i,j] + betap5*cwd[i,j] + betap6*tmn[i,j] + betap7*tmx[i,j] + # Climate effects
			  alpha[j] # Random effects
      } #j
    } #i
    }",fill = TRUE)
sink()

#########################################################################################

#### Specifications of JAGS run
# Initial values
# Specify initial values for mu.alpha, sigma.alpha, and beta1
inits <- function() list(mu.alpha = runif(1, -3, 3),
                         sigma.alpha = runif(1, 0, 5),
                         muq = runif(1, -3, 3),
                         betaq = runif(1, -3, 3),
                         betap1 = runif(1, -3, 3),
                         betap2 = runif(1, -3, 3),
                         betap3 = runif(1, -3, 3),
                         betap4 = runif(1, -3, 3),
                         betap5 = runif(1, -3, 3),
                         betap6 = runif(1, -3, 3),
                         betap7 = runif(1, -3, 3),
                         beta8 = runif(1, -3, 3))
# Monitored parameters
params <- c('muq', 'alpha',
            'betaq', 'betap1', 'betap2', 'betap3', 'betap4', 'betap5', 'betap6', 'betap7', 'betap8')
            #, 'beta8', 'alpha')
# MCMC specifications
ni=81000; nt=10; nc=3
# for jags.parfit(), burn-in iterations = n.adapt + n.update
n.adapt <- 500; n.update <- 500
nb <- n.adapt + n.update


#### Non-parallel jags()
# glmmOutput <- jags(data = jagsGLMMdata, 
#                 inits = inits, 
#                 parameters.to.save = params, 
#                 model.file = "climate_zib_glmm_model.jags", 
#                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
#                 working.directory = getwd())    
# print(glmmOutput)

#### Parallel JAGS
# Make a SOCK cluster using snow
cl <- makeCluster(3, type = "SOCK")
date()
# Call to jags.parfit
glmmOutput <- jags.parfit(cl, data = jagsGLMMdata,
                            params = params,
                            model = "climate_zib_glmm_model.jags",
                            inits = inits,
                            n.adapt = n.adapt, n.update = n.update,
                            n.iter = ni, thin = nt, n.chains = nc)
date()
stopCluster(cl) # Close the cluster
#### Compute statistics and save output
saveRDS(glmmOutput, file = "zib_glmm_jags_out_full.rds")
glmmdctab <- dctable(glmmOutput)
glmmResults <- data.frame(rbindlist(glmmdctab))
row.names(glmmResults) <- names(glmmdctab)
saveRDS(glmmResults, file = "zib_glmm_jags_out_params.rds")
print("SUCCESS!")

