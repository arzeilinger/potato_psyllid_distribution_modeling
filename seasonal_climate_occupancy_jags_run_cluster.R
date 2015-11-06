#### Static occupancy modeling of potato psyllid occurrences
#### Using 2 composite seasons as "seasons" and months as repeat surveys 
#### Following Van Strien et al. 2015 

#### Preliminaries
my_packages<-c('snow', 'rjags', 'R2jags', 'dclone', 
               'data.table')
lapply(my_packages, require, character.only=T)

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
b5~dnorm(0,0.01)
b6~dnorm(0,0.01)
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
    tmx[i,t]~dnorm(0,0.01)
  }
}

# Ecological process sub-model
for(i in 1:nsite){
  #z[i,1]~dbern(psi)
  for(t in 1:nseason.year){
     logit(muZ[i,t]) <- a[t] + b1*tmn[i,t] + b2*pow(tmn[i,t],2) + 
                        b3*tmx[i,t] + b4*pow(tmx[i,t],2) +
                        b5*aet[i,t] + b6*cwd[i,t] + 
                        eta[i]
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
      lp[i,j,t] <- alpha.p[t] + beta1.p*month[i,j,t] + beta1.p*pow(month[i,j,t], 2) + 
                   (delta1.p*all_lists[i,j,t])/(all_lists[i,j,t] + delta2.p)
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
occDataList <- readRDS("seasonal_climate_occupancy_model_data.rds")
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
b6zst <- runif(1, -3, 3)
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
#parameters <- c("p", "delta1.p", "delta2.p",
#		"beta1.p", "beta2.p","alpha.p")
parameters <- c("a", "b1", "b2", "b3", "b4", "b5", "b6", "psi.fs", "mean.p")


# MCMC parameters
ni=41000; nb=1000; nt=10; nc=3
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
row.names(occResults) <- names(occdctab)
saveRDS(occResults, file = "seasonal_climate_occ_jags_out_params.rds")
print(occResults)
print("SUCCESS!")
