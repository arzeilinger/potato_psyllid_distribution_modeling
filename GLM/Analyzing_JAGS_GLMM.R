#### Analyzing GLMM Bayesian output

#### Preliminaries
rm(list = ls())
my_packages<-c('data.table', 'snow', 'dclone', 'rjags', 'R2jags',
               'lattice', 'akima')
lapply(my_packages, require, character.only=T)

setwd("C:/Users/Adam/Documents/GitHub/potato_psyllid_distribution_modeling/GLM")

# load Bayesian GLMM results
glmmOutput <- readRDS("climate_glmm_jags_out_full.rds")
glmmdctab <- dctable(glmmOutput)
glmmResults <- data.frame(rbindlist(glmmdctab))
row.names(glmmResults) <- names(glmmdctab)

glmmResults <- readRDS("climate_glmm_jags_out_params.rds")
glmmResults <- glmmResults[,c("mean", "X2.5.", "X97.5.", "r.hat")]
names(glmmResults) <- c("mean", "cil", "ciu", "rhat")
glmmResults$params <- row.names(glmmResults)

# load detection data set
detectData <- readRDS("potato_psyllid_detection_dataset.rds")
str(detectData)

##################################################################################
#### Model predictions
# Use standardized old data
stdxmat <- detectData[,c("stdyear", "stdmonth", "stdlnlist_length", "stdaet", "stdcwd", "stdtmn", "stdtmx")]
beta <- glmmResults$mean[-9]
mu <- glmmResults$mean[9]
predFunc <- function(x){
  yv <- plogis(mu + beta[1]*x[1] + beta[2]*x[2] + beta[3]*x[2]^2 + # grand mean, year, and month covariates
                 beta[4]*x[3] + # list length
                 beta[5]*x[4] + beta[6]*x[5] + beta[7]*x[6] + beta[8]*x[7]) # climate covariates
  return(yv)
}

detectData$predocc <- apply(stdxmat, 1, predFunc)

plot(x = detectData$year, y = detectData$predocc)
plot(x = detectData$aet, y = detectData$predocc)
plot(x = detectData$tmn, y = detectData$predocc)
plot(x = detectData$month, y = detectData$predocc)

zz <- with(detectData, interp(x = year, y = month, z = predocc, duplicate = 'mean'))
pdf("year-month-occupancy_contourplot.pdf")
  filled.contour(zz, col = topo.colors(32),
                 xlab = "Year", ylab = "Month")
dev.off()