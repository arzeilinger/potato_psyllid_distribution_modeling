#### Analyzing GLMM Bayesian output

#### Preliminaries
my_packages<-c('data.table', 'snow', 'dclone', 'rjags', 'R2jags',
               'lattice', 'akima', 'tidyr', 'RCurl', 'foreign')
lapply(my_packages, require, character.only=T)

# Load model coefficient estimates
glmmResults <- readRDS("climate_glmm_jags_out_params.rds")
glmmResults <- glmmResults[,c("mean", "X2.5.", "X97.5.", "r.hat")]
names(glmmResults) <- c("mean", "cil", "ciu", "rhat")
glmmResults$params <- row.names(glmmResults)
# Covariate results
glmmResults[grep("beta", glmmResults$params),]

# load detection data set
detectData <- readRDS("potato_psyllid_detection_dataset.rds")
str(detectData)

##################################################################################
#### Model predictions
## Joining random effects (alphas) to detection dataset
# Make a data.frame that links the alpha estimate to the cellID, sorting the cellIDs is critical! 
siteRE <- data.frame(alpha = glmmResults[grep("alpha", glmmResults$params), "mean"],
                     cellID = sort(unique(detectData$cellID)))
# Add alphas to detection dataset for each cellID
detectData <- detectData[order(detectData$cellID),]
siteAlpha <- rep(NA, length(siteRE$cellID))
for(i in 1:length(siteRE$cellID)){
  cell.i <- siteRE$cellID[i]
  alpha.i <- siteRE$alpha[i]
  siteAlpha[which(detectData$cellID == cell.i)] <- alpha.i
}
detectData <- cbind(detectData, siteAlpha)

## Model coefficients
betas <- glmmResults[grep("beta", glmmResults$params), "mean"]
mu <- glmmResults[glmmResults$params == "mu", "mean"]

## Function to predict potato psyllid occupancy from model results
# includes site random effects
predFunc <- function(x){
  yv <- plogis(mu + beta[1]*x[1] + beta[2]*x[2] + beta[3]*x[2]^2 + # grand mean, year, and month covariates
                 beta[4]*x[3] + # list length
                 beta[5]*x[4] + beta[6]*x[5] + beta[7]*x[6] + beta[8]*x[7] + # climate covariates
                 x[8]) # site random effect (alpha)
  return(yv)
}

# Use standardized old data to get predicted occupancy
stdxmat <- detectData[,c("stdyear", "stdmonth", "stdlnlist_length", "stdaet", "stdcwd", "stdtmn", "stdtmx", "siteAlpha")]
# Predicted occupancy
detectData$predocc <- apply(stdxmat, 1, predFunc)

# plots
plot(x = detectData$year, y = detectData$predocc)
plot(x = detectData$aet, y = detectData$predocc)
plot(x = detectData$tmn, y = detectData$predocc)
plot(x = detectData$month, y = detectData$predocc)
plot(x = detectData$lnlist_length, y = detectData$predocc)

# trivariate plots with month and year
zz <- with(detectData, interp(x = year, y = month, z = predocc, duplicate = 'mean'))
pdf("year-month-occupancy_contourplot.pdf")
  filled.contour(zz, col = topo.colors(32),
                 xlab = "Year", ylab = "Month")
dev.off()