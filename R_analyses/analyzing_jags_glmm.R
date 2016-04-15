#### Analyzing GLMM JAGS output

#### Preliminaries
rm(list = ls())
my_packages<-c('data.table', 'lattice', 'akima', 'tidyr', 'RCurl', 'foreign')
lapply(my_packages, require, character.only=T)

source("R_functions/museum_specimen_analysis_functions.R")

## Load JAGS model coefficient estimates from GitHub
url <- "https://raw.githubusercontent.com/arzeilinger/potato_psyllid_distribution_modeling/master/output/climate_glmm_params.csv"
glmmResults <- getURL(url) %>% textConnection() %>% read.csv(., header = TRUE)
glmmResults[grep("beta", glmmResults$params),]

#### Analyzing ZIB GLMM model output

## Load JAGS coefficient estimates from local folder
glmmResults <- readRDS("output/climate_glmm_jags_out_params.rds")
glmmResults <- glmmResults[,c("mean", "X2.5.", "X97.5.", "r.hat")]
names(glmmResults) <- c("mean", "cil", "ciu", "rhat")
glmmResults$params <- row.names(glmmResults)
# Covariate results
glmmResults[grep("beta", glmmResults$params),]
write.csv(glmmResults, file = "output/climate_glmm_params.csv", row.names = FALSE)


## load detection dataset from GitHub
url <- "https://raw.githubusercontent.com/arzeilinger/potato_psyllid_distribution_modeling/master/output/potato_psyllid_detection_dataset.csv"
detectData <- getURL(url) %>% textConnection() %>% read.csv(., header = TRUE)
## Load detection dataset from local folder
detectData <- readRDS("output/potato_psyllid_detection_dataset.rds")
str(detectData)

##################################################################################
#### Model predictions
## Joining random effects (alphas) to detection dataset

# If the parameter output includes mu.alpha, need to change the name for the following code to work
glmmResults[glmmResults$params == "mu.alpha", "params"] <- "mu"

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

#### Predictions from full binomial model
# Setup all covariates, including quadratic and interaction
detectData$month2 <- detectData$stdmonth^2
detectData$llyr <- detectData$stdlnlist_length*detectData$stdyear
covars <- c("stdyear", "stdmonth", "month2", "stdlnlist_length", "llyr", "stdaet", "stdtmn", "stdtmx", "siteAlpha")

detectData$predOcc <- predFunc(betas = betas, covars = covars)


# plots
plot(x = detectData$year, y = detectData$detection)
lines(smooth.spline(detectData$year, detectData$predOcc))

plot(x = detectData$aet, y = detectData$predOcc)
plot(x = detectData$tmn, y = detectData$predOcc)
lines(smooth.spline(detectData$tmn, detectData$predOcc))

plot(x = detectData$month, y = detectData$predOcc)
lines(smooth.spline(detectData$month, detectData$predOcc))
plot(x = detectData$lnlist_length, y = detectData$predOcc)

# trivariate plots with month and year
zz <- with(detectData, interp(x = year, y = month, z = predOcc, duplicate = 'median'))
pdf("results/figures/year-month-occupancy_contourplot.pdf")
  filled.contour(zz, col = topo.colors(32), xlab = "Year", ylab = "Month")
dev.off()