# California WKT polygon
CApolygon <- 'POLYGON((-124.541015 41.902277, -120.058593 41.967659, -119.970703 39.095962, -114.082031 34.885930,-114.521484 32.546813,-120.410156 32.249974,  -125.419921 40.245991, -124.541015 41.902277))'
# WKT polygon of CA, AZ, UT, NV
CAUNpolygon <- 'POLYGON((-124.3828 41.994332, -110.979479 41.994332, -110.979479 41.073272, -109.045886 41.006981, -109.089831 31.309212, -110.979479 31.309212, -115.242175 32.651086, -117.263659 32.428807,  -121.130847 34.227845, -124.646472 40.340398, -124.3828 41.994332))'


# Function to generate species' lists by splitting a data.frame of records by their collectionIDs
make_lists <- function(records = NA, min.list.length = 3){
  records_list <- split(records, records$collectionID)
  min_length_lists <- as.numeric(which(unlist(lapply(records_list, function(x) nrow(x) >= min.list.length))))
  records_list <- records_list[unique(min_length_lists)]
  return(records_list)
}

# Function to standardize covariates for occupancy models
standardize <- function(covar){
  mean.covar <- mean(covar, na.rm = TRUE)
  sd.covar <- sd(covar[!is.na(covar)])
  return((covar - mean.covar)/sd.covar)
}

# Function to transform species lists into data frame
detectDataFunc <- function(spLists, focalSpecies = "Bactericera.cockerelli"){
  detectData <- data.frame(year = as.numeric(unlist(lapply(strsplit(names(spLists), ", "), function(x) x[1]))), 
                           cellID = unlist(lapply(strsplit(names(spLists), ", "), function(x) x[2])),
                           month = unlist(lapply(strsplit(names(spLists), ", "), function(x) x[3])),
                           ymo = unlist(lapply(1:length(spLists), function(x) unique(spLists[[x]]$ymo))),
                           season = unlist(lapply(1:length(spLists), function(x) unique(spLists[[x]]$Season))),
                           collectionID = unlist(lapply(1:length(spLists), function(x) unique(spLists[[x]]$collectionID))),
                           list_length = as.numeric(unlist(lapply(spLists, nrow))),
                           collectors_length = as.numeric(unlist(lapply(1:length(spLists), function(x) length(unique(spLists[[x]][,"Collector"]))))),
                           aet = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$aet, na.rm = TRUE))),
                           cwd = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$cwd, na.rm = TRUE))),
                           tmn = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$tmn, na.rm = TRUE))),
                           tmx = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$tmx, na.rm = TRUE))),
                           aetsd = unlist(lapply(1:length(spLists), function(x) sd(spLists[[x]]$aet, na.rm = TRUE))),
                           cwdsd = unlist(lapply(1:length(spLists), function(x) sd(spLists[[x]]$cwd, na.rm = TRUE))),
                           tmnsd = unlist(lapply(1:length(spLists), function(x) sd(spLists[[x]]$tmn, na.rm = TRUE))),
                           tmxsd = unlist(lapply(1:length(spLists), function(x) sd(spLists[[x]]$tmx, na.rm = TRUE))),
                           detection = unlist(lapply(spLists, function(x) as.numeric(any(x$Species == focalSpecies)))))
  detectData$month <- as.numeric(levels(detectData$month))[detectData$month]
  detectData$cellID <- as.numeric(levels(detectData$cellID))[detectData$cellID]
  return(detectData)
}


# Make data matrices for JAGS GLMM with spatial random effect
makeEcoDataMatrix <- function(var, data = detectData, fill = NA){
  dataf <- data[,c("ymo","cellID",var)]
  dataMatrix <- eval(substitute(spread(dataf, key = cellID, value = y, fill = fill), 
                                list(y = as.name(names(dataf)[3]))))
  return(as.matrix(dataMatrix[,-1]))
}


## Function to predict potato psyllid occupancy from model results
# includes site random effects
predFunc <- function(betas = betas, covars = covars, data = detectData){
  # betas are a vector of coefficient estimates from the models
  # covars are a list of column names of covariates from the data set
  # must square any quadratic terms and multiply any interactions first, an include as separate columns
  # covariate names must be in correct order
  # assuming that the model includes random site effects, which are included as the last column of the covariate matrix
  betasf <- c(betas, 1) # The additional "1" is for the random site effects 
  xmat <- data[,covars]
  yv <- plogis(betasf %*% t(xmat)) # inner product of betas and covariates
  return(t(yv))
}


################################################################################################################
#### Functions for analysis of count data with N-mixture model or Poisson GLM

# Function for Shannon's Diversity index
diversityFunc <- function(x){
  require(vegan)
  species <- unique(x[,"Species"])
  speciesCounts <- sapply(species, function(y) sum(x[,"Species"] == y), simplify = TRUE)
  H <- diversity(speciesCounts)
  return(H)
}

# Function to transform species counts lists into data frame
countDataFunc <- function(spLists, focalSpecies = "Bactericera.cockerelli"){
  countData <- data.frame(year = as.numeric(unlist(lapply(strsplit(names(spLists), ", "), function(x) x[1]))), 
                          cellID = unlist(lapply(strsplit(names(spLists), ", "), function(x) x[2])),
                          month = unlist(lapply(strsplit(names(spLists), ", "), function(x) x[3])),
                          season = unlist(lapply(1:length(spLists), function(x) unique(spLists[[x]]$Season))),
                          eventID = unlist(lapply(1:length(spLists), function(x) unique(spLists[[x]][,"eventID"]))),
                          list_length = as.numeric(unlist(lapply(spLists, function(x) length(unique(x$Species))))),
                          total_n = as.numeric(unlist(lapply(spLists, nrow))),
                          diversity = as.numeric(unlist(lapply(spLists, diversityFunc))),
                          collectors_length = as.numeric(unlist(lapply(1:length(spLists), function(x) length(unique(spLists[[x]][,"Collector"]))))),
                          aet = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$aet, na.rm = TRUE))),
                          cwd = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$cwd, na.rm = TRUE))),
                          tmn = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$tmn, na.rm = TRUE))),
                          tmx = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$tmx, na.rm = TRUE))),
                          count = unlist(lapply(spLists, function(x) sum(x$Species == focalSpecies))))
  countData$month <- as.numeric(levels(countData$month))[countData$month]
  countData$cellID <- as.numeric(levels(countData$cellID))[countData$cellID]
  return(countData)
}


#### Gio's wrapper for glmer(), to run multiple species
glmer_wrapper <- function(species_name = NA, detection_table = NA){
        options(na.action = "na.fail")
        # run full model
        glmer_full <- glmer(detection_table[,species_name] ~ 
                            year + list_length + jdn + I(jdn^2) + min_temp + total_precip + (1 | cellID), 
                            family = binomial, data = detection_table
                            )              
        # run model selection 
        glmer_model_summary <- dredge(glmer_full, subset = expression(dc(jdn, `I(jdn^2)`)))
        # extract all models
        glmer_models <- get.models(glmer_model_summary, subset = cumsum(weight) <= .95)
        # group output
        glmer_output <- list(glmer_full = glmer_full, glmer_model_summary = glmer_model_summary, glmer_models = glmer_models)
        # return output
        saveRDS(glmer_output, paste("output/", species_name, '_glmer_output', sep = ""))
        assign(paste(species_name, '_glmer_output', sep = ''), glmer_output, pos = 1)
}

