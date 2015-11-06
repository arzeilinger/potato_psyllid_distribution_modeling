#### Functions and polygons for potato psyllid occupancy analysis
# Function to generate species' lists by splitting a data.frame of records by their collectionIDs
make_lists <- function(records = NA, min.list.length = 4){
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

# Function to generate detection/non-detection data and associated covariates for potato psyllids
# Data need to be supplied as a list of data.frames, one for each "collection"
detectDataFunc <- function(spLists){
  detectData <- data.frame(year = as.numeric(unlist(lapply(strsplit(names(spLists), ", "), function(x) x[1]))), 
                           cellID = unlist(lapply(strsplit(names(spLists), ", "), function(x) x[2])),
                           month = unlist(lapply(strsplit(names(spLists), ", "), function(x) x[3])),
                           eventID = unlist(lapply(1:length(spLists), function(x) unique(spLists[[x]][,"eventID"]))),
                           list_length = as.numeric(unlist(lapply(spLists, nrow))),
                           cwd = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$cwd, na.rm = TRUE))),
                           aet = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$aet, na.rm = TRUE))),
                           tmn = unlist(lapply(1:length(spLists), function(x) mean(spLists[[x]]$tmn, na.rm = TRUE))),
                           detection = unlist(lapply(spLists, function(x) as.numeric(any(x$Species == "Bactericera.cockerelli")))))
  detectData$month <- as.numeric(levels(detectData$month))[detectData$month]
  detectData$cellID <- as.numeric(levels(detectData$cellID))[detectData$cellID]
  return(detectData)
}


#### Function to make a 2-dimensional dataframe for variables in the ecological submodel
makeDataEcoModel <- function(var, data){
  # var = the variable to reshape
  # data = data set which contains var
  varData <- aggregate(data[,var], by = list(data$season.year, data$cellID), mean, na.rm = TRUE)
  names(varData) <- c("season.year","cellID","var")
  varData <- as.matrix(spread(varData, key = season.year, value = var)[,-1])
  varData[is.nan(varData)] <- NA
  return(varData)
}

# California WKT polygon
CApolygon <- 'POLYGON((-124.541015 41.902277, -120.058593 41.967659, -119.970703 39.095962, -114.082031 34.885930,-114.521484 32.546813,-120.410156 32.249974,  -125.419921 40.245991, -124.541015 41.902277))'
# WKT polygon of CA, AZ, UT, NV
CAUNpolygon <- 'POLYGON((-124.3828 41.994332, -110.979479 41.994332, -110.979479 41.073272, -109.045886 41.006981, -109.089831 31.309212, -110.979479 31.309212, -115.242175 32.651086, -117.263659 32.428807,  -121.130847 34.227845, -124.646472 40.340398, -124.3828 41.994332))'
