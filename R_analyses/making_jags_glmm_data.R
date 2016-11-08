#### Bayesian binomial GLMM for list length analysis of potato psyllid occupancy
#### Inspired by Isaac et al. 2014
#### Specifying spatial random effect as in Van Strien et al. 2015 (I think) 

rm(list = ls())
#### Preliminaries
my_packages<-c('data.table', 'tidyr', 'lattice', 'dplyr', 'ggplot2', "lme4", "lmerTest")
lapply(my_packages, require, character.only=T)

## load functions
source("R_functions/museum_specimen_analysis_functions.R")

# Load species lists data set with climate data 
longLists <- readRDS("output/Hemip_Long_Lists_Climate_15km_Cells_2016-06-14.rds")

# Make into dataframe
longListsDF <- longLists %>% rbindlist() %>% as.data.frame()

# Collectors of potato psyllids, from RawRecords data set in making_species_lists
ppCollectors <- readRDS("output/Bactericera_cockerelli_collectors.rds")
# select only lists that contain collectors of potato psyllids
ppLists <- onlyCollectors(ppCollectors)

# Transform to data frame with pp detection
detectData <- detectDataFunc(ppLists) 
detectData <- dplyr::filter(detectData, !is.na(aet) & !is.na(cwd) & !is.na(tmn) & !is.na(tmx))
detectData$lnlist_length <- log(detectData$list_length)
detectData$seasonNum <- detectData$season %>% as.numeric() # Factor levels: 1 = autumn, 2 = spring, 3 = summer, 4 = winter
str(detectData)
table(detectData$detection)


# standardize numeric covariates, include as new variables in data frame
detectData$month2 <- detectData$month^2
covars <- c("year", "month", "lnlist_length", "aet", "cwd", "tmn", "tmx", "month2")
covars.i <- as.numeric(sapply(covars, function(x) which(names(detectData) == x), simplify = TRUE))
for(i in covars.i){
  var.i <- names(detectData)[i]
  stdname.i <- paste("std", var.i, sep = "")
  stdvar.i <- standardize(detectData[,var.i])
  detectData[,stdname.i] <- stdvar.i
}
str(detectData)


# Additional covariates for quadratic effects and interactions
detectData$llyr <- detectData$stdlnlist_length*detectData$stdyear
detectData$yearmonth <- detectData$stdyear*detectData$stdmonth
detectData$yearmonth2 <- detectData$stdyear*detectData$stdmonth*detectData$stdmonth


# Save detectData
saveRDS(detectData, file = "output/potato_psyllid_detection_dataset.rds")
write.csv(detectData, file = "output/potato_psyllid_detection_dataset.csv", row.names = FALSE)


#### Make list of "flat" data vectors for NIMBLE model
# make indices
N <- nrow(detectData)
siteID <- detectData$cellID %>% factor(., levels = unique(.)) %>% as.numeric()
nsite <- detectData$cellID %>% unique() %>% length()

nimbleData <- with(detectData, 
                   list(N = N,
                        nsite = nsite,
                        aet = stdaet,
                        tmn = stdtmn,
                        tmx = stdtmx,
                        year = stdyear,
                        month = stdmonth,
                        month2 = stdmonth2,
                        season = seasonNum,
                        list_length = stdlnlist_length,
                        year_list_length = llyr,
                        year_month = yearmonth,
                        year_month2 = yearmonth2,
                        y = detection,
                        siteID = siteID))
saveRDS(nimbleData, file = "output/data_nimble_zib.rds")

#saveRDS(nimbleData, file = "../zib_glmm/data/data_nimble_zib.rds")



#### Make JAGS data set, with lists as rows (i) and sites as columns (j)
detectionMatrix <- makeEcoDataMatrix("detection")
stdcovars <- paste("std", covars, sep = "")
jagsData <- lapply(stdcovars, function(x) makeEcoDataMatrix(x, fill = 0))
names(jagsData) <- stdcovars

jagsGLMMdata <- list(detectionMatrix = detectionMatrix,
                     year = jagsData$stdyear,
                     month = jagsData$stdmonth,
                     list_length = jagsData$stdlnlist_length,
                     aet = jagsData$stdaet,
                     cwd = jagsData$stdcwd,
                     tmn = jagsData$stdtmn,
                     tmx = jagsData$stdtmx,
                     nlist = nrow(detectionMatrix),
                     nsite = ncol(detectionMatrix))
saveRDS(jagsGLMMdata, file = "output/Data_JAGS_GLMM.rds")



#############################################################################################
#### Exploring intercorrelations among climate variables
climateVars <- detectData[,c("year", "aet", "cwd", "tmn", "tmx")]
tiff("results/figures/climate_variables_intercorrelations_plot.tif")
  pairs(climateVars)
dev.off()
# CWD is highly correlated with AET, probably should drop CWD


#### Year vs. climate variables
climateVars <- c("aet", "tmn", "tmx")
outdir <- "results/figures/"
for(i in 1:length(climateVars)){
  var.i <- climateVars[i]
  print(var.i)
  fileName <- paste(outdir,"plot_year_vs_", var.i, ".tif", sep="")
  tiff(fileName)
    plot(y = detectData[, var.i], x = detectData$year, 
         pch = 16, col = "black",
         ylab = var.i, xlab = "Year")
    lines(smooth.spline(detectData$year, detectData[, var.i], nknots = 4, tol = 1e-6), lwd = 2, col = "blue")
  dev.off()
}


#############################################################################################
#### List-level variation in climate variables

hist(detectData$aetsd)
hist(detectData$tmnsd)
hist(detectData$tmxsd)

climatesd <- detectData[,c("aetsd", "tmnsd", "tmxsd")]

for(i in 1:ncol(climatesd)){
  print(names(climatesd)[i])
  print(mean(climatesd[,i], na.rm = TRUE))
  print(median(climatesd[,i], na.rm = TRUE))
}


##############################################################################################
#### Variation in climate variables over whole data set
#### Using the data set imported at the top: output/Hemip_Long_Lists_Climate_15km_Cells_2016-06-14.rds

hist(longListsDF$aet)
hist(longListsDF$tmn)
hist(longListsDF$tmx)

for(i in 1:ncol(climateVars)){
  print(names(climateVars)[i])
  print(mean(climateVars[,i], na.rm = TRUE))
  print(sd(climateVars[,i]))
}


#### What proportion of data points had tmn < 7 deg C (developmental threshold of potato psyllid)?
hist(detectData$tmn)
nbt <- detectData$tmn[detectData$tmn < 7] %>% length()
propbt <- nbt/length(detectData$tmn)


detectData$belowThreshold <- ifelse(detectData$tmn < 7, 1, 0)
table(detectData$detection, detectData$belowThreshold)
with(detectData, chisq.test(detection, belowThreshold, simulate.p.value = TRUE))


#### Contributions of different collections
longListsDF[grep("CDFA", longListsDF$UID),] %>% nrow()
longListsDF[grep("AMNH", longListsDF$UID),] %>% nrow()


#### Prepare data set for archiving
#### Use longListsDF-- includes BCM data and duplicate records have been filtered out
archiveData <- longListsDF[,c("Species", "Family", "DecimalLongitude", "DecimalLatitude", "MaxErrorInMeters",
                              "County", "Locality", "Collector", "Associated_Taxon", "UID", "Date", "site", "cellID", "collectionID",
                              "aet", "cwd", "tmn", "tmx")]
# CDFA PPDC data
ppdc <- archiveData[grep("CDFA", archiveData$UID),]
write.csv(ppdc, file = "output/CDFA_records_for_archiving.csv", row.names = FALSE)

# AMNH data
amnh <- archiveData[grep("AMNH", archiveData$UID),]
write.csv(amnh, file = "output/AMNH_records_for_archiving.csv", row.names = FALSE)


################################################################################################################
#### Exploring other species in the lists
speciesFreq <- as.data.frame(table(longListsDF$Species))
names(speciesFreq) <- c("species", "count")
speciesFreq <- speciesFreq[order(speciesFreq$count),]

# Are there Rhopalosiphum sp. aphids in the data set?
speciesFreq[grep("Rhopalosiphum", speciesFreq$species),]