#### Bayesian binomial GLMM for list length analysis of potato psyllid occupancy
#### Inspired by Isaac et al. 2014
#### Specifying spatial random effect as in Van Strien et al. 2015 (I think) 

rm(list = ls())
#### Preliminaries
my_packages<-c('data.table', 'tidyr', 'lattice', 'dplyr', 'ggplot2')
lapply(my_packages, require, character.only=T)

## load functions
source("R_functions/museum_specimen_analysis_functions.R")

# Load species lists data set with climate data 
AllLists <- readRDS("output/All_Hemip_Lists_Climate_15km_Cells_2016-04-14.rds")

AllLists <- readRDS("output/All_Hemip_Lists_Climate_15km_Cells_2016-06-15.rds")


# Collectors of potato psyllids, from RawRecords data set in making_species_lists
ppCollectors <- readRDS("output/potato_psyllid_collectors.rds")

# Keep only long lists (ll >= 3)
longLists <- AllLists %>% rbindlist() %>% as.data.frame() %>% make_lists(., min.list.length = 3)

longLists <- readRDS("output/Hemip_Long_Lists_Climate_15km_Cells_2016-06-14.rds")
longListsDF <- longLists %>% rbindlist() %>% as.data.frame()

# select only lists that contain collectors of potato psyllids
ppcCollections <- longListsDF[longListsDF$Collector %in% ppCollectors, "collectionID"]
ppcData <- longListsDF[longListsDF$collectionID %in% ppcCollections,]
ppcLists <- ppcData %>% make_lists(., min.list.length = 3)

ppcLists <- onlyCollectors(ppCollectors)

# Transform to data frame with pp detection
detectData <- detectDataFunc(ppcLists) 
detectData <- dplyr::filter(detectData, !is.na(aet) & !is.na(cwd) & !is.na(tmn) & !is.na(tmx))
detectData$lnlist_length <- log(detectData$list_length)
detectData$seasonNum <- detectData$season %>% as.numeric() # Factor levels: 1 = autumn, 2 = spring, 3 = summer, 4 = winter
str(detectData)
table(detectData$detection)

#### Exploring intercorrelations among climate variables
climateVars <- detectData[,c("year", "aet", "cwd", "tmn", "tmx")]
tiff("results/figures/climate_variables_intercorrelations_plot.tif")
  pairs(climateVars)
dev.off()
# CWD is highly correlated with AET, probably should drop CWD

# standardize numeric covariates, include as new variables in data frame
covars <- c("year", "month", "lnlist_length", "aet", "cwd", "tmn", "tmx")
covars.i <- as.numeric(sapply(covars, function(x) which(names(detectData) == x), simplify = TRUE))
for(i in covars.i){
  var.i <- names(detectData)[i]
  stdname.i <- paste("std", var.i, sep = "")
  stdvar.i <- standardize(detectData[,var.i])
  detectData[,stdname.i] <- stdvar.i
}

# Additional covariates for quadratic effects and interactions
detectData$stdmonth2 <- detectData$stdmonth^2
detectData$stdllyr <- detectData$stdlnlist_length*detectData$stdyear
str(detectData)

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
                        year_list_length = stdllyr,
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


############################################################################################
#### Figures for lists

#### Histogram of list length
# Just potato psyllid occurrences
ppData <- detectData[detectData$detection == 1,]

list_length_histogram <- ggplot(detectData,aes(x=list_length)) + 
  geom_histogram(fill = "darkgrey", alpha = 1, binwidth = 1) +
  geom_histogram(data=subset(detectData,detection == 1),fill = "black", alpha = 1, binwidth = 1) +
  xlab("List length") + ylab("Frequency") + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) 

ggsave(filename = "results/figures/list_length_histogram.tiff", plot = list_length_histogram)


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