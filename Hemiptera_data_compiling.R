#### Compiling Hemiptera records from GBIF (including Essig), non-GBIF Essig specimens, CDFA, and American Museum
#### Keeping Heteroptera records in data set
#### For compiling only Sternorrhyncha and Auchenorrhyncha, see old data sets

rm(list = ls())
my.packages<-c('ecoengine', 'rgbif', 'lubridate', 
               'data.table', 'maps', 'taxize', 'ggplot2',
               'rgdal', 'rgeos', 'tidyr', 'ggmap', 'dplyr')
lapply(my.packages, require, character.only=T)

# Load functions
source("R_functions/museum_specimen_analysis_functions.R")
#setwd("C:/Users/Adam/Documents/UC Berkeley post doc/BIGCB/Pest Project/Potato psyllid")




##########################################################################################
#### Compiling GBIF and Essig records for all Hemiptera
##########################################################################################
# Getting GBIF records for Hemiptera
hemipKey <- head(name_lookup(query = "Hemiptera", rank = "order", return = "data"))$nubKey[1]
hemipCount <- occ_count(taxonKey = hemipKey, country = 'US', georeferenced = TRUE)
fieldNames <- c("name", "key", "catalogNumber", "order", "family", 
                "decimalLongitude", "decimalLatitude",
                "elevation", 
                "stateProvince", "county", "locality",
                "year", "month", "day",
                "recordedBy", "taxonRank")
hemipGBIF <- occ_search(taxonKey = hemipKey,
                        geometry = CAUNpolygon,
                        hasCoordinate = TRUE,
                        return = "data", limit = hemipCount,
                        fields = fieldNames)
hemipGBIF <- subset(hemipGBIF, taxonRank == "SPECIES" | taxonRank == "SUBSPECIES")
save(hemipGBIF, file = "data/GBIF_Records_Hemiptera.Rdata")

# Loading Essig Hemiptera records
hemipEssig <- read.csv("data/Essig_Hemiptera_Records_2015-08-31.csv", header = TRUE)
# Georeferenced Essig specimens not in database
essig.pp.etoh <- read.csv("data/Essig EtOH potato psyllid specimens georeferenced.csv", header = TRUE)

#### Combine 2 Essig datasets
# Select columns, families, and only states of interest for Hemiptera Essig records
SOI <- c("California", "Arizona", "Utah", "Nevada")
hemipEssig <- as.data.frame(rbindlist(lapply(SOI, function(x) hemipEssig[hemipEssig$StateProvince == x,])))
hemipEssig <- hemipEssig[!is.na(hemipEssig$SpecificEpithet),]
EssigcolumnNames <- c("ScientificName", "bnhm_id", "Ordr", "Family", 
                      "DecimalLongitude", "DecimalLatitude",
                      "MinElevationMeters", "MaxErrorInMeters",
                      "StateProvince", "County", "Locality",
                      "YearCollected", "MonthCollected", "DayCollected",
                      "Collector", "Associated_Taxon") 
hemipEssig <- hemipEssig[,EssigcolumnNames]
# Select columns and only States of Interest for Geolocate Psyllid records not in Essig, i.e., ETOH specimens
GeolocatecolumnNames <- c("ScientificName", "bnhm_id", "Ordr", "Family",
                          "Corrected.longitude", "Corrected.latitude",
                          "MinElevationMeters", "Corrected.uncertainty.radius",
                          "StateProvince", "County", "Locality",
                          "YearCollected", "MonthCollected", "DayCollected",
                          "Collector", "Associated_Taxon")
essig.pp.etoh <- as.data.frame(rbindlist(lapply(SOI, function(x) essig.pp.etoh[essig.pp.etoh$StateProvince == x,])))
essig.pp.etoh <- essig.pp.etoh[,GeolocatecolumnNames]

#### Combine GBIF, Essig, and non-Essig datasets
hemipGBIF$MaxErrorInMeters <- NA
hemipGBIF$Associated_Taxon <- NA
GBIFdat <- hemipGBIF[,c("name", "catalogNumber", "order", "family", 
                "decimalLongitude", "decimalLatitude",
                "elevation", "MaxErrorInMeters",
                "stateProvince", "county", "locality",
                "year", "month", "day",
                "recordedBy", "Associated_Taxon")]
# Removing Essig records from GBIF data, because error radius and associated taxon are missing from GBIF
GBIFdat <- GBIFdat[!GBIFdat$catalogNumber %in% hemipEssig$bnhm_id,]
# Essig records are duplicated within GBIF so number removed from GBIFdat > number of unique records in hemipEssig
names(essig.pp.etoh) <- names(hemipEssig)
hemipEssig <- rbind(hemipEssig, essig.pp.etoh)
names(GBIFdat) <- names(hemipEssig)
hemipRecords <- rbind(GBIFdat, hemipEssig)
hemipRecords <- hemipRecords[!is.na(hemipRecords$YearCollected),]
hemipRecords <- hemipRecords[!is.na(hemipRecords$MonthCollected),]
hemipRecords <- hemipRecords[!is.na(hemipRecords$DecimalLatitude),]

# Save Essig + GBIF Hemiptera records
save(hemipRecords, file = "output/Compiled_Essig-GBIF_Hemiptera_Records.Rdata")


#############################################################################
#### Combining Sterno + Aucheno records from Essig and GBIF with CDFA records
load("output/Compiled_Essig-GBIF_Hemiptera_Records.Rdata")
CDFArecords <- read.csv("data/CDFA_Hemiptera_Records.csv", header = TRUE)

#### Combine Essig, GBIF, and CDFA records
# Paratrioza is included in the database; change to Bactericera
CDFArecords[CDFArecords$GENUS == "Paratrioza",]
CDFArecords$GENUS[CDFArecords$GENUS == "Paratrioza"] <- "Bactericera"
# Make scientific names
CDFArecords$ScientificName <- paste(CDFArecords$GENUS, CDFArecords$SPECIES, sep = " ")
# Remove trap collected records -- these are from targeted surveys
CDFArecords <- CDFArecords[CDFArecords$TRAP.TYPE == "",]
# Split date collected to year, month, and day collected
CDFArecords$DATE.COLLECTED <- strptime(as.POSIXct(CDFArecords$DATE.COLLECTED, format = "%Y-%m-%d"), format = "%Y-%m-%d")
CDFArecords$YearCollected <- year(CDFArecords$DATE.COLLECTED)
CDFArecords$MonthCollected <- month(CDFArecords$DATE.COLLECTED)
CDFArecords$DayCollected <- day(CDFArecords$DATE.COLLECTED)
# Make unique record IDs from row numbers
CDFArecords$bnhm_id <- paste("CDFA", row.names(CDFArecords), sep = "")
CDFArecords$Ordr <- "Hemiptera"
CDFArecords <- CDFArecords[!is.na(CDFArecords$LATITUDE),]
CDFArecords$MinElevationMeters <- NA
CDFArecords$MaxErrorInMeters <- NA
CDFArecords <- CDFArecords[,c("ScientificName", "bnhm_id", "Ordr", "FAMILY", 
                                  "LONGITUDE", "LATITUDE", "MinElevationMeters", "MaxErrorInMeters",
                                  "STATE", "COUNTY", "MUNICIPALITY",
                                  "YearCollected", "MonthCollected", "DayCollected",
                                  "COLLECTOR", "HOST.PLANT")]
names(CDFArecords) <- names(hemipRecords)
# Combine GBIF-Essig records and CDFA records
hemipRecords <- rbind(hemipRecords, CDFArecords)
dim(hemipRecords)
# Save database as RDS file
saveRDS(hemipRecords, file = "output/Compiled_Hemiptera_records_CDFA-GBIF_2016-01-15.rds")


############################################################################################
#### Exploring and adding American Museum records
#### Data sets downloaded by Pete Oboyski for all Hemiptera in California 2016-01-05
amfiles <- list.files(path="data/American_museum_data_sets", full.names = TRUE)
AMrecords <- as.data.frame(rbindlist(lapply(amfiles, function(x) read.delim(x))))
saveRDS(AMrecords, file = "output/American_Museum_Full_Specimen_Data.rds")

#### Data munging
AMrecords <- readRDS("output/American_Museum_Full_Specimen_Data.rds")
# remove NAs from Lat and Genus
AMrecords$Lat <- as.numeric(levels(AMrecords$Lat))[AMrecords$Lat]
AMrecords <- AMrecords[!is.na(AMrecords$Lat),]
AMrecords <- AMrecords[!is.na(AMrecords$Genus),]
# remove NAs from date
AMrecords <- AMrecords[-grep("provided", AMrecords$Start_Date, ignore.case = TRUE),]
AMrecords <- AMrecords[!AMrecords$Start_Date == "",]
# Convert start dates to R-readable dates
AMrecords$startDate <- AMrecords$Start_Date %>% dmy()
goodDates <- AMrecords[!is.na(AMrecords$startDate),]
noday <- AMrecords[is.na(AMrecords$startDate),]
noday$startDate <- dmy(paste0("1 ",noday$Start_Date))
nonoDate <- noday[is.na(noday$startDate),]
nrow(nonoDate)
# 423 records that are still NA; a few have a good date but most are useless. Just give up and throw these away.
AMrecords <- rbind(goodDates, noday[!is.na(noday$startDate),])
nrow(AMrecords)
# Replace Paratrioza with Bactericera
AMrecords[AMrecords$Genus == "Paratrioza","Genus"] <- "Bactericera"
# Make scientific names
AMrecords$ScientificName <- paste(AMrecords$Genus, AMrecords$species, sep = " ")
# Split date collected to year, month, and day collected
AMrecords$YearCollected <- year(AMrecords$startDate)
AMrecords$MonthCollected <- month(AMrecords$startDate)
AMrecords$DayCollected <- day(AMrecords$startDate)
# Make unique record IDs from row numbers
AMrecords$bnhm_id <- AMrecords$PBIUSI
AMrecords$Ordr <- "Hemiptera"
AMrecords$Associated_Taxon <- paste(AMrecords$Host_Genus, AMrecords$Host_species, AMrecords$Host_Common_Name, sep = " ")
# What to do with the Spec_Count column? Any useful information there?
AMrecords <- AMrecords[,c("ScientificName", "bnhm_id", "Ordr", "Family", 
                          "Lon", "Lat", "Elev_m", "Lat_Lon_Accuracy",
                          "State_Prov", "Sec_Subdiv", "Locality",
                          "YearCollected", "MonthCollected", "DayCollected",
                          "Collector", "Associated_Taxon")]
# Load hemipRecords data set
hemipRecords <- readRDS("output/Compiled_Hemiptera_records_CDFA-GBIF_2016-01-15.rds")
names(AMrecords) <- names(hemipRecords)
saveRDS(AMrecords, "output/Compiled_AMNH_Hemiptera_records_2016-01-19.rds")
write.csv(AMrecords, "output/Compiled_AMNH_Hemiptera_records_2016-01-19.csv", row.names = FALSE)
hemipRecords <- rbind(hemipRecords, AMrecords)
table(hemipRecords[order(hemipRecords$Family),]$Family)
# Make Decade column
hemipRecords$Decade <- round(hemipRecords$YearCollected, -1)
hemipRecords <- hemipRecords[order(hemipRecords$Decade),]
hemipRecords$DecadeFactor <- as.numeric(factor(hemipRecords$Decade))

# Save database as RDS file
saveRDS(hemipRecords, file = "output/Compiled_Hemiptera_records_2016-01-19.rds")

nrow(hemipRecords[hemipRecords$ScientificName == "Bactericera cockerelli",])



############################################################################################
#### Exploring records and making maps
hemipRecords <- readRDS("output/Compiled_Hemiptera_records_2016-01-19.rds")
# Need to fix dates occurring after 2015 (!), for now just remove them
hemipRecords <- hemipRecords[hemipRecords$YearCollected >= 1900,]
hemipRecords <- hemipRecords[hemipRecords$YearCollected <= 2015,]
# Map Hemiptera records
# Potato psyllid points are in red
pdf("results/figures/Hemip_potato_psyllid_map.pdf")
  map("state", regions = c("california", "nevada", "utah", "arizona"))
  points(hemipRecords$DecimalLongitude, hemipRecords$DecimalLatitude, pch = 1, col = "blue")
  ppindex <- which(hemipRecords$ScientificName == "Bactericera cockerelli")
  points(hemipRecords$DecimalLongitude[ppindex], hemipRecords$DecimalLatitude[ppindex], pch = 1, col = "red")
dev.off()

#### Exploring potato psyllid records
ppRecords <- hemipRecords[hemipRecords$ScientificName == "Bactericera cockerelli",]

map("county", regions = c("california"))
points(ppRecords[ppRecords$StateProvince == "California",]$DecimalLongitude, 
       ppRecords[ppRecords$StateProvince == "California",]$DecimalLatitude, 
       pch = 1, col = "blue")

ppCounties <- table(ppRecords[ppRecords$StateProvince == "California",]$County)
write.csv(ppCounties, file = "Potato psyllid data/CA_counties_potato_psyllid_records_2014-04-20.csv", row.names = FALSE)

map("county", regions = c("california"))
points(ppRecords[ppRecords$StateProvince == "California" & ppRecords$County == "",]$DecimalLongitude, 
       ppRecords[ppRecords$StateProvince == "California" & ppRecords$County == "",]$DecimalLatitude, 
       pch = 1, col = "blue")

#### Make a map of California with PP record points colored by decade
ppRecords$Decade <- round(ppRecords$YearCollected, -1)
ppRecords <- ppRecords[order(ppRecords$Decade),]
ppRecords$DecadeFactor <- as.numeric(factor(ppRecords$Decade))

tiff("results/figures/Potato_psyllid_decadal_map.tif")
  map("county", regions = c("california"))
  points(ppRecords[ppRecords$StateProvince == "California",]$DecimalLongitude, 
         ppRecords[ppRecords$StateProvince == "California",]$DecimalLatitude,
         pch = 1, col = ppRecords$DecadeFactor)
dev.off()



#####################################################################
#### Map with ggmap
caMap <- map_data("county", regions = "california")
qplot(long, lat, data = caMap)


caMap <- get_map(location = "California", zoom = 6, maptype = "toner-hybrid")
cappRecords <- ppRecords[ppRecords$StateProvince == "California",]
cappMap <-   ggmap(caMap) + geom_point(aes(x = DecimalLongitude, y = DecimalLatitude, col = Decade),
                  data = cappRecords)
ggsave(filename = "results/figures/Potato_psyllid_decadal_map3.tiff", cappMap,
       #height = 14, width = 12, units = "in", 
       dpi = 600)


#### Exploring records of other pest species
#### To assess potential of using for PA/PO data comparison
# First, explore dataset without CDFA records. Want to know what's in the museum records
hemipMuseum <- hemipRecords[-grep("CDFA", hemipRecords$bnhm_id),]

#### sharpshooters
ss <- c("Graphocephala atropunctata", "Draeculacephala minerva", "Xyphon fulgida", "Homalodisca vitripennis")
sapply(ss, function(x) nrow(SAMuseum[SAMuseum$ScientificName == x,]), simplify = TRUE)
# Very few specimens
sstable <- as.data.frame(table(SAMuseum[SAMuseum$Family == "Cicadellidae",]$ScientificName))
sstable[order(sstable$Freq),]


#### aphids
aphidtable <- as.data.frame(table(SAMuseum[SAMuseum$Family == "Aphididae",]$ScientificName))
aphidtable[order(aphidtable$Freq),]


#### lygus hesperus
lygus <- hemipMuseum[hemipMuseum$ScientificName == "Lygus hesperus",]
nrow(lygus)

## Map lygus occurrences with ggmap
# caMap <- map_data("county", regions = "california")
# qplot(long, lat, data = caMap)
caMap <- get_map(location = "California", zoom = 6, maptype = "toner-hybrid")
calhRecords <- lygus[lygus$StateProvince == "California",]
calhMap <-   ggmap(caMap) + geom_point(aes(x = DecimalLongitude, y = DecimalLatitude, col = Decade),
                  data = calhRecords)
ggsave(filename = "results/figures/lygus_decadal_map.pdf", calhMap,
       height = 14, width = 12, units = "in", 
       dpi = 600)

# Histogram of year collected
pdf("results/figures/lygus_year_collected_histogram.pdf")
  hist(calhRecords$YearCollected)
dev.off()

# What plants are the Lygus associated with?
table(calhRecords$Associated_Taxon)
# How many associated with cotton?
cotton.lygus <- as.data.frame(rbindlist(lapply(c("Gossypium", "cotton"), function(x) calhRecords[grep(x,calhRecords$Associated_Taxon, ignore.case = TRUE),])))
# Looks like only two Lygus records on cotton