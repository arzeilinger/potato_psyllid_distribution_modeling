#### Making lists of species for list length analysis
##### SETTING THINGS UP #####
rm(list = ls())
## Set working directory
#setwd("C:/Users/Adam/Documents/UC Berkeley post doc/BIGCB/Pest Project")

## Load libraries
my_packages<-c('rgdal', 'sp', 'maptools', 'rgeos', 'ggplot2', 
               'data.table', 'maps', 'tidyr', 'dplyr', 
               'raster', 'cowplot', 'ggmap', 'ggsn')
lapply(my_packages, require, character.only=T)

## Load functions and polygons
source("R_functions/museum_specimen_analysis_functions.R")

#### Reference raster for California
#### Using BCM 2014 raster for Minimum temperature for 11/2000
#### Resolution of BCM 2014 data = 270 m x 270 m
#### Rdata file loaded is a list with raster as 1st item, metadata as 2nd item
Ref_raster_data <- readRDS("data/CA_BCM_HST_Monthly_tmn_2000_11.Rdata")
Ref_raster <- Ref_raster_data[[1]]
nkm <- 15 # Cell dimension, as number of km on a side
# To get cells 10 km x 10 km, aggregate by 37
ncells <- round(3.703704*nkm, digits = 0) # Number of cells to aggregate into 1 new cell to get desired cell size
Ref_raster <- aggregate(Ref_raster, ncells) 
# Need to redefine extent of BCM raster in terms of lat and long 
dem <- raster("data/ca_270m_t6.asc")
extent(Ref_raster) <- extent(dem)
# Save and plot reference raster 
saveRDS(Ref_raster, file = "output/reference_raster.rds")
tiff("output/reference_raster.tif")
  plot(Ref_raster)
dev.off()

#### Load data
# Species records
hemipRecords <- readRDS("output/Compiled_Hemiptera_records_2016-01-19.rds")

# #### Check for records of various species from presence-absenced data sets
# rglpsyllids <- hemipRecords[hemipRecords$ScientificName == "Glycaspis brimblecombei",] # Red gum lerp psyllid; 39 records
# gss <- hemipRecords[hemipRecords$ScientificName == "Draeculacephala minerva",] # green sharpshooter; 267 records
# rss <- hemipRecords[grep("Carneocephala", hemipRecords$ScientificName),] # red-headed sharpshooter; 0 records
# bgss <- hemipRecords[hemipRecords$ScientificName == "Graphocephala atropunctata",] # blue-green sharpshooter; 96 records 
# hb <- hemipRecords[hemipRecords$ScientificName == "Murgantia histrionica",] # harlequin bug; 360 records
# blackscale <- hemipRecords[hemipRecords$ScientificName == "Saissetia oleae",] # black scale; 93 records
# redscale <- hemipRecords[hemipRecords$ScientificName == "Aonidiella aurantii",]

#### Potato psyllids and associated host plants
pp <- hemipRecords[hemipRecords$ScientificName == "Bactericera cockerelli",]
hostPlants <- unique(pp$Associated_Taxon)
hostPlants[grep("Solanum", hostPlants)] # 9 associated Solanum species
sd <- pp[grep("dulcamara", pp$Associated_Taxon),] # No records of Solanum dulcamara

# #### Spittlebugs; superfamily Cercopoidea
# spittleFamilies <- c("Cercopidae", "Aphrophoridae", "Clastopteridae", "Epipyidae", "Machaerotidae")
# spittlebug <- as.data.frame(rbindlist(lapply(spittleFamilies, function(x) hemipRecords[hemipRecords$Family == x,])))
# write.csv(spittlebug, file = "output/spittlebug_museum_records_2016-04-15.csv", row.names = FALSE)

# # Remove carnivorous and aquatic Heteropteran families
# List of non-herbivorous families (predaceous, paristic, and aquatic)
noherbFamilies <- c("Reduviidae", "Geocoridae", "Nabidae", "Cimicidae", "Pieidae",
                    "Phymatidae", "Anthocoridae", "Polyctenidae", "Mesoveliidae",
                    "Gerridae", "Nepidae", "Belostomatidae", "Gelastocoridae", "Hydrometridae",
                    "Naucoridae", "Notonectidae", "Corixidae", "Ochteridae", "Hebridae",
                    "Macroveliidae", "Veliidae", "Saldidae", "Leptopodidae")
# Remove acquatic and predaceous families
noherbFamilies.i <- which(hemipRecords$Family %in% noherbFamilies)
herbFamilies <- hemipRecords$Family[-noherbFamilies.i]
RawRecords <- hemipRecords[which(hemipRecords$Family %in% herbFamilies),]
RawRecords$UID <- as.character(RawRecords$bnhm_id)
RawRecords$Species <- as.character(RawRecords$ScientificName)
RawRecords$Date <- as.character(paste(RawRecords$YearCollected, RawRecords$MonthCollected, RawRecords$DayCollected, sep = "-"))
RawRecords$County <- as.character(RawRecords$County)
RawRecords$Locality <- as.character(RawRecords$Locality)
#RawRecords$Source <- as.character(RawRecords$Source)
RawRecords$Collector <- as.character(RawRecords$Collector)
RawRecords$Error <- as.numeric(as.character(RawRecords$MaxErrorInMeters))


##### DATA PROCESSING #####
## Process species records
# Remove records with a georeferencing error > 5000; need to keep points without georeference error measurements as 90% of records are NA
Records <- RawRecords[RawRecords$Error < 5000 | is.na(RawRecords$Error),]
# Records <- RawRecords
# Edit species names
# Merge Myzus persicae names
Records$ScientificName[Records$ScientificName == "Myzus (Nectarosiphon) persicae"] <- "Myzus persicae"
Records$Species <- gsub(" ", ".", Records$Species)

# Make lists of only collectors of specific focal species 
ppCollectors <- findCollectors("Bactericera.cockerelli")
lygusCollectors <- findCollectors("Lygus.hesperus")
myzusCollectors <- findCollectors("Myzus.persicae")
rpadiCollectors <- findCollectors("Rhopalosiphum.padi")

# Remove bad species names: ones without epithet or "undetermined"
speciesNames <- unique(Records$ScientificName)
badNames <- c(speciesNames[grep(" sp.", speciesNames, fixed = TRUE)],
              speciesNames[grep(" spp.", speciesNames, fixed = TRUE)],
              speciesNames[grep("undetermined", speciesNames)])
badNames.i <- which(speciesNames %in% badNames)
goodNames <- speciesNames[-badNames.i]
Records <- Records[which(Records$ScientificName %in% goodNames),]

# Add important fields for filtering and combining observations
Records$jdn <- strptime(as.POSIXct(Records$Date, format = "%Y-%m-%d"), format = "%Y-%m-%d")$yday+1 # Calculating julian date number from individual dates
Records$Season <- factor(ifelse(Records$MonthCollected == 12 | Records$MonthCollected <= 2, "winter",
                            ifelse(Records$MonthCollected > 2 & Records$MonthCollected <= 5, "spring",
                                   ifelse(Records$MonthCollected > 5 & Records$MonthCollected <= 8, "summer",
                                          "autumn"))))
#### Constructing cells, event IDs, and collection IDs
Records$site <- apply(Records[c("DecimalLongitude", "DecimalLatitude")], 1, function(x) paste(x[1], x[2], sep = ", "))
## When using BCM 2014 raster as Reference raster, need to transform lat, long.
RecordSP <- Records
coordinates(RecordSP) <- RecordSP[,c("DecimalLongitude", "DecimalLatitude")]
projection(RecordSP) <- CRS('+proj=longlat +datum=WGS84')
RecordSP <- spTransform(RecordSP,CRS('+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'))
Records$cellID <- factor(cellFromXY(Ref_raster, RecordSP)) # Note: this line of code requires that a reference raster has been previously generated
# Add a filtering "event" (i.e. year*cellID) field
Records$eventID <- apply(Records[c("YearCollected", "cellID")], 1, function(x) paste(x[1], x[2], sep = ", "))
# Add a collection ID (i.e. event*MonthCollected)
Records$collectionID <- apply(Records[c("eventID", "MonthCollected")], 1, function(x) paste(x[1], x[2], sep = ", "))
# Select only records after 1900
Records <- subset(Records, YearCollected >= 1900 & YearCollected < 2015)
# Subset to only California
Records <- Records[Records$StateProvince == "California",]
# Remove NAs
Records <- Records[which(!is.na(Records$cellID) & !is.na(Records$Species)),]

# Keep dataset of replicated Records
duplRecords <- Records
# This is the data set submitted to datadryad.org
dryadRecords <- duplRecords %>% dplyr::select(-Season, -Error, -DecadeFactor, -Date)
saveRDS(dryadRecords, file = "output/hemiptera_records_for_datadryad.rds")
write.csv(dryadRecords, file = "output/hemiptera_records_for_datadryad.csv", row.names = FALSE)

### Remove duplicate observations that include identical occupancy information 
Records <- Records[!duplicated(Records[c("Species", "collectionID")]),]


############################################################################
#### Add climate data from BCM 2014 rasters
############################################################################
## Climate extraction functions
source("R_functions/Extract_climate_data_functions.R")
cdir <- "C:/Users/Adam/Documents/UC Berkeley post doc/BIGCB/Pest Project/Climate data/BCM2014/Rasters"

#### Note: For annual Tmin and Tmax, I'm using water year summaries.
# Water years run from Oct 1 of the previous year to Sept 30 of the current year.

#### For Tmin: 
# Assume that the minimum occurred in Jan or Feb.
# I should extract the Tmin of the previous winter from when the collection was made.
# If the record was collected between Jan and Sept (month = 1:9), extract the annual min of the year collected.
# If the record was collected between Oct and Dec (month = 10:12), extract the annual min from the year previous to the year collected.
Records$tmnYear <- Records %>% with(., ifelse(MonthCollected >= 1 & MonthCollected <= 9, YearCollected, YearCollected - 1))

#### For Tmax:
# Assume that the max occurred in July, August, or September
# Make a column of t-1 where t = YearCollected
# This ensures that the temperature occurred prior to when the insects were collected.
Records$tmxYear <- Records$YearCollected - 1

Records <- subset(Records, tmxYear >= 1900)
# First, extract monthly aet and cwd data
Records <- Records %>% extractClimateMonthly(., colNames = c("DecimalLongitude", "DecimalLatitude", "YearCollected", "MonthCollected"),
                                             fac = c("aet", "cwd"), cdir = cdir) %>%
                   # Second, extract annual min summary of tmn data, for year previous to YearCollected
                   extractClimateAnnually(., colNames = c("DecimalLongitude", "DecimalLatitude", "tmnYear"),
                                          fac = c("tmn"), summaryType = "min", cdir = cdir) %>%
                   # Third, extract annual max summary of tmx data, for year previous to YearCollected
                   extractClimateAnnually(., colNames = c("DecimalLongitude", "DecimalLatitude", "tmxYear"),
                                          fac = c("tmx"), summaryType = "max", cdir = cdir)


#################################################################################
## Make species lists from Essig, GBIF, and CDFA datasets
longLists <- make_lists(Records, min.list.length = 3)
#AllLists <- make_lists(Records, 1)
speciesNames <- unique(Records$Species)
#saveRDS(AllLists, file = "output/All_Hemip_Lists_Climate_15km_Cells_2016-06-15.rds")
saveRDS(longLists, file = "output/Hemip_Long_Lists_Climate_15km_Cells_2016-06-14.rds")

##############################################################################################
#### Exploring the cleaned data set of duplicated records
## map of all Hemiptera records, with potato psyllid records in different color
duplRecords <- readRDS("output/hemiptera_records_for_datadryad.rds")
ppRecords <- duplRecords[duplRecords$Species == "Bactericera.cockerelli",]


HemipteraCAMap <- ggplot(duplRecords, aes(x = DecimalLongitude, y = DecimalLatitude)) +
  borders(database = "state", regions = "california") +
  geom_point(colour = "darkgrey", pch = 1, size = 1) +
  geom_point(data = subset(duplRecords,Species == "Bactericera.cockerelli"), colour = "black", size = 1) +
  scale_x_continuous(name = "Longitude") +
  scale_y_continuous(name = "Latitude") +
  ggsn::north(location = "bottomleft", symbol = 10, 
              x.min = min(duplRecords$DecimalLongitude), x.max = max(duplRecords$DecimalLongitude),
              y.min = min(duplRecords$DecimalLatitude), y.max = max(duplRecords$DecimalLatitude)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 
  
HemipteraCAMap

## Hemiptera CA map using map function
# tiff("results/figures/all_hemip_records_CAmap.tif")
#   map("state", regions = c("california"))
#   map.axes(cex.axis = 1.4)
#   points(duplRecords$DecimalLongitude,
#          duplRecords$DecimalLatitude,
#          pch = 1, col = "grey", cex = 1.8)
#   points(ppRecords$DecimalLongitude,
#          ppRecords$DecimalLatitude,
#          pch = 16, col = "black", cex = 1.8)
# dev.off()


#### Histogram of year collected for all specimens
all_hemip_histogram <- ggplot(duplRecords,aes(x=YearCollected)) + 
  geom_histogram(fill = "darkgrey", alpha = 1, binwidth = 1) +
  geom_histogram(data=subset(duplRecords,Species == "Bactericera.cockerelli"),fill = "black", alpha = 1, binwidth = 1) +
  xlab("Year collected") + ylab("Frequency") + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 
all_hemip_histogram

ggsave(filename = "results/figures/all_hemip_records_year_histogram.tiff", 
       plot = all_hemip_histogram)

fig1 <- plot_grid(HemipteraCAMap, all_hemip_histogram, align = "v", nrow = 2, ncol = 1, labels = "auto")
ggsave(filename = "results/figures/occupancy_figures/final_figures/figure1.tiff",
       plot = fig1,
       width = 3, height = 6, units = "in", dpi = 600)





## Map of CA raster cells as polygons
#### For California state boundary polygon
## Download States boundaries (might take time)
out <- getData('GADM', country='United States', level=1)
## Extract California state
California <- out[out$NAME_1 %in% 'California',]
## Reproject California boundary
California <- spTransform(California, projection(Ref_raster))

## Convert BCM reference raster to polygons
CApoly <- rasterToPolygons(Ref_raster, dissolve = TRUE) 

## Load species lists data set
# Load species lists data set with climate data 
longLists <- readRDS("output/Hemip_Long_Lists_Climate_15km_Cells_2016-06-14.rds")
# Make into dataframe
longListsDF <- longLists %>% rbindlist() %>% as.data.frame()
# select potato psyllid records in the long lists data set
ppFromLists <- longListsDF %>% dplyr::filter(Species == "Bactericera.cockerelli")
ppFromLists <- SpatialPointsDataFrame(ppFromLists[,c("decimalLatitude","decimalLongitude")])
coordinates(ppFromLists) <- c("DecimalLatitude","DecimalLongitude")
ppFromLists <- spTransform(ppFromLists, CRS(projection(Ref_raster)))
#### Can't get lat/long of records transformed properly

tiff("results/figures/CA_polygon.tif")
  plot(California)
  plot(CApoly, add = TRUE, border = "grey", col = "white")
  # points(ppFromLists$DecimalLongitude,
  #        ppFromLists$DecmialLatitude,
  #        pch = 16, col = "black", cex = 0.8)
dev.off()



##############################################################################################
#### Exploring the lists
## Combine lists back into data set and save for extracting climate data
listData <- as.data.frame(rbindlist(AllLists))
listData[listData$Species == "Bactericera.cockerelli",]
nrow(listData[listData$Species == "Bactericera.cockerelli",])
saveRDS(listData, "output/All_Hemip_species_lists_2016-02-18.rds")

# How are the list lengths distributed?
ListLength <- as.numeric(unlist(lapply(PALists, function(x) nrow(x))))
hist(ListLength, breaks = seq(4,40,1))

# Where are the lists distributed?
xyLists <- data.frame(xyFromCell(Ref_raster, cell = as.numeric(levels(listData$cellID))[listData$cellID]))
tiff("results/figures/pp_list_map.tif")
  map("state", regions = c("california"))
  points(xyLists$x, xyLists$y, pch = 1, col = "grey35", cex = 1.8)
  points(listData[listData$Species == "Bactericera.cockerelli",]$DecimalLongitude,
         listData[listData$Species == "Bactericera.cockerelli",]$DecimalLatitude,
         pch = 16, col = "darkgreen", cex = 1.8)
dev.off()

# What years and months are represented in the lists?
listID <- unique(listData$collectionID)
listdate <- data.frame(Year = as.numeric(unlist(lapply(strsplit(listID, ", "), function(x) x[1]))),
                       cellID = as.numeric(unlist(lapply(strsplit(listID, ", "), function(x) x[2]))),
                       Season = unlist(lapply(strsplit(listID, ", "), function(x) x[3])))
listdate <- listdate[!duplicated(listdate[c("Year", "Season")]),]
#write.csv(listdate, file = "Potato psyllid data/Species_list_collections_2015-04-21.csv", row.names = FALSE)



########################################################################################################
## Make figures for lab meeting
ppRecords <- SArecords[SArecords$Species == "Bactericera cockerelli" & SArecords$StateProvince == "California",]
# map of all Stern + Auchen records
tiff("results/figures/all_potato_psyllid_records_CAmap.tif")
  map("state", regions = c("california"))
  points(ppRecords$DecimalLongitude,
         ppRecords$DecimalLatitude,
         pch = 1, col = "darkgreen", cex = 1.8)
dev.off()

tiff("results/figures/potato_psyllid_histogram.tif")
  hist(ppRecords$YearCollected, ylab = "Frequency", xlab = "Year Collected",
       main = "", xlim = c(1900,2020))
dev.off()

# map of all Stern + Auchen records
CArecords <- SArecords[SArecords$StateProvince == "California",]
tiff("results/figures/all_SA_records_CAmap.tif")
  map("state", regions = c("california"))
  points(CArecords$DecimalLongitude,
         CArecords$DecimalLatitude,
         pch = 1, col = "darkgreen", cex = 1.8)
dev.off()

tiff("results/figures/SA_histogram.tif")
  hist(CArecords$YearCollected, ylab = "Frequency", xlab = "Year Collected",
       main = "", xlim = c(1900,2020))
dev.off()


# Map of CA raster cells as polygons
CAcrop <- readWKT(CApolygon)
CAraster <- crop(Ref_raster, CAcrop)
plot(CAraster)
CApoly <- rasterToPolygons(CAraster) 
listpp <- listData[listData$Species == "Bactericera cockerelli",]
tiff("results/figures/CA_polygon.tif")
  map("state", regions="california")
  lines(CApoly, lty = 2, col = "darkgrey")
  points(listpp$DecimalLongitude,
         listpp$DecimalLatitude,
         pch = 1, col = "darkgreen", cex = 1.8)
dev.off()



#### Make figures for Essig talk 2016-09-23
tiff("results/figures/potato_psyllids_records_CAmap.tif")
  map("state", regions = c("california"))
  map.axes(cex.axis = 1.4)
  # points(duplRecords$DecimalLongitude,
  #        duplRecords$DecimalLatitude,
  #        pch = 1, col = "grey", cex = 1.8)
  points(ppRecords$DecimalLongitude,
         ppRecords$DecimalLatitude,
         pch = 16, col = "black", cex = 1.8)
dev.off()


#### Histogram of year collected for potato psyllids specimens
pp_histogram <- ggplot(ppRecords,aes(x=YearCollected)) + 
  geom_histogram(fill = "black", alpha = 1, binwidth = 1) +
  #geom_histogram(data=subset(duplRecords,Species == "Bactericera.cockerelli"),fill = "black", alpha = 1, binwidth = 1) +
  xlab("Year collected") + ylab("Frequency") + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

ggsave(filename = "results/figures/potato_psyllid_records_year_histogram.tiff", 
       plot = pp_histogram)
