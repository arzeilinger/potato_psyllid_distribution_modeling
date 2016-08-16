library(sp)
library(fields)
#### Load reference raster
Ref_raster <- readRDS("output/reference_raster.rds")
### Create empty output raster
empty_raster_df <- as.data.frame(Ref_raster, xy = TRUE)
empty_raster_df$layer <- NA
#### For California state boundary polygon
## Download States boundaries (might take time)
out <- getData('GADM', country='United States', level=1)
## Extract California state
California <- out[out$NAME_1 %in% 'California',]
## Reproject California boundary
California <- spTransform(California, projection(Ref_raster))
#### Replace cell values with P(occupancy) for cells with data
#### P(occupancy) values are averaged over 20 years for each cell
yearsForMaps <- c(1920, 1950, 1990) # The years for which each raster map will begin
nyears <- 20 # Number of years combined in each raster map
i <- 2
for(i in 1:length(yearsForMaps)){
        year.i <- yearsForMaps[i]
        # Select only years of interest
        rasterData <- detectData[detectData$year >= year.i & detectData$year <= (year.i+nyears), c("year", "cellID", "pocc")]
        print(year.i)
        print(dim(table(rasterData$cellID, rasterData$year)))
        # Average P(occupancy) over years for each cell
        rasterSummary <- rasterData %>% group_by(cellID) %>% summarise(meanOcc = mean(pocc)) %>% as.data.frame()
        # Use the empty raster to create an output raster
        poccMap <- empty_raster_df
        poccMap$layer[rasterSummary$cellID] <- rasterSummary$meanOcc
        # Create raster from poccMap data.frame 
        poccMap <- rasterFromXYZ(poccMap)
        # Set extent as lat/long coordinates and plot
        # extent(poccMap) <- extent(California)
        fileName <- paste("occupancy_raster_map_", year.i, ".tif", sep="")
        tiff(fileName)
        print(rasterVis::levelplot(poccMap, margin = FALSE, par.settings = GrTheme(region = brewer.pal(9, 'Greys'))) +
                      latticeExtra::layer(sp.polygons(California)))
        dev.off()
        # # Alternative method of plotting both raster and California state border
        # tiff(fileName)
        #   plot(poccMap)
        #   map("state", regions = c("california"), add = TRUE)
        # dev.off()
}

