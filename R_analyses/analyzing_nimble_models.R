##################################################################
#### Analyzing NIMBLE models

my.packages <- c("coda", "lattice", "akima", "raster",
                 "tidyr", "dplyr", "maps", "rasterVis",
                 "sp", "fields")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/museum_specimen_analysis_functions.R")

#### FOR OCCUPANCY MODEL
#### Loading saved MCMC run, saved as list, "samplesList"
load(file = 'output/MCMC_month_list2.RData')
# Directory for figures from occupancy model
outdir <- "results/figures/occupancy_figures/"

#### FOR GLMM MODEL
#### Loading saved MCMC run, saved as list, "samplesList"
load(file = 'output/MCMC_glmm_list2.RData')
# Directory for figures from glmm model
outdir <- "results/figures/glmm_figures/"


#######################################################################
#### Assessing convergence and summarizing/plotting results
#### Assessing convergence of only covariates and mu.alpha. Memory requirements too great to assess for all p_occ[i]

# Make mcmc.list with only covariates and mu.alpha
mcmcs <- mcmc.list(lapply(1:length(samplesList), function(x) as.mcmc(samplesList[[x]][,c(1:15,ncol(samplesList[[x]]))])))

## Convergence diagnostics
## Rhat
rhat <- coda::gelman.diag(mcmcs, autoburnin = FALSE)
## Effective sample size
pd <- effectiveSize(mcmcs)
# Combine into one table
rhatSummary <- rhat[[1]] %>% round(., digits = 2) 
rhatSummary <- paste(rhatSummary[,1], " (", rhatSummary[,2], ")", sep="")
diagnostics <- cbind(rhatSummary, round(pd, digits = 2))
diagnostics
#write.csv(diagnostics, file = "results/occupancy_MCMC_diagnostics.csv", row.names = FALSE)

write.csv(diagnostics, file = "results/glmm_MCMC_diagnostics.csv", row.names = FALSE)

## Posterior Density Plots
pdf(paste(outdir, "trace_and_posterior_density_plots.pdf", sep=""))
  plot(mcmcs[[2]], ask = FALSE)
dev.off()


###############################################################################################
#### Posterior Inferences
#### Mean and 95% Credible intervals
results <- as.data.frame(cbind(apply(samplesList[[2]], 2, mean),
                               apply(samplesList[[2]], 2, function(x) quantile(x, 0.025)),
                               apply(samplesList[[2]], 2, function(x) quantile(x, 0.975))))
names(results) <- c("mean", "cil", "ciu")
results$params <- row.names(results)
results[1:15,] # Coefficient results

# Make results table for ms
resultsTable <- rbind(results[-grep("p_occ", results$params), c("mean", "cil", "ciu")]) %>% round(., digits = 2)
resultsTable$summary <- with(resultsTable, paste(mean, " [", cil, ", ", ciu, "]", sep = ""))
# # Save results for occupancy model
# saveRDS(results, "results/occupancy_model_results.rds")
# write.csv(results, file = "results/occupancy_model_results.csv", row.names = FALSE)

# Save results for glmm model
saveRDS(results, "results/glmm_model_results.rds")
write.csv(results, file = "results/glmm_model_results.csv", row.names = FALSE)

##############################################################################################################
#### Plots

# # Load occupancy MCMC results
results <- readRDS("results/occupancy_model_results.rds")

# Load GLMM MCMC results
#results <- readRDS("results/glmm_model_results.rds")

resultsPars <- results[-grep("p_occ", results$params), c("mean", "cil", "ciu", "params")]
# Occupancy model covars
resultsPars$covar <- c("det_intercept", "list_length", "year_list_length", "aet", "tmn", "tmx", "year", "month", "month2", NA, NA)
# GLMM model covars
#resultsPars$covar <- c("list_length", "year_list_length", "aet", "tmn", "tmx", "year", "month", "month2", NA, NA)

resultsPars


#### Plotting P(occupancy) against covariates
pocc <- results[grep("p_occ", results$params),]
# Load detection data set
detectData <- readRDS("output/potato_psyllid_detection_dataset.rds")
detectData$pocc <- pocc$mean

#### Year vs P(occupancy)
yearline <- coefline("year", "stdyear")
tiff(paste(outdir,"year_vs_pocc.tif",sep=""))
  plot(x = detectData$year, y = detectData$pocc,
       xlab = list("Year Collected", cex = 1.4), 
       ylab = list("Probability of occupancy", cex = 1.4),
       cex.axis = 1.3)
  lines(smooth.spline(detectData$year, detectData$pocc, nknots = 4, tol = 1e-6, df = 3), lwd = 2)
  #lines(yearline$covar, yearline$predOcc, lwd = 2, lty = 1)
  #abline(lm(detectData$pocc ~ detectData$year), lty = 1, lwd = 2)
dev.off()
  

#### List length vs P(occupancy)
llline <- coefline("list_length", "stdlnlist_length")
tiff(paste(outdir,"list_length_vs_pocc.tif",sep=""))
  plot(x = detectData$list_length, y = detectData$pocc,
       xlab = list("Length of species lists", cex = 1.4), 
       ylab = list("Probability of occupancy", cex = 1.4),
       cex.axis = 1.3)
  lines(smooth.spline(detectData$list_length, detectData$pocc, nknots = 4, tol = 1e-20), lwd = 2)
  #lines(llline$covar, llline$predOcc, lwd = 2, lty = 1)
dev.off()


#### Month vs. P(occupancy)
# Doesn't inlcude coefficient line
tiff(paste(outdir,"month_vs_pocc.tif",sep=""))
  plot(x = detectData$month, y = detectData$pocc,
       xlab = list("Month collected", cex = 1.4), 
       ylab = list("Probability of occupancy", cex = 1.4),
       cex.axis = 1.3)
  lines(smooth.spline(detectData$month, detectData$pocc, nknots = 8, tol = 1e-6), lwd = 2)
dev.off()

# AET vs. P(occupancy)
tiff(paste(outdir,"aet_vs_pocc.tif",sep=""))
  plot(x = detectData$aet, y = detectData$pocc,
       xlab = list("Actual evapotranspiration", cex = 1.4), 
       ylab = list("Probability of occupancy", cex = 1.4),
       cex.axis = 1.3)
  lines(smooth.spline(detectData$aet, detectData$pocc, nknots = 8, tol = 1e-6), lwd = 2)
dev.off()

# tmn vs. P(occupancy)
tiff(paste(outdir,"tmn_vs_pocc.tif",sep=""))
  plot(x = detectData$tmn, y = detectData$pocc,
       xlab = list("Annual minimum temperature (deg C)", cex = 1.4), 
       ylab = list("Probability of occupancy", cex = 1.4),
       cex.axis = 1.3)
  lines(smooth.spline(detectData$tmn, detectData$pocc, nknots = 8, tol = 1e-6), lwd = 2)
dev.off()



#### trivariate plots with month and year
zz <- with(detectData, interp(x = year, y = month, z = pocc, duplicate = 'mean'))
tiff(paste(outdir,"year-month-occupancy_contourplot_nimble_occupancy_greyscale.tif",sep=""))
  filled.contour(zz, col = gray(seq(0,1,by=0.05)), 
                 xlab = list("Year collected", cex = 1.4), 
                 ylab = list("Month collected", cex = 1.4),
                 cex.axis = 1.3)
dev.off()

tiff(paste(outdir,"year-month-occupancy_contourplot_nimble_occupancy_colormap.tif",sep=""))
  filled.contour(zz, col = topo.colors(32), 
                 xlab = list("Year collected", cex = 1.4), 
                 ylab = list("Month collected", cex = 1.4),
                 cex.axis = 1.3)
dev.off()


#### Making raster maps of P(occupancy)
# Read in reference raster
# This reference raster was used to generate the cellIDs, in "Making_species_list.R" script
# Each cell is 15km x 15km
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
nyears <- 25 # Number of years combined in each raster map
#i <- 2
#### for loop to make maps with different-sized points instead of raster cells
# The size of symbols are relative to the P(occupancy) value
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
        poccMap_points <- poccMap %>% dplyr::filter(complete.cases(.))
        #### Add points for potato psyllid detections
        ppData <- detectData[detectData$year >= year.i & detectData$year <= (year.i+nyears) & detectData$detection == 1, c("year", "cellID")]
        # CellIDs where potato psyllids were collected
        ppData <- unique(ppData$cellID)
        ppMap <- empty_raster_df
        ppMap$layer[ppData] <- 1
        ppMap <- ppMap %>% dplyr::filter(complete.cases(.))
        ppMap <- semi_join(x = poccMap, y = ppMap, by = c("x", "y"))
        fileName <- paste(outdir,"occupancy_raster_map_", year.i, "_points.tif", sep="")
        tiff(fileName)
          plot(California, border = "darkgrey")
          points(poccMap_points[, 1:2], pch = 1, cex = poccMap_points$layer * 3)
          points(ppMap[,1:2], pch = 16, cex = ppMap$layer * 3)
        dev.off()
}


#### for loop to create P(occupancy) raster maps
# for(i in 1:length(yearsForMaps)){
#         year.i <- yearsForMaps[i]
#         # Select only years of interest
#         rasterData <- detectData[detectData$year >= year.i & detectData$year <= (year.i+nyears), c("year", "cellID", "pocc")]
#         print(year.i)
#         print(dim(table(rasterData$cellID, rasterData$year)))
#         # Average P(occupancy) over years for each cell
#         rasterSummary <- rasterData %>% group_by(cellID) %>% summarise(meanOcc = mean(pocc)) %>% as.data.frame()
#         # Use the empty raster to create an output raster
#         poccMap <- empty_raster_df
#         poccMap$layer[rasterSummary$cellID] <- rasterSummary$meanOcc
#         # Create raster from poccMap data.frame 
#         poccMap <- rasterFromXYZ(poccMap)
#         # Set extent as lat/long coordinates and plot
#         # extent(poccMap) <- extent(California)
#         fileName <- paste(outdir, "occupancy_raster_map_", year.i, ".tif", sep="")
#         tiff(fileName)
#         print(rasterVis::levelplot(poccMap, margin = FALSE, par.settings = GrTheme(region = brewer.pal(9, 'Greys'))) +
#                       latticeExtra::layer(sp.polygons(California)))
#         dev.off()
#         # # Alternative method of plotting both raster and California state border
#         # tiff(fileName)
#         #   plot(poccMap)
#         #   map("state", regions = c("california"), add = TRUE)
#         # dev.off()
# }
