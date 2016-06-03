##################################################################
#### Analyzing NIMBLE models

my.packages <- c("coda", "lattice", "akima", "raster",
                 "tidyr", "dplyr", "maps", "rasterVis")
lapply(my.packages, require, character.only = TRUE)

#### FOR OCCUPANCY MODEL
#### Loading saved MCMC run, saved as list, "samplesList"
load(file = 'output/MCMC_month_list.RData')
# Directory for figures from occupancy model
outdir <- "results/figures/occupancy_figures/"

#### FOR GLMM MODEL
#### Loading saved MCMC run, saved as list, "samplesList"
# load(file = 'output/MCMC_glmm_list.RData')
# # Directory for figures from glmm model
# outdir <- "results/figures/glmm_figures/"


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
cbind(rhatSummary, round(pd, digits = 2))

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
# Save results for occupancy model
saveRDS(results, "results/occupancy_model_results.rds")

# # Save results for glmm model
# saveRDS(results, "results/glmm_model_results.rds")

##############################################################################################################
#### Plots

# Load occupancy MCMC results
results <- readRDS("results/occupancy_model_results.rds")

# # Load GLMM MCMC results
# results <- readRDS("results/glmm_model_results.rds")

#### Plotting P(occupancy) against covariates
pocc <- results[grep("p_occ", results$params),]
# Load detection data set
detectData <- readRDS("output/potato_psyllid_detection_dataset.rds")
detectData$pocc <- pocc$mean

#### Year vs P(occupancy)
# Making a line based on coefficient estimates
# Is not very insightful
# yearv <- as.data.frame(cbind(unique(detectData$year), unique(detectData$stdyear)))
# names(yearv) <- c("year", "stdyear")
# yearv <- yearv[order(yearv$stdyear),]
# mu_alpha <- results[results$params == "mu_alpha", "mean"]
# betaYear <- results[results$params == "beta[7]", "mean"]
# yearv$predOcc <- plogis(mu_alpha + betaYear*yearv$stdyear)

tiff(paste(outdir,"year_vs_pocc.tif",sep=""))
  plot(x = detectData$year, y = detectData$pocc,
       xlab = list("Year Collected", cex = 1.4), 
       ylab = list("Probability of occupancy", cex = 1.4),
       cex.axis = 1.3)
  lines(smooth.spline(detectData$year, detectData$pocc, nknots = 4, tol = 1e-6, df = 4), lwd = 2)
dev.off()
  
# List length vs P(occupancy)
tiff(paste(outdir,"list_length_vs_pocc.tif",sep=""))
  plot(x = detectData$list_length, y = detectData$pocc,
       xlab = list("Length of species lists", cex = 1.4), 
       ylab = list("Probability of occupancy", cex = 1.4),
       cex.axis = 1.3)
  lines(smooth.spline(detectData$list_length, detectData$pocc, nknots = 4, tol = 1e-20), lwd = 2)
dev.off()

# Month vs. P(occupancy)
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
Ref_raster <- readRDS("output/reference_raster.rds")
rasterCells <- data.frame(cellID = 1:ncell(Ref_raster))

#### For California state boundary polygon
# Download States boundaries (might take time)
out <- getData('GADM', country='United States', level=1)
# Extract California state
California <- out[out$NAME_1 %in% 'California',]

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
  rasterSummary <- rasterData %>% group_by(cellID) %>% summarise(meanOcc = mean(pocc))
  # Make a data.frame with a meanOcc value for each cell and set NAs to 0
  rasterValues <- left_join(rasterCells, rasterSummary, by = "cellID")
  rasterValues[is.na(rasterValues$meanOcc), "meanOcc"] <- 0
  poccMap <- Ref_raster
  # Replace raster values with P(occupancy)
  poccMap[] <- rasterValues$meanOcc
  # Set extent as lat/long coordinates and plot
  extent(poccMap) <- extent(California)
  fileName <- paste(outdir,"occupancy_raster_map_", year.i, ".tif", sep="")
  tiff(fileName)
    print(rasterVis::levelplot(poccMap, margin = FALSE, par.settings = GrTheme(region = brewer.pal(9, 'Greys'))) +
      layer(sp.polygons(California)))
  dev.off()
  # # Alternative method of plotting both raster and California state border
  # tiff(fileName)
  #   plot(poccMap)
  #   map("state", regions = c("california"), add = TRUE)
  # dev.off()
}

