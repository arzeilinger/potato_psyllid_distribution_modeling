##################################################################
#### Analyzing NIMBLE models

my.packages <- c("coda", "lattice", "akima", "raster",
                 "tidyr", "dplyr", "maps", "rasterVis")
lapply(my.packages, require, character.only = TRUE)

#### FOR OCCUPANCY MODEL
#### Loading saved MCMC run, sames as list, "samplesList"
load(file = 'output/MCMC_month_list.RData')

#### FOR GLMM MODEL
#### Loading saved MCMC run, sames as list, "samplesList"
# load(file = 'output/MCMC_glmm_list.RData')


#######################################################################
#### Assessing convergence and summarizing/plotting results
#### Assessing convergence of only covariates and mu.alpha. Memory requirements too great to assess for all p_occ[i]

# Make mcmc.list with only covariates and mu.alpha
mcmcs <- mcmc.list(lapply(1:length(samplesList), function(x) as.mcmc(samplesList[[x]][,1:15])))

## Rhat
coda::gelman.diag(mcmcs, autoburnin = FALSE)
## Effective sample size
effectiveSize(mcmcs)

## Posterior Density Plots
pdf("results/figures/trace_and_posterior_density_plots_occupancy.pdf")
  plot(mcmcs[[2]], ask = FALSE)
dev.off()


#### Posterior Inferences
#### Mean and 95% Credible intervals
results <- as.data.frame(cbind(apply(samplesList[[2]], 2, mean),
                               apply(samplesList[[2]], 2, function(x) quantile(x, 0.025)),
                               apply(samplesList[[2]], 2, function(x) quantile(x, 0.975))))
names(results) <- c("mean", "cil", "ciu")
results$params <- row.names(results)
results[1:15,] # Coefficient results


##############################################################################################################
#### Plots

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

tiff("results/figures/occupancy_year_vs_pocc.tif")
  plot(x = detectData$year, y = detectData$pocc,
       xlab = list("Year Collected", cex = 1.4), 
       ylab = list("Probability of occupancy", cex = 1.4),
       cex.axis = 1.3)
  lines(smooth.spline(detectData$year, detectData$pocc, nknots = 4, tol = 1e-6, df = 4), lwd = 2)
dev.off()
  
# List length vs P(occupancy)
tiff("results/figures/occupancy_list_length_vs_pocc.tif")
  plot(x = detectData$list_length, y = detectData$pocc,
       xlab = list("Length of species lists", cex = 1.4), 
       ylab = list("Probability of occupancy", cex = 1.4),
       cex.axis = 1.3)
  lines(smooth.spline(detectData$list_length, detectData$pocc, nknots = 4, tol = 1e-20), lwd = 2)
dev.off()

# Month vs. P(occupancy)
tiff("results/figures/occupancy_month_vs_pocc.tif")
  plot(x = detectData$month, y = detectData$pocc,
       xlab = list("Month collected", cex = 1.4), 
       ylab = list("Probability of occupancy", cex = 1.4),
       cex.axis = 1.3)
  lines(smooth.spline(detectData$month, detectData$pocc, nknots = 8, tol = 1e-6), lwd = 2)
dev.off()


#### trivariate plots with month and year
zz <- with(detectData, interp(x = year, y = month, z = pocc, duplicate = 'mean'))
tiff("results/figures/year-month-occupancy_contourplot_nimble_occupancy.tif")
  filled.contour(zz, col = gray(seq(0,1,by=0.05)), 
                 xlab = list("Year collected", cex = 1.4), 
                 ylab = list("Month collected", cex = 1.4),
                 cex.axis = 1.3)
dev.off()

tiff("results/figures/year-month-occupancy_contourplot_nimble_occupancy_colormap.tif")
  filled.contour(zz, col = topo.colors(32), 
                 xlab = list("Year collected", cex = 1.4), 
                 ylab = list("Month collected", cex = 1.4),
                 cex.axis = 1.3)
dev.off()


#### Making raster maps of P(occupancy)
# Read in reference raster
Ref_raster <- readRDS("output/reference_raster.rds")
rasterCells <- data.frame(cellID = 1:ncell(Ref_raster))
extent(Ref_raster) <- spTransform(Ref_raster,CRS('+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'))


# Download States boundaries (might take time)
out <- getData('GADM', country='United States', level=1)
# Extract California state
California <- out[out$NAME_1 %in% 'California',]

#### Replace cell values
#### 1930 - 1940 map
yearsForMaps <- c(1920, 1950, 1990)
for(i in 1:length(yearsForMaps)){
  year.i <- yearsForMaps[i]
  rasterData <- detectData[detectData$year >= year.i & detectData$year <= (year.i+21), c("year", "cellID", "pocc")]
  print(year.i)
  print(table(rasterData$cellID, rasterData$year))
  rasterSummary <- rasterData %>% group_by(cellID) %>% summarise(meanOcc = mean(pocc))
  rasterValues <- left_join(rasterCells, rasterSummary, by = "cellID")
  rasterValues[is.na(rasterValues$meanOcc), "meanOcc"] <- 1
  poccMap <- Ref_raster
  poccMap[] <- rasterValues$meanOcc
  extent(poccMap) <- extent(California)
  fileName <- paste("results/figures/occupancy_raster_map_", year.i, ".tif", sep="")
  tiff(fileName)
    levelplot(poccMap, margin = FALSE, par.settings = GrTheme()) +
      layer(sp.polygons(California))
  dev.off()
}

