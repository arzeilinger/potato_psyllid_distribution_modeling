#### Script to extract climate data from rasters

extractClimateMonthly <- function(ll, colNames = c("long", "lat", "year_nm", "month_nm"), fac=c('aet', 'cwd'),cdir=getwd()) {
    ## ll is a dataframe that includes columns with longitude, latitude, year, and month, with column names specified by "colNames" statement
    require(sp)
    require(data.table)
    require(raster)
    long <- colNames[1]
    lat <- colNames[2]
    year_nm <- colNames[3]
    month_nm <- colNames[4]
    # Get good extent from other raster
    dem <- raster('C:/Users/Adam/Documents/UC Berkeley post doc/BIGCB/Pest Project/Climate data/ca_270m_t6.asc')
    good.extent <- extent(dem)
    # Load a BCM raster to use in case of error
    #errorRaster <- readRDS(paste(cdir,'/tmx/CA_BCM_HST_Monthly_tmx_1920_01.Rdata',sep=''))[[1]]
    # Prepend single-digit months with 0s
    ll[,month_nm] <- formatC(ll[,month_nm], width = 2, format = "d", flag = "0") 
    ll$ymo <- paste(ll[,year_nm], ll[,month_nm], sep = "_")
    ymoLists <- split(ll, ll$ymo)
#     if (FALSE) plot(llSP.aea) # I'm not clear what this is?
    lfac <- length(fac)
    #climvals <- matrix(NA,nrow=length(ymoLists),ncol=lfac)
    #i=1
    extractFunc <- function(x){
      listName <- names(ymoLists)[x]
      listData <- ymoLists[[x]]
      listDataSP <- listData
      coordinates(listDataSP) <- listDataSP[,c(long,lat)]
      projection(listDataSP) <- CRS('+proj=longlat +datum=WGS84')
      #class(llSP)
      #listDataSP <- SpatialPoints(listDataSP[,c(long, lat)], proj4string = CRS('+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'))
      listDataSP <- spTransform(listDataSP,CRS('+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'))
      fileNames <- lapply(fac, function(y) paste(cdir,"/",y,'/CA_BCM_HST_Monthly_',y,'_',listName,'.Rdata',sep=''))
      rass <- lapply(fileNames, function(z) readRDS(z)[[1]])
#       rass <- lapply(fileNames, function(z) tryCatch(readRDS(z)[[1]], 
#                                                      error = function(e) 
#                                                        setValues(errorRaster, rep(-99, ncell(errorRaster)))))
      rasStack <- stack(rass)
      extent(rasStack) <- good.extent
      vals <- raster::extract(rasStack, listDataSP)
      results <- cbind(listData, vals)
      names(results) <- c(names(listData), fac)
      return(results)
    }
    climvalsList <- lapply(1:length(ymoLists), extractFunc)
    climvals <- as.data.frame(rbindlist(climvalsList))
    #ll <- data.frame(ll,climvals)
    return(climvals)
}



#### Extract climate values from yearly 
extractClimateAnnually <- function(ll, colNames = c("long", "lat", "year_nm"), fac=c('tmn'), summaryType = "ave", cdir=getwd()) {
    ## ll is a dataframe that includes columns with longitude, latitude, year, and month, with column names specified by "colNames" statement
    ## summaryType is the kind of annual summary calculated from the monthly rasters, likely either 'min', 'max', or 'mean'
    require(sp)
    require(data.table)
    require(raster)
    long <- colNames[1]
    lat <- colNames[2]
    year_nm <- colNames[3]
    #month_nm <- colNames[4]
    dem <- raster('C:/Users/Adam/Documents/UC Berkeley post doc/BIGCB/Pest Project/Climate data/ca_270m_t6.asc')
    good.extent <- extent(dem)
    # Load a BCM raster to use in case of error
    # errorRaster <- readRDS(paste(cdir,'/tmx/CA_BCM_HST_Monthly_tmx_1920_01.Rdata',sep=''))[[1]]
    # Prepend single-digit months with 0s
    # ll[,month_nm] <- formatC(ll[,month_nm], width = 2, format = "d", flag = "0") 
    #ll$ymo <- paste(ll[,year_nm], ll[,month_nm], sep = "_")
    yrLists <- split(ll, ll[,year_nm])
#     if (FALSE) plot(llSP.aea) # I'm not clear what this is?
    lfac <- length(fac)
    #climvals <- matrix(NA,nrow=length(ymoLists),ncol=lfac)
    #i=1
    extractFunc <- function(x){
      listName <- names(yrLists)[x]
      listData <- yrLists[[x]]
      listDataSP <- listData
      coordinates(listDataSP) <- listDataSP[,c(long,lat)]
      projection(listDataSP) <- CRS('+proj=longlat +datum=WGS84')
      #class(llSP)
      #listDataSP <- SpatialPoints(listDataSP[,c(long, lat)], proj4string = CRS('+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'))
      listDataSP <- spTransform(listDataSP,CRS('+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'))
      fileNames <- lapply(fac, function(fac) paste(cdir,"/",fac,'/Annual summaries/BCM2014_',fac,listName,'_wy_',summaryType,'_HST.Rdata',sep=''))
      # fileNames <- lapply(fac, function(fac) paste(cdir,'/BCM2014_',fac,listName,'_wy_',summaryType,'_HST.Rdata',sep='')) # For files directly from Ackerly lab (yearly average)
      rass <- lapply(fileNames, function(z) readRDS(z)[[1]])
#       rass <- lapply(fileNames, function(z) tryCatch(readRDS(z)[[1]], 
#                                                      error = function(e) 
#                                                        setValues(errorRaster, rep(-99, ncell(errorRaster)))))
      rasStack <- stack(rass)
      extent(rasStack) <- good.extent
      vals <- tryCatch(raster::extract(rasStack, listDataSP),
                       error = function(e)
                         rep(NA, nrow(listDataSP)))
      results <- cbind(listData, vals)
      names(results) <- c(names(listData), fac)
      return(results)
    }
    climvalsList <- lapply(1:length(yrLists), extractFunc)
    climvals <- as.data.frame(rbindlist(climvalsList))
    #ll <- data.frame(ll,climvals)
    return(climvals)
}


#### Function to generate randomly sampled background values
xMonthlyRand <- function(ll, nRand=10,fac=c('ppt','tmn','tmx'),cdir=getwd()) {
    ## ll is a dataframe that includes 'long', 'lat', and 'year_nm', 'month_nm'
    require(raster)
    dem <- raster(CAsp)
    good.extent <- extent(dem)
    # Prepend months with a "0" if single digit
    ll$month_nm <- formatC(ll$month_nm, width = 2, format = "d", flag = "0") 
    ll$ymo <- paste(ll$year_nm, ll$month_nm, sep="_")
    # Make a list of data sets split by year and month
    ymoLists <- split(ll,ll$ymo)
    # Number of climate variables
    lfac <- length(fac)
    # Vectorized function to loop through year-month lists, load relevant rasters, and extract values
    randExtractFunc <- function(x){
      # Select and load relevant raster
      listName <- names(ymoLists)[x]
      fname <- paste(cdir,"/",fac.i,'/CA_BCM_HST_Monthly_',fac.i,'_',listName,'.Rdata',sep='')
      ras <- readRDS(fname)[[1]]
      extent(ras) <- good.extent
      rass <- stack(ras)
      # if multiple climate variables, select additional rasters
      if (lfac>1) {
        for (i in 2:lfac) {
          fname <- paste(cdir,"/",fac[i],'/CA_BCM_HST_Monthly_',fac[i],'_',listName,'.Rdata',sep='')
          ras.more <- readRDS(fname)[[1]]
          extent(ras.more) <- good.extent
          rass <- addLayer(rass,ras.more)
        }
      }
      # Get values from raster stack
      vals <- getValues(rass)
      rsel <- which(complete.cases(vals))
      vals <- vals[rsel,]
      # Randomly select points and extract values from each raster
      rSamp <- sample(nrow(vals),nRand)
      climvals <- vals[rSamp,]
      # Combine coordinate points and climate values into a data.frame
      xy <- xyFromCell(rass,1:length(rass[[1]]))
      xy <- xy[rsel,]
      xy <- xy[rSamp,]
      climvals <- data.frame(xy,climvals)
      names(climvals) <- c('long','lat',fac)
      return(climvals)
    }
    # Loop through all year-month lists
    climvals <- lapply(1:length(ymoLists), randExtractFunc)
    names(climvals) <- names(ymoLists)
    return(climvals)
}
