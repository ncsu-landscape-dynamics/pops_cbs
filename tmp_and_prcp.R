library(terra)
library(sf)
library(raster)
library(doParallel)
library(foreach)

# Temperature germination everity function
tsv <- function(tmean) {
  tsv_values <- ifelse(tmean <= 12, 0, 
                            ifelse(tmean > 12 & tmean < 32, 100*(-0.07+0.88*(exp(-0.5*(log(tmean/23.08)/0.28)^2))),
                                        ifelse(tmean >= 32, 0,0)))
  return(tsv_values)
}

# Precipitation severity indicator function
psv <- function(prcp) {
  prcp_values <- ifelse(prcp < 0.25, 0, 1)
  return(prcp_values)
}

# Read in florida plss data
florida <- vect("/Users/evandadson/Downloads/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")

# outpath <- "/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/citrus_black_spot/"

# Testing Raster creation
test_out <- "/Users/evandadson/Desktop/test/"

files = list.files(path= '/Volumes/cmjone25/Data/Original/Daymet/precip', pattern="\\.nc4$", full.names = TRUE)
precip <- lapply(files, rast)

files = list.files(path= '/Volumes/cmjone25/Data/Original/Daymet/tmin', pattern="\\.nc4$", full.names = TRUE)
temp_min <- lapply(files, rast)

files = list.files(path= '/Volumes/cmjone25/Data/Original/Daymet/tmax', pattern="\\.nc4$", full.names = TRUE)
temp_max <- lapply(files, rast)
 
for (i in 32:38) {
    # prcp <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/precip/daymet_v3_prcp_", years, "_na.nc4"))[[i]]
    # tmax <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/tmin/daymet_v3_tmin_", years, "_na.nc4"))[[i]]
    # tmin <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/tmax/daymet_v3_tmax_", years, "_na.nc4"))[[i]]
  
  # Crop precip raster to florida extent
  prcp <- terra::crop(precip[[i]], florida)
  
  # Project prcp to the crs of florida
  prcp <- project(prcp, crs(florida))
  
  # Reclassify precipitation values of raster
  psv_values <- apply(values(prcp), MARGIN = c(1,2), FUN = psv)
  
  # Assign values of raster to reclassified values
  values(prcp) <- psv_values
  
  ## Temperature rasters
  tmin <- terra::crop(temp_min[[i]], florida)
  tmax <- terra::crop(temp_max[[i]], florida)

  tmin <- project(tmin, crs(florida))
  tmax <- project(tmax, crs(florida))
  
  tmean <- (values(tmin) + values(tmax)) / 2
  
  values(tmin) <- tmean
  
  tsv_values <- apply(tmean, MARGIN = c(1,2), FUN = tsv)
  
  values(tmin) <- tsv_values
  
  writeRaster(prcp, paste0(test_out, "temp_coeff_", i+1979, "_.tif"), overwrite = TRUE)
  writeRaster(tmin, paste0(test_out, "prcp_coeff_", i+1979, "_.tif"), overwrite = TRUE)
}
