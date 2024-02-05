library(terra)



# Path to file output location
outpath <- "/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/citrus_black_spot/"
start_year <- 2011
end_year <- 2017
days_per_year <- 365



# Temperature germination function
temp_fun <- function(tmean) {
  temp_values <- ifelse(tmean > 12 & tmean < 32, 100*(-0.07+0.88*(exp(-0.5*(log(tmean/23.08)/0.28)^2))), 0)
  return(temp_values)
}



# Precipitation indicator function
prcp_fun <- function(prcp) {
  prcp_values <- ifelse(prcp < 0.25, 0, 1)
  return(prcp_values)
}



# Florida boundary
florida <- vect("/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/citrus_black_spot/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
florida <- florida[florida$NAME == "Florida"]



# Create temp and prcip raster files for each year and day
for (year in seq(start_year, end_year)) {
  
  for (day in 1:days_per_year) {
    
    # Read rasters
    prcp <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/precip/daymet_v3_prcp_", 2011, "_na.nc4"))[[18]]
    tmax <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/tmin/daymet_v3_tmin_", year, "_na.nc4"))[[day]]
    tmin <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/tmax/daymet_v3_tmax_", year, "_na.nc4"))[[day]]
    
    # Project florida onto prcp crs
    florida <- terra::project(florida, prcp)
    
    # Crop prcp raster to florida extent
    prcp <- terra::crop(prcp, florida)
    
    # Apply precipitation indicator function to raster
    prcp_values <- app(prcp, fun = prcp_fun)
    
    # Crop temp rasters to florida extent
    tmin <- terra::crop(tmin, florida)
    tmax <- terra::crop(tmax, florida)
    
    # Calculate mean temperature
    tmean <- (tmin + tmax) / 2
    
    # Apply temperature severity function to raster
    temp_values <- app(tmean, fun = temp_fun)
    
    # Write rasters
    writeRaster(prcp_values, paste0(outpath, "precip/","prcp_coeff_", year, "_", day, "_.tif"), overwrite = TRUE)
    writeRaster(temp_values, paste0(outtpath, "temp/","temp_coeff_", year, "_", day, "_.tif"), overwrite = TRUE)
    
  }
  
}


