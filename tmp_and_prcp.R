library(terra)



# Path to file output location
outpath <- "Z:/Data/Raster/USA/pops_casestudies/citrus_black_spot/"
start_year <- 2010
end_year <- 2022



# Temperature germination function
temp_fun <- function(tmean) {
  temp_values <- terra::ifel(tmean > 12 & tmean < 32, 100*(-0.07+0.88*(exp(-0.5*(log(tmean/23.08)/0.28)^2))), 0)
  return(temp_values)
}



# Precipitation indicator function
prcp_fun <- function(prcp) {
  prcp_values <- terra::ifel(prcp < 0.25, 0, 1)
  return(prcp_values)
}



# US Census. 2019. “TIGER/Line Shapefiles.” Accessed on Feb 22, 2024.
florida <- vect(paste0(outpath, "cb_2018_us_state_20m/cb_2018_us_state_20m.shp"))
florida <- florida[florida$NAME == "Florida"]



# Create temp and prcip raster files for each year and day
for (year in seq(start_year, end_year)) {
  
  # Read raster files conditionally based on year
  if (year < 2018) {
    
    prcp <- rast(paste0("Z:/Data/Original/Daymet/precip/daymet_v3_prcp_", year, "_na.nc4"))
    tmin <- rast(paste0("Z:/Data/Original/Daymet/tmin/daymet_v3_tmin_", year, "_na.nc4"))
    tmax <- rast(paste0("Z:/Data/Original/Daymet/tmax/daymet_v3_tmax_", year, "_na.nc4"))
    
    # Project florida onto prcp crs
    florida <- terra::project(florida, prcp)
    
    # Crop prcp raster to florida extent
    prcp <- terra::crop(prcp, florida)
    
    # Apply precipitation indicator function to raster
    prcp_values <- prcp_fun(prcp)
    
    # Project florida onto prcp crs
    florida <- terra::project(florida, tmin)
    
    # Crop temp rasters to florida extent
    tmin <- terra::crop(tmin, florida)
    tmax <- terra::crop(tmax, florida)
    
    # Calculate mean temperature
    tmean <- (tmin + tmax) / 2
    
    # Apply temperature germination function to raster
    temp_values <- temp_fun(tmean)
    
    # Write raster files to outpath location
    writeRaster(prcp_values, paste0(outpath, "precip/","prcp_coeff_", year, "_.tif"), overwrite = TRUE)
    writeRaster(temp_values, paste0(outpath, "temp/","temp_coeff_", year, "_.tif"), overwrite = TRUE)
    
  }
  
  else {
    
    prcp <- rast(paste0("Z:/Data/Original/Daymet/precip/daymet_v4_daily_na_prcp_", year, ".nc"))
    tmax <- rast(paste0("Z:/Data/Original/Daymet/tmin/daymet_v4_daily_na_tmin_", year, ".nc"))
    tmin <- rast(paste0("Z:/Data/Original/Daymet/tmax/daymet_v4_daily_na_tmax_", year, ".nc"))
    
    # Project florida onto prcp crs
    florida <- terra::project(florida, prcp)
    
    # Crop prcp raster to florida extent
    prcp <- terra::crop(prcp, florida)
    
    # Apply precipitation indicator function to raster
    prcp_values <- prcp_fun(prcp)
    
    # Crop temp rasters to florida extent
    tmin <- terra::crop(tmin, florida)
    tmax <- terra::crop(tmax, florida)
    
    # Calculate mean temperature
    tmean <- (tmin + tmax) / 2
    
    # Apply temperature germination function to raster
    temp_values <- temp_fun(tmean)
    
    # Write raster files to outpath location
    writeRaster(prcp_values, paste0(outpath, "precip/","prcp_coeff_", year, "_.tif"), overwrite = TRUE)
    writeRaster(temp_values, paste0(outpath, "temp/","temp_coeff_", year, "_.tif"), overwrite = TRUE)

    
  }
  
}

start_year_2 = 2018

for (year in seq(start_year_2, end_year)) {
  prcp <- rast(paste0("Z:/Data/Original/Daymet/precip/daymet_v4_daily_na_prcp_", year, ".nc"))
  tmax <- rast(paste0("Z:/Data/Original/Daymet/tmin/daymet_v4_daily_na_tmin_", year, ".nc"))
  tmin <- rast(paste0("Z:/Data/Original/Daymet/tmax/daymet_v4_daily_na_tmax_", year, ".nc"))
  
  # Project florida onto prcp crs
  florida <- terra::project(florida, prcp)
  
  # Crop prcp raster to florida extent
  prcp <- terra::crop(prcp, florida)
  
  # Apply precipitation indicator function to raster
  prcp_values <- prcp_fun(prcp)
  
  # Crop temp rasters to florida extent
  tmin <- terra::crop(tmin, florida)
  tmax <- terra::crop(tmax, florida)
  
  # Calculate mean temperature
  tmean <- (tmin + tmax) / 2
  
  # Apply temperature germination function to raster
  temp_values <- temp_fun(tmean)
  
  # Write raster files to outpath location
  writeRaster(prcp_values, paste0(outpath, "precip/","prcp_coeff_", year, ".tif"), overwrite = TRUE)
  writeRaster(temp_values, paste0(outpath, "temp/","temp_coeff_", year, ".tif"), overwrite = TRUE)
}


# Temperature values between 0 and 1
for (year in seq(start_year, end_year)) {
  temp_file <- rast(paste0(outpath, "temp/","temp_coeff_", year, ".tif"))
  temp_file <- abs(temp_file/100)
  writeRaster(temp_file, paste0(outpath, "temp/","temp_coeff_", year, ".tif"), overwrite = TRUE)
}




