library(terra)
# library(lubridate)

# Temperature germination severity function
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

# Florida shapefile
florida <- vect("/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/citrus_black_spot/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
florida <- florida[florida$NAME == "Florida"]

# Path to file output location
outpath <- "/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/citrus_black_spot/"

# Create temp and prcp rasters for each year and day
for (year in 2011:2017) {
  for (day in 1:365) {
    prcp <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/precip/daymet_v3_prcp_", year, "_na.nc4"))[[day]]
    tmax <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/tmin/daymet_v3_tmin_", year, "_na.nc4"))[[day]]
    tmin <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/tmax/daymet_v3_tmax_", year, "_na.nc4"))[[day]]
    
    # Project florida onto prcp crs
    florida <- terra::project(florida, prcp)
    
    # Crop prcp raster to florida extent
    prcp <- terra::crop(prcp, florida)
    
    # Apply precipitation indicator function to raster
    psv_values <- app(prcp, fun = psv)
    
    # Crop temp rasters to florida extent
    tmin <- terra::crop(tmin, florida)
    tmax <- terra::crop(tmax, florida)
    
    # Calculate mean temperature
    tmean <- (tmin + tmax) / 2
    
    # Apply temperature severity function to raster
    tsv_values <- app(tmean, fun = tsv)
    
    writeRaster(psv_values, paste0(test_out, "temp_coeffs_", year, "/", "temp_coeff_", year, "_", day, "_.tif"), overwrite = TRUE)
    writeRaster(tsv_values, paste0(test_out, "prcp_coeffs_", year, "/", "prcp_coeff_", year, "_", day, "_.tif"), overwrite = TRUE)
  }
}
