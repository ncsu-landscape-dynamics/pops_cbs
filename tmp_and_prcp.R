library(terra)
source("temp_fun.R")


# Path to file output location
outpath <- "Z:/Data/Raster/USA/pops_casestudies/citrus_black_spot/"
start_year <- 2010
end_year <- 2022

wang_temp <- function(tmean){return(-0.07+0.88*(exp(-0.5*(log(tmean/23.08)/0.28)^2)))}

# Temperature germination function
temp_fun <- function(tmean) {
  temp_values <- terra::ifel(tmean > 12.3 & tmean < 32,(1/wang_temp(23.08))*wang_temp(tmean), 0)
  return(temp_values)
}



# Precipitation indicator function
prcp_fun <- function(prcp) {
  prcp_values <- terra::ifel(prcp < 2.5, 0, 1)
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
    writeRaster(prcp_values, paste0(outpath, "precip/","prcp_coeff_", year, ".tif"), overwrite = TRUE)
    writeRaster(temp_values, paste0(outpath, "temp/","temp_coeff_", year, ".tif"), overwrite = TRUE)
    
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
    writeRaster(prcp_values, paste0(outpath, "precip/","prcp_coeff_", year, ".tif"), overwrite = TRUE)
    writeRaster(temp_values, paste0(outpath, "temp/","temp_coeff_", year, ".tif"), overwrite = TRUE)

    
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

cbs_path = "Z:/Data/Raster/USA/pops_casestudies/citrus_black_spot/"
cbs_out = "Z:/Data/Raster/USA/pops_casestudies/citrus_black_spot/outputs/"

# Weather forecasting
for (year in seq(2010, 2022)) {
  temp <- rast(paste0(cbs_path, "temp/temp_coeff_", year, ".tif"))
  precip <- rast(paste0(cbs_path, "precip/prcp_coeff_", year, ".tif"))
  assign(paste0("temp", year), temp)
  assign(paste0("precip", year), precip)
}
# Average temperature files
crs(temp2010)<-crs(temp2011)<-crs(temp2012)<-crs(temp2013)<-crs(temp2014)<-crs(temp2015)<-crs(temp2016)<-crs(temp2017)<-crs(temp2018)<-crs(temp2019)<-crs(temp2020)<-crs(temp2021)<-crs(temp2022)
names(temp2010)<-names(temp2011)<-names(temp2012)<-names(temp2013)<-names(temp2014)<-names(temp2015)<-names(temp2016)<-names(temp2017)<-names(temp2018)<-names(temp2019)<-names(temp2020)<-names(temp2021)<-names(temp2022)
temp_mean <- mosaic(temp2010, temp2011, temp2012, temp2013, temp2014, temp2015, temp2016, 
                    temp2017, temp2018, temp2019, temp2020, temp2021, temp2022, fun="mean")
temp_sd <- stdev(temp2010, temp2011, temp2012, temp2013, temp2014, temp2015, temp2016, 
                 temp2017, temp2018, temp2019, temp2020, temp2021, temp2022)
names(temp_mean) <- paste0("tmean_", seq(1,365))
names(temp_sd) <- paste0("tstdev_", seq(1,365))
writeRaster(temp_mean, paste0(cbs_path, "temp/average_temp.tif"), overwrite=F)
writeRaster(temp_sd, paste0(cbs_path, "temp/sd_temp.tif"), overwrite=F)

# Average precip files
crs(precip2010)<-crs(precip2011)<-crs(precip2012)<-crs(precip2013)<-crs(precip2014)<-crs(precip2015)<-crs(precip2016)<-crs(precip2017)<-crs(precip2018)<-crs(precip2019)<-crs(precip2020)<-crs(precip2021)<-crs(precip2022)
names(precip2010)<-names(precip2011)<-names(precip2012)<-names(precip2013)<-names(precip2014)<-names(precip2015)<-names(precip2016)<-names(precip2017)<-names(precip2018)<-names(precip2019)<-names(precip2020)<-names(precip2021)<-names(precip2022)
precip_mean <- mosaic(precip2010, precip2011, precip2012, precip2013, precip2014, precip2015, precip2016, 
                    precip2017, precip2018, precip2019, precip2020, precip2021, precip2022, fun="mean")
precip_sd <- stdev(precip2010, precip2011, precip2012, precip2013, precip2014, precip2015, precip2016, 
                 precip2017, precip2018, precip2019, precip2020, precip2021, precip2022)
names(precip_mean) <- paste0("pmean_", seq(1,365))
names(precip_sd) <- paste0("pstdev_", seq(1,365))
time(precip_mean) <- NULL
time(precip_sd) <- NULL
writeRaster(precip_mean, paste0(cbs_path, "precip/average_precip.tif"), overwrite=F)
writeRaster(precip_sd, paste0(cbs_path, "precip/sd_precip.tif"), overwrite=F)

