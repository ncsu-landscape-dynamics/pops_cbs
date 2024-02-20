library(terra)



# Path to file output location
outpath <- "/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/citrus_black_spot/"
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



# Florida boundary
florida <- vect("/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/citrus_black_spot/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
florida <- florida[florida$NAME == "Florida"]



# Create temp and prcip raster files for each year and day
for (year in seq(start_year, end_year)) {
  
  # Read raster files conditionally based on year
  if (year < 2018) {
    
    prcp <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/precip/daymet_v3_prcp_", year, "_na.nc4"))
    tmax <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/tmin/daymet_v3_tmin_", year, "_na.nc4"))
    tmin <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/tmax/daymet_v3_tmax_", year, "_na.nc4"))
    
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
  
  else {
      
      tmax <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/tmin/daymet_v4_daily_na_tmin_", year, ".nc"))
      tmin <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/tmax/daymet_v4_daily_na_tmax_", year, ".nc"))
      
      # Project florida onto prcp crs
      florida <- terra::project(florida, tmax)
      
      # Crop temp rasters to florida extent
      tmin <- terra::crop(tmin, florida)
      tmax <- terra::crop(tmax, florida)
      
      # Calculate mean temperature
      tmean <- (tmin + tmax) / 2
      
      # Apply temperature germination function to raster
      temp_values <- temp_fun(tmean)
      
      # Write raster files to outpath location
      writeRaster(temp_values, paste0(outpath, "temp/","temp_coeff_", year, "_.tif"), overwrite = TRUE)
    
  }
  
}



# Function to take create raster of host data
infection_rast <- function(x, cbs_host) {
  #' host_ext
  #'
  #' Description
  #' @param x spatRaster containing CBS occurrence data
  #' @param cbs_data spatVector containing spatial polygons of host species
  #' @return spatRaster mapped to extent of cbs hosts
  #ext(x) <- ext(cbs_host)
  host_data <- rasterize(cbs_host, x, field = "area", cover = TRUE)
  host_data <- host_data*100
  #host_data <- resample(host_data, x)
  return(host_data)
}



# Host raster
cbs <- read.csv("/Users/evandadson/Library/CloudStorage/GoogleDrive-edadson@ncsu.edu/Shared drives/Data/Original/pest-occurrence/Citrus_Black_Spot/MasterCBS2022.csv")
cbs_host <- vect("/Users/evandadson/Library/CloudStorage/GoogleDrive-edadson@ncsu.edu/Shared drives/Data/Original/pest-occurrence/Citrus_Black_Spot/CHRPMultiBlocks.shp")

cbs$Long = abs(cbs$Long)
cbs$count = 1
cbs_host$area = 1

for (year in seq(start_year, end_year)) {
  
  # Split based on year
  cbs_year <- cbs[grep(year, as.Date(cbs$Receive.Date, format = "%m/%d/%y")),]
  
  # Convert to raster with longlat referenced coordinates
  cbs_rast <- rast(cbs_year, type = "xzy", crs = crs(cbs_host), 
                   extent = ext(cbs_host), digits = 6)

  names(cbs_rast) = "count"
  values(cbs_rast) = 1
  cbs_host_rast <- infection_rast(cbs_rast, cbs_host)
  #cbs_inf <- infection_rast(cbs_rast, cbs_host)
  assign(paste0("cbs_", year), cbs_host_rast)
}