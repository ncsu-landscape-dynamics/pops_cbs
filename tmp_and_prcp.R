library(terra)
library(sf)

florida <- vect("/Users/evandadson/Downloads/plss_tr_2004/plss_tr_2004.shp", crs = "+proj=lcc +lat_0=42.5 +lon_0=-100 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")

# outpath <- "/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/citrus_black_spot/"

# Testing Raster creation
test_out <- "/Users/evandadson/Desktop/test"

for (year in 2011:2021) {
  prcp <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/precip/daymet_v3_prcp_", year, "_na.nc4"))
  tmin <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/tmin/daymet_v3_tmin_", year, "_na.nc4"))
  tmax <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/tmax/daymet_v3_tmax_", year, "_na.nc4"))
  
  prcp <- project(florida, prcp)
  tmin <- project(florida, tmin)
  tmax <- project(florida, tmax)
  
  prcp <- terra::crop(prcp, florida)
  tmin <- terra::crop(tmin, florida)
  tmax <- terra::crop(tmax, florida)
  
  tsv_values <- tsv(tmin, tmax)
  prcp_values <- psv(prcp)
  
  writeRaster(tsv_values, paste0(test_out, "temp_", year, "_", day, ".tif"))
  writeRaster(prcp_values, paste0(test_out, "prcp_", year, "_", day, ".tif"))
}
