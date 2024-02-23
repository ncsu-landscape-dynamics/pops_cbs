library(terra)

# Set up parameters
start_year <- 2010
end_year <- 2022

# Outpath
outpath <- "/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/citrus_black_spot/"

# Read occurrence and host data
cbs <- read.csv("/Volumes/cmjone25/Data/Original/pest-occurrence/Citrus_Black_Spot/MasterCBS2022.csv")
cbs_host <- vect("/Volumes/cmjone25/Data/Original/pest-occurrence/Citrus_Black_Spot/CHRPMultiBlocks.shp")

# Correct out of place long value
cbs$Long = -1*abs(cbs$Long)

# Initialize count column for unique occurrences
cbs$positive = 1
cbs_host$area = 1

# Plot overall infection over collection period against host data
plot(cbs_host, xlim=c(-88,-79), ylim=c(25,31), col="light yellow", border="light gray")
points(cbs$Long, cbs$Lat, cex = 0.5, col = 'red')

host_rast <- function(z, cbs_host) {
  #' host_ext
  #'
  #' Description
  #' @param z spatRaster
  #' @param cbs_data spatVector containing spatial polygons of host species
  #' @return spatRaster mapped to extent of cbs hosts
  cbs_host <- project(cbs_host, crs(z))
  host_data <- rasterize(cbs_host, z, field = "area", cover = TRUE)
  host_data <- host_data*100
  return(host_data)
}

# Set up input objects
weather <- rast("/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/citrus_black_spot/precip/prcp_coeff_2016_.tif")
host <- weather[[1]]
host_cbs <- host_rast(weather, cbs_host)
writeRaster(host_cbs, paste0(outpath, "host/host.tiff"))

# Yearly infection
for (year in seq(start_year, end_year)) {
  sp <- cbs[grep(year, as.Date(cbs$Receive.Date, format = "%m/%d/%y")),]
  sp1 <- vect(sp, geom = c("Long", "Lat"), crs=crs(cbs_host))
  host<-project(host,crs(sp1))
  infections_year <- terra::rasterize(sp1, host, field = "positive", fun = sum)
  infections_year <- project(infections_year, crs(weather))
  writeRaster(infections_year, paste0(outpath, "infection/cbs_", year, ".tiff"))
}

# Seasonal infection
for (season in 1:13) {
  sp <- cbs[cbs$Season == season,]
  sp1 <- vect(sp, geom = c("Long", "Lat"), crs = crs(cbs_host))
  host<-project(host,crs(sp1))
  infections_season <- terra::rasterize(sp1, host, field = "positive", fun = sum)
  infections_season <- project(infections_season, crs(weather))
  writeRaster(infections_year, paste0(outpath, "infection/cbs_", season, ".tiff"))
}
