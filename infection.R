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
weather <- rast(paste0(outpath, "precip/prcp_coeff_2016.tif"))
host <- weather[[1]]
host_cbs <- host_rast(weather, cbs_host)
writeRaster(host_cbs, paste0(outpath, "host/host.tif"), overwrite = T)

# Test abandoned acreage assuming 70% remains in infected counties
counties <- vect("/Volumes/cmjone25/Data/Vector/USA/us_lower_48_counties.gpkg")
counties <- counties[counties$STATE_NAME == "Florida"]
counties <- project(counties, crs(host_cbs))
counties <- counties[counties$NAME == "Brevard" | counties$NAME == "Charlotte" | counties$NAME == "Collier" | counties$NAME == "De Soto" | counties$NAME == "Glades" | counties$NAME == "Hardee" | counties$NAME == "Hendry" | counties$NAME == "Hernando" | counties$NAME == "Highlands" | counties$NAME == "Hillsborough" | counties$NAME == "Indian River" | counties$NAME == "Lake" | counties$NAME == "Lee" | counties$NAME == "Manatee" | counties$NAME == "Marion" | counties$NAME == "Martin" | counties$NAME == "Okeechobee" | counties$NAME == "Orange" | counties$NAME == "Osceola" | counties$NAME == "Pasco" | counties$NAME == "Polk" | counties$NAME == "Putnam" | counties$NAME == "St. Lucie" | counties$NAME == "Sarasota" | counties$NAME == "Seminole" | counties$NAME == "Volusia"]
mask <- rasterize(counties, host_cbs, field=1)
adjusted_raster <- ifel(mask == 1, host_cbs * 0.649967327, host_cbs)
writeRaster(adjusted_raster, paste0(outpath, "host/abandoned.tif"), overwrite = T)
percent_host <- c(1.126334132, 1.096928774, 1.078414289, 1.066434328, 1.037464605, 1, 0.948159442, 0.894576345, 0.873230233, 0.843171422, 0.828795469, 0.804835548, 0.741015029, 0.649967327)

# Assign abandoned grove and total populations files
for (year in seq(start_year, end_year)) {
  # abandoned grove data
  adjusted_raster <- ifel(mask == 1, host_cbs*percent_host[year-2009], host_cbs)
  adjusted_raster <- ifel(mask == 1 & adjusted_raster > 100, 100, host_cbs)
  writeRaster(adjusted_raster, paste0(outpath, "host/abandoned", year, ".tif"), overwrite = T)
  
  total_pops_file = adjusted_raster
  total_pops_file = 100*total_pops_file*(1/total_pops_file)
  writeRaster(total_pops_file, paste0(outpath, "total_pops_", year, ".tif"), overwrite = T)
}

months_to_check <- c("September", "October", "November", "December")
pattern <- paste(months_to_check, collapse = "|")

for (year in seq(start_year, end_year)) {
  # infections from september to december
  sp <- cbs[grep(year, as.Date(cbs$Receive.Date, format = "%m/%d/%y")),]
  sp1 <- sp[grep(pattern, months(as.Date(sp$Receive.Date, format = "%m/%d/%y"))), ]
  sp2 <- vect(sp1, geom = c("Long", "Lat"), crs=crs(cbs_host))
  sp2 <- project(sp2, crs(host))
  infections_year <- terra::rasterize(sp2, host, field = "positive", fun = sum)
  writeRaster(infections_year, paste0(outpath, "infection/inf_after_sep_", year, ".tif"), overwrite = T)
}

# Yearly infection
for (year in seq(start_year, end_year)) {
  # yearly infections
  sp <- cbs[grep(year, as.Date(cbs$Receive.Date, format = "%m/%d/%y")),]
  sp1 <- vect(sp, geom = c("Long", "Lat"), crs=crs(cbs_host))
  sp1 <- project(sp1, crs(host))
  infections_year <- terra::rasterize(sp1, host, field = "positive", fun = sum)
  writeRaster(infections_year, paste0(outpath, "infection/cbs_", year, ".tif"), overwrite = T)
  
  # Infections after treatment is applied
  inf_after_trt <- round(infections_year*0.171)
  inf_after_trt <- terra::ifel(inf_after_trt == 0, NA, inf_after_trt)
  writeRaster(inf_after_trt, paste0(outpath, "infection/inf_after_", year, ".tif"), overwrite = T)
}

# Seasonal infection
for (season in 1:13) {
  sp <- cbs[cbs$Season == season,]
  sp1 <- vect(sp, geom = c("Long", "Lat"), crs = crs(cbs_host))
  sp1 <- project(sp1, crs(host))
  infections_season <- terra::rasterize(sp1, host, field = "positive", fun = sum)
  writeRaster(infections_year, paste0(outpath, "infection/cbs_", season, ".tiff"),
              overwrite = T)
  
  # Infections after treatment is applied
  inf_after_trt <- round(infections_season*0.171)
  inf_after_trt <- terra::ifel(inf_after_trt == 0, NA, inf_after_trt)
  writeRaster(inf_after_trt, paste0(outpath, "infection/season/inf_after_", season, ".tif"), overwrite = T)
}

