library(terra)
library(lubridate)
library(dynamicSDM)

# Set up parameters
start_year <- 2010
end_year <- 2022

# Outpath
outpath <- "Z:/Data/Raster/USA/pops_casestudies/citrus_black_spot/"

# Read occurrence and host data
cbs <- read.csv("Z:/Data/Original/pest-occurrence/Citrus_Black_Spot/MasterCBS2022.csv")
cbs_host <- vect("Z:/Data/Original/pest-occurrence/Citrus_Black_Spot/CHRPMultiBlocks.shp")

# Correct out of place long value
cbs$Long = -1*abs(cbs$Long)
cbs$Receive.Date <- gsub('//', '/', cbs$Receive.Date)
# Initialize count column for unique occurrences
cbs$positive = 1
cbs_host$area = 1

cbs$year = year(as.Date(cbs$Receive.Date, format = "%m/%d/%y"))
cbs$month = month(as.Date(cbs$Receive.Date, format = "%m/%d/%y"))
cbs$day = day(as.Date(cbs$Receive.Date, format = "%m/%d/%y"))

# Plot overall infection over collection period against host data
plot(cbs_host, xlim=c(-88,-79), ylim=c(25,31), col="light yellow", border="light gray")
zoom(cbs_host)
points(positives$Long, positives$Lat, pch = 16, cex = 1, col = 'red')
points(ps_abs$y, ps_abs$x, pch = 16, cex = 1.5, col = 'green')

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
counties <- vect("Z:/Data/Vector/USA/us_lower_48_counties.gpkg")
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
  # infections from September to December
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

library(lubridate)

ps_abs <- dynamicSDM::spatiotemp_pseudoabs(spatial.method = "buffer", 
                                           temporal.method = "buffer", 
                                           occ.data = data.frame(cbind(x = cbs$Long, y = cbs$Lat, 
                                                                       year = cbs$year, 
                                                                       month = cbs$month, 
                                                                       day = cbs$day)),
                                           spatial.buffer = 5000, 
                                           n.pseudoabs = 497,
                                           temporal.ext = c("2010-01-01", "2022-12-31"), 
                                           temporal.buffer = 0,
                                           prj = "+proj=longlat +datum=WGS84")

# Plot overall infection over collection period against host data
plot(cbs_host, xlim=c(-88,-79), ylim=c(25,31), col="light yellow", border="light gray")
plot(cbs$Long, cbs$Lat, pch = 16, cex = 1, col = 'red')
points(ps_abs$x, ps_abs$y, pch = 16, cex = 1, col = 'green')

names(ps_abs)[1:2] <- c("Long", "Lat")
ps_abs$positive <- 0

cbs_true <- cbs[,c("Long", "Lat", "year", "month", "day", "positive")]
all_points <- rbind(cbs_true, ps_abs)

# Yearly infection
for (year in seq(start_year, end_year)) {
  # yearly infections
  sp <- all_points[all_points$year == year,]
  sp1 <- vect(sp, geom = c("Long", "Lat"), crs=crs(cbs_host))
  sp1 <- project(sp1, crs(host))
  infections_year <- terra::rasterize(sp1, host, field = "positive", fun = "sum")
  writeRaster(infections_year, paste0(outpath, "infection/cbs_pa_", year, ".tif"), overwrite = T)
  
  # Infections after treatment is applied
  inf_after_trt <- round(infections_year*0.171)
  inf_after_trt <- terra::ifel(inf_after_trt == 0, NA, inf_after_trt)
  writeRaster(inf_after_trt, paste0(outpath, "infection/inf_after_pa_", year, ".tif"), overwrite = T)
}
