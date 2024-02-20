library(terra)

# Set up parameters
start_year <- 2010
end_year <- 2022

# Read occurrence and host data
cbs <- read.csv("/Users/evandadson/Library/CloudStorage/GoogleDrive-edadson@ncsu.edu/Shared drives/Data/Original/pest-occurrence/Citrus_Black_Spot/MasterCBS2022.csv")
cbs_host <- vect("/Users/evandadson/Library/CloudStorage/GoogleDrive-edadson@ncsu.edu/Shared drives/Data/Original/pest-occurrence/Citrus_Black_Spot/CHRPMultiBlocks.shp")

# Correct out of place long value
cbs$Long = -1*abs(cbs$Long)

# Initialize count column for unique occurrences
cbs$count = 1

# Plot overall infection over collection period against host data
plot(cbs_host, xlim=c(-88,-79), ylim=c(25,31), col="light yellow", border="light gray")
points(cbs$Long, cbs$Lat, cex = 0.5, col = 'red')

# Set up input objects
weather <- rast("-path")
host <- weather[1]
host_cbs <- host_rast(host, cbs_host)
writeRaster(host_cbs, "host.tiff")
# r <- rast(cbs_host, ncol = 100, nrow = 100)
# sp <- vect(cbs, geom = c("Long", "Lat"), crs = crs(cbs_host))
# new <- rasterize(sp, host, "count", fun = sum)

# Yearly infection
for (year in seq(start_year, end_year)) {
  sp <- cbs[grep(year, as.Date(cbs$Receive.Date, format = "%m/%d/%y")),]
  sp1 <- vect(sp, geom = c("Long", "Lat"), crs = crs(cbs_host))
  infections_year <- rasterize(sp1, host, "count", fun = sum)
  #assign(paste0("cbs_", year), infections_year)
  wrtieRaster(infections_year, paste0("cbs_", year, ".tiff"))
}

# Seasonal infection
for (season in 1:13) {
  sp <- cbs[cbs$Season == season,]
  sp1 <- vect(sp, geom = c("Long", "Lat"), crs = crs(cbs_host))
  infections_year <- rasterize(sp1, r, "count", fun = sum)
  assign(paste0("cbs_season_", season), infections_year)
}
