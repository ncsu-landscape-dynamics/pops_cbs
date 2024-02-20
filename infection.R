library(terra)

start_year <- 2010
end_year <- 2022

cbs <- read.csv("/Users/evandadson/Library/CloudStorage/GoogleDrive-edadson@ncsu.edu/Shared drives/Data/Original/pest-occurrence/Citrus_Black_Spot/MasterCBS2022.csv")
cbs_host <- vect("/Users/evandadson/Library/CloudStorage/GoogleDrive-edadson@ncsu.edu/Shared drives/Data/Original/pest-occurrence/Citrus_Black_Spot/CHRPMultiBlocks.shp")

cbs$Long = -1*abs(cbs$Long)
cbs$count = 1
cbs_host$area = 1

plot(cbs_host, xlim=c(-88,-79), ylim=c(25,31), col="light yellow", border="light gray")
points(cbs$Long, cbs$Lat, cex = 0.5, col = 'red')


r <- rast(cbs_host, ncol = 100, nrow = 100)
sp <- vect(cbs, geom = c("Long", "Lat"), crs = crs(cbs_host))
new <- rasterize(sp, r, "count", fun = sum)


r <- rast(cbs_host, ncol = 100, nrow = 100)
for (year in seq(start_year, end_year)) {
  sp <- cbs[grep(year, as.Date(cbs$Receive.Date, format = "%m/%d/%y")),]
  sp1 <- vect(sp, geom = c("Long", "Lat"), crs = crs(cbs_host))
  infections_year <- rasterize(sp1, r, "count", fun = sum)
  assign(paste0("cbs_", year), infections_year)
}
