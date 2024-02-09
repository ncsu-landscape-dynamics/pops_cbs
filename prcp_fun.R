prcp_fun <- function(prcp) {
  #' prcp_fun
  #' 
  #' 
  #' @param prcp spatRaster of precipitation values in mm/day
  #' @return spatRaster of indicator (1 = conductive for infection, 0 = not)
  prcp_values <- terra::ifel(prcp < 0.25, 0, 1)
  return(prcp_values)
}