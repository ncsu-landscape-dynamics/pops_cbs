host_rast <- function(z, cbs_data) {
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