host_ext <- function(z, cbs_host) {
  #' host_ext
  #'
  #' Description
  #' @param z spatRaster
  #' @param host_data spatVector containing spatial polygons of host species
  #' @return spatRaster mapped to extent of cbs hosts
  cbs_host <- terra::project(cbs_host, crs(z))
  z <- terra::crop(z, ext(cbs_host))
  z <- terra::mask(z, cbs_host)
  return(z)
}