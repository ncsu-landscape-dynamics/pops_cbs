temp_fun <- function(tmean) {
  #' temp_fun
  #' 
  #' 
  #' @param tmean spatRaster of mean temperature values
  #' @return spatRaster of temperature germination results (Wang et al., 2014)
  temp_values <- terra::ifel(tmean > 12 & tmean < 32, 100*(-0.07+0.88*(exp(-0.5*(log(tmean/23.08)/0.28)^2))), 0)
  return(temp_values)
}