temp_fun <- function(tmean) {
  #' temp_fun
  #' 
  #' 
  #' @param tmean spatRaster of mean temperature values
  #' @return spatRaster of temperature germination results (Er et al., 2013)
  temp_values <- terra::ifel(tmean > 4 & tmean < 37, round(floor((625000/14473))*(0.00000009868*(tmean-4)^2.9212*(37-tmean)^1.3909),2), 0)
  return(temp_values)
}