tsv <- function(tmin, tmax) {
  tmean = (tmin + tmax) / 2
  tsv_values <- terra::ifel(tmean <= 12, 0, 
                           terra::ifel(tmean > 12 & tmean < 32, 100*(-0.07+0.88*(exp(-0.5*(log(tmean/23.08)/0.28)^2))),
                                       terra::ifel(tmean >= 32, 0,0)))
}