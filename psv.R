psv <- function(prcp) {
  prcp_values <- terra::ifel(prcp < 0.25, 0, 1)
}