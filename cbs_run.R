# Jones, C., Jones, S., Petrasova, A., Petras, V., Gaydos, D., Skrip, M.,
# Takeuchi, Y., Bigsby, K., and Meentemeyer, R., 2021. Iteratively forecasting
# biological invasions with PoPS and a little help from our friends. Frontiers
# in Ecology and the Environment DOI: 10.1002/fee.2357

install.packages("remotes")
remotes::install_github("ncsu-landscape-dynamics/rpops")
library(PoPS)
library(terra)

cbs_path = "Z:/Data/Raster/USA/pops_casestudies/citrus_black_spot/"
cbs_out = "Z:/Data/Raster/USA/pops_casestudies/citrus_black_spot/outputs/"

# calibrated means and covariance matrices
for (year in seq(2010, 2021)) {
 cal_means <- read.csv(paste0(cbs_out, "posterior_means_", year, ".csv"))
 cal_cov <- read.csv(paste0(cbs_out, "posterior_cov_matrix_", year, ".csv"))
 assign(paste0("means", year), cal_means[[1]])
 assign(paste0("cov", year), cal_cov)
}

bayesian_mnn_checks <- function(prior_means,
                                prior_cov_matrix,
                                calibrated_means,
                                calibrated_cov_matrix,
                                prior_weight, weight) {
  checks_passed <- TRUE
  if (length(prior_means) == length(calibrated_means) && prior_weight > 0) {
    posterior_means <- prior_means * prior_weight + calibrated_means * weight
  } else if (prior_weight == 0) {
    posterior_means <- calibrated_means
  } else {
    checks_passed <- FALSE
    failed_check <- prior_means_error
  }
  
  if (nrow(prior_cov_matrix) == nrow(calibrated_cov_matrix) &&
      ncol(prior_cov_matrix) == ncol(calibrated_cov_matrix) &&
      prior_weight > 0) {
    posterior_cov_matrix <- prior_cov_matrix * prior_weight +
      calibrated_cov_matrix * weight
  } else if (prior_weight == 0) {
    posterior_cov_matrix <- calibrated_cov_matrix
  } else {
    checks_passed <- FALSE
    failed_check <- prior_cov_matrix_error
  }
  
  if (checks_passed) {
    outs <- list(checks_passed, posterior_means, posterior_cov_matrix)
    names(outs) <- c("checks_passed", "posterior_means", "posterior_cov_matrix")
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- failed_check_list
    return(outs)
  }
}

# Combined parameter calibration
cal10_11 <- bayesian_mnn_checks(means2010,
                                cov2010,
                                means2011,
                                cov2011,
                                95/(95+28), 28/(95+28))

cal11_12 <- bayesian_mnn_checks(cal10_11$posterior_means,
                                cal10_11$posterior_cov_matrix,
                                means2012,
                                cov2012,
                                123/(123+26), 26)

cal12_13 <- bayesian_mnn_checks(cal11_12$posterior_means,
                                cal11_12$posterior_cov_matrix,
                                means2013,
                                cov2013,
                                149, 39)

cal13_14 <- bayesian_mnn_checks(cal12_13$posterior_means,
                                cal12_13$posterior_cov_matrix,
                                means2014,
                                cov2014,
                                188, 33)

cal14_15 <- bayesian_mnn_checks(cal13_14$posterior_means,
                                cal13_14$posterior_cov_matrix,
                                means2015,
                                cov2015,
                                221, 25)

cal15_16 <- bayesian_mnn_checks(cal14_15$posterior_means,
                                cal14_15$posterior_cov_matrix,
                                means2016,
                                cov2016,
                                246, 9)

cal16_17 <- bayesian_mnn_checks(cal15_16$posterior_means,
                                cal15_16$posterior_cov_matrix,
                                means2017,
                                cov2017,
                                255, 13)

cal17_18 <- bayesian_mnn_checks(cal16_17$posterior_means,
                                cal16_17$posterior_cov_matrix,
                                means2018,
                                cov2018,
                                268, 36)

cal18_19 <- bayesian_mnn_checks(cal17_18$posterior_means,
                                cal17_18$posterior_cov_matrix,
                                means2019,
                                cov2019,
                                304, 83)

cal19_20 <- bayesian_mnn_checks(cal18_19$posterior_means,
                                cal18_19$posterior_cov_matrix,
                                means2020,
                                cov2020,
                                387, 57)

cal20_21 <- bayesian_mnn_checks(cal19_20$posterior_means,
                                cal19_20$posterior_cov_matrix,
                                means2021,
                                cov2021,
                                444, 18)

parameter_means = cal10_11$posterior_means
parameter_cov_matrix = cal10_11$posterior_cov_matrix

start_time <- Sys.time()

run_2011 <- pops_multirun(
  infected_file_list = paste0(cbs_path, "infection/cbs_2011.tif"),
  host_file_list = paste0(cbs_path, "host/host.tif"),
  total_populations_file = paste0(cbs_path, "total_pops_file.tif"),
  parameter_means,
  parameter_cov_matrix,
  pest_host_table = paste0(cbs_path, "pest_host_table_cbs.csv"),
  competency_table = paste0(cbs_path, "competency_table_cbs.csv"),
  temp = TRUE,
  temperature_coefficient_file = paste0(cbs_path, "temp/temp_coeff_2012.tif"),
  precip = TRUE,
  precipitation_coefficient_file = paste0(cbs_path, "precip/prcp_coeff_2012_.tif"),
  model_type = "SI",
  latency_period = 0,
  time_step = "day",
  season_month_start = 4,
  season_month_end = 9,
  start_date = "2012-01-01",
  end_date = "2012-12-31",
  use_survival_rates = FALSE,
  survival_rate_month = 3,
  survival_rate_day = 15,
  survival_rates_file = "",
  use_lethal_temperature = FALSE,
  temperature_file = "",
  lethal_temperature = -12.87,
  lethal_temperature_month = 1,
  mortality_frequency = "day",
  mortality_frequency_n = 1,
  management = TRUE,
  treatment_dates = "2012-04-01",
  treatments_file = paste0(cbs_path, "trt.tif"),
  treatment_method = "ratio",
  natural_kernel_type = "cauchy",
  anthropogenic_kernel_type = "cauchy",
  natural_dir = "NONE",
  anthropogenic_dir = "NONE",
  number_of_iterations = 10,
  number_of_cores = 7,
  pesticide_duration = 180,
  pesticide_efficacy = 0.829,
  random_seed = NULL,
  output_frequency = "year",
  output_frequency_n = 1,
  movements_file = "",
  use_movements = FALSE,
  start_exposed = FALSE,
  generate_stochasticity = TRUE,
  establishment_stochasticity = TRUE,
  movement_stochasticity = TRUE,
  dispersal_stochasticity = TRUE,
  establishment_probability = 0.5,
  dispersal_percentage = 0.99,
  quarantine_areas_file = "",
  use_quarantine = FALSE,
  use_spreadrates = FALSE,
  use_overpopulation_movements = FALSE,
  overpopulation_percentage = 0,
  leaving_percentage = 0,
  leaving_scale_coefficient = 1,
  exposed_file_list = "",
  mask = NULL,
  write_outputs = "None",
  output_folder_path = cbs_out,
  network_filename = "",
  network_movement = "walk",
  use_initial_condition_uncertainty = FALSE,
  use_host_uncertainty = FALSE,
  weather_type = "deterministic",
  temperature_coefficient_sd_file = "",
  precipitation_coefficient_sd_file = "",
  dispersers_to_soils_percentage = 0,
  quarantine_directions = "",
  multiple_random_seeds = FALSE,
  file_random_seeds = NULL,
  use_soils = FALSE,
  soil_starting_pest_file = "",
  start_with_soil_populations = FALSE,
  county_level_infection_data = FALSE
)

end_time <- Sys.time()
time_taken <- round(end_time-start_time, 2)
time_taken

file_name <- paste(cbs_out, "multirun_outputs.rdata", sep = "")
save(run_2011, file = file_name)
