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
cal11_12 <- bayesian_mnn_checks(means2010,
                                cov2010,
                                means2011,
                                cov2011,
                                97/(97+28), 28/(97+28))

cal11_13 <- bayesian_mnn_checks(cal11_12$posterior_means,
                                cal11_12$posterior_cov_matrix,
                                means2012,
                                cov2012,
                                125/(125+26), 26/(125+26))

cal11_14 <- bayesian_mnn_checks(cal11_13$posterior_means,
                                cal11_13$posterior_cov_matrix,
                                means2013,
                                cov2013,
                                151/(151+39), 39/(151+39))

cal11_15 <- bayesian_mnn_checks(cal11_14$posterior_means,
                                cal11_14$posterior_cov_matrix,
                                means2014,
                                cov2014,
                                190/(190+33), 33/(190+33))

cal11_16 <- bayesian_mnn_checks(cal11_15$posterior_means,
                                cal11_15$posterior_cov_matrix,
                                means2015,
                                cov2015,
                                223/(223+25), 25/(223+25))

cal11_17 <- bayesian_mnn_checks(cal11_16$posterior_means,
                                cal11_16$posterior_cov_matrix,
                                means2016,
                                cov2016,
                                248/(248+9), 9/(248+9))

cal11_18 <- bayesian_mnn_checks(cal11_17$posterior_means,
                                cal11_17$posterior_cov_matrix,
                                means2017,
                                cov2017,
                                257/(257+13), 13/(257+13))

cal11_19 <- bayesian_mnn_checks(cal11_18$posterior_means,
                                cal11_18$posterior_cov_matrix,
                                means2018,
                                cov2018,
                                270/(270+36), 36/(270+36))

cal11_20 <- bayesian_mnn_checks(cal11_19$posterior_means,
                                cal11_19$posterior_cov_matrix,
                                means2019,
                                cov2019,
                                306/(306+83), 83/(306+83))

cal11_21 <- bayesian_mnn_checks(cal11_20$posterior_means,
                                cal11_20$posterior_cov_matrix,
                                means2020,
                                cov2020,
                                389/(389+57), 57/(389+57))

cal11_22 <- bayesian_mnn_checks(cal11_21$posterior_means,
                                cal11_21$posterior_cov_matrix,
                                means2021,
                                cov2021,
                                446/(446+18), 18/(446+18))

parameter_means = cal11_22$posterior_means
parameter_cov_matrix = cal11_22$posterior_cov_matrix

# Forecasting without preventative management
start_time <- Sys.time()

cbs_nm <- pops_multirun(
  infected_file_list = paste0(cbs_path, "infection/cbs_2022.tif"),
  host_file_list = paste0(cbs_path, "host/host.tif"),
  total_populations_file = paste0(cbs_path, "total_pops_file.tif"),
  parameter_means,
  parameter_cov_matrix,
  pest_host_table = paste0(cbs_path, "pest_host_table_cbs.csv"),
  competency_table = paste0(cbs_path, "competency_table_cbs.csv"),
  temp = TRUE,
  temperature_coefficient_file = paste0(cbs_path, "temp/average_temp.tif"),
  precip = TRUE,
  precipitation_coefficient_file = paste0(cbs_path, "precip/average_precip.tif"),
  model_type = "SI",
  latency_period = 0,
  time_step = "day",
  season_month_start = 4,
  season_month_end = 9,
  start_date = "2023-01-01",
  end_date = "2050-12-31",
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
  management = FALSE,
  treatment_dates = "2023-04-01",
  treatments_file = paste0(cbs_path, "trt.tif"),
  treatment_method = "ratio",
  natural_kernel_type = "cauchy",
  anthropogenic_kernel_type = "cauchy",
  natural_dir = "NONE",
  anthropogenic_dir = "NONE",
  number_of_iterations = 100,
  number_of_cores = 26,
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
  weather_type = "probabilistic",
  temperature_coefficient_sd_file = paste0(cbs_path, "temp/sd_temp.tif"),
  precipitation_coefficient_sd_file = paste0(cbs_path, "precip/sd_precip.tif"),
  dispersers_to_soils_percentage = 0,
  quarantine_directions = "",
  multiple_random_seeds = FALSE,
  file_random_seeds = NULL,
  use_soils = FALSE,
  soil_starting_pest_file = "",
  start_with_soil_populations = FALSE,
  county_level_infection_data = FALSE
)

file_name <- paste0(cbs_out, "pops_runs/no_management/", "simulation_probability.tif")
writeRaster(cbs_nm$probability, file_name, overwrite=T)
file_name <- paste0(cbs_out, "pops_runs/no_management/", "simulation_mean.tif")
writeRaster(cbs_nm$simulation_mean, file_name, overwrite=T)
file_name <- paste0(cbs_out, "pops_runs/no_management/", "simulation_sd.tif")
writeRaster(cbs_nm$simulation_sd, file_name, overwrite=T)
file_name <- paste0(cbs_out, "pops_runs/no_management/", "simulation_min.tif")
writeRaster(cbs_nm$simulation_min, file_name, overwrite=T)
file_name <- paste0(cbs_out, "pops_runs/no_management/", "simulation_max.tif")
writeRaster(cbs_nm$simulation_max, file_name, overwrite=T)
file_name <- paste0(cbs_out, "pops_runs/no_management/", "multirun_outputs.rdata")
save(cbs_nm, file = file_name)

for (yr in seq(22, 49)) {
  if (yr == 22) {
    infected_file_list = paste0(cbs_path, "infection/cbs_2022.tif")
    use_initial_condition_uncertainty = FALSE
  }
  else {
    infected_file_list = paste0(cbs_out, "pops_runs/manage/", "simulation_uncert", yr, ".tif")
    use_initial_condition_uncertainty = TRUE
  }
  cbs_manage <- pops_multirun(
    infected_file_list,
    host_file_list = paste0(cbs_path, "host/host.tif"),
    total_populations_file = paste0(cbs_path, "total_pops_file.tif"),
    parameter_means,
    parameter_cov_matrix,
    pest_host_table = paste0(cbs_path, "pest_host_table_cbs.csv"),
    competency_table = paste0(cbs_path, "competency_table_cbs.csv"),
    temp = TRUE,
    temperature_coefficient_file = paste0(cbs_path, "temp/average_temp.tif"),
    precip = TRUE,
    precipitation_coefficient_file = paste0(cbs_path, "precip/average_precip.tif"),
    model_type = "SI",
    latency_period = 0,
    time_step = "day",
    season_month_start = 4,
    season_month_end = 9,
    start_date = paste0("20", yr+1, "-01-01"),
    end_date = paste0("20", yr+1, "-12-31"),
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
    treatment_dates = paste0("20", yr+1, "-04-01"),
    treatments_file = paste0(cbs_path, "trt.tif"),
    treatment_method = "ratio",
    natural_kernel_type = "cauchy",
    anthropogenic_kernel_type = "cauchy",
    natural_dir = "NONE",
    anthropogenic_dir = "NONE",
    number_of_iterations = 100,
    number_of_cores = 14,
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
    use_initial_condition_uncertainty,
    use_host_uncertainty = FALSE,
    weather_type = "probabilistic",
    temperature_coefficient_sd_file = paste0(cbs_path, "temp/sd_temp.tif"),
    precipitation_coefficient_sd_file = paste0(cbs_path, "precip/sd_precip.tif"),
    dispersers_to_soils_percentage = 0,
    quarantine_directions = "",
    multiple_random_seeds = FALSE,
    file_random_seeds = NULL,
    use_soils = FALSE,
    soil_starting_pest_file = "",
    start_with_soil_populations = FALSE,
    county_level_infection_data = FALSE
  )
  simulation_uncert <- c(cbs_manage$simulation_mean, cbs_manage$simulation_sd)
  names(simulation_uncert) <- c("sim_mean", "sim_sd")
  file_name <- paste0(cbs_out, "pops_runs/manage/", "simulation_uncert", yr+1, ".tif")
  writeRaster(simulation_uncert, file_name, overwrite=T)
  file_name <- paste0(cbs_out, "pops_runs/manage/", "simulation_probability", yr+1, ".tif")
  writeRaster(cbs_manage$probability, file_name, overwrite=T)
  file_name <- paste0(cbs_out, "pops_runs/manage/", "simulation_mean", yr+1, ".tif")
  writeRaster(cbs_manage$simulation_mean, file_name, overwrite=T)
  file_name <- paste0(cbs_out, "pops_runs/manage/", "simulation_sd", yr+1, ".tif")
  writeRaster(cbs_manage$simulation_sd, file_name, overwrite=T)
  file_name <- paste0(cbs_out, "pops_runs/manage/", "simulation_min", yr+1, ".tif")
  writeRaster(cbs_manage$simulation_min, file_name, overwrite=T)
  file_name <- paste0(cbs_out, "pops_runs/manage/", "simulation_max", yr+1, ".tif")
  writeRaster(cbs_manage$simulation_max, file_name, overwrite=T)
  file_name <- paste0(cbs_out, "pops_runs/manage/", "multirun_outputs", yr+1, ".rdata")
  save(cbs_manage, file = file_name)
}

end_time <- Sys.time()
time_taken <- round(end_time-start_time, 2)
time_taken
