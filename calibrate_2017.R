# Jones, C., Jones, S., Petrasova, A., Petras, V., Gaydos, D., Skrip, M.,
# Takeuchi, Y., Bigsby, K., and Meentemeyer, R., 2021. Iteratively forecasting
# biological invasions with PoPS and a little help from our friends. Frontiers
# in Ecology and the Environment DOI: 10.1002/fee.2357

# install.packages("remotes")
# remotes::install_github("ncsu-landscape-dynamics/rpops")
library(PoPS)
library(terra)

cbs_path = "Z:/Data/Raster/USA/pops_casestudies/citrus_black_spot/"
cbs_out = "Z:/Data/Raster/USA/pops_casestudies/citrus_black_spot/outputs/"

# load calibration outputs
load(paste0(cbs_out, "calibration_outputs_2016.rdata"))
prior_means <- cal_2016$posterior_means
prior_cov_matrix <- cal_2016$posterior_cov_matrix

start_time <- Sys.time()

# Calibration for PoPS model
cal_2017 <- PoPS::calibrate(
  infected_years_file = paste0(cbs_path, "infection/cbs_2018.tif"),
  number_of_observations = 13,
  prior_number_of_observations = 257,
  prior_means = prior_means,
  prior_cov_matrix = prior_cov_matrix,
  params_to_estimate = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  number_of_generations = 7,
  generation_size = 1000,
  pest_host_table = paste0(cbs_path, "pest_host_table_cbs.csv"),
  competency_table = paste0(cbs_path, "competency_table_cbs.csv"),
  infected_file_list = paste0(cbs_path, "infection/cbs_2017.tif"),
  host_file_list = paste0(cbs_path, "host/host.tif"),
  total_populations_file = paste0(cbs_path, "total_pops_file.tif"),
  temp = TRUE,
  temperature_coefficient_file = paste0(cbs_path, "temp/temp_coeff_2018.tif"),
  precip = TRUE,
  precipitation_coefficient_file = paste0(cbs_path, "precip/prcp_coeff_2018.tif"),
  model_type = "SI",
  latency_period = 0,
  time_step = 'day',
  season_month_start = 4,
  season_month_end = 9,
  start_date = "2018-01-01",
  end_date = "2018-12-31",
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
  treatment_dates = "2018-04-01",
  treatments_file = paste0(cbs_path, "trt.tif"),
  treatment_method = "ratio",
  natural_kernel_type = "cauchy",
  anthropogenic_kernel_type = "cauchy",
  natural_dir = "NONE",
  natural_kappa = 0,
  anthropogenic_dir = "NONE",
  anthropogenic_kappa = 0,
  pesticide_duration = 180,
  pesticide_efficacy = 0.829,
  mask = NULL,
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
  calibration_method = "ABC",
  number_of_iterations = 1e+06,
  exposed_file_list = "",
  verbose = TRUE,
  write_outputs = c("summary_outputs"),
  output_folder_path = cbs_out,
  network_filename = "",
  network_movement = "walk",
  success_metric = "rmse",
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
time_taken <- round(end_time - start_time, 2)
time_taken

file_name <- paste0(cbs_out, "posterior_means_2017.csv")
write.csv(cal_2017$posterior_means, file_name, row.names = FALSE)

file_name <- paste0(cbs_out, "posterior_cov_matrix_2017.csv")
write.csv(cal_2017$posterior_cov_matrix, file_name, row.names = FALSE)

file_name <- paste(cbs_out, "calibration_outputs_2017.rdata", sep = "")
save(cal_2017, file = file_name)
