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

# posterior means from previous calibration
prior_means <- read.csv(paste0(cbs_out, "posterior_means_2021.csv"))
prior_means <- prior_means$x
prior_means <- prior_means[1:6]

# posterior covariance matrix from previous calibration
prior_cov_matrix <- read.csv(paste0(cbs_out, "posterior_cov_matrix_2021.csv"))
prior_cov_matrix <- as.matrix(prior_cov_matrix)

# Calibration for PoPS model
PoPS::calibrate(
  infected_years_file = paste0(cbs_path, "inf_years_file13.tif"),
  number_of_observations = 18,
  prior_number_of_observations = 478,
  prior_means = prior_means,
  prior_cov_matrix = prior_cov_matrix,
  params_to_estimate = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  number_of_generations = 7,
  generation_size = 1000,
  pest_host_table = paste0(cbs_path, "pest_host_table_cbs.csv"),
  competency_table = paste0(cbs_path, "competency_table_cbs.csv"),
  infected_file_list = paste0(cbs_path, "infection/cbs_2013.tif"),
  host_file_list = paste0(cbs_path, "host/host.tif"),
  total_populations_file = paste0(cbs_path, "total_pops_file.tif"),
  temp = TRUE,
  temperature_coefficient_file = paste0(cbs_path, "temp/temp_coeff_2013_.tif"),
  precip = TRUE,
  precipitation_coefficient_file = paste0(cbs_path, "precip/prcp_coeff_2013_.tif"),
  model_type = "SI",
  latency_period = 0,
  time_step = 'day',
  season_month_start = 4,
  season_month_end = 9,
  start_date = "2022-01-01",
  end_date = "2022-12-31",
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
  treatment_dates = c('2022_04_01',
                      '2022_05_01',
                      '2022_06_01',
                      '2022_07_01',
                      '2022_08_01',
                      '2022_09_01'),
  treatments_file = c(paste0(cbs_path, "host/host.tiff"),
                      paste0(cbs_path, "host/host.tiff"),
                      paste0(cbs_path, "host/host.tiff"),
                      paste0(cbs_path, "host/host.tiff"),
                      paste0(cbs_path, "host/host.tiff"),
                      paste0(cbs_path, "host/host.tiff")),
  treatment_method = "ratio",
  natural_kernel_type = "cauchy",
  anthropogenic_kernel_type = "cauchy",
  natural_dir = "NONE",
  natural_kappa = 0,
  anthropogenic_dir = "NONE",
  anthropogenic_kappa = 0,
  pesticide_duration = c(30,30,30,30,30,30),
  pesticide_efficacy = 0.829,
  mask = NULL,
  output_frequency = "day",
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
  number_of_iterations = 1e+05,
  exposed_file_list = "",
  verbose = TRUE,
  write_outputs = "summary_outputs",
  output_folder_path = paste0(cbs_path, "outputs/"),
  network_filename = "",
  network_movement = "walk",
  success_metric = "mcc",
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

file_name <- paste0(cbs_out, "posterior_means_2022.csv")
write.csv(cal_2021$posterior_means, file_name, row.names = FALSE)

file_name <- paste0(cbs_out, "posterior_cov_matrix_2022.csv")
write.csv(cal_2021$posterior_cov_matrix, file_name, row.names = FALSE)