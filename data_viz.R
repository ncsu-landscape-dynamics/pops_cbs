# Data Visualizations for Calibration Data
library(ggplot2)

cbs_path = "Z:/Data/Raster/USA/pops_casestudies/citrus_black_spot/"
cbs_out = "Z:/Data/Raster/USA/pops_casestudies/citrus_black_spot/outputs/"

load(paste0(cbs_out, "calibration_outputs_2010.rdata"))

posterior_means <- cal_2010$posterior_means
raw_calibration_data <- as.data.frame(cal_2010$raw_calibration_data)

# Histogram with density plot
ggplot(raw_calibration_data, aes(x=V1)) + 
  geom_density(alpha=.2, fill="blue") +
  geom_vline(aes(xintercept=posterior_means[1])) +
  ggtitle("Distribution of Reproductive Rate")

# Histogram with density plot
ggplot(raw_calibration_data, aes(x=V2)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white") +
  geom_density(alpha=.2, fill="blue") +
  geom_vline(aes(xintercept=posterior_means[2])) +
  ggtitle("Distribution of Natural Dispersal Distance")

# Histogram with density plot
ggplot(raw_calibration_data, aes(x=V3)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white") +
  geom_density(alpha=.2, fill="blue") +
  geom_vline(aes(xintercept=posterior_means[3])) +
  ggtitle("Distribution of Percent Natural Dispersal")

# Histogram with density plot
ggplot(raw_calibration_data, aes(x=V4)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white") +
  geom_density(alpha=.2, fill="blue") +
  geom_vline(aes(xintercept=posterior_means[4])) +
  ggtitle("Distribution of Anthropogenic Dispersal Distance")
