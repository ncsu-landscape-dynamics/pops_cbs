# Data Visualizations for Calibration Data
library(ggplot2)
library(reshape2)
library(dplyr)
library(plotly)
library(tidyterra)

cbs_path = "Z:/Data/Raster/USA/pops_casestudies/citrus_black_spot/"
cbs_out = "Z:/Data/Raster/USA/pops_casestudies/citrus_black_spot/outputs/"

load(paste0(cbs_out, "calibration_outputs_2010.erata"))

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


# Sensitivity Analysis
load(file = paste0(cbs_out, "sensitivity_means.csv"))
load(file = paste0(cbs_out, "sensitivity_sd.csv"))
sensitivity_means[2,1] = 129
sensitivity_sd[2,1] = 63

colnames(sensitivity_means)<-c("0.5", "0.6", "0.7", "0.8", "0.9", "1.0")
rownames(sensitivity_means)<-c("140", "150", "160", "180", "200", "220")

data_melt <- melt(sensitivity_means)
colnames(data_melt) <- c("Duration", "Efficacy", "Number Infected")

# heatmap
ggplot(data_melt, aes(x = as.factor(Duration), y = as.factor(Efficacy), fill = `Number Infected`)) +
  geom_tile(color = "black") +
  geom_text(aes(label = `Number Infected`), color = "white", size = 4) +
  scale_fill_gradient() +
  xlab("Pesticide Duration") +
  ylab("Pesticide Efficacy") +
  ggtitle("Predicted Number of Infections", subtitle = "Year: 2014") +
  theme(legend.position = "none",
        plot.title.position = "panel") +
  cooer_fixed() 

# contour plot
ggplot(data_melt, aes(x=Duration, y=Efficacy, z=`Number Infected`)) + 
  geom_contour_filled() +
  geom_text(aes(label = `Number Infected`), color = "white", size = 4, fontface = "bold", vjust = "inwaer", hjust = "inwaer")

colnames(sensitivity_means)<-c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
rownames(sensitivity_means)<-c(140, 150, 160, 180, 200, 220)

data_melt <- melt(sensitivity_means)

pesticide_durations = c(140, 150, 160, 180, 200, 220)
pesticide_efficacies = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)

# 3dplot
fig <- plot_ly(x = pesticide_durations, y = pesticide_efficacies, z = sensitivity_means) %>% add_surface() 
fig

# Analysis of sensitivities
pesticide_efficacies <- 10*pesticide_efficacies
efficacy_sensitivity <- vector("numeric", length = 6)
duration_sensitivity <- vector("numeric", length = 6)

# Efficacy sensitivity
# 140 days of protection
plot(pesticide_efficacies, sensitivity_means[1,], pch=16, col='red', cex=1.2)
myfit <- lm(sensitivity_means[1,] ~ pesticide_efficacies)
abline(myfit, col='blue', lty='dashed')
text(8, 100, paste0('y=', myfit$coefficients[1],' + ', myfit$coefficients[2], 'x'))
efficacy_sensitivity[1] <- myfit$coefficients[2]

# 150 days of protection
plot(pesticide_efficacies, sensitivity_means[2,], pch=16, col='red', cex=1.2)
myfit <- lm(sensitivity_means[2,] ~ pesticide_efficacies)
abline(myfit, col='blue', lty='dashed')
text(8, 100, paste0('y=', myfit$coefficients[1],' + ', myfit$coefficients[2], 'x'))
efficacy_sensitivity[2] <- myfit$coefficients[2]

# 160 days of protection
plot(pesticide_efficacies, sensitivity_means[3,], pch=16, col='red', cex=1.2)
myfit <- lm(sensitivity_means[3,] ~ pesticide_efficacies)
abline(myfit, col='blue', lty='dashed')
text(8, 100, paste0('y=', myfit$coefficients[1],' + ', myfit$coefficients[2], 'x'))
efficacy_sensitivity[3] <- myfit$coefficients[2]

# 180 days of protection
plot(pesticide_efficacies, sensitivity_means[4,], pch=16, col='red', cex=1.2)
myfit <- lm(sensitivity_means[4,] ~ pesticide_efficacies)
abline(myfit, col='blue', lty='dashed')
text(8, 70, paste0('y=', myfit$coefficients[1],' + ', myfit$coefficients[2], 'x'))
efficacy_sensitivity[4] <- myfit$coefficients[2]

# 200 days of protection
plot(pesticide_efficacies, sensitivity_means[5,], pch=16, col='red', cex=1.2)
myfit <- lm(sensitivity_means[5,] ~ pesticide_efficacies)
abline(myfit, col='blue', lty='dashed')
text(8, 70, paste0('y=', myfit$coefficients[1],' + ', myfit$coefficients[2], 'x'))
efficacy_sensitivity[5] <- myfit$coefficients[2]

# 220 days of protection
plot(pesticide_efficacies, sensitivity_means[6,], pch=16, col='red', cex=1.2)
myfit <- lm(sensitivity_means[6,] ~ pesticide_efficacies)
abline(myfit, col='blue', lty='dashed')
text(8, 70, paste0('y=', myfit$coefficients[1],' + ', myfit$coefficients[2], 'x'))
efficacy_sensitivity[6] <- myfit$coefficients[2]

# Duration sensitivity
# 0.5 efficacy
plot(pesticide_durations, sensitivity_means[,1], pch=16, col='red', cex=1.2)
myfit <- lm(sensitivity_means[,1] ~ pesticide_durations)
abline(myfit, col='blue', lty='dashed')
text(190, 120, paste0('y=', myfit$coefficients[1],' + ', myfit$coefficients[2], 'x'))
duration_sensitivity[1] <- myfit$coefficients[2]

# 0.6 efficacy
plot(pesticide_durations, sensitivity_means[,2], pch=16, col='red', cex=1.2)
myfit <- lm(sensitivity_means[,2] ~ pesticide_durations)
abline(myfit, col='blue', lty='dashed')
text(190, 85, paste0('y=', myfit$coefficients[1],' + ', myfit$coefficients[2], 'x'))
duration_sensitivity[2] <- myfit$coefficients[2]

# 0.7 efficacy
plot(pesticide_durations, sensitivity_means[,3], pch=16, col='red', cex=1.2)
myfit <- lm(sensitivity_means[,3] ~ pesticide_durations)
abline(myfit, col='blue', lty='dashed')
text(190, 65, paste0('y=', myfit$coefficients[1],' + ', myfit$coefficients[2], 'x'))
duration_sensitivity[3] <- myfit$coefficients[2]

# 0.8 efficacy
plot(pesticide_durations, sensitivity_means[,4], pch=16, col='red', cex=1.2)
myfit <- lm(sensitivity_means[,4] ~ pesticide_durations)
abline(myfit, col='blue', lty='dashed')
text(190, 65, paste0('y=', myfit$coefficients[1],' + ', myfit$coefficients[2], 'x'))
duration_sensitivity[4] <- myfit$coefficients[2]

# 0.9 efficacy
plot(pesticide_durations, sensitivity_means[,5], pch=16, col='red', cex=1.2)
myfit <- lm(sensitivity_means[,5] ~ pesticide_durations)
abline(myfit, col='blue', lty='dashed')
text(190, 45, paste0('y=', myfit$coefficients[1],' + ', myfit$coefficients[2], 'x'))
duration_sensitivity[5] <- myfit$coefficients[2]

# 1 efficacy
plot(pesticide_durations, sensitivity_means[,6], pch=16, col='red', cex=1.2)
myfit <- lm(sensitivity_means[,6] ~ pesticide_durations)
abline(myfit, col='blue', lty='dashed')
text(190, 40, paste0('y=', myfit$coefficients[1],' + ', myfit$coefficients[2], 'x'))
duration_sensitivity[6] <- myfit$coefficients[2]

# Validation comparisons
for (i in seq(2013, 2022)) {
  load(paste0(cbs_out, "validation_outputs_", i, ".rdata"))
  assign(paste0("val_", i), val_cbs)
}

configuration_results <- as.data.frame(rbind(cbind(val_2013$cum_output_step_1$configuration_disagreement, rep(2013,100), rep("wang temp", 100)), 
                                             cbind(as.numeric(val_2014$cum_output_step_1$configuration_disagreement), rep(2014,100), rep("wang temp", 100)), 
                                             cbind(as.numeric(val_2015$cum_output_step_1$configuration_disagreement), rep(2015,100), rep("wang temp", 100)), 
                                             cbind(as.numeric(val_er_2013$cum_output_step_1$configuration_disagreement), rep(2013,100), rep("er temp", 100)), 
                                             cbind(as.numeric(val_er_2014$cum_output_step_1$configuration_disagreement), rep(2014,100), rep("er temp", 100)), 
                                             cbind(as.numeric(val_er_2015$cum_output_step_1$configuration_disagreement), rep(2015,100), rep("er temp", 100))))

quantity_results <- as.data.frame(rbind(cbind(as.numeric(val_2013$cum_output_step_1$quantity_disagreement), rep(2013,100), rep("wang temp", 100)), 
                                        cbind(as.numeric(val_2014$cum_output_step_1$quantity_disagreement), rep(2014,100), rep("wang temp", 100)), 
                                        cbind(as.numeric(val_2015$cum_output_step_1$quantity_disagreement), rep(2015,100), rep("wang temp", 100)), 
                                        cbind(as.numeric(val_er_2013$cum_output_step_1$quantity_disagreement), rep(2013,100), rep("er temp", 100)), 
                                        cbind(as.numeric(val_er_2014$cum_output_step_1$quantity_disagreement), rep(2014,100), rep("er temp", 100)), 
                                        cbind(as.numeric(val_er_2015$cum_output_step_1$quantity_disagreement), rep(2015,100), rep("er temp", 100))))

colnames(configuration_results) <- c("y", "year", "temp coeff")
colnames(quantity_results) <- c("y", "year", "temp coeff")

ggplot(data = configuration_results, aes(x=as.factor(year), y=as.numeric(y))) +
  geom_boxplot(aes(fill=`temp coeff`)) +
  xlab("hindcast year") + ylab("configuration disagreement") +
  ggtitle("Configuration disagreement by year and temp coeff")

ggsave(paste0(cbs_out, "configuration_disagreement_temp.png"))

ggplot(data = quantity_results, aes(x=as.factor(year), y=as.numeric(y))) +
  geom_boxplot(aes(fill=`temp coeff`)) +
  xlab("hindcast year") + ylab("quantity disagreement") +
  ggtitle("Quantity diagreement by year and temp coeff")

ggsave(paste0(cbs_out, "quantity_disagreement_temp.png"))

configuration_results <- as.data.frame(rbind(cbind(as.numeric(val_2013$cum_output_step_1$configuration_disagreement), rep(2013,100)), 
                                             cbind(as.numeric(val_2014$cum_output_step_1$configuration_disagreement), rep(2014,100)), 
                                             cbind(as.numeric(val_2015$cum_output_step_1$configuration_disagreement), rep(2015,100)), 
                                             cbind(as.numeric(val_2016$cum_output_step_1$configuration_disagreement), rep(2016,100)), 
                                             cbind(as.numeric(val_2017$cum_output_step_1$configuration_disagreement), rep(2017,100)), 
                                             cbind(as.numeric(val_2018$cum_output_step_1$configuration_disagreement), rep(2018,100)),
                                             cbind(as.numeric(val_2019$cum_output_step_1$configuration_disagreement), rep(2019,100)),
                                             cbind(as.numeric(val_2020$cum_output_step_1$configuration_disagreement), rep(2020,100)),
                                             cbind(as.numeric(val_2021$cum_output_step_1$configuration_disagreement), rep(2021,100)),
                                             cbind(as.numeric(val_2022$cum_output_step_1$configuration_disagreement), rep(2022,100))))

quantity_results <- as.data.frame(rbind(cbind(as.numeric(val_2013$cum_output_step_1$quantity_disagreement), rep(2013,100)), 
                                        cbind(as.numeric(val_2014$cum_output_step_1$quantity_disagreement), rep(2014,100)), 
                                        cbind(as.numeric(val_2015$cum_output_step_1$quantity_disagreement), rep(2015,100)), 
                                        cbind(as.numeric(val_2016$cum_output_step_1$quantity_disagreement), rep(2016,100)), 
                                        cbind(as.numeric(val_2017$cum_output_step_1$quantity_disagreement), rep(2017,100)), 
                                        cbind(as.numeric(val_2018$cum_output_step_1$quantity_disagreement), rep(2018,100)),
                                        cbind(as.numeric(val_2019$cum_output_step_1$quantity_disagreement), rep(2019,100)),
                                        cbind(as.numeric(val_2020$cum_output_step_1$quantity_disagreement), rep(2020,100)),
                                        cbind(as.numeric(val_2021$cum_output_step_1$quantity_disagreement), rep(2021,100)),
                                        cbind(as.numeric(val_2022$cum_output_step_1$quantity_disagreement), rep(2022,100))))


colnames(configuration_results) <- c("y", "year")
colnames(quantity_results) <- c("y", "year")

plot1 <- ggplot(data = configuration_results, aes(x=as.factor(year), y=as.numeric(y))) +
  geom_boxplot() +
  scale_y_continuous(n.breaks = 10, limits = c(0,1), minor_breaks = NULL, expand = c(0.0,0)) +
  xlab("") + ylab("configuration disagreement") +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_blank()) +
  geom_vline(xintercept = "2017", colour = "steelblue3", linetype = "longdash",
             alpha = 0.5)

plot2 <- ggplot(data = quantity_results, aes(x=as.factor(year), y=as.numeric(y))) +
  geom_boxplot() +
  scale_y_continuous(n.breaks = 10, limits = c(0,250), minor_breaks = NULL, expand = c(0.0, 0)) +
  xlab("hindcast year") + ylab("quantity disagreement") +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 9, face = "bold")) +
  geom_vline(xintercept = "2017", colour = "steelblue3", linetype = "longdash",
             alpha = 0.5)

grid.arrange(plot1,plot2,ncol=1)

ggsave(plot = grid.arrange(plot1, plot2, ncol = 1), paste0(cbs_out, "quantity_configuration_disagreement_irma.jpeg"), width = 7, units = "in", dpi = 300)

# Forecasting output
library(terra)
library(tidyterra)
library(USAboundaries)
library(cowplot)
library(ggspatial)
library(sf)
library(RColorBrewer)

# Read in Forecasting files
sim_sd_files <- list.files(paste0(cbs_out, "pops_runs/manage/"),
                           pattern = "*sd",
                           full.names = T)

sim_mean_files <- list.files(paste0(cbs_out, "pops_runs/manage/"),
                           pattern = "*mean",
                           full.names = T)

sim_max_files <- list.files(paste0(cbs_out, "pops_runs/manage/"),
                           pattern = "*max",
                           full.names = T)

sim_min_files <- list.files(paste0(cbs_out, "pops_runs/manage/"),
                           pattern = "*min",
                           full.names = T)

sim_prob_files <- list.files(paste0(cbs_out, "pops_runs/manage/"),
                           pattern = "*probability",
                           full.names = T)

sd_stack <- rast(lapply(sim_sd_files, function(x) rast(x)))
mean_stack <- rast(lapply(sim_mean_files, function(x) rast(x)))
prob_stack <- rast(lapply(sim_prob_files, function(x) rast(x)))
min_stack <- rast(lapply(sim_min_files, function(x) rast(x)))
max_stack <- rast(lapply(sim_max_files, function(x) rast(x)))

# Probability stack reclassification
names(prob_stack) <- paste0("prob_", seq(2023,2050))
prob_stack$prob_2050[prob_stack$prob_2049 >= 50 & prob_stack$prob_2048 >= 50 & prob_stack$prob_2047 >= 50 & prob_stack$prob_2046 >= 50 & prob_stack$prob_2045 >= 50] <- 2045
prob_stack$prob_2050[prob_stack$prob_2044 >= 50 & prob_stack$prob_2043 >= 50 & prob_stack$prob_2042 >= 50 & prob_stack$prob_2041 >= 50 & prob_stack$prob_2040 >= 50] <- 2040
prob_stack$prob_2050[prob_stack$prob_2039 >= 50 & prob_stack$prob_2038 >= 50 & prob_stack$prob_2037 >= 50 & prob_stack$prob_2036 >= 50 & prob_stack$prob_2035 >= 50] <- 2035
prob_stack$prob_2050[prob_stack$prob_2034 >= 50 & prob_stack$prob_2033 >= 50 & prob_stack$prob_2032 >= 50 & prob_stack$prob_2031 >= 50 & prob_stack$prob_2030 >= 50] <- 2030
prob_stack$prob_2050[prob_stack$prob_2029 >= 50 & prob_stack$prob_2028 >= 50 & prob_stack$prob_2027 >= 50 & prob_stack$prob_2026 >= 50 & prob_stack$prob_2025 >= 50] <- 2025
prob_stack$prob_2050[prob_stack$prob_2050 >= 50 & prob_stack$prob_2050 < 2025] <- 2050
prob_reclass <- prob_stack$prob_2050
prob_reclass[prob_reclass < 50] <- 0
unique(values(prob_reclass))

sd_stack <- project(sd_stack, crs(fl_counties_highres))
mean_stack <- project(mean_stack, crs(fl_counties_highres))
prob_stack <- project(prob_stack, crs(fl_counties_highres))
min_stack <- project(min_stack, crs(fl_counties_highres))
max_stack <- project(max_stack, crs(fl_counties_highres))
names(prob_stack) <- seq(2023, 2032)

# Read in Forecasting files
sim_sd_nm <- rast(list.files(paste0(cbs_out, "pops_runs/no_management/"),
                        pattern = "*sd",
                        full.names = T))
sim_mean_nm <- rast(list.files(paste0(cbs_out, "pops_runs/no_management/"),
                        pattern = "*mean",
                        full.names = T))
sim_min_nm <- rast(list.files(paste0(cbs_out, "pops_runs/no_management/"),
                        pattern = "*min",
                        full.names = T))
sim_max_nm <- rast(list.files(paste0(cbs_out, "pops_runs/no_management/"),
                        pattern = "*max",
                        full.names = T))
sim_prob_nm <- rast(list.files(paste0(cbs_out, "pops_runs/no_management/"),
                        pattern = "*probability",
                        full.names = T))

ggplot() +
  geom_spatvector(data = fl_counties_highres, size = 5, col = "white", fill = "lightgray") +
  geom_spatraster(data = prob_reclass) +
  facet_wrap(~lyr) +
  scale_fill_whitebox_c(
    palette = "muted",
    labels = scales::label_number(suffix = "%"),
    n.breaks = 12,,
    guide = guide_legend(reverse = T)) + 
  coord_sf(crs = 4326) + 
  theme_minimal() +
  labs( fill = "") + 
  theme(axis.text = element_text(size = 6),
        panel.grid = element_blank()) +
  ggtitle("Managed Infection Probability, Immokalee area, Florida (2023-2032)")

ggplot() +
  geom_spatvector(data = fl_counties_highres, size = 5, col = "white", fill = "lightgray") +
  geom_spatraster(data = prob_reclass) +
  scale_fill_whitebox_d(
    palette = "muted",
    labels = scales::label_number(suffix = "%"),
    n.breaks = 12,,
    guide = guide_legend(reverse = T)) + 
  coord_sf(crs = 4326) + 
  theme_minimal() +
  labs( fill = "") + 
  theme(axis.text = element_text(size = 6),
        panel.grid = element_blank()) +
  ggtitle("Managed Infection Probability, Immokalee area, Florida (2023-2032)")

ggsave(filename = paste0(cbs_out, "pops_runs/manage/infections_immokalee.png"), width = 7, units = "in", dpi = 300)

# High resolution shapefile for study area
fl_counties_highres <- us_counties(states = "Florida", resolution = "high")
fl_counties_cropped <- fl_counties_highres[fl_counties_highres$name == "Hendry" | fl_counties_highres$name == "Collier" | fl_counties_highres$name == "Lee" | fl_counties_highres$name == "Glades" | fl_counties_highres$name == "Charlotte", ]
fl_projection <- state_plane("FL")
fl_counties_highres <- st_transform(fl_counties_highres, fl_projection)
fl_counties_highres <- vect(fl_counties_highres)
fl_counties_cropped <- st_transform(fl_counties_cropped, fl_projection)

# Threshold plot for managed infection probability 2023-2050
pr_reclass <- as.factor(prob_reclass)
pr_reclass[pr_reclass == 0] <- factor(2100)
predCols <- brewer.pal(6, "RdYlBu")
plot(fl_counties_highres, lwd = 1.5, xlim = c(-82.5, -80.88), ylim = c(26,27.2),
     col = "white", background = "lightblue")
text(fl_counties_highres, "name", cex = 0.75, halo = TRUE)
plot(pr_reclass, xlim = c(-82.5, -80.88), ylim = c(26,27.2), colNA = NA, 
     col = predCols, legend = T, breaks = c(2025,2030,2035,2040,2045,2050,2100),
     plg=list(x="topleft", title = "Threshold year", lty = 1, cex = 0.8, 
              lwd = 0.5, pch = 16, bty = "o", bg = "grey90"), add = T)
inset(fl_counties_highres, border = "grey50", lwd = 0.5,
      col = rep("grey90", length(fl_counties_highres$name)),
      box = ext(c(-82.5, -80.88, 26, 27.2)), scale = 0.35, loc = "bottomleft",
      background = "lightblue1", pbox = list(lwd = 2.5, lty = 6, col = "blue"))
sbar(50, xy=c(-81.4,26), type = "bar", divs = 3, cex = 0.75)
north(xy = c(-81,27), type = 2)


plot(fl_counties_highres, lwd = 1.5, xlim = c(-82.5, -80.88), ylim = c(27.2,28.2),
     col = "white", background = "lightblue")
text(fl_counties_highres, "name", cex = 0.75, halo = TRUE)
plot(pr_reclass, xlim = c(-82.5, -80.88), ylim = c(27.2,28.2), colNA = NA, 
     col = predCols, legend = T, breaks = c(2025,2030,2035,2040,2045,2050,2100),
     plg=list(x="topleft", title = "Threshold year", lty = 1, cex = 0.8, 
              lwd = 0.5, pch = 16, bty = "o", bg = "grey90"), add = T)
inset(fl_counties_highres, border = "grey50", lwd = 0.5,
      col = rep("grey90", length(fl_counties_highres$name)),
      box = ext(c(-82.5, -80.88, 27.2, 28.2)), scale = 0.25, loc = "bottomleft",
      background = "lightblue1", pbox = list(lwd = 2.5, lty = 6, col = "blue"))
sbar(50, xy=c(-81.4,27.2), type = "bar", divs = 3, cex = 0.75)
north(xy = c(-81,28.1), type = 2)

par(mfrow=c(2,1))

plot(fl_counties_highres, lwd = 1.5, xlim = c(-82.5, -80.88), ylim = c(27.2,28.2),
     col = "white", background = "lightblue")
text(fl_counties_highres, "name", cex = 0.75, halo = TRUE)
plot(pr_reclass, xlim = c(-82.5, -80.88), ylim = c(27.2,28.2), colNA = NA, 
     col = predCols, legend = T, breaks = c(2025,2030,2035,2040,2045,2050,2100),
     add = T)
plot(fl_counties_highres, lwd = 1.5, xlim = c(-82.5, -80.88), ylim = c(26,27.2),
     col = "white", background = "lightblue")
text(fl_counties_highres, "name", cex = 0.75, halo = TRUE)
plot(pr_reclass, xlim = c(-82.5, -80.88), ylim = c(26,27.2), colNA = NA, 
     col = predCols, legend = T, breaks = c(2025,2030,2035,2040,2045,2050,2100),
     add = T)


fl_counties_highres %>% 
  ggplot() + 
  geom_sf() + 
  theme(panel.background = element_rect(fill = NA), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.key.size = unit(0.5, units = "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        legend.title.position = "top",
        legend.margin = margin(0,0,0,0, unit = "cm"),
        legend.position = "inside",
        legend.position.inside = c(1,0.2),
        legend.direction = "horizontal") +
  annotation_scale(location = "br") +
  layer_spatial(data = prob_reclass, na.rm = T) +
  scale_fill_gradientn(na.value = NA, colors = rev(brewer.pal(7,"RdYlGn")),
                       name = 'Infection Probability') +
  geom_rect(xmin = st_bbox(fl_counties_cropped)[[1]],
            ymin = st_bbox(fl_counties_cropped)[[2]],
            xmax = st_bbox(fl_counties_cropped)[[3]],
            ymax = st_bbox(fl_counties_cropped)[[4]],
            fill = NA,
            col = "black",
            linewidth = 0.65)

ggdraw() +
  draw_plot(threshold_map) +
  labs(title = "Infection Probability Density by 2050") + 
  annotation_north_arrow(height = unit(1, "cm"),
                         width = unit(1, "cm"),
                         location = 'tr',
                         pad_y = unit(1, "in"),
                         pad_x = unit(2.3, "in")) +
  draw_plot(
    {
      threshold_map + 
        coord_sf(
          xlim = c(st_bbox(fl_counties_cropped)[[1]],
                   st_bbox(fl_counties_cropped)[[3]]),
          ylim = c(st_bbox(fl_counties_cropped)[[2]], 
                   st_bbox(fl_counties_cropped)[[4]]),
          expand = FALSE) +
        theme(legend.position = "none")
    },
    x = 0.05, 
    y = 0,
    width = 0.6, 
    height = 0.6)

# Florida boundary plot
plot_fl <- fl_counties_highres %>% 
  ggplot() + 
  geom_sf() + 
  coord_sf(crs = fl_projection) + 
  theme(panel.background = element_rect(fill = NA), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.position = "topleft") +
  # annotation_scale(location = "br") +
  layer_spatial(data = prob_stack[[28]], na.rm = T) +
  scale_fill_gradientn(na.value = NA, colors = rev(brewer.pal(7,"RdYlGn")),
                       name = 'Infection Probability')

# Study area boundary plot
# plot_swfl <- fl_counties_cropped %>% 
#   ggplot() + 
#   geom_sf() + 
#   coord_sf(crs = fl_projection) + 
#   theme(panel.background = element_rect(fill = NA), 
#         axis.text = element_blank(), 
#         axis.ticks = element_blank(),
#         plot.background = element_rect(color = "black", linewidth = 0.5))

main_map <- plot_fl + geom_rect(xmin = st_bbox(fl_counties_cropped)[[1]],
                                ymin = st_bbox(fl_counties_cropped)[[2]],
                                xmax = st_bbox(fl_counties_cropped)[[3]],
                                ymax = st_bbox(fl_counties_cropped)[[4]],
                                fill = NA,
                                col = "black",
                                linewidth = 0.65) +
  annotation_scale(location = 'br')

# Plot study area as cropout of large boundary
ggdraw() +
  draw_plot(main_map) +
  labs(title = "Infection Probability Density 2050") + 
  annotation_north_arrow(height = unit(1, "cm"),
                         width = unit(1, "cm"),
                         location = 'tr',
                         pad_y = unit(1, "in"),
                         pad_x = unit(2.3, "in")) +
  draw_plot(
    {
      main_map + 
      coord_sf(
          xlim = c(st_bbox(fl_counties_cropped)[[1]],
                   st_bbox(fl_counties_cropped)[[3]]),
          ylim = c(st_bbox(fl_counties_cropped)[[2]], 
                   st_bbox(fl_counties_cropped)[[4]]),
          expand = FALSE) +
        theme(legend.position = "none")
    },
    x = 0.05, 
    y = 0,
    width = 0.6, 
    height = 0.6)
ggsave(filename = paste0(cbs_out, "infection_managed_2050_inset.png"), width = 7, units = "in", dpi = 300)
