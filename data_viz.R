# Data Visualizations for Calibration Data
library(ggplot2)
library(reshape2)
library(dplyr)
library(plotly)

cbs_path = "/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/citrus_black_spot/"
cbs_out = "/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/citrus_black_spot/outputs/"

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
  coord_fixed() 

# contour plot
ggplot(data_melt, aes(x=Duration, y=Efficacy, z=`Number Infected`)) + 
  geom_contour_filled() +
  geom_text(aes(label = `Number Infected`), color = "white", size = 4, fontface = "bold", vjust = "inward", hjust = "inward")

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
