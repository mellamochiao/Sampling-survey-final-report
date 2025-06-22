# mean
mean_y <- mean(population$value)

# Calculate mean within each stratum (mean_yi)
mean_yi_gender <- prop_samples_gender %>%
  group_by(stratum_gender) %>%
  summarise(mean_value = mean(value)) %>%
  pull(mean_value)
mean_yi_location <- prop_samples_location %>%
  group_by(stratum_location) %>%
  summarise(mean_value = mean(value)) %>%
  pull(mean_value)

# Calculate SSB
SSB_gender <- sum(prop_sample_sizes_gender * (mean_yi_gender - mean_y)^2)
SSB_location <- sum(prop_sample_sizes_location * (mean_yi_location - mean_y)^2)

# Print results
cat("SSB based on Gender Stratification:", SSB_gender, "\n")
cat("SSB based on Location Stratification:", SSB_location, "\n")
