library(dplyr)

# Define population characteristics
set.seed(42)
population_size <- 2000

# 設定母體為 Urban 的 500 對夫妻和 Suburban 的 500 對夫妻
genders <- rep(c('male', 'female'), population_size / 2)
locations <- rep(c(rep('Urban', population_size / 4), rep('Suburban', population_size / 4)), 2)

# 設定 population data frame
population <- data.frame(
  gender = genders,
  location = locations,
  value = rpois(population_size, 50)  # variable of interest
)

# 設定不同的分層方式（性別、地區）
stratify_by_gender <- function(gender) {
  return(ifelse(gender == 'male', 'male', 'female'))
}

stratify_by_location <- function(location) {
  return(ifelse(location == 'Urban', 'Urban', 'Suburban'))
}

# Add stratum columns to the population data
population <- population %>%
  mutate(
    stratum_gender = stratify_by_gender(gender),
    stratum_location = stratify_by_location(location)
  )

# 2 種分層方式，分別對計算各層的樣本數
total_sample_size <- 50
gender_strata_sizes <- table(population$stratum_gender)
location_strata_sizes <- table(population$stratum_location)

prop_sample_sizes_gender <- round(prop.table(gender_strata_sizes) * total_sample_size)
prop_sample_sizes_location <- round(prop.table(location_strata_sizes) * total_sample_size)

# 根據不同的分層會拿到對樣本的 function
sample_stratified <- function(population, strat_column, sample_sizes) {
  samples <- lapply(names(sample_sizes), function(stratum) {
    stratum_population <- population %>%
      filter(([[strat_column]]) == stratum)
    stratum_sample <- stratum_population %>%
      sample_n(sample_sizes[[stratum]], replace = TRUE)
    return(stratum_sample)
  })
  return(bind_rows(samples))
}

# Proportional allocation sampling
prop_samples_gender <- sample_stratified(population, "stratum_gender", prop_sample_sizes_gender)
prop_samples_location <- sample_stratified(population, "stratum_location", prop_sample_sizes_location)

# 母體 mean
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
