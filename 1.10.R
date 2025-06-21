population <- bivariate_poisson_data$y1 + bivariate_poisson_data$y2
N <- length(population)

# Sample sizes to evaluate
sample_sizes <- c(10, 30, 50, 100, 200)
num_simulations <- 1000

# Create storage for results
results <- data.frame(
  Sample_Size = sample_sizes,
  RE_Theoretical = NA,
  RE_Empirical = NA
)

# Loop over sample sizes
for (i in seq_along(sample_sizes)) {
  n <- sample_sizes[i]
  var_wor <- numeric(num_simulations)
  var_wr <- numeric(num_simulations)
  
  for (j in 1:num_simulations) {
    # SRSWOR (without replacement)
    sample_wor <- sample(population, n, replace = FALSE)
    fpc <- (N - n) / N
    var_wor[j] <- var(sample_wor) * fpc / n
    
    # SRSWR (with replacement)
    sample_wr <- sample(population, n, replace = TRUE)
    var_wr[j] <- var(sample_wr) / n
  }
  
  # Compute RE
  re_theoretical <- (N - n) / (N - 1)
  re_empirical <- mean(var_wor) / mean(var_wr)
  
  results$RE_Theoretical[i] <- re_theoretical
  results$RE_Empirical[i] <- re_empirical
}

# Print the final result table
print(results)

library(ggplot2)

ggplot(results, aes(x = RE_Theoretical, y = RE_Empirical)) +
  geom_point(size = 3, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Relative Efficiency: Theoretical vs. Empirical",
    x = "Theoretical RE = (N - n)/(N - 1)",
    y = "Empirical RE"
  ) +
  theme_minimal()

