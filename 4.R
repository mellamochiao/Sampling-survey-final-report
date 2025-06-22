#R code - analysis - ( Question 4 )
set.seed(200)
library(dplyr)

# setting
M <- 60		# 60 clusters
m <- 5			# choose 5 cluster from M
K <- 20			# 20 samples from each cluster
n_total <- m * K  
B <- 1000         	# conduct 1000 times
rho <- 0.85       	# ICC

# cluster 
cluster_means <- rnorm(M, mean = 100, sd = 10)  
residual_sd <- sqrt((1 - rho) * 20)             
population <- do.call(rbind, lapply(1:M, function(cl) {
  data.frame(
    cluster = cl,
    value = rnorm(K, mean = cluster_means[cl], sd = residual_sd))
}))

# True population mean
true_mean <- mean(population$value)  

# result
results <- replicate(B, {
  sampled_cl <- sample(1:M, m)
  cluster_sample <- filter(population, cluster %in% sampled_cl)
  mean_cluster <- mean(cluster_sample$value)
  
  # SRSWOR
  srs_sample <- slice_sample(population, n = n_total)
  mean_srs <- mean(srs_sample$value)
  
  c(cluster = mean_cluster, srs = mean_srs)
})
results_df <- as.data.frame(t(results))
summary_df <- data.frame(
  Method = c("Cluster Sampling", "SRSWOR"),
  Sample_Mean = c(mean(results_df$cluster), mean(results_df$srs)),
  True_Mean = true_mean,
  Bias = c(mean(results_df$cluster) - true_mean,mean(results_df$srs) - true_mean),
  Variance = c(var(results_df$cluster),var(results_df$srs)),
  MSE = c(mean((results_df$cluster - true_mean)^2),mean((results_df$srs - true_mean)^2))
)

# RE
summary_df$Relative_Efficiency_vs_Cluster <- summary_df$MSE / summary_df$MSE[1]
summary_df$Relative_Efficiency_vs_Cluster[1] <- "-" 
print(summary_df, row.names = FALSE)

#R code - conduct plot - ( Question 4 )
library(ggplot2)
library(tidyr)
library(dplyr)
results_long <- results_df %>%
  pivot_longer(cols = c(cluster, srs), names_to = "Method", values_to = "SampleMean") %>%
  mutate(Method = recode(Method,
                         "cluster" = "Cluster Sampling",
                         "srs" = "SRSWOR"))

ggplot(results_long, aes(x = SampleMean, fill = Method, color = Method)) +
  geom_density(alpha = 0.4) +
  labs(title = "Sampling Distribution of Sample Means",
       x = "Sample Mean", y = "Density") +
  theme_minimal()
