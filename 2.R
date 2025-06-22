#(1)
set.seed(200)
N <- 1000
iter <- 1000
sample_sizes <- c(10, 30, 50, 100)
lambda1 <- 2
lambda2 <- 3
lambda3_list <- c(0.27273, 1.05629, 2.47481, 5.79815, 22.45491)
cor_target <- c(0.1, 0.3, 0.5, 0.7, 0.9)
bias_matrix <- matrix(0, nrow = length(sample_sizes), ncol = length(lambda3_list))
rownames(bias_matrix) <- paste0("n=", sample_sizes)
colnames(bias_matrix) <- paste0("Cor≈", cor_target)
for (c in 1:length(lambda3_list)) {
  lambda3 <- lambda3_list[c]
  x1 <- rpois(N, lambda1)
  x2 <- rpois(N, lambda2)
  x3 <- rpois(N, lambda3)
  y1 <- x1 + x3 
  y2 <- x2 + x3   
  mu_y <- mean(y1)
  mu_x <- mean(y2)
  for (s in 1:length(sample_sizes)) {
    n <- sample_sizes[s]
    mu_r_hat_all <- numeric(iter)
    for (i in 1:iter) {
      index <- sample(1:N, n, replace = FALSE)
      y_sample <- y1[index]
      x_sample <- y2[index]
      mu_r_hat <- mean(y_sample) / mean(x_sample) * mu_x
      mu_r_hat_all[i] <- mu_r_hat
    }
    bias_matrix[s, c] <- mean(mu_r_hat_all) - mu_y
  }
}
bias_df <- as.data.frame(bias_matrix)
bias_df$SampleSize <- rownames(bias_matrix)
print(bias_df)
library(ggplot2)
library(reshape2)
bias_long <- melt(bias_df, id.vars = "SampleSize", variable.name = "Correlation", value.name = "Bias")
bias_long$SampleSize <- factor(bias_long$SampleSize, levels = paste0("n=", sample_sizes))
ggplot(bias_long, aes(x = Correlation, y = Bias, group = SampleSize, color = SampleSize)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(title = "Bias of μ̂_r vs Cor(x, y)",
       x = "Correlation (Cor(x, y))",
       y = "Bias of μ̂_r",
       color = "Sample Size n") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

#(2)
set.seed(200)
N <- 1000        
iter <- 1000     
sample_sizes <- c(10, 30, 50, 100)
lambda1 <- 2
lambda2 <- 3
lambda3_list <- c(0.27273, 1.05629, 2.47481, 5.79815, 22.45491)
cor_target <- c(0.1, 0.3, 0.5, 0.7, 0.9)
mse_matrix <- matrix(0, nrow = length(sample_sizes), ncol = length(lambda3_list))
rownames(mse_matrix) <- paste0("n=", sample_sizes)
colnames(mse_matrix) <- paste0("Cor≈", cor_target)
for (c in 1:length(lambda3_list)) {
  lambda3 <- lambda3_list[c]
  x1 <- rpois(N, lambda1)
  x2 <- rpois(N, lambda2)
  x3 <- rpois(N, lambda3)
  y1 <- x1 + x3 
  y2 <- x2 + x3  
  mu_y <- mean(y1)
  mu_x <- mean(y2)
  for (s in 1:length(sample_sizes)) {
    n <- sample_sizes[s]
    mu_r_hat_all <- numeric(iter)
    for (i in 1:iter) {
      index <- sample(1:N, n, replace = FALSE)
      y_sample <- y1[index]
      x_sample <- y2[index]
      mu_r_hat <- mean(y_sample) / mean(x_sample) * mu_x
      mu_r_hat_all[i] <- mu_r_hat
    }
    mse_matrix[s, c] <- mean((mu_r_hat_all - mu_y)^2)
  }
}
mse_df <- as.data.frame(mse_matrix)
mse_df$SampleSize <- rownames(mse_matrix)
print(mse_df)
mse_long <- melt(mse_df, id.vars = "SampleSize", variable.name = "Correlation", value.name = "MSE")
mse_long$SampleSize <- factor(mse_long$SampleSize, levels = paste0("n=", sample_sizes))
mse_long$Correlation <- as.numeric(gsub("Cor≈", "", mse_long$Correlation))
ggplot(mse_long, aes(x = Correlation, y = MSE, group = SampleSize, color = SampleSize)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(title = "Empirical MSE of μ̂_r vs Cor(x, y)",
       x = "Correlation (Cor(x, y))",
       y = "Empirical MSE of μ̂_r",
       color = "Sample Size n") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))
#(3)
set.seed(200)
N <- 1000
iter <- 1000
sample_sizes <- c(10, 30, 50, 100)
lambda1 <- 2
lambda2 <- 3
lambda3 <- 10 
x1 <- rpois(N, lambda1)
x2 <- rpois(N, lambda2)
x3 <- rpois(N, lambda3)
y1 <- x1 + x3   
y2 <- x2 + x3   
mu_y <- mean(y1)
mu_x <- mean(y2)
coverage_table <- data.frame()
for (n in sample_sizes) {
  coverage_hat <- 0  
  coverage_tilde <- 0  
  for (i in 1:iter) {
    index <- sample(1:N, n, replace = FALSE)
    x_s <- y2[index]
    y_s <- y1[index]
    x_bar <- mean(x_s)
    y_bar <- mean(y_s)
    s_x2 <- var(x_s)
    s_y2 <- var(y_s)
    s_xy <- cov(x_s, y_s)
    mu_r_hat <- y_bar / x_bar * mu_x
    var_hat <- (1/n) * (s_y2 / x_bar^2 + (y_bar^2 * s_x2) / x_bar^4 - 2 * y_bar * s_xy / x_bar^3) * mu_x^2
    se_hat <- sqrt(var_hat)
    mu_r_vec <- y_s / x_s * mu_x
    se_tilde <- sd(mu_r_vec) / sqrt(n)
    ci_hat <- c(mu_r_hat - 1.96 * se_hat, mu_r_hat + 1.96 * se_hat)
    ci_tilde <- c(mu_r_hat - 1.96 * se_tilde, mu_r_hat + 1.96 * se_tilde)
    if (mu_y >= ci_hat[1] && mu_y <= ci_hat[2]) coverage_hat <- coverage_hat + 1
    if (mu_y >= ci_tilde[1] && mu_y <= ci_tilde[2]) coverage_tilde <- coverage_tilde + 1
  }
  coverage_table <- rbind(
    coverage_table,
    data.frame(
      SampleSize = n,
      Coverage_Hat = round(coverage_hat / iter, 4),
      Coverage_Tilde = round(coverage_tilde / iter, 4)
    )
  )
}
print(coverage_table)
coverage_long <- melt(coverage_table, id.vars = "SampleSize", variable.name = "Estimator", value.name = "Coverage")
coverage_long$Estimator <- factor(coverage_long$Estimator,
                                  levels = c("Coverage_Hat", "Coverage_Tilde"),
                                  labels = c("Var_hat (Linearized)", "Var_tilde (Empirical)"))
ggplot(coverage_long, aes(x = SampleSize, y = Coverage, color = Estimator, group = Estimator)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  labs(title = "95% CI Coverage Probability vs Sample Size",
       x = "Sample Size (n)",
       y = "Coverage Probability",
       color = "Variance Estimator") +
  theme_minimal() +
  theme(text = element_text(size = 14))
