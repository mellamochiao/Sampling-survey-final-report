library(ggplot2)

set.seed(200)
N <- 1000
x <- rpois(N, 10)
mux <- mean(x)
n_vec <- c(20, 50, 100, 200)
iter <- 1000
df <- data.frame()  

### (1) SRSWOR 
cat("(1) Unbiasedness of ȳ under SRSWOR:\n")
for (n in n_vec) {
  ybar <- numeric(iter)
  for (i in 1:iter) {
    s <- sample(x, n, replace = FALSE)
    ybar[i] <- mean(s)
  }
  cat(sprintf("n = %d, E[ȳ] = %.4f (True μ = %.4f)\n", n, mean(ybar), mux))
  
  df <- rbind(df, data.frame(
    sample_mean = ybar,
    n = as.factor(n),
    type = "SRSWOR (ȳ)"
  ))
}

### (2) SRSWR 
cat("\n(2) Unbiasedness of ȳₙ under SRSWR:\n")
for (n in n_vec) {
  ybarn <- numeric(iter)
  for (i in 1:iter) {
    s <- sample(x, n, replace = TRUE)
    ybarn[i] <- mean(s)
  }
  cat(sprintf("n = %d, E[ȳₙ] = %.4f (True μ = %.4f)\n", n, mean(ybarn), mux))
  
  df <- rbind(df, data.frame(
    sample_mean = ybarn,
    n = as.factor(n),
    type = "SRSWR (ȳₙ)"
  ))
}

### (3) SRSWR (weighted average) 
cat("\n(3) Unbiasedness of ȳᵥ under SRSWR (weighted):\n")
for (n in n_vec) {
  ybarv <- numeric(iter)
  for (i in 1:iter) {
    s <- sample(x, n, replace = TRUE)
    tab <- table(s)
    yv <- sum(as.numeric(names(tab)) * (tab / n))
    ybarv[i] <- yv
  }
  cat(sprintf("n = %d, E[ȳᵥ] = %.4f (True μ = %.4f)\n", n, mean(ybarv), mux))
  
  df <- rbind(df, data.frame(
    sample_mean = ybarv,
    n = as.factor(n),
    type = "SRSWR (ȳᵥ)"
  ))
}

###Boxplot
ggplot(df, aes(x = n, y = sample_mean, fill = type)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.75) +
  geom_hline(yintercept = mux, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Comparison of Sample Means under SRSWOR and SRSWR",
       subtitle = "Dashed red line = True μ from population",
       y = "Sample Mean", x = "Sample Size (n)",
       fill = "Estimator Type") +
  theme_minimal(base_size = 14)

#(4)
df_s2 <- data.frame()

cat("(4) E[S²] vs σ² under SRSWOR:\n")
for (n in n_vec) {
  s2 <- numeric(iter)
  for (i in 1:iter) {
    s <- sample(x, n, replace = FALSE)
    s2[i] <- var(s)
  }
  avg_s2 <- mean(s2)
  cat(sprintf("n = %d, E[S²] = %.4f (True σ² = %.4f)\n", n, avg_s2, sigma2))
  
  df_s2 <- rbind(df_s2, data.frame(
    n = factor(n),
    value = avg_s2,
    type = "E[S²]",
    truth = sigma2
  ))
}

#(5)Compare Var(ȳ) and empirical variance of ȳ under SRSWOR
for (n in n_vec) {
  ybar <- numeric(iter)
  for (i in 1:iter) {
    s <- sample(x, n, replace = FALSE)
    ybar[i] <- mean(s)
  }
  var_theory <- (N - n)/N * sigma2 / n
  var_emp <- var(ybar)
  cat(sprintf("\n(5) n=%d, Var(ȳ)=%.5f (empirical), %.5f (theoretical)", n, var_emp, var_theory))
}

#(6) Please compare Var(¯yn) and the empirical variance of y¯n under SRSWR.
for (n in n_vec) {
  ybarn <- numeric(iter)
  for (i in 1:iter) {
    s <- sample(x, n, replace = TRUE)
    ybarn[i] <- mean(s)
  }
  var_emp <- var(ybarn)
  cat(sprintf("\n(6) n=%d, Var(ȳₙ)=%.5f (empirical)", n, var_emp))
}
var_df <- data.frame()

### (5) SRSWOR：Var(ȳ) empirical vs theoretical
cat("(5) Var(ȳ) under SRSWOR:\n")
for (n in n_vec) {
  ybar <- numeric(iter)
  for (i in 1:iter) {
    s <- sample(x, n, replace = FALSE)
    ybar[i] <- mean(s)
  }
  var_emp <- var(ybar)
  var_theory <- (N - n)/N * sigma2 / n
  cat(sprintf("n = %d, Var(ȳ) = %.5f (empirical), %.5f (theoretical)\n", n, var_emp, var_theory))
  
  var_df <- rbind(var_df, data.frame(
    n = factor(n),
    value = var_emp,
    type = "SRSWOR (empirical)"
  ))
  var_df <- rbind(var_df, data.frame(
    n = factor(n),
    value = var_theory,
    type = "SRSWOR (theoretical)"
  ))
}

### (6) SRSWR：Var(ȳₙ) empirical vs theoretical
cat("\n(6) Var(ȳₙ) under SRSWR:\n")
for (n in n_vec) {
  ybarn <- numeric(iter)
  for (i in 1:iter) {
    s <- sample(x, n, replace = TRUE)
    ybarn[i] <- mean(s)
  }
  var_emp <- var(ybarn)
  var_theory <- sigma2 / n
  cat(sprintf("n = %d, Var(ȳₙ) = %.5f (empirical), %.5f (theoretical)\n", n, var_emp, var_theory))
  
  var_df <- rbind(var_df, data.frame(
    n = factor(n),
    value = var_emp,
    type = "SRSWR (empirical)"
  ))
  var_df <- rbind(var_df, data.frame(
    n = factor(n),
    value = var_theory,
    type = "SRSWR (theoretical)"
  ))
}


ggplot(var_df, aes(x = n, y = value, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Empirical vs Theoretical Variance of Sample Mean",
       subtitle = "Comparison under SRSWOR and SRSWR",
       x = "Sample Size (n)",
       y = "Variance",
       fill = "Method") +
  theme_minimal(base_size = 14)


# (7) (8)
df_7 <- data.frame()
df_8 <- data.frame()

# ---------- (7) ----------
cat(" (7) Var̂(ȳ) under SRSWOR\n")
for (n in n_vec) {
  ybar <- numeric(iter)
  var_hat <- numeric(iter)
  for (i in 1:iter) {
    s <- sample(x, n, replace = FALSE)
    ybar[i] <- mean(s)
    S2 <- var(s)
    var_hat[i] <- (N - n)/N * S2 / n
  }
  theo <- (N - n)/N * sigma2 / n
  avg_hat <- mean(var_hat)
  cat(sprintf("n = %d, E[ȳ] = %.4f, E[Var̂(ȳ)] = %.5f (Theory = %.5f)\n",
              n, mean(ybar), avg_hat, theo))
  df_7 <- rbind(df_7, data.frame(
    n = factor(n),
    value = avg_hat,
    theory = theo,
    type = "SRSWOR: E[Var̂(ȳ)]"
  ))
}

# ---------- (8) ----------
cat("(8) Var̂(ȳₙ) under SRSWR\n")
for (n in n_vec) {
  ybarn <- numeric(iter)
  var_hat <- numeric(iter)
  for (i in 1:iter) {
    s <- sample(x, n, replace = TRUE)
    ybarn[i] <- mean(s)
    S2 <- var(s)
    var_hat[i] <- S2 / n
  }
  theo <- sigma2 / n
  avg_hat <- mean(var_hat)
  cat(sprintf("n = %d, E[ȳₙ] = %.4f, E[Var̂(ȳₙ)] = %.5f (Theory = %.5f)\n",
              n, mean(ybarn), avg_hat, theo))
  df_8 <- rbind(df_8, data.frame(
    n = factor(n),
    value = avg_hat,
    theory = theo,
    type = "SRSWR: E[Var̂(ȳₙ)]"
  ))
}
#plot
ggplot(df_all, aes(x = n, y = value, fill = type)) +
  geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.8), color = "black") +
  geom_point(aes(y = theory), color = "red", shape = 17, size = 3,
             position = position_dodge(0.8)) +
  labs(title = "Variance Estimators for Sample Means (7 & 8)",
       subtitle = "Red triangles = theoretical variance from formulas",
       x = "Sample Size (n)", y = "Variance Estimate",
       fill = "Estimator Type") +
  theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
