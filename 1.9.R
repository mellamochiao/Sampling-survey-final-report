#install.packages("bivpois")
library(bivpois)
#install.packages("MASS")
library(MASS)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("reshape2")
library(reshape2)

set.seed(200)

N <- 1000
# 指定泊松分佈的參數
lambda1 <- 2
lambda2 <- 3
lambda3 <- 10
# 生成三個獨立的泊松隨機變量
x1 <- rpois(N, lambda1)
x2 <- rpois(N, lambda2)
x3 <- rpois(N, lambda3)
# 生成雙變量泊松分佈
y1 <- x1 + x3
y2 <- x2 + x3
# 將生成的數據存儲在一個數據框中
bivariate_poisson_data <- data.frame(y1, y2)

population <- bivariate_poisson_data$y1 + bivariate_poisson_data$y2
mu <- mean(population)  # 真實母體平均

sample_sizes <- c(10, 30, 50, 100)
num_replicates <- 1000

coverage_results <- data.frame(
  Sample_Size = sample_sizes,
  Coverage_Probability_z = NA,
  Coverage_Probability_t = NA
)

for (i in seq_along(sample_sizes)) {
  n <- sample_sizes[i]
  z_hits <- 0
  t_hits <- 0
  
  for (j in 1:num_replicates) {
    # SRSWOR 抽樣（不放回）
    sample_data <- sample(population, n, replace = FALSE)
    ybar <- mean(sample_data)
    s <- sd(sample_data)
    
    # 信賴區間（z與t）
    z_crit <- qnorm(0.975)
    t_crit <- qt(0.975, df = n - 1)
    
    z_interval <- c(ybar - z_crit * s / sqrt(n), ybar + z_crit * s / sqrt(n))
    t_interval <- c(ybar - t_crit * s / sqrt(n), ybar + t_crit * s / sqrt(n))
    
    if (mu >= z_interval[1] && mu <= z_interval[2]) z_hits <- z_hits + 1
    if (mu >= t_interval[1] && mu <= t_interval[2]) t_hits <- t_hits + 1
  }
  
  # 儲存結果
  coverage_results$Coverage_Probability_z[i] <- z_hits / num_replicates
  coverage_results$Coverage_Probability_t[i] <- t_hits / num_replicates
}

# 顯示結果
print(coverage_results)

library(ggplot2)
library(reshape2)  # 如果你要轉換 data frame 長格式

# 先轉為長格式以便 ggplot 使用
coverage_long <- melt(coverage_results, id.vars = "Sample_Size",
                      variable.name = "Method", value.name = "Coverage")

# 畫圖
ggplot(coverage_long, aes(x = Sample_Size, y = Coverage, color = Method)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("red", "blue"),
                     labels = c("Z distribution", "T distribution")) +
  labs(
    title = "95% CI Coverage Probability",
    x = "n",
    y = "Coverage Probability",
    color = "methods"
  ) +
  theme_minimal() +
  ylim(0.9, 0.97)


# 畫差值圖
# Add difference column
coverage_results$Coverage_Diff_t_minus_z <- coverage_results$Coverage_Probability_t - coverage_results$Coverage_Probability_z

# Plot the difference
ggplot(coverage_results, aes(x = Sample_Size, y = Coverage_Diff_t_minus_z)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_line(linewidth = 1, color = "darkgreen") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Difference in Coverage: T - Z",
    x = "Sample Size (n)",
    y = "Coverage Difference"
  ) +
  theme_minimal()







