# Task 2
library(ggplot2)
library(patchwork)

# Prepare the data
url <- "https://users.ssc.wisc.edu/~behansen/econometrics/cps09mar.txt"
data <- read.table(url, header = FALSE)

# Extract the variables
x1 <- data[, 5]  # hours
x2 <- data[, 6]  # education
x3 <- data[, 4]  # earnings


data <- cbind(1, x1, x2, x3)
data <- as.matrix(data)

# Simulation function
sim <- function(data, s, sigma2, beta) {
  n <- nrow(data)
  results <- matrix(NA, nrow = s, ncol = 3) 
  
  for(i in 1:s) { 
    # Generate y 
    e <- rnorm(n, 0, sqrt(sigma2))
    y <- data %*% beta + e
    
    # Data for the misspecified model without x2
    X <- data[, c(1, 2, 4)]
    X_transpose <- t(X) # X'
    XtX_inv <- solve(X_transpose %*% X) # (X'X)^-1
    
    simulated_beta <- XtX_inv %*% X_transpose %*% y 
    
    # Store the results
    results[i, 1] <- simulated_beta[1]  
    results[i, 2] <- simulated_beta[2]  
    results[i, 3] <- simulated_beta[3]
    
  }
  return(results)
}

# Parameters and replication size that are given in the task
beta <- c(20, 0.03, 0.7, 1.4)
var <- 6
s <- 300

# Simulate
sim_results <- sim(data = data, s = s, sigma2 = var, beta = beta)
sim_results <- as.data.frame(sim_results)
colnames(sim_results) <- c("Beta 0", "Beta 1", "Beta 3")

# Histogram for beta0 and beta3
create_histogram <- function(data, beta_name, beta_label, fill_color = "steelblue") {
  ggplot(data, aes(x = .data[[beta_name]])) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 30,
      fill = fill_color,
      color = "white",
      alpha = 0.8
    ) +
    labs(
      x = beta_label,
      y = "Frequency",
      title = " "
    ) +
    theme_minimal()
}

hist_beta0 <- create_histogram(
  data = sim_results,
  beta_name = "Beta 0",
  beta_label = expression(beta[0]^s)
)

hist_beta3 <- create_histogram(
  data = sim_results,
  beta_name = "Beta 3", 
  beta_label = expression(beta[3]^s)
)

# Combine and plot the histograms
hist_beta0 + hist_beta3



#Deep dive into Beta1, extract the mean, sample standard errors for the CI
means <- apply(sim_results[2], 2, mean)
standard_errors <- apply(sim_results[2], 2, sd) / sqrt(s)

ci_lower <- means - qt(0.975, s - 1) * standard_errors
ci_upper <- means + qt(0.975, s - 1) * standard_errors

ci_results <- data.frame(
  beta = "beta_1",
  mean = means,
  ci_lower = ci_lower,
  ci_upper = ci_upper
)

print(ci_results)

ggplot(sim_results, aes(x = `Beta 1`)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "steelblue",
    color = "white"
  ) +
  geom_vline(
    xintercept = c(ci_lower, ci_upper),
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(
    x = expression(beta[1]^s),
    y = "Frequency",
    title = " "
  ) +
  theme_minimal()

