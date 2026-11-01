library(ggplot2)
# ------------------------------------------------------------------------------------------------

# Exercise 1a and 1b

# Load the data
x <- c(4.06, 3.63, 0.41, 2.45, 0.22, 3.39, 0.27, 0.69, 0.32, 0.09, 4.00, 
       3.63, 2.65, 3.06, 3.77, 2.29, 2.31, 1.51)

y <- c(168.44, 178.66, 58.81, 189.29, 14.79, 178.70, 103.69, 103.30, 57.30, 
       15.10, 172.28, 172.30, 194.32, 193.20, 167.50, 191.50, 173.98, 174.60)

# Fit the non-linear regression model using log transformation
log_y <- log(y)
# We use cbind to create the design matrix
X_mat <- cbind(1, x, log(x))

# Solve beta_hat using theorem from lecture 3 slide 18 (theorem 3.1)
beta_hat <- solve(t(X_mat) %*% X_mat) %*% t(X_mat) %*% log_y

# Create a data frame for plotting
df <- data.frame(x = x, y = y)

# Set coefficients
b0 <- beta_hat[1]
b1 <- beta_hat[2]
b2 <- beta_hat[3]

# Fitted regression model 
model_curve <- function(x) {
  exp(b0) * exp(b1 * x) * (x ^ b2)
}

# Plot the data and the fitted curve
ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(color = "Observed data")) +
  stat_function(fun = model_curve, 
                size = 1, 
                aes(color = "Fitted curve")) +
  labs(x = "x",
       y = "y",
       color = "Data") +
  scale_color_manual(values = c("Observed data" = "black", "Fitted curve" = "red")) +
  theme_minimal()

# Exercise 1.c

x_critical <- - (b2 / b1)
x_critical

# See if second derivative is negative at critical point

-b2/(x_critical^2)

#-----------------------------------------------------------------------------------------------