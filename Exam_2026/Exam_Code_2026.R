#PREPARATIONS

set.seed(2026)       

# DATA GENERATING PROCESS -------------------------------------------------------------

n <- 100 #sample size
k <- 2 # number of observed predictors

# Step 1 - Generate predictors
# We include a one for the intercept
x1 <- runif(n, 0, 3)
X <- cbind (1, x1) # design matrix for observed predictors (including intercept)

# True parameters
beta <- matrix(c(3, 0.9), nrow = k)

# We introduce heteroskedasticity through an exponential form

# We make the variance dependent on x1
lambda <- 2.5
sigma2 <- exp(X[, 2] * lambda)

# Generate error terms 
e <- rnorm(n, mean = 0, sd = sqrt(sigma2))

# Calculate Y 
Y <- X %*% beta + e

# Plot the heteroskedasticity
plot(X[, 2], e,
     xlab = "x1", 
     ylab = "Residuals (e)")
# Add a horizontal line at 0
abline(h = 0, lty = 2, col = "red")
# --------------------------------------------------------------------------------------

# Step 2 - OLS estimation
XtX <- t(X) %*% X
XtY <- t(X) %*% Y
beta_hat <- solve(XtX) %*% XtY

# Calculate residuals to estimate the variance of the error term
e_hat <- Y - X %*% beta_hat

Xtxinv <- solve(t(X) %*% X) # bread of the sandwich estimator

# Estimate the variance of the error term with HC1 adjustment 
# using corresponding formula
n <- nrow(X)
k <- ncol(X)
adjustment <- n / (n - k) 

Xe <- X * as.vector(e_hat)
meat <- t(Xe) %*% Xe # meat of the sandwich estimator

# Robust covariance matrix

var_ols_hc1 <- adjustment * (Xtxinv %*% meat %*% Xtxinv)

# Extract standard errors

se_ols_hc1 <- sqrt(diag(var_ols_hc1))

# Table for results
ols_results <- cbind(
  True_Beta = as.vector(beta),
  Estimate = as.vector(beta_hat),
  Std_Error = as.vector(se_ols_hc1)
)

print(ols_results)

# -------------------------------------------------------------------------------------

# Step 3: Feasible GLS Implementation

# Fetch residuals from OLS step
u_hat <- Y - X %*% beta_hat

ln_e_hat_sq <- log(u_hat^2)

# We use OLS to estimate lambda hat using OLS
lambda_hat <- solve(t(X) %*% X) %*% t(X) %*% ln_e_hat_sq 

# Construct the diagonal weight matrix

variance_i <- exp(X[, 2] * lambda_hat[2])

W <- diag(as.vector(1 / abs(variance_i))) 

# We calculate FGLS with the given formula 
XtWX <- t(X) %*% W %*% X
XtWX_inv <- solve(XtWX)
XtWy <- t(X) %*% W %*% Y

#Beta FGLS
beta_fgls <- XtWX_inv %*% XtWy

# Standard errors for FGLS
se_fgls <- sqrt(diag(XtWX_inv)) 

fgls_results <- cbind(
  True_Beta = as.vector(beta),
  Estimate = as.vector(beta_fgls),
  Std_Error = as.vector(se_fgls)
)

print(fgls_results)

# -----------------------------------------------------------------------------------------------------------

# Step 4: Monte Carlo Simulation 

# Parameters for the simulation

set.seed(2026)
reps <- 1000
n <- 1000 # This parameter is changed for the comparison 100 vs 1000
beta_true <- c(3, 0.9)
lambda <- 2.5

# Save results
results_ols <- matrix(NA, nrow = reps, ncol = 2)
results_ols_se <- matrix(NA, nrow = reps, ncol = 2)
results_fgls <- matrix(NA, nrow = reps, ncol = 2)
results_fgls_se <- matrix(NA, nrow = reps, ncol = 2)

# Simualtion step
for (i in 1:reps) {
  x1 <- runif(n, 0, 3)
  X <- cbind(1, x1)
  sigma2 <- exp(X[, 2] * lambda)
  e <- rnorm(n, mean = 0, sd = sqrt(sigma2))
  Y <- X %*% beta_true + e
  
  # OLS estimation
  XtX <- t(X) %*% X
  XtY <- t(X) %*% Y
  Xtxinv <- solve(XtX)
  beta_hat_sim <- solve(XtX) %*% XtY
  
  # Save results
  results_ols[i, ] <- beta_hat_sim
  
  # Robust HC1
  
  e_hat <- Y - X %*% beta_hat_sim
  n_obs <- nrow(X)
  k_params <- ncol(X)
  adjustment <- n_obs / (n_obs - k_params) # HC1 adjustment
  
  Xe <- X * as.vector(e_hat)
  meat <- t(Xe) %*% Xe # meat of the sandwich estimator)
  
  var_ols_hc1 <- adjustment * (Xtxinv %*% meat %*% Xtxinv)
  
  # Results OLS
  results_ols_se[i, ] <- sqrt(diag(var_ols_hc1))
  
  # FGLS estimation
  # Fetch residuals from OLS step
  
 ln_e_hat_sq <- log(e_hat^2)
 # Estimate lambda hat using OLS
 lambda_hat_sim <- solve(t(X) %*% X) %*% t(X) %*% ln_e_hat_sq
 
 # Construct the diagonal weight matri
 variance_i <- exp(X[, 2] * lambda_hat_sim[2])
 
 W <- diag(as.vector(1 / variance_i)) 
 
 # Finally calculate FGLS estimates using the formula
 XtWX <- t(X) %*% W %*% X
 XtWX_inv <- solve(XtWX)
 XtWy <- t(X) %*% W %*% Y
 
 #Beta FGLS
 beta_fgls <- XtWX_inv %*% XtWy
 
 # Save results
 results_fgls[i, ] <- beta_fgls
 results_fgls_se[i, ] <- sqrt(diag(XtWX_inv))
}

# Calculate Bias, Variance and Standard Error
bias_ols <- colMeans(results_ols) - beta_true
var_ols <- apply(results_ols, 2, var)
std_error_ols <- sqrt(var_ols)

bias_fgls <- colMeans(results_fgls) - beta_true
var_fgls <- apply(results_fgls, 2, var)
std_error_fgls <- sqrt(var_fgls)

mse_ols <- bias_ols^2 + var_ols
mse_fgls <- bias_fgls^2 + var_fgls

# Collect results in table
results_sim <- data.frame(
  Method = c("OLS_b0", "OLS_b1", "FGLS_b0", "FGLS_b1"),
  Bias = c(bias_ols, bias_fgls),
  Variance = c(var_ols, var_fgls),
  Std_Error = c(std_error_ols, std_error_fgls),
  MSE = c(mse_ols, mse_fgls)
)

# Print results
print(results_sim)

# -----------------------------------------------------------------------------------





