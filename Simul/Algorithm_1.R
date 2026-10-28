# Global quasi-likelihood ratio test (Algorithm 1).
# Parameter order: theta = (psi_x, psi_y, psi_z, phi_xy, phi_yz, phi_xz).
# Observed LR: 2 * n * (loss under null - loss under alternative).
# Sample W and V are evaluated at the null estimate.
# No additional R packages are required.

# 1. Generate one observation from a separate Gibbs sampler
draw_one_sample = function(burn_in, theta) {
  X = Y = Z = 1
  
  for (b in seq_len(burn_in + 1)) {
    X = rpois(1, exp(theta[1] - theta[4] * Y - theta[6] * Z))
    Y = rpois(1, exp(theta[2] - theta[4] * X - theta[5] * Z))
    Z = rpois(1, exp(theta[3] - theta[6] * X - theta[5] * Y))
  }
  
  c(X = X, Y = Y, Z = Z)
}

# 2. Compute the mean of negative log quasi-likelihood.
loss = function(theta, samples) {
  X = samples[, "X"]
  Y = samples[, "Y"]
  Z = samples[, "Z"]
  
  eta_x = theta[1] - theta[4] * Y - theta[6] * Z
  eta_y = theta[2] - theta[4] * X - theta[5] * Z
  eta_z = theta[3] - theta[6] * X - theta[5] * Y
  
  mean(exp(eta_x) - X * eta_x +
         exp(eta_y) - Y * eta_y +
         exp(eta_z) - Z * eta_z)
}

# 3. Compute the loss Hessian W, gradient outer product V, and gradient.
WV = function(theta, samples) {
  n = nrow(samples)
  X = samples[, "X"]
  Y = samples[, "Y"]
  Z = samples[, "Z"]
  
  D = list(
    cbind(1, 0, 0, -Y, 0, -Z),
    cbind(0, 1, 0, -X, -Z, 0),
    cbind(0, 0, 1, 0, -Y, -X)
  )
  
  response = cbind(X, Y, Z)
  g = matrix(0, n, 6)
  W = matrix(0, 6, 6)
  
  for (j in 1:3) {
    A = D[[j]]
    lambda = exp(drop(A %*% theta))
    g = g + A * (lambda - response[, j])
    W = W + crossprod(A, A * lambda) / n
  }
  
  list(W = W, V = crossprod(g) / n, gradient = colMeans(g))
}

gradient = function(theta, samples) WV(theta, samples)$gradient

# 4. Estimate theta subject to phi >= 0; return only the estimate.
fit_for_test = function(start, samples) {
  
  fit = optim(
    par = start,
    fn = loss,
    gr = gradient,
    samples = samples,
    method = "L-BFGS-B",
    lower = c(rep(-Inf, 3), rep(0, 3)),
    control = list(maxit = 1000, factr = 1e3, pgtol = 1e-9)
  )
  
  # A convergence warning if required
  if (fit$convergence != 0) warning(fit$message)
  fit$par
}

# 5. Simulate M marginal null LR values using Algorithm 1.
# The quadratic expression is used here to obtain null critical values.
simulate_marginal_null = function(Sigma, W, n, M = 1000, B=c(4:6), U=c(1:3)) {

  W = (W + t(W)) / 2
  Sigma = (Sigma + t(Sigma)) / 2
  invisible(chol(W))
  
  # Generate the power set of B, including the empty set.
  C_set = list(integer(0))
  for (j in B) {
    C_set = c(C_set, lapply(C_set, function(C) c(C, j)))
  }
  
  # Compute the matrices for each C once, before drawing samples.
  W_UU_inv = solve(W[U, U])
  faces = lapply(C_set, function(C) {
    A = setdiff(B, C)
    F = c(C, U)
    
    K = matrix(0, length(F), length(A))
    if (length(A)) {
      K = solve(W[F, F, drop = FALSE], W[F, A, drop = FALSE]) # K=W_{FF}^{-1} W_{FA}
    }
    
    W_A_given_F = W[A, A, drop = FALSE] - W[A, F, drop = FALSE] %*% K
    W_C_given_U = W[C, C, drop = FALSE] -
      W[C, U, drop = FALSE] %*% W_UU_inv %*% W[U, C, drop = FALSE]
    
    list(C = C, A = A, F = F, K = K,
         W_A_given_F = W_A_given_F, W_C_given_U = W_C_given_U)
  })
  
  # Each row is one draw from N(0, Sigma / n).
  eig = eigen(Sigma, symmetric = TRUE)
  if (min(eig$values) < -1e-10 * max(abs(eig$values))) {
    stop("Sigma has a materially negative eigenvalue.")
  }
  Sigma_root = eig$vectors %*% diag(sqrt(pmax(eig$values, 0) / n), 6)
  delta_tilde = matrix(rnorm(M * 6), M, 6) %*% t(Sigma_root) # delta_tilde ~ N(0, Sigma/n)
  LR = rep(NA_real_, M)
  
  for (m in seq_len(M)) {
    d = delta_tilde[m, ]
    
    for (face in faces) {
      delta_F = d[face$F] + drop(face$K %*% d[face$A])
      delta_C = delta_F[seq_along(face$C)]
      dual = drop(face$W_A_given_F %*% d[face$A])
      
      # Empty C or A makes its corresponding condition automatically true.
      if (all(delta_C > 0) && all(dual <= 0)) {
        LR[m] = if (length(face$C) == 0) 0 else
          n * drop(crossprod(delta_C, face$W_C_given_U %*% delta_C))
        break
      }
    }
  }
  
  if (anyNA(LR)) stop("A Gaussian draw could not be assigned to a face.")
  LR
}

# 6. Compute the observed LR and simulated critical values for one dataset.
# alpha can be a scalar or a vector of significance levels.
test_one_dataset = function(samples, M = 1000, alpha = 0.05, 
                            B=4:6, U=1:3) {
  n = nrow(samples)
  means = colMeans(samples[, c("X", "Y", "Z"), drop = FALSE])
  
  if (any(means <= 0)) stop("An all-zero column has no finite null intercept.")
  
  # Fit the null (phi = 0) and the alternative parameter space (phi >= 0).
  theta_null = c(unname(log(means)), 0, 0, 0)
  theta_hat = fit_for_test(start = theta_null, samples = samples)
  
  # Evaluate sample W and V at the null estimate in every scenario.
  wv = WV(theta_null, samples)
  W_hat = wv$W
  W_inv = solve(W_hat)
  Sigma_hat = W_inv %*% wv$V %*% W_inv
  
  LR_obs = 2 * n * (loss(theta_null, samples) - loss(theta_hat, samples))
  if (!is.finite(LR_obs) || LR_obs < -1e-7) {
    stop("Invalid observed quasi-likelihood LR.")
  }
  LR_obs = max(LR_obs, 0) 
  
  LR_null = simulate_marginal_null(Sigma = Sigma_hat, W = W_hat, n = n, M = M, B = B, U = U)
  critical_LR = quantile(LR_null, probs = 1 - alpha, type = 1, names = FALSE)
  
  list(LR_obs = LR_obs, critical_LR = critical_LR)
}

# 7. Repeat the test R times; return the rejection rates and LR values.
# Numerical errors stop the run instead of silently excluding repetitions.
run_rejection_study = function(theta_true, n = 100, burn_in = 50,
                               R = 1000, M = 1000, alpha = 0.05, B=4:6, U=1:3,
                               seed = 1234, label = "scenario") {
  set.seed(seed)
  data_seeds = sample.int(.Machine$integer.max, R)
  null_seeds = sample.int(.Machine$integer.max, R)
  
  LR_obs = numeric(R)
  critical_LR = matrix(NA_real_, R, length(alpha))
  colnames(critical_LR) = paste0("alpha_", alpha)
  
  for (r in seq_len(R)) {
    set.seed(data_seeds[r])
    samples = t(replicate(n, draw_one_sample(burn_in = burn_in, theta = theta_true)))
    
    set.seed(null_seeds[r])
    test = test_one_dataset(samples, M = M, alpha = alpha, B = B, U = U)
    LR_obs[r] = test$LR_obs
    critical_LR[r, ] = test$critical_LR
    
    if (r %% 100 == 0 || r == R) {
      message(sprintf("%s, n = %d: %d / %d", label, n, r, R))
    }
  }
  
  # Compare the R observed values with each column of critical values.
  rate = vapply(seq_along(alpha), function(j) {
    mean(LR_obs > critical_LR[, j])
  }, numeric(1))
  
  summary = data.frame(scenario = label, n = n, alpha = alpha, rate = rate)
  list(summary = summary, LR_obs = LR_obs, critical_LR = critical_LR)
}

# 8. Run the simulation grid.
# Reuse each dataset and its null draws across all three significance levels.
theta_settings = list(
  null = c(1, 1, 1, 0, 0, 0),
  alternative = c(1, 1, 1, 0.1, 0.1, 0.1)
)
scenario_seeds = c(null = 1234, alternative = 1234)
n_values = c(10, 50, 100, 200)
alpha_values = c(0.01, 0.05, 0.10)

results = list()
for (n_current in n_values) {
  for (scenario in names(theta_settings)) {
    key = paste(scenario, n_current, sep = "_")
    results[[key]] = run_rejection_study(
      theta_true = theta_settings[[scenario]],
      n = n_current, burn_in = 50, R = 1000, M = 1000, B= 4:6, U=1:3,
      alpha = alpha_values, seed = scenario_seeds[[scenario]], label = scenario
    )
  }
}

# 9. Show Type I error and power side by side, grouped by candidate alpha.
global_summary = do.call(rbind, lapply(results, function(x) x$summary))
rownames(global_summary) = NULL

null_table = global_summary[global_summary$scenario == "null", c("n", "alpha", "rate")]
power_table = global_summary[global_summary$scenario == "alternative", c("n", "alpha", "rate")]
names(null_table)[3] = "type1_error"
names(power_table)[3] = "power"

summary_table = merge(null_table, power_table, by = c("n", "alpha"))
summary_table = summary_table[order(summary_table$alpha, summary_table$n), ]
rownames(summary_table) = NULL

for (a in alpha_values) {
  tab = summary_table[summary_table$alpha == a, c("n", "type1_error", "power")]
  tab$type1_error = sprintf("%.3f", tab$type1_error)
  tab$power = sprintf("%.3f", tab$power)
  names(tab) = c("n", "Type I error", "Power")
  
  cat(sprintf("\nalpha = %.2f\n", a))
  print(tab, row.names = FALSE, right = TRUE)
}
