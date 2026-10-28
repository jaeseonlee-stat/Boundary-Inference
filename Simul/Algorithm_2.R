# Conditional quasi-likelihood ratio test (Algorithm 2).
# Parameter order: theta = (psi_x, psi_y, psi_z, phi_xy, phi_yz, phi_xz).
# Observed LR: 2 * n * (loss under null - loss under alternative).
# Sample W and V are evaluated at the null estimate.
# No additional R packages are required.

# 1. Generate one observation from a separate Gibbs sampler.
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

# 4. Estimate theta subject to constraint phi >= 0; return only the estimate.
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

# 5. Simulate M conditional null LR values using Algorithm 2.
# C is fixed at the observed set of positive interaction estimates.
simulate_conditional_LR = function(Sigma, W, n, C, M = 1000,
                                   batch_size = 4096, max_draws = 10000000, B=4:6, U=1:3) {
  if (length(C) == 0) return(numeric(M))
  
  A = setdiff(B, C)
  F = c(C, U)
  W = (W + t(W)) / 2
  Sigma = (Sigma + t(Sigma)) / 2
  invisible(chol(W))
  
  K = matrix(0, length(F), length(A))
  if (length(A)) {
    K = solve(W[F, F, drop = FALSE], W[F, A, drop = FALSE])
  }
  W_A_given_F = W[A, A, drop = FALSE] - W[A, F, drop = FALSE] %*% K
  W_C_given_U = W[C, C, drop = FALSE] -
    W[C, U, drop = FALSE] %*% solve(W[U, U], W[U, C, drop = FALSE])
  W_C_given_U = (W_C_given_U + t(W_C_given_U)) / 2
  
  eig = eigen(Sigma, symmetric = TRUE)
  if (min(eig$values) < -1e-10 * max(abs(eig$values))) {
    stop("Sigma has a materially negative eigenvalue.")
  }
  Sigma_root = eig$vectors %*% diag(sqrt(pmax(eig$values, 0) / n), 6)
  
  LR = numeric(M)
  m = proposals = 0
  
  # Draw proposals in batches and retain only those satisfying the fixed C.
  while (m < M && proposals < max_draws) {
    batch_n = min(batch_size, max_draws - proposals)
    d = matrix(rnorm(batch_n * 6), batch_n, 6, byrow = TRUE) %*% t(Sigma_root)
    proposals = proposals + batch_n
    
    delta_F = d[, F, drop = FALSE] + d[, A, drop = FALSE] %*% t(K)
    delta_C = delta_F[, seq_along(C), drop = FALSE]
    accepted = rowSums(delta_C > 0) == length(C)
    
    if (length(A)) {
      dual = d[, A, drop = FALSE] %*% t(W_A_given_F)
      accepted = accepted & rowSums(dual <= 0) == length(A)
    }
    
    keep = head(which(accepted), M - m)
    if (length(keep)) {
      z = delta_C[keep, , drop = FALSE]
      LR[m + seq_along(keep)] = n * rowSums((z %*% W_C_given_U) * z)
      m = m + length(keep)
    }
  }
  
  if (m < M) {
    stop("Fewer than M draws were accepted; increase max_draws and retry the same dataset.")
  }
  if (any(!is.finite(LR)) || any(LR < 0)) stop("Invalid simulated conditional LR.")
  LR
}

# 6. Compute the observed LR and conditional critical values for one dataset.
# alpha can be a scalar or a vector of significance levels.
test_one_dataset_A2 = function(samples, M = 1000, alpha = 0.05, B=4:6, U=1:3,
                               max_draws = 10000000) {
  n = nrow(samples)
  means = colMeans(samples[, c("X", "Y", "Z"), drop = FALSE])
  if (any(means <= 0)) stop("An all-zero column has no finite null intercept.")
  
  # Use the same observed statistic as the marginal test.
  theta_null = c(unname(log(means)), 0, 0, 0)
  theta_hat = fit_for_test(start = theta_null, samples = samples)
  LR_obs = 2 * n * (loss(theta_null, samples) - loss(theta_hat, samples))
  if (!is.finite(LR_obs) || LR_obs < -1e-7) {
    stop("The observed LR is nonfinite or materially negative.")
  }
  LR_obs = max(LR_obs, 0)
  
  # Treat phi estimates at or below 1e-8 as numerical zeros for face selection.
  C_hat = B[theta_hat[B] > 1e-8]
  
  # An empty C has a degenerate conditional distribution and never rejects.
  if (length(C_hat) == 0) {
    if (LR_obs > 1e-6) {
      warning("Empty C_hat but LR_obs is positive; check optimizer convergence.")
    }
    return(list(LR_obs = LR_obs, critical_LR = rep(0, length(alpha)),
                reject = rep(FALSE, length(alpha)), C_hat = C_hat))
  }
  
  # Evaluate sample W and V at the null estimate.
  wv = WV(theta_null, samples)
  W_inv = solve(wv$W)
  Sigma_hat = W_inv %*% wv$V %*% W_inv
  
  LR_null = simulate_conditional_LR(
    Sigma = Sigma_hat, W = wv$W, n = n, C = C_hat, M = M, B=B, U=U,
    max_draws = max_draws
  )
  critical_LR = quantile(LR_null, probs = 1 - alpha, type = 1, names = FALSE)
  
  list(LR_obs = LR_obs, critical_LR = critical_LR,
       reject = LR_obs > critical_LR, C_hat = C_hat)
}

# 7. Repeat the conditional test R times.
# Numerical errors stop the run instead of silently excluding repetitions.
run_rejection_study_A2 = function(theta_true, n = 100, burn_in = 50, B = 4:6, U = 1:3,
                                  R = 1000, M = 1000, alpha = 0.05,
                                  seed = 1234, label = "scenario",
                                  max_draws = 10000000) {
  set.seed(seed)
  data_seeds = sample.int(.Machine$integer.max, R)
  null_seeds = sample.int(.Machine$integer.max, R)
  
  LR_obs = numeric(R)
  critical_LR = matrix(NA_real_, R, length(alpha))
  colnames(critical_LR) = paste0("alpha_", alpha)
  reject = matrix(FALSE, R, length(alpha))
  C_hat = vector("list", R)
  
  for (r in seq_len(R)) {
    set.seed(data_seeds[r])
    samples = t(replicate(n, draw_one_sample(burn_in = burn_in, theta = theta_true)))
    
    set.seed(null_seeds[r])
    test = test_one_dataset_A2(samples, M = M, alpha = alpha, max_draws = max_draws, B=B, U=U)
    LR_obs[r] = test$LR_obs
    critical_LR[r, ] = test$critical_LR
    reject[r, ] = test$reject
    C_hat[r] = list(test$C_hat)
    
    if (r %% 100 == 0 || r == R) {
      message(sprintf("%s, n = %d: %d / %d", label, n, r, R))
    }
  }
  
  # All R datasets count, including those with an empty C (non-rejections).
  summary = data.frame(scenario = label, n = n, alpha = alpha, rate = colMeans(reject))
  list(summary = summary, LR_obs = LR_obs, critical_LR = critical_LR, C_hat = C_hat)
}

# 8. Run the simulation grid.
# Reuse each dataset and its conditional null draws across all alpha values.
theta_settings = list(
  null = c(1, 1, 1, 0, 0, 0),
  alternative = c(1, 1, 1, 0.1, 0.1, 0.1)
)
scenario_seeds = c(null = 1234, alternative = 1234)
n_values = c(10, 50, 100, 200)
alpha_values = c(0.01, 0.05, 0.10)

results_A2 = list()
for (n_current in n_values) {
  for (scenario in names(theta_settings)) {
    key = paste(scenario, n_current, sep = "_")
    results_A2[[key]] = run_rejection_study_A2(
      theta_true = theta_settings[[scenario]],
      n = n_current, burn_in = 50, R = 1000, M = 1000, B = 4:6, U = 1:3,
      alpha = alpha_values, seed = scenario_seeds[[scenario]], label = scenario
    )
  }
}

# 9. Show Type I error and power side by side, grouped by candidate alpha.
global_summary_A2 = do.call(rbind, lapply(results_A2, function(x) x$summary))
rownames(global_summary_A2) = NULL

null_table = global_summary_A2[global_summary_A2$scenario == "null", c("n", "alpha", "rate")]
power_table = global_summary_A2[global_summary_A2$scenario == "alternative", c("n", "alpha", "rate")]
names(null_table)[3] = "type1_error"
names(power_table)[3] = "power"

summary_table_A2 = merge(null_table, power_table, by = c("n", "alpha"))
summary_table_A2 = summary_table_A2[order(summary_table_A2$alpha, summary_table_A2$n), ]
rownames(summary_table_A2) = NULL

for (a in alpha_values) {
  tab = summary_table_A2[summary_table_A2$alpha == a, c("n", "type1_error", "power")]
  tab$type1_error = sprintf("%.3f", tab$type1_error)
  tab$power = sprintf("%.3f", tab$power)
  names(tab) = c("n", "Type I error", "Power")
  
  cat(sprintf("\nalpha = %.2f\n", a))
  print(tab, row.names = FALSE, right = TRUE)
}
