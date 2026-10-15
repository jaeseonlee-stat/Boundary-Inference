# Simulate the boundary-constrained estimator for p = 2.
#
# Settings:
# B = {1}, U = {2}.
# The first coordinate satisfies delta_1 >= xi_1,
# while the second coordinate is unrestricted.
#
# delta_tilde ~ N(0, Sigma/n).

sim_delta_U = function(nsim, n, W, Sigma, xi, seed = 123){
  set.seed(seed)
  
  # check W \succ 0, \xi <= 0
  
  # draw delta_tilde ~ N(0, Sigma/n)
  Z = matrix( rnorm(2 * nsim), ncol = 2)
  delta_tilde = Z %*% chol(Sigma/n)
  
  delta_hat = matrix(0, nrow = nsim, ncol = 2)
  face = character(nsim)
  
  for (r in 1:nsim){
    
    d1 = delta_tilde[r, 1]; d2 = delta_tilde[r, 2]
    
    # determine C
    if ( d1 > xi[1]){ # C={1}
      
      delta_hat[r, ] = c(d1, d2)
      face[r] = "{1}"
      
      
    } else {
      
      delta_hat[r, ] = c(xi[1], d2 + (W[2,1]/W[2,2]) * (d1 - xi[1]))
      face[r] = "empty"
    }
  }
  
  data.frame(
    tilde1 = delta_tilde[, 1],
    tilde2 = delta_tilde[, 2],
    hat1 = delta_hat[, 1],
    hat2 = delta_hat[, 2],
    face = face
  )
}

# Calculating the skewness
skewness = function(x) {
  z = x - mean(x)
  m2 = mean(z^2)
  
  if (m2 == 0) return(NA_real_)
  
  mean(z^3) / m2^(3/2)
}


# Run simulations
n = 1
nsim = 2000000

# Please choose appropriate W and Sigma for each Example from the below
# Example 1
W = matrix(c(1, -0.99, -0.99, 1), nrow = 2)
Sigma = matrix(c(1, 0.99, 0.99, 1), nrow = 2)

# Example 2
W     = matrix(c(10, -5, -5, 3), nrow = 2)
Sigma = matrix(c(2, 7, 7, 26), nrow = 2)

# Example 3
W = matrix(c(1, -0.9, -0.9, 1), nrow = 2)
Sigma = matrix(c(5, -2, -2, 1), nrow = 2)

# Check if W and Sigma are positive definite
eigen(W); eigen(Sigma)

# xi1 values
xi_values = c(0, -0.1, -0.2, -0.5, -1, -2)

sims = lapply(xi_values, function(xi) {
  sim_delta_U(nsim, n, W, Sigma, xi = c(xi, 0), seed = 123)
})

# Simulated distributions for each C and marginal
u0 = lapply(sims, function(d) d$hat2[d$face == "empty"])
u1 = lapply(sims, function(d) d$hat2[d$face == "{1}"])
u  = lapply(sims, function(d) d$hat2)

# simulated mixture frequencies (%).
pct0 = 100 * lengths(u0) / nsim
pct1 = 100 * lengths(u1) / nsim

data.frame(
  xi = xi_values,
  C_empty = round(pct0, 3),
  C_1 = round(pct1, 3)
)

# Collect the three distributions.
u_groups = list(
  "C = empty" = u0,
  "C = {1}" = u1,
  "Marginal" = u
)

# Estimate densities using a common bandwidth and x-axis
# across different xi values within each distribution.
estimate_densities = function(u) {
  
  bw_each = vapply(u, bw.nrd0, numeric(1))
  bw = 1.1 * max(bw_each)
  
  xlim = c(
    min(vapply(u, min, numeric(1))),
    max(vapply(u, max, numeric(1)))
  ) + c(-3, 3) * bw
  
  dens = lapply(
    u, density,
    bw = bw,
    n = 2048,
    from = xlim[1],
    to = xlim[2]
  )
  
  return(dens)
}

# Calculate densities for all three distributions.
dens_groups = lapply(u_groups, estimate_densities)

# Plot the conditional and marginal distributions.
cols = c("black", "blue", "darkgreen", "red", "magenta", "grey")
ltys = c(1, 2, 6, 4, 5, 1)

par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
for (j in seq_along(dens_groups)) {
  
  dens = dens_groups[[j]]
  
  ys = vapply(dens, function(d) d$y, numeric(2048))
  
  # Conditional distributions: include mixture probabilities.
  if (j == 1) {
    pct = pct0
  } else if (j == 2) {
    pct = pct1
  }
  
  if (j <= 2) {
    
    details = sprintf("(%.3f%%)", pct)
    
    legend_labels = as.expression(
      lapply(seq_along(xi_values), function(i) {
        bquote(xi[1] == .(xi_values[i]) ~ .(details[i]))
      })
    )
    
    ylab = "Conditional density"
    
  } else {
    
    legend_labels = as.expression(
      lapply(xi_values, function(xi) {
        bquote(xi[1] == .(xi))
      })
    )
    
    ylab = "Marginal density"
  }
  
  matplot(
    dens[[1]]$x, ys,
    type = "l",
    col = cols,
    lty = ltys,
    lwd = 3,
    ylim = c(0, 1.2 * max(ys)),
    bty = "l",
    xlab = expression(hat(delta)[2]),
    ylab = ylab,
    main = names(dens_groups)[j]
  )
  
  legend(
    "topright",
    legend = legend_labels,
    col = cols,
    lty = ltys,
    lwd = 2.5,
    cex = 0.75,
    bty = "n"
  )
}

# dev.off()

# Calculate skewness for all distributions.

skewness_table = sapply(
  u_groups,
  function(u) {
    vapply(u, skewness, numeric(1))
  }
)

skewness_table = data.frame(
  xi = xi_values,
  round(skewness_table, 4),
  check.names = FALSE
)

print(skewness_table, row.names = FALSE)
