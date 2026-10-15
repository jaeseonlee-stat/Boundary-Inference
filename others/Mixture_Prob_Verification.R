# Sanity check for mixture probabilities for p=3.
#
# We compare the simulated frequencies that the
# boundary-constrained estimator lies on each C
# with the corresponding calculated probabilities.
#
# Settings:
# B = {1,2,3}, U = empty
# delta_tilde ~ N(0, Sigma/n)
# delta_hat = argmin Q(delta) subject to delta >= 0
#
# The quadratic objective is
# Q(delta) = a^\top delta + 0.5 * delta^\top W delta
# where a = - W %*% delta_tilde.

set.seed(123) 

n = 1
nsim = 10000
p = 3

W = matrix(c(1, .5, .3,
             .5, 1, .4,
             .3, .4, 1), nrow = 3, byrow = T)

Sigma = matrix(c(
  1, .6, .2,
  .6, 1, .5,
  .2, .5, 1), nrow = 3, byrow = TRUE)

eigen(W)$values
eigen(Sigma)$values

# 1. Generate delta_tilde ~ N_3(0, Sigma/n)
Z = matrix(rnorm(nsim * p), nrow = nsim, ncol = p)
delta_tilde = Z %*% chol(Sigma / n)

# 2. Minimise Q(delta) directly
# Here, we set B={1,2,3}
optim_delta_hat = function(delta_tilde, W){
  
  p = length(delta_tilde)
  
  # delta_tilde = - W^{-1} a
  # therefore a = - W delta_tilde
  a = -as.vector(W %*% delta_tilde)
  
  # objective function
  Q = function(delta){
    sum(a * delta) +
      0.5 * drop(t(delta) %*% W %*% delta)
  }
  
  # gradient of Q
  grad_Q = function(delta){
    as.vector(a + W %*% delta)
  }
  
  # constrained minimization
  fit = optim(
    par    = pmax(delta_tilde, 0),
    fn     = Q,
    gr     = grad_Q,
    method = "L-BFGS-B",
    lower  = rep(0, p),
    upper  = rep(Inf, p)
  )
  
  if (fit$convergence != 0) {warning(fit$message)}
  
  return(fit$par)
}

delta_hat = matrix(
  NA,
  nrow = nsim,
  ncol = p
)

for(r in 1:nsim){
  delta_hat[r, ] =  optim_delta_hat(delta_tilde = delta_tilde[r, ], W = W)
}


get_C = function(delta_hat, tol = 1e-8){
  
  C = which(delta_hat > tol)
  
  if(length(C) == 0){
    return("empty")
  }
  
  paste0("{", paste(C, collapse = ","), "}")
}

face_hat = apply(delta_hat, 1, get_C)

C_names = c(
  "empty",
  "{1}",
  "{2}",
  "{3}",
  "{1,2}",
  "{1,3}",
  "{2,3}",
  "{1,2,3}"
)

frequencies = table(factor(face_hat, levels = C_names)) / nsim
frequencies # empirical proportion from nsim = 10000

# calculated mixture probabilities
library(mvtnorm)

mixture_prob = function(C, W, Sigma){
  
  p = nrow(W)
  B = 1:p
  A = setdiff(B, C)
  
  # T matrix:
  # Z_C = T %*% delta_tilde
  T = matrix(0, nrow = p, ncol = p)
  
  row_now = 1
  
  # --------------------------------------------------
  # First part:
  # delta_C^F
  # = delta_tilde_C
  #   + W_CC^{-1} W_CA delta_tilde_A
  # --------------------------------------------------
  
  if(length(C) > 0){
    
    nc = length(C)
    
    T[row_now:(row_now + nc - 1), C] =
      diag(nc)
    
    if(length(A) > 0){
      
      W_CC = W[C, C, drop = FALSE]
      W_CA = W[C, A, drop = FALSE]
      
      T[row_now:(row_now + nc - 1), A] =
        solve(W_CC, W_CA)
    }
    
    row_now = row_now + nc
  }
  
  
  # --------------------------------------------------
  # Second part:
  # - W_{A|C} delta_tilde_A
  # --------------------------------------------------
  
  if(length(A) > 0){
    
    na = length(A)
    
    W_AA = W[A, A, drop = FALSE]
    
    if(length(C) == 0){
      
      W_A_C = W_AA
      
    } else {
      
      W_AC = W[A, C, drop = FALSE]
      W_CC = W[C, C, drop = FALSE]
      W_CA = W[C, A, drop = FALSE]
      
      W_A_C =
        W_AA -
        W_AC %*% solve(W_CC, W_CA)
    }
    
    T[row_now:(row_now + na - 1), A] =
      -W_A_C
  }
  
  
  # covariance of Z_C
  V = T %*% (Sigma/n) %*% t(T)
  
  # P(Z_C > 0)
  prob = pmvnorm(
    lower = rep(0, p),
    upper = rep(Inf, p),
    mean  = rep(0, p),
    sigma = V
  )
  
  return(as.numeric(prob))
}

C_list = list(
  integer(0),
  c(1),
  c(2),
  c(3),
  c(1,2),
  c(1,3),
  c(2,3),
  c(1,2,3)
)

C_names = c(
  "empty",
  "{1}",
  "{2}",
  "{3}",
  "{1,2}",
  "{1,3}",
  "{2,3}",
  "{1,2,3}"
)

calculated_probs = sapply(
  C_list,
  function(C){
    mixture_prob(
      C = C,
      W = W,
      Sigma = Sigma
    )
  }
)

names(calculated_probs) = C_names

comparison = data.frame(
  C = C_names,
  Simulation = as.numeric(frequencies),
  Probability = as.numeric(calculated_probs),
  Difference = as.numeric(frequencies) - as.numeric(calculated_probs)
)
comparison[, 2:4] =
  round(comparison[, 2:4], 4)

comparison
