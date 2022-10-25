###########################################################################
###########################################################################
###                                                                     ###
###                       LIBRARIES AND FUNCTIONS                       ###
###                                                                     ###
###########################################################################
###########################################################################



##### Authors 
# Hung Tong, University of Alabama
# Cristina Tortora, San Jose State University


library(mvtnorm)
library(cluster)
library(ContaminatedMixt)
library(MixGHD)
library(mice)
library(rootSolve)
library(gtools)
library(cubature)
library(MixtureMissing)
library(Matrix)

#------------------------------------------------------------------------
# Select some observations within a data set and generate missing values
#------------------------------------------------------------------------
hide.values <- function(data, n_miss = 10) {
  mis_indices <- sample(1:nrow(data), n_miss)
  p <- ncol(data)
  
  for (i in mis_indices) {
    new_xi <- data[i, ]
    new_xi[sample(1:p, sample(1:(p - 1), 1))] <- NA
    data[i, ] <- new_xi
  }
  
  return(data)
}

#-----------------------------------------------------
# Probability density function of a MSCN distribution
#-----------------------------------------------------
dMSCN <- function(x, mu = rep(0, d), Gamma = NULL, Lambda = NULL, Sigma = NULL, alpha = rep(0.99, d), eta = rep(1.01, d)) {
  if (is.vector(mu)) {
    d <- length(mu)
  } else {
    d <- 1
  }
  
  if (!is.null(Sigma)) {
    eigen.decomp <- eigen(Sigma)
    Gamma <- eigen.decomp$vectors
    Lambda <- diag(eigen.decomp$values)
  } else {
    if (nrow(Gamma) != ncol(Gamma)) {
      stop('Gamma must be a square matrix')
    }
    if (nrow(Lambda) != ncol(Lambda)) {
      stop('Lambda must be a square matrix')
    }
    if (any(Lambda[upper.tri(Lambda)] != 0, Lambda[lower.tri(Lambda)] != 0)) {
      stop('Lambda must be a diagonal matrix')
    }
    if (nrow(Gamma) != nrow(Lambda)) {
      stop('Gamma and Lambda must have the same dimensions')
    }
  }
  
  if (any(alpha < 0 | alpha > 1)) {
    stop('alpha must be a vector with elements in (0,1)')
  }
  if (any(eta < 1)) {
    stop('eta must be a vector with elements greater than 1')
  }
  
  if (is.matrix(x)) {
    n <- nrow(x)
  }
  if (is.vector(x) & d > 1) {
    n <- 1
  }
  if (is.vector(x) & d == 1) {
    n <- length(x)
  }
  
  PDF <- rep(1, n)
  for (i in 1:n) {
    if (is.matrix(x)) {
      y <- t(Gamma) %*% (x[i, ] - mu)
    }
    if (is.vector(x) & d > 1) {
      y <- t(Gamma) %*% (x - mu)
    }
    if (is.vector(x) & d == 1) {
      y <- x - mu
    }
    
    lambda <- diag(Lambda)
    for (h in 1:d) {
      # res[h] <- dCN(y[h], 0, lambda[h], alpha[h], eta[h])
      good.norm <- dnorm(y[h], mean = 0, sd = sqrt(lambda[h]))
      bad.norm <- dnorm(y[h], mean = 0, sd = sqrt(eta[h] * lambda[h]))
      PDF[i] <- PDF[i] * (alpha[h] * good.norm + (1 - alpha[h]) * bad.norm)
    }
  }
  
  PDF <- (PDF < 10^-323) * 10^-323 + (PDF > 10^-323) * PDF
  return(PDF)
}

#-------------------------------------------------
# Generate realizations of a MSCN random variable
#-------------------------------------------------
rMSCN <- function(n, mu = rep(0, d), Gamma = NULL, Lambda = NULL, Sigma = NULL, alpha = rep(0.99, d), eta = rep(1.01, d)) {
  if (!is.vector(mu)) {
    stop('mu must be a vector')
  }
  d <- length(mu)
  
  if (!is.null(Sigma)) {
    eigen.decomp <- eigen(Sigma)
    Gamma <- eigen.decomp$vectors
    Lambda <- diag(eigen.decomp$values)
  } else {
    if (nrow(Gamma) != ncol(Gamma)) {
      stop('Gamma must be a square matrix')
    }
    if (nrow(Lambda) != ncol(Lambda)) {
      stop('Lambda must be a square matrix')
    }
    if (any(Lambda[upper.tri(Lambda)] != 0, Lambda[lower.tri(Lambda)] != 0)) {
      stop('Lambda must be a diagonal matrix')
    }
    if (nrow(Gamma) != nrow(Lambda)) {
      stop('Gamma and Lambda must have the same dimensions')
    }
  }
  
  if (any(alpha < 0 | alpha > 1)) {
    stop('alpha must be a vector with elements in (0,1)')
  }
  if (any(eta < 1)) {
    stop('eta must be a vector with elements greater than 1')
  }
  
  X <- matrix(NA, nrow = d, ncol = n)
  V <- matrix(NA, nrow = n, ncol = d)
  for (i in 1:n) {
    Y <- replicate(1, rnorm(d))
    V[i, ] <- sapply(1:d, function(h) rbinom(1, 1, prob = alpha[h]))
    W <- (V[i, ] + (1 - V[i, ]) / eta)^-1
    W <- diag(W)
    X[, i] <- mu + Gamma %*% Lambda^(0.5) %*% W^(0.5) %*% Y
  }
  return(t(X))
}

#-----------------------------------
# Find the trace of a square matrix
#-----------------------------------
tr <- function(x) {
  if (!is.matrix(x)) {
    stop('x must be a matrix.')
  }
  if (nrow(x) != ncol(x)) {
    stop('x must be a square matrix.')
  }
  return(sum(diag(x)))
}

#-----------------------------------------------------------------------------
# Marginal density function of a MSCN random variable, e.g. f(xi^o | Z_ig = 1)
#-----------------------------------------------------------------------------

find.z <- function(x, mu = rep(0, d), Gamma = NULL, Lambda = NULL,alpha = rep(0.99, d), eta = rep(1.01, d)) {
  if (!is.vector(mu) | !is.vector(x)) {
    stop('x and mu must be a vector')
  }
  
  if (length(x) != length(mu)) {
    stop('x and mu must have the same dimensions')
  }
  
  if (nrow(Gamma) != ncol(Gamma)) {
    stop('Gamma must be a square matrix')
  }
  if (nrow(Lambda) != ncol(Lambda)) {
    stop('Lambda must be a square matrix')
  }
  if (any(Lambda[upper.tri(Lambda)] != 0, Lambda[lower.tri(Lambda)] != 0)) {
    stop('Lambda must be a diagonal matrix')
  }
  if (nrow(Gamma) != nrow(Lambda)) {
    stop('Gamma and Lambda must have the same dimensions')
  }
  
  if (any(alpha < 0 | alpha > 1)) {
    stop('alpha must be a vector with elements in (0,1)')
  }
  if (any(eta < 1)) {
    stop('eta must be a vector with elements greater than 1')
  }
  
  m <- is.na(x)
  o <- !m
  d <- length(x)
  q <- sum(o)
  A <- Gamma %*% Lambda^(0.5)
  A <- A[o, ]
  if (is.vector(A)) {
    A <- t(as.matrix(A))
  }
  x_o <- x[o]
  mu_o <- mu[o]
  
  func <- function(x) {
    
    vals <- apply(x, 2, function(t) {
      
      chrc_prod <- Reduce('*', lapply(1:d, function(k) {
        good <- alpha[k] * exp(-0.5 * (crossprod(t, A[, k]))^2)
        bad <- (1 - alpha[k]) * exp(-0.5 * eta[k] * (crossprod(t, A[, k]))^2)
        return(good + bad)
      }))
      
      return(
        exp(-1i * crossprod(t, x_o)) * exp(1i * crossprod(t, mu_o)) * chrc_prod
      )
      
    })
    
    return(
      matrix(vals, ncol = ncol(x))
    )
    
  }
  
  hcuba <- hcubature(
    f               = func,
    lowerLimit      = rep(-Inf, q),
    upperLimit      = rep(Inf, q),
    tol             = 1e-10, 
    maxEval         = 500,   
    vectorInterface = TRUE
  )
  
  res <- (2 * pi)^-q * hcuba$integral
  
  if (res < 0) {
    return(0)
  }
  return(res)
}

#---------------------------------------------------------------------------------
# Find f(xi^o | Z_ig = 1, V_ihg = v_ihg), where v_ihg = 1 means xh_is_good = TRUE
# v_ihg = 1 means xh_is_good = FALSE
#---------------------------------------------------------------------------------

find.v <- function(x, h, xh_is_good = TRUE, mu = rep(0, d), Gamma = NULL, Lambda = NULL,alpha = rep(0.99, d), eta = rep(1.01, d)) {
  if (!is.vector(mu) | !is.vector(x)) {
    stop('x and mu must be a vector')
  }
  
  if (length(x) != length(mu)) {
    stop('x and mu must have the same dimensions')
  }
  
  if (nrow(Gamma) != ncol(Gamma)) {
    stop('Gamma must be a square matrix')
  }
  if (nrow(Lambda) != ncol(Lambda)) {
    stop('Lambda must be a square matrix')
  }
  if (any(Lambda[upper.tri(Lambda)] != 0, Lambda[lower.tri(Lambda)] != 0)) {
    stop('Lambda must be a diagonal matrix')
  }
  if (nrow(Gamma) != nrow(Lambda)) {
    stop('Gamma and Lambda must have the same dimensions')
  }
  
  if (any(alpha < 0 | alpha > 1)) {
    stop('alpha must be a vector with elements in (0,1)')
  }
  if (any(eta < 1)) {
    stop('eta must be a vector with elements greater than 1')
  }
  
  m <- is.na(x)
  o <- !m
  d <- length(x)
  q <- sum(o)
  
  if (!(h %in% 1:d)) {
    stop('h must be a whole number from 1 to d')
  }
  
  A <- Gamma %*% Lambda^(0.5)
  A <- A[o, ]
  if (is.vector(A)) {
    A <- t(as.matrix(A))
  }
  x_o <- x[o]
  mu_o <- mu[o]
  
  func <- function(x) {
    
    vals <- apply(x, 2, function(t) {
      
      chrc_prod <- Reduce('*', lapply(1:d, function(k) {
        if (h == k) {
          if (xh_is_good) {
            return(exp(-0.5 * (crossprod(t, A[, k]))^2))
          } else {
            return(exp(-0.5 * eta[k] * (crossprod(t, A[, k]))^2))
          }
        } else {
          good <- alpha[k] * exp(-0.5 * (crossprod(t, A[, k]))^2)
          bad <- (1 - alpha[k]) * exp(-0.5 * eta[k] * (crossprod(t, A[, k]))^2)
          return(good + bad)
        }
      }))
      
      return(
        exp(-1i * crossprod(t, x_o)) * exp(1i * crossprod(t, mu_o)) * chrc_prod
      )
      
    })
    
    return(
      matrix(vals, ncol = ncol(x))
    )
    
  }
  
  hcuba <- hcubature(
    f               = func,
    lowerLimit      = rep(-Inf, q),
    upperLimit      = rep(Inf, q),
    tol             = 1e-10, 
    maxEval         = 500,   
    vectorInterface = TRUE
  )
  
  res <- (2 * pi)^-q * hcuba$integral
  
  if (res < 0) {
    return(0)
  }
  return(res)
}

#---------------------------------
# Find x tilde and x double tilde
#---------------------------------
find.x.tilde <- function(x, h, v_ihg, xh_is_good = TRUE, mu = rep(0, d), Gamma = NULL, Lambda = NULL, alpha = rep(0.99, d), eta = rep(1.01, d)) {
  if (!is.vector(mu) | !is.vector(x)) {
    stop('x and mu must be a vector')
  }
  
  if (length(x) != length(mu)) {
    stop('x and mu must have the same dimensions')
  }
  
  if (nrow(Gamma) != ncol(Gamma)) {
    stop('Gamma must be a square matrix')
  }
  if (nrow(Lambda) != ncol(Lambda)) {
    stop('Lambda must be a square matrix')
  }
  if (any(Lambda[upper.tri(Lambda)] != 0, Lambda[lower.tri(Lambda)] != 0)) {
    stop('Lambda must be a diagonal matrix')
  }
  if (nrow(Gamma) != nrow(Lambda)) {
    stop('Gamma and Lambda must have the same dimensions')
  }
  
  if (any(alpha < 0 | alpha > 1)) {
    stop('alpha must be a vector with elements in (0,1)')
  }
  if (any(eta < 1)) {
    stop('eta must be a vector with elements greater than 1')
  }
  
  m <- is.na(x)
  o <- !m
  d <- length(x)
  q <- sum(o)
  
  if (!(h %in% 1:d)) {
    stop('h must be a whole number from 1 to d')
  }
  
  x_o <- x[o]
  mu_o <- mu[o]
  mu_m <- mu[m]
  
  V <- gtools::permutations(n = 2, r = d, v = 0:1, repeats.allowed = TRUE)
  if (xh_is_good) {
    V <- V[V[, h] == 1, ]
  } else {
    V <- V[V[, h] == 0, ]
  }
  
  x_tilde <- Reduce('+', lapply(1:nrow(V), function(k) {
    W_v <- diag(d)
    diag(W_v) <- sapply(1:d, function(l) {
      ifelse(V[k, l] == 1, 1, eta[l])
    })
    
    Z <- Gamma %*% W_v %*% Lambda %*% t(Gamma)
    Z_mo <- Z[m, o]
    Z_oo_inv <- solve(Z[o, o])
    
    num <- density.v(V[k, ], alpha) * mvtnorm::dmvnorm(x_o, mu_o, as.matrix(Z[o, o]))
    den <- v_ihg
    
    temp <- num / den
    if (is.infinite(temp)) {
      temp <- num
    }
    if (is.nan(temp)) {
      temp <- 0
    }
    
    (mu_m + Z_mo %*% Z_oo_inv %*% (x_o - mu_o)) * temp
  }))
  
  var_x_m <- Reduce('+', lapply(1:nrow(V), function(k) {
    W_v <- diag(d)
    diag(W_v) <- sapply(1:d, function(l) {
      ifelse(V[k, l] == 1, 1, eta[l])
    })
    
    Z <- Gamma %*% W_v %*% Lambda %*% t(Gamma)
    Z_mm <- Z[m, m]
    Z_mo <- Z[m, o]
    Z_om <- Z[o, m]
    Z_oo_inv <- solve(Z[o, o])
    
    num <- density.v(V[k, ], alpha) * mvtnorm::dmvnorm(x_o, mu_o, as.matrix(Z[o, o]))
    den <- v_ihg
    
    temp <- num / den
    if (is.infinite(temp)) {
      temp <- num
    }
    if (is.nan(temp)) {
      temp <- 0
    }
    
    (Z_mm - Z_mo %*% Z_oo_inv %*% Z_om) * temp
  }))
  
  list(x_tilde = x_tilde, var_x_m = var_x_m)
}

#--------------------------------------------
# Find the probablity that V = (v1, ..., vd)
# For instance V = (0, 1, 0, 0, 1) 
#--------------------------------------------
density.v <- function(v, alpha) {
  if (!is.vector(v) | !is.vector(alpha)) {
    stop('v and alpha must be a vector')
  }
  
  if (length(v) != length(alpha)) {
    stop('v and alpha must have the same dimensions')
  }
  d <- length(v)
  
  Reduce('*', lapply(1:d, function(k) {
    alpha[k]^v[k] * (1 - alpha[k])^(1 - v[k])
  }))
}

#--------------------------------------------
# Find
#--------------------------------------------
V.to.W <- function(V, eta) {
  if (!is.vector(eta)) {
    stop('eta must be a vector')
  }
  
  if (!is.matrix(V)) {
    stop('V must be matrix')
  }
  
  d <- length(eta)
  if (d != ncol(V) | 2^d != nrow(V)) {
    stop('Indicator matrix V does not match the number of variates')
  }
  
  inv_weights <- matrix(NA, nrow = 2^d, ncol = d)
  for (h in 1:d) {
    inv_weights[, h] <- ifelse(V[, h] == 1, 1, eta[h])
  }
  
  return(inv_weights)
}

getall <- function(loglik) {
  if (length(loglik) < 3) {
    return(Inf)
  }
  n = length(loglik)
  lm1 = loglik[n]
  lm  = loglik[(n-1)]
  lm_1  = loglik[(n-2)]
  am = (lm1 - lm)/(lm - lm_1)
  lm1.Inf = lm + (lm1 - lm)/(1-am)
  val = lm1.Inf - lm
  if (is.nan(val)) val=0
  if (val < 0) val= 1
  return( val )
}

###########################################################################
###########################################################################
###                                                                     ###
###                          PLR FACTORIZATION                          ###
###                                                                     ###
###########################################################################
###########################################################################

#-----------------------------#
# Strictly elimination matrix #
#-----------------------------#

Lbar <- function(n){
  
  Lb      <- matrix(0,(n*(n-1)/2),n^2)
  indcol  <- sapply(0:(n-2), function(i) 
    sapply((i*n+(i+2)):(n+i*n), function(j) j )
  )
  
  c1      <- 1:(n*(n-1)/2)
  c2      <- unlist(indcol)
  ind     <- cbind(c1,c2)
  Lb[ind] <- 1
  
  return(Lb)
  
}

#-------------------------------------#
# Stricly half-vectorization operator #
#-------------------------------------#

bvech <- function(A){   # A: square marrix
  
  n   <- ncol(A)
  res <- Lbar(n)%*%c(A)
  
  return(res)
  
}

#--------------#
# From Q to Pl #
#--------------#

QPl <- function(Q   # Orthogonal matrix
){
  
  luQ    <- Matrix::lu(Q)
  eluQ   <- Matrix::expand(luQ)
  l      <- bvech(as.matrix(eluQ$L))
  U      <- eluQ$U
  P      <- eluQ$P
  return(list(      # It returns:
    P=P,  # Permutation matrix
    l=l   # Vector of n(n-1)/2 elements
  )
  )
  
}


#--------------#
# From Q to PL #
#--------------#

QPL <- function(Q   # Orthogonal matrix
){
  
  luQ    <- lu(Q)
  eluQ   <- expand(luQ)
  L      <- eluQ$L
  U      <- eluQ$U
  P      <- eluQ$P
  return(list(      # It returns:
    P=P,  # Permutation matrix
    L=L   # Unit lower triangular matrix
  )
  )
  
}

#--------------#
# From Pl to Q #
#--------------#

PlQ <- function(P,  # Permutation matrix
                l   # Vector of n(n-1)/2 elements
){
  
  n  <- (1+sqrt(1+8*length(l)))/2
  L  <- diag(n) + matrix(t(Lbar(n))%*%l,n,n) 
  
  Q <- P %*% L %*% solve(qr.R(qr(P%*%L)), tol = 1e-50) 
  
  return(Q)         # It returns an orthogonal matrix
  
}

#--------------#
# From PL to Q #
#--------------#

PLQ <- function(P,  # Permutation matrix
                L   # Unit lower triangular matrix
){
  
  Q <- P %*% L %*% solve(qr.R(qr(P%*%L)), tol = 1e-50) 
  
  return(Q)         # It returns an orthogonal matrix
  
}

############################################################################
############################################################################
###                                                                      ###
###                    EM ALGORITHM FOR MIXTURE MODEL                    ###
###                                                                      ###
############################################################################
############################################################################

mixture_incomplete_mscn <- function(x, G, max_iter = 30, epsilon = 0.01) {
  if (G < 2) {
    stop("G must be at least 2.")
  }
  
  if(is.null(ncol(x)) | ncol(x) < 2) {
    stop("x for mixture_mvcn must be at least bivariate.")
  }
  
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  
  # Initialization
  d <- ncol(x)
  n <- nrow(x)
  py <- rep(1/G, G)
  # alpha <- matrix(0.6, nrow = d, ncol = G)
  # eta <- matrix(1.4, nrow = d, ncol = G)
  alpha <- matrix(0.999, nrow = d, ncol = G)
  eta   <- matrix(1.001, nrow = d, ncol = G)
  x_o <- x[complete.cases(x), ]
  init <- pam(x_o, G)
  mu <- init$medoids
  Sigma <- array(NA, dim = c(d, d, G))
  Lambda <- array(NA, dim = c(d, d, G))
  Gama <- array(NA, dim = c(d, d, G))
  for (g in 1:G) {
    Sigma[, , g] <- var(x_o[init$clustering == g, ])
    spec.decomp <- eigen(Sigma[, , g])
    Gama[, , g] <- spec.decomp$vectors
    Lambda[, , g] <- diag(spec.decomp$values)
  }
  z <- matrix(NA, nrow = n, ncol = G)
  z_tilde <- matrix(NA, nrow = n, ncol = G)
  v_good <- array(NA, dim = c(n, d, G))
  v_tilde <- array(NA, dim = c(n, d, G))
  v_hat <- array(NA, dim = c(n, 2^d, G))
  V <- permutations(n = 2, r = d, v = 0:1, repeats.allowed = TRUE)
  W <- array(NA, dim = c(2^d, d, G))
  A <- array(0, dim = c(d, d, n, G))
  B <- array(0, dim = c(n, d, G))
  E <- array(0, dim = c(d, d, n, G))
  sA <- array(NA, dim = c(d, d, G))
  sB <- array(NA, dim = c(1, d, G))
  sE <- array(NA, dim = c(d, d, G))
  sigma_tilde <- array(NA, dim = c(d, d, n, 2^d, G))
  
  x_tilde <- array(rep(x, d * G), dim = c(n, d, 2^d, G))
  s <- array(NA, dim = c(n, d, G))
  
  iter <- 0
  l <- NULL
  
  while (iter < max_iter &  getall(l) > epsilon) {
    
    tim <- system.time({
      
      ### E-step
      for (g in 1:G) {
        for (i in 1:n) {
          # cat('\n', iter + 1, g, i)
          m <- is.na(x[i, ])
          if (any(m)) {
            z[i, g] <- find.z(x[i, ], mu = mu[g, ], Gamma = Gama[, , g], Lambda = Lambda[, , g],
                              alpha = alpha[, g], eta = eta[, g])
            
            
            for (h in 1:d) {
              v_good[i, h, g] <- find.v(x[i, ], h = h, xh_is_good = TRUE, mu = mu[g, ],
                                        Gamma = Gama[, , g], Lambda = Lambda[, , g],
                                        alpha = alpha[, g], eta = eta[, g])
              v_tilde[i, h, g] <- alpha[h, g] * v_good[i, h, g] / z[i, g]
            }
          } else {
            z[i, g] <- dMSCN(x[i, ], mu = mu[g, ], Gamma = Gama[, , g], Lambda = Lambda[, , g],
                             alpha = alpha[, g], eta = eta[, g])
            
            y <- t(Gama[, , g]) %*% (x[i, ] - mu[g, ])
            
            for (h in 1:d) {
              good_norm <- alpha[h, g] * dnorm(y[h], mean = 0, sd = sqrt(Lambda[h, h, g]))
              bad_norm <- (1 - alpha[h, g]) * dnorm(y[h], mean = 0, sd = sqrt(eta[h, g] * Lambda[h, h, g]))
              v_tilde[i, h, g] <- good_norm / (good_norm + bad_norm)
            }
          }
        }
      }
      
    })[3]
    
    # cat(tim, '\n')
    
    v_tilde[is.nan(v_tilde)] <- 0
    v_tilde[is.infinite(v_tilde)] <- 1
    
    for (i in 1:n) {
      ss <- py * z[i, ]
      if (sum(ss) == 0) {
        ss <- rep(1, G)
      }
      z_tilde[i, ] <- ss / sum(ss)
    }
    N <- colSums(z_tilde)
    
    for (g in 1:G) {
      W[, , g] <- V.to.W(V, eta = eta[, g])
      
      for (i in 1:n) {
        m <- is.na(x[i, ])
        
        if (any(m)) {
          o <- !m
          
          x_o <- x[i, o]
          x_m <- x[i, m]
          mu_o <- mu[g, o]
          mu_m <- mu[g, m]
          
          for (k in 1:(2^d)) {
            M <- Gama[, , g] %*% diag(W[k, , g]) %*% Lambda[, , g] %*% t(Gama[, , g])
            M_mo <- M[m, o]
            M_om <- t(M_mo)
            M_oo_inv <- solve(M[o, o])
            M_mm <- M[m, m]
            
            v_hat[i, k, g] <- density.v(V[k, ], alpha[, g]) * mvtnorm::dmvnorm(x_o, mu_o, as.matrix(M[o, o])) / z[i, g]
            x_tilde[i, m, k, g] <- mu_m + M_mo %*% M_oo_inv %*% (x_o - mu_o)
            
            sigma_tilde[o, o, i, k, g] <- tcrossprod(x_o - mu_o)
            sigma_tilde[o, m, i, k, g] <- tcrossprod(x_o - mu_o, x_tilde[i, m, k, g] - mu_m)
            sigma_tilde[m, o, i, k, g] <- t(sigma_tilde[o, m, i, k, g])
            
            sigma_tilde[m, m, i, k, g] <- tcrossprod(x_tilde[i, m, k, g] - mu_m)
            sigma_tilde[m, m, i, k, g] <- sigma_tilde[m, m, i, k, g] + M_mm - M_mo %*% M_oo_inv %*% M_om
          }
        } else {
          for (k in 1:(2^d)) {
            Z <- Gama[, , g] %*% diag(W[k, , g]) %*% Lambda[, , g] %*% t(Gama[, , g])
            v_hat[i, k, g] <- density.v(V[k, ], alpha[, g]) * mvtnorm::dmvnorm(x[i, ], mu[g, ], as.matrix(Z)) / z[i, g]
          }
          
          sigma_tilde[, , i, , g] <- rep(tcrossprod(x[i, ] - mu[g, ]), 2^d)
        }
        
        v_hat[i, is.nan(v_hat[i, , g]), g] <- 0 
        v_hat[i, is.infinite(v_hat[i, , g]), g] <- 1 
        
        A[, , i, g] <- Reduce('+', lapply(1:(2^d), function(k) {
          inv_W <- diag(W[k, , g]^-1)
          return(v_hat[i, k, g] * inv_W)
        }))
        
        B[i, , g] <- Reduce('+', lapply(1:(2^d), function(k) {
          inv_W <- diag(W[k, , g]^-1)
          return(v_hat[i, k, g] * inv_W %*% t(Gama[, , g]) %*%  x_tilde[i, , k, g])
        }))
        
        E[, , i, g] <- Reduce('+', lapply(1:(2^d), function(k) {
          inv_W <- diag(W[k, , g]^-1)
          return(v_hat[i, k, g] * t(Gama[, , g]) %*% sigma_tilde[, , i, k, g] %*% Gama[, , g] %*% inv_W)
        }))
        
      }
      
      sA[, , g] <- Reduce('+', lapply(1:n, function(i) {
        return(z_tilde[i, g] * A[, , i, g])
      }))
      
      sB[, , g] <- Reduce('+', lapply(1:n, function(i) {
        return(z_tilde[i, g] * B[i, , g])
      }))
      
      sE[, , g] <- Reduce('+', lapply(1:n, function(i) {
        return(z_tilde[i, g] * E[, , i, g])
      }))
    }
    
    ### M-step
    for (g in 1:G) {
      py[g] <- N[g] / n
      
      for (h in 1:d) {
        alpha[h, g] <- sum(z_tilde[, g] * v_tilde[, h, g]) / N[g]
        alpha[h, g] <- max(alpha[h, g], 0.5)
        alpha[h, g] <- min(alpha[h, g], 1)
      }
      
      mu[g, ] <- Gama[, , g] %*% solve(sA[, , g]) %*% sB[, , g]
      Lambda[, , g] <- sE[, , g]
      Lambda[, , g] <- diag(diag(Lambda[, , g])) / N[g]
      
      eta[, g] <- multiroot(function(x) {
        W[, , g] <- V.to.W(V, eta = x)
        
        H <- array(0, dim = c(d, d, n))
        sH <- array(NA, dim = c(d, d))
        
        EE <- array(0, dim = c(d, d, n))
        sEE <- array(NA, dim = c(d, d))
        
        for (i in 1:n) {
          H[, , i] <- Reduce('+', lapply(1:(2^d), function(k) {
            inv_W <- diag(W[k, , g]^-1)
            dW <- diag(ifelse(V[k, ] == 1, 0, 1))
            return(v_hat[i, k, g] * inv_W %*% dW)
          }))
          
          EE[, , i] <- Reduce('+', lapply(1:(2^d), function(k) {
            inv_W <- diag(W[k, , g]^-1)
            dW <- diag(ifelse(V[k, ] == 1, 0, 1))
            res <- v_hat[i, k, g] * inv_W %*% solve(Lambda[, , g]) %*% t(Gama[, , g])
            res <- res %*% sigma_tilde[, , i, k, g] %*% Gama[, , g] %*% inv_W %*% dW
            return(res)
          }))
        }
        
        sH <- Reduce('+', lapply(1:n, function(i) {
          return(z_tilde[i, g] * H[, , i])
        }))
        
        sEE <- Reduce('+', lapply(1:n, function(i) {
          return(z_tilde[i, g] * EE[, , i])
        }))
        
        sEE <- diag(diag(sEE))
        diag(sH - sEE)
      }, start = eta[, g], maxiter = 100)$root
      
      for (h in 1:d) {
        if (eta[h, g] < 1.001) {
          eta[h, g] <- 1.001
          alpha[h, g] <- 0.999
        }
      }
      # eta[eta[, g] < 1.001, g] <- 1.001
      
      Gama0 <- Gama[, , g]
      tempfac <- QPl(Gama0)
      P       <- tempfac$P
      initial.values <- tempfac$l
      
      obj_f <- function(par, P, z_tilde, v_hat, W, Lambda, sigma_tilde) {
        # Gama <- matrix(par, nrow = d, ncol = d)
        Gama <- PlQ(P = P, l = par)
        
        obj <- 0
        for (i in 1:n) {
          temp <- Reduce('+', lapply(1:(2^d), function(k) {
            inv_W <- diag(W[k, ]^-1)
            Trce <- tr(Gama %*% inv_W %*% solve(Lambda) %*% t(Gama) %*% sigma_tilde[, , i, k])
            return(v_hat[i, k] * Trce)
          }))
          
          obj <- obj + z_tilde[i] * temp
        }
        
        return(obj)
      }
      
      resl <- optim(par = initial.values, fn = obj_f, method = 'BFGS',
                    P = P, z_tilde = z_tilde[, g], v_hat = v_hat[, , g], W = W[, , g],
                    Lambda = Lambda[, , g], sigma_tilde = sigma_tilde[, , , , g], control = list(maxit = 500))
      
      f.obj <- resl$value
      est   <- resl$par
      
      Gama[, , g] <- PlQ(P = P, l = est)
      
    }
    
    # Calculate loglikelihood
    fx <- matrix(0, nrow = nrow(x), ncol = G)
    for (g in 1:G) {
      for (i in 1:n) {
        m <- is.na(x[i, ])
        if (any(m)) {
          fx[i, g] <- find.z(x[i, ], mu = mu[g, ], Gamma = Gama[, , g], Lambda = Lambda[, , g],
                             alpha = alpha[, g], eta = eta[, g])
        } else {
          fx[i, g] <- dMSCN(x[i, ], mu = mu[g, ], Gamma = Gama[, , g], Lambda = Lambda[, , g],
                            alpha = alpha[, g], eta = eta[, g])
        }
      }
    }
    lik <- apply(fx, 1, function(xi, wt) {
      return(sum(wt * xi))
    }, wt = py)
    lik[lik <= 10^(-323)] <- 10^(-323)
    loglik <- sum(log(lik))
    l <- c(l, loglik)
    plot(l, type = 'b', pch = 16)
    # print(l)
    # print(sort(l))
    # print(l == sort(l))
    iter <- iter + 1
    # cat('Iteration: ', iter, all(l == sort(l)), '\n') 
  }
  
  # update cluster membership
  clust <- apply(z_tilde, 1, which.max)
  
  # identify bad observations
  
  detect1 <- array(1/n,c(n, d, G), dimnames = list(1:n, paste("dim_",1:d, sep=""), paste("cluster_", 1:G, sep = "")))
  detect2 <- array(1/n,c(n, d), dimnames = list(1:n, paste("dim_", 1:d, sep="")))
  for(h in 1:d){
    for(i in 1:n){
      detect2[i, h] <- ifelse(v_tilde[i, h, clust[i]] > 1/2, "*", "bad")
    }
    for(g in 1:G){
      detect1[, h, g] <- ifelse(v_tilde[, h, g] > 1/2, "*", "bad")
    }
  }
  
  outlier <- rep(0, n)
  for(i in 1:n) {
    if (any(v_tilde[i, , clust[i]] < 0.5)) {
      outlier[i] <- 1
    }
  }
  
  out <- list(mu = mu, Gamma = Gama, Lambda = Lambda,pi = py, alpha = alpha, 
              eta = eta, cluster = clust, loglikelihood = l, 
              iteration = iter, detect1 = detect1, detect2 = detect2, outlier = outlier, v_tilde = v_tilde)
  return(out)
}
