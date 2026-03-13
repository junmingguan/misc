fastica_r1update_1iter = function(X,w){
  w= w/sqrt(sum(w^2))
  P = t(X) %*% w
  G = tanh(P)
  G2 = 1-tanh(P)^2
  w = X %*% G - sum(G2) * w
  w = w/sqrt(sum(w^2))
  return(w)
}

fastica_r1update = function(init_w, X1, n_iter = 100) {
  w = init_w
  for(i in 1:n_iter)
  {
    w = fastica_r1update_1iter(X1,w)
  }
  return(w)
}
preprocess = function(X, n.comp=10, center = TRUE){
  n <- nrow(X)
  p <- ncol(X)
  if (center) X <- scale(X, scale = FALSE)
  X <- t(X)

  ## This appears to be equivalant to X1 = t(svd(X)$v[,1:n.comp])
  V <- X %*% t(X)/n
  s <- La.svd(V)
  D <- diag(c(1/sqrt(s$d)))
  K <- D %*% t(s$u)
  K <- matrix(K[1:n.comp, ], n.comp, p)
  X1 <- K %*% X
  return(X1)
}

preprocess2 = function(X, n.comp=10, center = FALSE){
  n <- nrow(X)
  p <- ncol(X)
  if (center) X <- scale(X, scale = FALSE)
  # X <- t(X)

  svd.X <- svd(X)
  X1 <- svd.X$u[,1:n.comp] %*% diag(rep(1, n.comp)) %*% t(svd.X$v[,1:n.comp])
  return(X1)
}


ebica_parallel <- function(X, K, d, n, max_iter, tol, s, S_init) {
  prior_probs <- numeric(K)

  if (!is.null(S_init)) {
    W <- X %*% S_init
  } else {
    W <- matrix(rnorm(d * K), nrow = d, ncol = K)
  }

  if (K > 1) {
    svd_res <- svd(W)
    W <- svd_res$u %*% t(svd_res$v)
  } else {
    W <- W / sqrt(sum(W^2))
  }

  for (iter in 1:max_iter) {
    W_old <- W
    Y <- t(W) %*% X

    S_plus <- matrix(0, nrow = n, ncol = K)
    for (k in 1:K) {
      y_k <- if (K == 1) as.numeric(Y) else Y[k, ]

      ebnm_fit <- ebnm_rademacher(x = y_k, s = s)
      S_plus[, k] <- ebnm_fit$posterior$mean

      prior_probs[k] <- ebnm_fit$fitted_g$pi[1]
    }

    W_plus <- X %*% S_plus

    is_converged <- FALSE
    if (K > 1) {
      svd_res <- svd(W_plus)
      W <- svd_res$u %*% t(svd_res$v)
      if (max(1 - abs(colSums(W * W_old))) < tol) {
        is_converged <- TRUE
      }
    } else {
      W <- W_plus / sqrt(sum(W_plus^2))
      if (1 - abs(sum(W * W_old)) < tol) {
        is_converged <- TRUE
      }
    }

    if (is_converged || iter == max_iter) {
      if (is_converged) cat("Converged in", iter, "iterations.\n")
      break
    }
  }

  return(list(
    W = W,
    S = t(X) %*% W,
    S_plus = S_plus,
    prior_probs = prior_probs
  ))
}

ebica_deflation <- function(X, K, d, n, max_iter, tol, s, S_init) {
  prior_probs <- numeric(K)
  S_plus <- matrix(0, nrow = n, ncol = K)

  W_found <- matrix(0, nrow = d, ncol = K)

  for (k in 1:K) {
    if (!is.null(S_init)) {
      # w <- as.numeric(X %*% S_init)
      w <- as.numeric(X %*% S_init[, k])
    } else {
      w <- rnorm(d)
    }
    w <- w / sqrt(sum(w^2))

    for (iter in 1:max_iter) {
      w_old <- w
      y <- as.numeric(t(w) %*% X)

      ebnm_fit <- ebnm_rademacher(x = y, s = s)
      s_plus <- ebnm_fit$posterior$mean

      current_pi <- ebnm_fit$fitted_g$pi[1]

      w_plus <- as.numeric(X %*% s_plus)

      if (k > 1) {
        A <- as.matrix(W_found[, 1:(k-1)])
        w_plus <- w_plus - A %*% (t(A) %*% w_plus)
      }

      w <- as.numeric(w_plus / sqrt(sum(w_plus^2)))

      if (1 - abs(sum(w * w_old)) < tol || iter == max_iter) {
        if (1 - abs(sum(w * w_old)) < tol) {
          cat("Component", k, "converged in", iter, "iterations.\n")
        }

        prior_probs[k] <- current_pi
        S_plus[, k] <- s_plus
        break
      }
    }
    W_found[, k] <- w
  }

  return(list(
    W = W_found,
    S = t(X) %*% W_found,
    S_plus = S_plus,
    prior_probs = prior_probs
  ))
}

ebica <- function(X, K = 1, alg_typ = "parallel", max_iter = 1000,
                  tol = 1e-6, s = 1.0, S_init = NULL) {
  d <- nrow(X)
  n <- ncol(X)

  if (!alg_typ %in% c("parallel", "deflation")) {
    stop("alg_typ must be either 'parallel' or 'deflation'.")
  }
  if (K > d) {
    stop("K cannot be greater than the number of features (d).")
  }

  if (!is.null(S_init)) {
    if (is.vector(S_init)) S_init <- matrix(S_init, ncol = 1)

    if (K == 1) {
      if (nrow(S_init) != n || ncol(S_init) != K) {
        stop("S_init must be a matrix with n rows and K columns.")
      }
    }
  }

  if (alg_typ == "parallel" || K == 1) {
    return(ebica_parallel(X, K, d, n, max_iter, tol, s, S_init))
  } else {
    return(ebica_deflation(X, K, d, n, max_iter, tol, s, S_init))
  }
}



# new

#' Compute the Normal means log likelihood
#'
#' @param x A vector of observations in an EBNM problem
#' @param s A vector of standard errors in an EBNM problem
#' @param Et A vector of posterior means
#' @param Et2 A vector of posterior mean of second moments
#'
#' @references This function is copied from flashier:::normal.means.loglik
#' (https://github.com/willwerscheid/flashier/blob/master/R/objective.R)

normal_means_loglik <- function(x, s, Et, Et2) {
  idx <- is.finite(s) & s > 0
  x <- x[idx]
  s <- s[idx]
  Et <- Et[idx]
  Et2 <- Et2[idx]

  return(-0.5 * sum(log(2 * pi * s^2) + (1 / s^2) * (Et2 - 2 * x * Et + x^2)))
}

ebica_parallel <- function(X, K, d, n, max_iter, tol, s, S_init, conv_crit, symmetric_rademacher) {
  prior_probs <- numeric(K)
  obj_history <- numeric(max_iter)

  if (!is.null(S_init)) {
    W <- X %*% S_init
  } else {
    W <- matrix(rnorm(d * K), nrow = d, ncol = K)
  }

  if (K > 1) {
    svd_res <- svd(W)
    W <- svd_res$u %*% t(svd_res$v)
  } else {
    W <- matrix(W / sqrt(sum(W^2)), ncol = 1)
  }

  for (iter in 1:max_iter) {
    W_old <- W
    Y <- t(W) %*% X

    S_plus <- matrix(0, nrow = n, ncol = K)
    KL_total <- 0
    sum_E_S2 <- 0

    for (k in 1:K) {
      y_k <- Y[k, ]

      if (symmetric_rademacher) {
        ebnm_fit <- ebnm_symm_rademacher(x = y_k, s = s)
      } else {
        ebnm_fit <- ebnm_rademacher(x = y_k, s = s)
      }

      post_mean <- ebnm_fit$posterior$mean
      post_sec_moment <- ebnm_fit$posterior$sd^2 + post_mean^2

      S_plus[, k] <- post_mean
      prior_probs[k] <- ebnm_fit$fitted_g$pi[1]

      Eloglik <- normal_means_loglik(x = y_k, s = s,
                                            Et = post_mean,
                                            Et2 = post_sec_moment)

      KL_k <- Eloglik - ebnm_fit$log_likelihood

      KL_total <- KL_total + KL_k
      sum_E_S2 <- sum_E_S2 + sum(post_sec_moment)
    }

    W_plus <- X %*% S_plus

    R2 <- sum(X^2) - 2 * sum(as.vector(Y) * as.vector(t(S_plus))) + sum_E_S2
    global_Eloglik <- -0.5 * (d * n) * log(2 * pi * s^2) - R2 / (2 * s^2)
    elbo <- global_Eloglik - KL_total
    obj_history[iter] <- elbo


    is_converged <- FALSE
    if (conv_crit == "elbo") {
      if (iter > 1 && (elbo - obj_history[iter - 1]) < tol) is_converged <- TRUE
    } else {
      if (K > 1) {
        svd_res <- svd(W_plus)
        W <- svd_res$u %*% t(svd_res$v)
        if (max(1 - abs(colSums(W * W_old))) < tol) is_converged <- TRUE
      } else {
        W <- matrix(W_plus / sqrt(sum(W_plus^2)), ncol = 1)
        if (1 - abs(sum(W * W_old)) < tol) is_converged <- TRUE
      }
    }

    if (conv_crit == "elbo") {
      if (K > 1) {
        svd_res <- svd(W_plus)
        W <- svd_res$u %*% t(svd_res$v)
      } else {
        W <- matrix(W_plus / sqrt(sum(W_plus^2)), ncol = 1)
      }
    }

    if (is_converged || iter == max_iter) {
      if (is_converged) cat("Converged in", iter, "iterations. Final ELBO:", elbo, "\n")
      break
    }
  }

  return(list(
    W = W,
    S = t(X) %*% W,
    S_plus = S_plus,
    prior_probs = prior_probs,
    obj = obj_history[1:iter]
  ))
}

ebica_deflation <- function(X, K, d, n, max_iter, tol, s, S_init, conv_crit, symmetric_rademacher) {
  prior_probs <- numeric(K)
  S_plus <- matrix(0, nrow = n, ncol = K)
  W_found <- matrix(0, nrow = d, ncol = K)
  obj_history <- list()

  KL_accumulated <- 0
  sum_E_S2_accumulated <- 0
  cross_term_accumulated <- 0

  sum_X2 <- sum(X^2)

  for (k in 1:K) {
    if (!is.null(S_init)) {
      w <- as.numeric(X %*% S_init[, k])
    } else {
      w <- rnorm(d)
    }
    w <- w / sqrt(sum(w^2))

    comp_obj <- numeric(max_iter)

    for (iter in 1:max_iter) {
      w_old <- w
      y <- as.numeric(t(w) %*% X)

      if (symmetric_rademacher) {
        ebnm_fit <- ebnm_symm_rademacher(x = y, s = s)
      } else {
        ebnm_fit <- ebnm_rademacher(x = y, s = s)
      }

      post_mean <- ebnm_fit$posterior$mean
      post_sec_moment <- ebnm_fit$posterior$sd^2 + post_mean^2

      Eloglik_1D <- normal_means_loglik(x = y, s = s,
                                            Et = post_mean,
                                            Et2 = post_sec_moment)
      KL_k <- Eloglik_1D - ebnm_fit$log_likelihood

      current_total_KL <- KL_accumulated + KL_k
      current_total_E_S2 <- sum_E_S2_accumulated + sum(post_sec_moment)
      current_cross_term <- cross_term_accumulated + 2 * sum(y * post_mean)

      R2 <- sum_X2 - current_cross_term + current_total_E_S2
      global_Eloglik <- -0.5 * (d * n) * log(2 * pi * s^2) - R2 / (2 * s^2)

      global_elbo <- global_Eloglik - current_total_KL
      comp_obj[iter] <- global_elbo

      s_plus <- post_mean
      current_pi <- ebnm_fit$fitted_g$pi[1]
      w_plus <- as.numeric(X %*% s_plus)

      if (k > 1) {
        A <- as.matrix(W_found[, 1:(k-1)])
        w_plus <- w_plus - A %*% (t(A) %*% w_plus)
      }

      w <- as.numeric(w_plus / sqrt(sum(w_plus^2)))

      is_converged <- FALSE
      if (conv_crit == "elbo") {
        if (iter > 1 && (global_elbo - comp_obj[iter - 1]) < tol) is_converged <- TRUE
      } else {
        if (1 - abs(sum(w * w_old)) < tol) is_converged <- TRUE
      }

      if (is_converged || iter == max_iter) {
        if (is_converged) cat("Component", k, "converged in", iter, "iterations. Global ELBO:", global_elbo, "\n")

        KL_accumulated <- KL_accumulated + KL_k
        sum_E_S2_accumulated <- sum_E_S2_accumulated + sum(post_sec_moment)
        cross_term_accumulated <- cross_term_accumulated + 2 * sum(y * post_mean)

        prior_probs[k] <- current_pi
        S_plus[, k] <- s_plus
        break
      }
    }
    W_found[, k] <- w
    obj_history[[k]] <- comp_obj[1:iter]
  }

  return(list(
    W = W_found,
    S = t(X) %*% W_found,
    S_plus = S_plus,
    prior_probs = prior_probs,
    obj = obj_history
  ))
}

ebica <- function(X, K = 1, alg_typ = "parallel", symmetric_rademacher = FALSE,
                  conv_crit = "weights", max_iter = 1000, tol = 1e-6, s = 1.0, S_init = NULL) {
  d <- nrow(X)
  n <- ncol(X)

  if (!alg_typ %in% c("parallel", "deflation")) {
    stop("alg_typ must be either 'parallel' or 'deflation'.")
  }
  if (!conv_crit %in% c("weights", "elbo")) {
    stop("conv_crit must be either 'weights' or 'elbo'.")
  }
  if (!is.logical(symmetric_rademacher)) {
    stop("symmetric_rademacher must be TRUE or FALSE.")
  }
  if (K > d) {
    stop("K cannot be greater than the number of features (d).")
  }

  if (!is.null(S_init)) {
    if (is.vector(S_init)) S_init <- matrix(S_init, ncol = 1)
    if (K == 1) {
      if (nrow(S_init) != n || ncol(S_init) != K) {
        stop("S_init must be a matrix with n rows and K columns.")
      }
    }
  }

  if (alg_typ == "parallel" || K == 1) {
    return(ebica_parallel(X, K, d, n, max_iter, tol, s, S_init, conv_crit, symmetric_rademacher))
  } else {
    return(ebica_deflation(X, K, d, n, max_iter, tol, s, S_init, conv_crit, symmetric_rademacher))
  }
}
