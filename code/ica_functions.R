fastica_rank_one_update_1iter = function(X, w, include_G2=TRUE){
  w= w/sqrt(sum(w^2))
  P = t(X) %*% w
  G = tanh(P)
  G2 = 0
  if(include_G2) # Note: if you don't include G2, this is like the EBCD updates rather than ICA and I believe it will always try to maximize logcosh
    {G2 = 1-tanh(P)^2}
  w = X %*% G - sum(G2) * w
  w = w/sqrt(sum(w^2))
  return(w)
}

fastica_rank_one_update = function(X1, init_w, include_G2 = TRUE, n_iter = 100) {
  w = init_w
  for(i in 1:n_iter)
  {
    w = fastica_rank_one_update_1iter(X1, w, include_G2)
  }
  return(w)
}

fastica_rank_K_update_1iter = function(X, W, include_G2 = TRUE) {
  W <- orthonormalize_W(W)
  P <- t(X) %*% W
  G <- tanh(P)
  G2_sum <- rep(0, ncol(W))
  if (include_G2) {
    G2_sum <- colSums(1 - tanh(P)^2)
  }
  W_plus <- X %*% G - W %*% diag(G2_sum, nrow = length(G2_sum))
  orthonormalize_W(W_plus)
}

fastica_rank_K_update = function(X1, init_W = NULL, K = NULL, include_G2 = TRUE, n_iter = 100) {
  if (is.null(init_W)) {
    if (is.null(K)) {
      stop("Either init_W or K must be supplied.")
    }
    init_W <- matrix(rnorm(nrow(X1) * K), nrow = nrow(X1), ncol = K)
  } else {
    if (is.vector(init_W)) init_W <- matrix(init_W, ncol = 1)
    if (!is.null(K) && ncol(init_W) != K) {
      stop("If supplied, K must match ncol(init_W).")
    }
    K <- ncol(init_W)
  }

  if (K > nrow(X1)) {
    stop("K cannot be greater than nrow(X1).")
  }
  if (any(!is.finite(init_W))) {
    stop("init_W must contain only finite values.")
  }
  if (any(colSums(init_W^2) == 0)) {
    stop("Each column of init_W must have nonzero norm.")
  }

  W <- orthonormalize_W(init_W)
  for (i in 1:n_iter) {
    W <- fastica_rank_K_update_1iter(X1, W, include_G2)
  }
  W
}

fastica_rank_one_update_noncentered = function(init_w, X, n.comp = 10, n_iter = 100) {
  w = init_w
  X1 <- preprocess(X)
  X1 = rbind(X1, rep(1,ncol(X1)))
  for(i in 1:n_iter)
  {
    w = fastica_rank_one_update_1iter(X1,w)
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
  X1 <- svd.X$u[,1:n.comp] %*% t(svd.X$v[,1:n.comp])
  return(X1)
}


compute_objective = function(X,w){
  P = t(X) %*% w
  return(mean(log(cosh(P))))
}


compute_objective_s = function(s){
  return(mean(log(cosh(s))))
}



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

check_ebica_inits <- function(S_init, W_init, n, d, K) {
  if (!is.null(S_init) && !is.null(W_init)) {
    stop("Only one of S_init and W_init can be supplied.")
  }

  if (!is.null(S_init)) {
    if (is.vector(S_init)) S_init <- matrix(S_init, ncol = 1)
    if (any(!is.finite(S_init))) {
      stop("S_init must contain only finite values.")
    }
    if (nrow(S_init) != n || ncol(S_init) != K) {
      stop("S_init must be a matrix with n rows and K columns.")
    }
  }

  if (!is.null(W_init)) {
    if (is.vector(W_init)) W_init <- matrix(W_init, ncol = 1)
    if (any(!is.finite(W_init))) {
      stop("W_init must contain only finite values.")
    }
    if (nrow(W_init) != d || ncol(W_init) != K) {
      stop("W_init must be a matrix with d rows and K columns.")
    }
    if (any(colSums(W_init^2) == 0)) {
      stop("Each column of W_init must have nonzero norm.")
    }
  }

  list(S_init = S_init, W_init = W_init)
}

orthonormalize_W <- function(W) {
  if (ncol(W) > 1) {
    svd_res <- svd(W)
    return(svd_res$u %*% t(svd_res$v))
  }
  matrix(W / sqrt(sum(W^2)), ncol = 1)
}

init_residual_sd <- function(X, s) {
  if (is.null(s)) {
    return(max(sqrt(mean(X^2)), sqrt(.Machine$double.eps)))
  }
  if (length(s) != 1 || !is.finite(s) || s <= 0) {
    stop("s must be NULL or a single positive finite value.")
  }
  s
}

estimate_residual_sd <- function(R2, d, n) {
  max(sqrt(R2 / (d * n)), sqrt(.Machine$double.eps))
}

residual_variance_path <- function(residual_sd) {
  if (is.list(residual_sd)) {
    return(lapply(residual_sd, function(x) x^2))
  }
  residual_sd^2
}

get_prior_prob <- function(ebnm_fit) {
  if (!is.null(ebnm_fit$fitted_g$pi)) {
    return(ebnm_fit$fitted_g$pi[1])
  }
  NA_real_
}

check_ebnm_fn_list <- function(ebnm_fn, K) {
  if (is.function(ebnm_fn)) {
    return(rep(list(ebnm_fn), K))
  }
  if (!is.list(ebnm_fn) || length(ebnm_fn) != K || !all(vapply(ebnm_fn, is.function, logical(1)))) {
    stop("ebnm_fn must be a function or a list of K functions.")
  }
  ebnm_fn
}

ebica_parallel <- function(X, K, d, n, max_iter, tol, s, estimate_s, S_init, W_init,
                           conv_crit, symmetric_rademacher, verbose) {
  prior_probs <- numeric(K)
  obj_history <- numeric(max_iter)
  residual_sd <- numeric(max_iter)
  s_current <- init_residual_sd(X, s)

  if (!is.null(W_init)) {
    W <- W_init
  } else if (!is.null(S_init)) {
    W <- X %*% S_init
  } else {
    W <- matrix(rnorm(d * K), nrow = d, ncol = K)
  }
  W <- orthonormalize_W(W)

  for (iter in 1:max_iter) {
    W_old <- W
    Y <- t(W) %*% X

    S_plus <- matrix(0, nrow = n, ncol = K)
    KL_total <- 0
    sum_E_S2 <- 0

    for (k in 1:K) {
      y_k <- Y[k, ]

      if (symmetric_rademacher) {
        ebnm_fit <- ebnm_rademacher(x = y_k, s = s_current)
      } else {
        ebnm_fit <- ebnm_generalized_rademacher(x = y_k, s = s_current)
      }

      post_mean <- ebnm_fit$posterior$mean
      post_sec_moment <- ebnm_fit$posterior$sd^2 + post_mean^2

      S_plus[, k] <- post_mean
      prior_probs[k] <- get_prior_prob(ebnm_fit)

      Eloglik <- normal_means_loglik(x = y_k, s = s_current,
                                            Et = post_mean,
                                            Et2 = post_sec_moment)

      KL_k <- Eloglik - ebnm_fit$log_likelihood

      KL_total <- KL_total + KL_k
      sum_E_S2 <- sum_E_S2 + sum(post_sec_moment)
    }

    W_plus <- X %*% S_plus

    R2 <- sum(X^2) - 2 * sum(as.vector(Y) * as.vector(t(S_plus))) + sum_E_S2
    if (estimate_s) {
      s_current <- estimate_residual_sd(R2, d, n)
    }
    residual_sd[iter] <- s_current
    global_Eloglik <- -0.5 * (d * n) * log(2 * pi * s_current^2) - R2 / (2 * s_current^2)
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
      if (is_converged && verbose > 0) cat("Converged in", iter, "iterations. Final ELBO:", elbo, "\n")
      break
    }
  }

  return(list(
    W = W,
    S = t(X) %*% W,
    S_plus = S_plus,
    prior_probs = prior_probs,
    residual_sd = s_current,
    residual_variance = s_current^2,
    residual_sd_path = residual_sd[1:iter],
    residual_variance_path = residual_variance_path(residual_sd[1:iter]),
    obj = obj_history[1:iter]
  ))
}

ebica_deflation <- function(X, K, d, n, max_iter, tol, s, estimate_s, S_init, W_init,
                            conv_crit, symmetric_rademacher, verbose) {
  prior_probs <- numeric(K)
  S_plus <- matrix(0, nrow = n, ncol = K)
  W_found <- matrix(0, nrow = d, ncol = K)
  obj_history <- list()
  residual_sd <- list()
  s_current <- init_residual_sd(X, s)

  KL_accumulated <- 0
  sum_E_S2_accumulated <- 0
  cross_term_accumulated <- 0

  sum_X2 <- sum(X^2)

  for (k in 1:K) {
    if (!is.null(W_init)) {
      w <- as.numeric(W_init[, k])
    } else if (!is.null(S_init)) {
      w <- as.numeric(X %*% S_init[, k])
    } else {
      w <- rnorm(d)
    }
    w <- w / sqrt(sum(w^2))

    comp_obj <- numeric(max_iter)
    comp_residual_sd <- numeric(max_iter)

    for (iter in 1:max_iter) {
      w_old <- w
      y <- as.numeric(t(w) %*% X)

      if (symmetric_rademacher) {
        ebnm_fit <- ebnm_rademacher(x = y, s = s_current)
      } else {
        ebnm_fit <- ebnm_generalized_rademacher(x = y, s = s_current)
      }

      post_mean <- ebnm_fit$posterior$mean
      post_sec_moment <- ebnm_fit$posterior$sd^2 + post_mean^2

      Eloglik_1D <- normal_means_loglik(x = y, s = s_current,
                                            Et = post_mean,
                                            Et2 = post_sec_moment)
      KL_k <- Eloglik_1D - ebnm_fit$log_likelihood

      current_total_KL <- KL_accumulated + KL_k
      current_total_E_S2 <- sum_E_S2_accumulated + sum(post_sec_moment)
      current_cross_term <- cross_term_accumulated + 2 * sum(y * post_mean)

      R2 <- sum_X2 - current_cross_term + current_total_E_S2
      if (estimate_s) {
        s_current <- estimate_residual_sd(R2, d, n)
      }
      comp_residual_sd[iter] <- s_current
      global_Eloglik <- -0.5 * (d * n) * log(2 * pi * s_current^2) - R2 / (2 * s_current^2)

      global_elbo <- global_Eloglik - current_total_KL
      comp_obj[iter] <- global_elbo

      s_plus <- post_mean
      current_pi <- get_prior_prob(ebnm_fit)
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
        if (is_converged && verbose > 0) cat("Component", k, "converged in", iter, "iterations. Global ELBO:", global_elbo, "\n")

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
    residual_sd[[k]] <- comp_residual_sd[1:iter]
  }

  return(list(
    W = W_found,
    S = t(X) %*% W_found,
    S_plus = S_plus,
    prior_probs = prior_probs,
    residual_sd = s_current,
    residual_variance = s_current^2,
    residual_sd_path = residual_sd,
    residual_variance_path = residual_variance_path(residual_sd),
    obj = obj_history
  ))
}

ebica <- function(X, K = 1, alg_typ = "parallel", symmetric_rademacher = FALSE,
                  conv_crit = "weights", max_iter = 1000, tol = 1e-6, s = NULL,
                  S_init = NULL, W_init = NULL, verbose = 0) {
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
  if (length(verbose) != 1 || !is.numeric(verbose) || !is.finite(verbose)) {
    stop("verbose must be a single finite numeric value.")
  }
  if (K > d) {
    stop("K cannot be greater than the number of features (d).")
  }

  inits <- check_ebica_inits(S_init, W_init, n, d, K)
  S_init <- inits$S_init
  W_init <- inits$W_init
  estimate_s <- is.null(s)

  if (alg_typ == "parallel" || K == 1) {
    return(ebica_parallel(X, K, d, n, max_iter, tol, s, estimate_s, S_init, W_init,
                          conv_crit, symmetric_rademacher, verbose))
  } else {
    return(ebica_deflation(X, K, d, n, max_iter, tol, s, estimate_s, S_init, W_init,
                           conv_crit, symmetric_rademacher, verbose))
  }
}



# parallel update for general source prior
ebica_generalized_parallel <- function(X, K, ebnm_fn, max_iter = 1000, tol = 1e-6, s = NULL,
                                   S_init = NULL, W_init = NULL, conv_crit = 'elbo',
                                   verbose = 0) {
  d <- nrow(X)
  n <- ncol(X)
  prior_probs <- numeric(K)
  obj_history <- numeric(max_iter)
  residual_sd <- numeric(max_iter)
  estimate_s <- is.null(s)
  s_current <- init_residual_sd(X, s)
  ebnm_fns <- check_ebnm_fn_list(ebnm_fn, K)

  inits <- check_ebica_inits(S_init, W_init, n, d, K)
  S_init <- inits$S_init
  W_init <- inits$W_init

  if (!is.null(W_init)) {
    W <- W_init
  } else if (!is.null(S_init)) {
    W <- X %*% S_init
  } else {
    W <- matrix(rnorm(d * K), nrow = d, ncol = K)
  }
  W <- orthonormalize_W(W)

  for (iter in 1:max_iter) {
    W_old <- W
    Y <- t(W) %*% X

    S_plus <- matrix(0, nrow = n, ncol = K)
    KL_total <- 0
    sum_E_S2 <- 0

    for (k in 1:K) {
      y_k <- Y[k, ]

      ebnm_fit <- ebnm_fns[[k]](x = y_k, s = s_current)

      post_mean <- ebnm_fit$posterior$mean
      post_sec_moment <- ebnm_fit$posterior$sd^2 + post_mean^2

      S_plus[, k] <- post_mean
      prior_probs[k] <- get_prior_prob(ebnm_fit)

      Eloglik <- normal_means_loglik(x = y_k, s = s_current,
                                     Et = post_mean,
                                     Et2 = post_sec_moment)

      KL_k <- Eloglik - ebnm_fit$log_likelihood

      KL_total <- KL_total + KL_k
      sum_E_S2 <- sum_E_S2 + sum(post_sec_moment)
    }

    W_plus <- X %*% S_plus

    R2 <- sum(X^2) - 2 * sum(as.vector(Y) * as.vector(t(S_plus))) + sum_E_S2
    if (estimate_s) {
      s_current <- estimate_residual_sd(R2, d, n)
    }
    residual_sd[iter] <- s_current
    global_Eloglik <- -0.5 * (d * n) * log(2 * pi * s_current^2) - R2 / (2 * s_current^2)
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
      if (is_converged && verbose > 0) cat("Converged in", iter, "iterations. Final ELBO:", elbo, "\n")
      break
    }
  }

  return(list(
    W = W,
    S = t(X) %*% W,
    S_plus = S_plus,
    prior_probs = prior_probs,
    residual_sd = s_current,
    residual_variance = s_current^2,
    residual_sd_path = residual_sd[1:iter],
    residual_variance_path = residual_variance_path(residual_sd[1:iter]),
    obj = obj_history[1:iter]
  ))
}
