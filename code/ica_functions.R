fastica_rank_one_update_1iter = function(X, w, include_G2=TRUE, pi = NULL){
  w= w/sqrt(sum(w^2))
  P = t(X) %*% w
  if (!is.null(pi)) {
    if (length(pi) != 1 || !is.finite(pi) || pi <= 0 || pi >= 1) {
      stop("pi must be NULL or a single value between 0 and 1.")
    }
    P = P + 0.5 * log(pi / (1 - pi))
  }
  G = tanh(P)
  G2 = 0
  if(include_G2)
    {G2 = 1-tanh(P)^2}
  w = X %*% G - sum(G2) * w
  w = w/sqrt(sum(w^2))
  return(w)
}

fastica_rank_one_update = function(X1, init_w, include_G2 = TRUE, n_iter = 100, pi = NULL) {
  w = init_w
  for(i in 1:n_iter)
  {
    w = fastica_rank_one_update_1iter(X1, w, include_G2, pi)
  }
  return(w)
}

fastica_vectorized_rank_one_update = function(X1, init_W = NULL, seeds = 1:100,
                                              include_G2 = TRUE, n_iter = 100,
                                              pi = NULL) {

  if (is.null(init_W)) {
    init_W = sapply(seeds, function(seed) {
      set.seed(seed)
      rnorm(nrow(X1))
    })
  }

  W = sweep(init_W, 2, sqrt(colSums(init_W^2)), "/")

  for(i in 1:n_iter) {
    P = t(X1) %*% W
    if (!is.null(pi)) {
      if (length(pi) != 1 || !is.finite(pi) || pi <= 0 || pi >= 1) {
        stop("pi must be NULL or a single value between 0 and 1.")
      }
      P = P + 0.5 * log(pi / (1 - pi))
    }
    G = tanh(P)
    if (include_G2) {
      G2 = 1 - tanh(P)^2
      W = X1 %*% G - sweep(W, 2, colSums(G2), "*")
    } else {
      W = X1 %*% G
    }
    W = sweep(W, 2, sqrt(colSums(W^2)), "/")
  }
  W
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

  ## This appears to be equivalant to X1 = sqrt(n) * t(svd(X)$v[,1:n.comp])
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


compute_whitened_residual_sd <- function(X, n.comp = 10, sigma, center = TRUE) {
  n <- nrow(X)
  p <- ncol(X)
  sigma2 <- sigma^2

  if (center) {
    X <- scale(X, scale = FALSE)
  }

  Xt <- t(X)

  V <- Xt %*% t(Xt) / n

  s <- La.svd(V)

  # Finite-sample correction due to centering
  sigma2_eff <- if (center) sigma2 * (1 - 1 / n) else sigma2

  se <- sqrt(sigma2_eff / s$d[1:n.comp])

  S <- matrix(
    rep(se, times = n),
    nrow = n.comp,
    ncol = n
  )

  return(S)
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

get_prior <- function(ebnm_fit) {
  return(ebnm_fit$fitted_g)
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


    if (K > 1) {
      svd_res <- svd(W_plus)
      W <- svd_res$u %*% t(svd_res$v)
    } else {
      W <- matrix(W_plus / sqrt(sum(W_plus^2)), ncol = 1)
    }

    is_converged <- FALSE
    if (conv_crit == "elbo") {
      if (iter > 1 && (elbo - obj_history[iter - 1]) < tol) is_converged <- TRUE
    } else if (conv_crit == "weights") {
      if (K > 1) {
        if (max(1 - abs(colSums(W * W_old))) < tol) is_converged <- TRUE
      } else {
        if (1 - abs(sum(W * W_old)) < tol) is_converged <- TRUE
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
      } else if (conv_crit == "weights") {
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
                  conv_crit = "maxiter", max_iter = 1000, tol = 1e-6, s = NULL,
                  S_init = NULL, W_init = NULL, verbose = 0) {
  d <- nrow(X)
  n <- ncol(X)

  if (!alg_typ %in% c("parallel", "deflation")) {
    stop("alg_typ must be either 'parallel' or 'deflation'.")
  }
  if (!conv_crit %in% c("maxiter", "weights", "elbo")) {
    stop("conv_crit must be one of 'maxiter', 'weights', or 'elbo'.")
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

ebica_rank_one_unconstrained <- function(X, symmetric_rademacher = FALSE,
                                         conv_crit = "maxiter", max_iter = 1000,
                                         tol = 1e-6, s = NULL, S_init = NULL,
                                         W_init = NULL, verbose = 0) {
  d <- nrow(X)
  n <- ncol(X)

  if (!conv_crit %in% c("maxiter", "weights", "elbo")) {
    stop("conv_crit must be one of 'maxiter', 'weights', or 'elbo'.")
  }
  if (!is.logical(symmetric_rademacher)) {
    stop("symmetric_rademacher must be TRUE or FALSE.")
  }
  if (length(verbose) != 1 || !is.numeric(verbose) || !is.finite(verbose)) {
    stop("verbose must be a single finite numeric value.")
  }
  if (length(max_iter) != 1 || !is.finite(max_iter) || max_iter <= 0 || max_iter != as.integer(max_iter)) {
    stop("max_iter must be a single positive integer.")
  }

  inits <- check_ebica_inits(S_init, W_init, n, d, 1)
  S_init <- inits$S_init
  W_init <- inits$W_init
  estimate_s <- is.null(s)
  s_current <- init_residual_sd(X, s)

  if (!is.null(W_init)) {
    w <- as.numeric(W_init[, 1])
  } else if (!is.null(S_init)) {
    S_init2 <- sum(S_init[, 1]^2)
    if (S_init2 == 0) {
      stop("S_init must have nonzero norm.")
    }
    w <- as.numeric(X %*% S_init[, 1]) / S_init2
  } else {
    w <- rnorm(d)
    # w <- w / sqrt(sum(w^2))
  }
  if (sum(w^2) == 0) {
    stop("Initial w must have nonzero norm.")
  }

  obj_history <- numeric(max_iter)
  residual_sd <- numeric(max_iter)

  for (iter in 1:max_iter) {
    w_old <- w
    w_norm2 <- sum(w^2)
    y <- as.numeric(t(w) %*% X) / w_norm2
    y_se <- s_current / sqrt(w_norm2)

    if (symmetric_rademacher) {
      ebnm_fit <- ebnm_rademacher(x = y, s = y_se)
    } else {
      ebnm_fit <- ebnm_generalized_rademacher(x = y, s = y_se)
    }

    post_mean <- ebnm_fit$posterior$mean
    post_sec_moment <- ebnm_fit$posterior$sd^2 + post_mean^2
    sum_E_S2 <- sum(post_sec_moment)
    if (sum_E_S2 <= 0) {
      stop("The posterior second moment is zero; cannot estimate w by MLE.")
    }

    Eloglik <- normal_means_loglik(x = y, s = y_se,
                                   Et = post_mean,
                                   Et2 = post_sec_moment)
    KL <- Eloglik - ebnm_fit$log_likelihood

    w <- as.numeric(X %*% post_mean) / sum_E_S2
    fitted <- as.numeric(t(w) %*% X)
    R2 <- sum(X^2) - 2 * sum(fitted * post_mean) + sum(w^2) * sum_E_S2
    if (estimate_s) {
      s_current <- estimate_residual_sd(R2, d, n)
    }
    residual_sd[iter] <- s_current
    global_Eloglik <- -0.5 * (d * n) * log(2 * pi * s_current^2) - R2 / (2 * s_current^2)
    elbo <- global_Eloglik - KL
    obj_history[iter] <- elbo

    is_converged <- FALSE
    if (conv_crit == "elbo") {
      if (iter > 1 && (elbo - obj_history[iter - 1]) < tol) is_converged <- TRUE
    } else if (conv_crit == "weights") {
      rel_w_change <- sqrt(sum((w - w_old)^2)) / max(sqrt(sum(w_old^2)), sqrt(.Machine$double.eps))
      if (rel_w_change < tol) is_converged <- TRUE
    }

    if (is_converged || iter == max_iter) {
      if (is_converged && verbose > 0) cat("Converged in", iter, "iterations. Final ELBO:", elbo, "\n")
      break
    }
  }

  w_norm2 <- sum(w^2)
  return(list(
    W = matrix(w, ncol = 1),
    S = as.numeric(t(w) %*% X) / w_norm2,
    S_plus = matrix(post_mean, ncol = 1),
    prior_probs = get_prior_prob(ebnm_fit),
    residual_sd = s_current,
    residual_variance = s_current^2,
    residual_sd_path = residual_sd[1:iter],
    residual_variance_path = residual_variance_path(residual_sd[1:iter]),
    obj = obj_history[1:iter]
  ))
}



# parallel update for general source prior
ebica_generalized_parallel <- function(X, K, ebnm_fn, max_iter = 1000, tol = 1e-6, s = NULL,
                                   S_init = NULL, W_init = NULL, conv_crit = 'maxiter',
                                   verbose = 0) {
  d <- nrow(X)
  n <- ncol(X)
  # prior_probs <- numeric(K)
  fitted_priors <- vector("list", K)
  obj_history <- numeric(max_iter)
  residual_sd <- numeric(max_iter)
  estimate_s <- is.null(s)
  s_current <- init_residual_sd(X, s)
  ebnm_fns <- check_ebnm_fn_list(ebnm_fn, K)

  if (!conv_crit %in% c("maxiter", "weights", "elbo")) {
    stop("conv_crit must be one of 'maxiter', 'weights', or 'elbo'.")
  }

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
      # prior_probs[k] <- get_prior_prob(ebnm_fit)
      fitted_priors[[k]] <- get_prior(ebnm_fit)

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


    if (K > 1) {
      svd_res <- svd(W_plus)
      W <- svd_res$u %*% t(svd_res$v)
    } else {
      W <- matrix(W_plus / sqrt(sum(W_plus^2)), ncol = 1)
    }

    is_converged <- FALSE
    if (conv_crit == "elbo") {
      if (iter > 1 && (elbo - obj_history[iter - 1]) < tol) is_converged <- TRUE
    } else if (conv_crit == "weights") {
      if (K > 1) {
        if (max(1 - abs(colSums(W * W_old))) < tol) is_converged <- TRUE
      } else {
        if (1 - abs(sum(W * W_old)) < tol) is_converged <- TRUE
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
    # prior_probs = prior_probs,
    fitted_priors = fitted_priors,
    residual_sd = s_current,
    residual_variance = s_current^2,
    residual_sd_path = residual_sd[1:iter],
    residual_variance_path = residual_variance_path(residual_sd[1:iter]),
    obj = obj_history[1:iter]
  ))
}

## ----rank-one-ica-solvers----
# Independent rank-1 ICA fits. X is whitened PCs x samples; each column of
# starts is a separate unit-sphere direction. No intercept or deflation.
rank1_ica_contrast <- function(y, contrast, alpha = 1) {
  contrast <- match.arg(contrast, c("logcosh", "quartic", "cubic", "abs", "logphi"))
  answer <- switch(contrast,
    logcosh = list(G = abs(y) + log1p(exp(-2 * abs(y))) - log(2),
                   g = tanh(y), gprime = 1 - tanh(y)^2),
    quartic = list(G = y^4, g = 4 * y^3, gprime = 12 * y^2),
    cubic = list(G = y^3 / 3, g = y^2, gprime = 2 * y),
    abs = list(G = y * abs(y) / 2, g = abs(y), gprime = sign(y)),
    logphi = {
      u <- alpha * y
      logcdf <- pnorm(u, log.p = TRUE)
      h <- exp(dnorm(u, log = TRUE) - logcdf)
      list(G = logcdf, g = alpha * h, gprime = -alpha^2 * h * (u + h))
    })
  lapply(answer, function(z) { dim(z) <- dim(y); z })
}

.rank1_prepare <- function(X, starts, maxiter, tol) {
  starts <- as.matrix(starts)
  stopifnot(is.matrix(X), is.numeric(X), all(is.finite(X)),
    nrow(starts) == nrow(X), ncol(starts) > 0L,
    all(is.finite(starts)), all(colSums(starts^2) > 0),
    length(maxiter) == 1L, maxiter >= 1L, maxiter == as.integer(maxiter))
  if (!missing(tol)) stopifnot(length(tol) == 1L, is.finite(tol), tol > 0)
  sweep(starts, 2, sqrt(colSums(starts^2)), "/")
}

.rank1_gap <- function(a, b, symmetric = FALSE) {
  gap <- sqrt(colSums((a - b)^2))
  if (symmetric) gap <- pmin(gap, sqrt(colSums((a + b)^2)))
  gap
}

.rank1_fastica_step <- function(X, W, contrast, alpha = 1) {
  c <- rank1_ica_contrast(t(X) %*% W, contrast, alpha)
  # Newton correction: E[X g(X'w)] - E[g'(X'w)] w.
  raw <- X %*% c$g / ncol(X) - sweep(W, 2, colMeans(c$gprime), "*")
  norms <- sqrt(colSums(raw^2))
  good <- is.finite(norms) & norms > 1e-14
  candidate <- W
  candidate[, good] <- sweep(raw[, good, drop = FALSE], 2, norms[good], "/")
  list(W = candidate, good = good)
}

# Run a fixed number of ordinary FastICA updates for each starting direction.
rank1_fastica_updates <- function(X, starts, contrast, n_iter = 50L, alpha = 1) {
  W <- .rank1_prepare(X, starts, n_iter)
  for (i in seq_len(n_iter)) {
    W <- .rank1_fastica_step(X, W, contrast, alpha)$W
  }
  W
}

fit_rank1_fastica <- function(X, starts, contrast, maxiter = 5000L,
                              tol = 1e-7, progress = function(...) NULL, alpha = 1) {
  W <- .rank1_prepare(X, starts, maxiter, tol)
  symmetric <- contrast %in% c("logcosh", "quartic")
  nr <- ncol(W)
  active <- seq_len(nr)
  lagged <- W
  iterations <- integer(nr)
  status <- rep("iteration_limit", nr)
  trace <- vector("list", maxiter)
  for (iter in seq_len(maxiter)) {
    old <- W[, active, drop = FALSE]
    step <- .rank1_fastica_step(X, old, contrast, alpha)
    gap <- .rank1_gap(step$W, old, symmetric)
    two_gap <- .rank1_gap(step$W, lagged[, active, drop = FALSE], symmetric)
    lagged[, active] <- old
    W[, active] <- step$W
    iterations[active] <- iter
    done <- step$good & gap < tol
    cycle <- step$good & !done & iter > 1L & two_gap < tol
    # A small backward gap can occur before an unstable fixed point/cycle.
    # Confirm the residual from the returned candidate before stopping it.
    check <- which(done | cycle)
    if (length(check)) {
      candidate <- step$W[, check, drop = FALSE]
      forward1 <- .rank1_fastica_step(X, candidate, contrast, alpha)
      forward2 <- .rank1_fastica_step(X, forward1$W, contrast, alpha)
      valid <- forward1$good & forward2$good
      one <- .rank1_gap(forward1$W, candidate, symmetric)
      two <- .rank1_gap(forward2$W, candidate, symmetric)
      step$good[check] <- valid
      done[check] <- valid & one < tol
      cycle[check] <- valid & !done[check] & two < tol
    }
    status[active[done]] <- "converged"
    status[active[cycle]] <- "two_cycle"
    status[active[!step$good]] <- "numerical_failure"
    active <- active[step$good & !done & !cycle]
    trace[[iter]] <- data.frame(iteration = iter, active = length(active),
      max_step = max(gap), median_step = median(gap))
    if (iter == 1L || iter %% 100L == 0L || !length(active))
      progress(sprintf("FastICA %s: iteration %d, %d active", contrast, iter, length(active)))
    if (!length(active)) break
  }
  next1 <- .rank1_fastica_step(X, W, contrast, alpha)
  next2 <- .rank1_fastica_step(X, next1$W, contrast, alpha)
  residual <- .rank1_gap(next1$W, W, symmetric)
  two_residual <- .rank1_gap(next2$W, W, symmetric)
  good <- next1$good & next2$good & status != "numerical_failure"
  status[good] <- ifelse(iterations[good] >= maxiter,
    "iteration_limit", "tolerance_not_met")
  status[good & residual < tol] <- "converged"
  status[good & residual >= tol & two_residual < tol] <- "two_cycle"
  status[!good] <- "numerical_failure"
  c <- rank1_ica_contrast(t(X) %*% W, contrast, alpha)
  gradient <- X %*% c$g / ncol(X)
  gradient <- gradient - sweep(W, 2, colSums(W * gradient), "*")
  diagnostics <- data.frame(start = seq_len(nr), iterations = iterations,
    status = status, objective = colMeans(c$G),
    gradient_norm = sqrt(colSums(gradient^2)),
    one_step_residual = residual, two_step_residual = two_residual)
  list(W = W, diagnostics = diagnostics, trace = do.call(rbind, trace[seq_len(iter)]))
}

fit_rank1_pg <- function(X, starts, contrast, maxiter = 5000L,
                         tol = 1e-7, progress = function(...) NULL, alpha = 1) {
  W <- .rank1_prepare(X, starts, maxiter, tol)
  nr <- ncol(W)
  Xt <- t(X)
  objective <- function(w) colMeans(rank1_ica_contrast(Xt %*% w, contrast, alpha)$G)
  values <- initial <- objective(W)
  eta <- rep(1, nr)
  iterations <- backtracks <- integer(nr)
  status <- rep("iteration_limit", nr)
  active <- seq_len(nr)
  trace <- vector("list", maxiter + 1L)
  trace[[1]] <- data.frame(iteration = 0L, active = nr,
    objective_mean = mean(values), min_objective_gain = 0, max_gradient = NA_real_)
  for (iter in seq_len(maxiter)) {
    old <- W[, active, drop = FALSE]
    c <- rank1_ica_contrast(Xt %*% old, contrast, alpha)
    gradient <- X %*% c$g / ncol(X)
    tangent <- gradient - sweep(old, 2, colSums(old * gradient), "*")
    norm2 <- colSums(tangent^2)
    done <- is.finite(norm2) & norm2 < tol^2
    bad <- !is.finite(norm2)
    status[active[done]] <- "converged"
    status[active[bad]] <- "numerical_failure"
    pending <- which(!done & !bad)
    gain <- numeric()
    accepted_ids <- integer()
    for (bt in seq_len(60L)) {
      if (!length(pending)) break
      ids <- active[pending]
      candidate <- old[, pending, drop = FALSE] +
        sweep(tangent[, pending, drop = FALSE], 2, eta[ids], "*")
      candidate <- sweep(candidate, 2, sqrt(colSums(candidate^2)), "/")
      proposed <- objective(candidate)
      # Require monotone Armijo ascent after normalization to the sphere.
      accept <- is.finite(proposed) & proposed >= values[ids] +
        1e-4 * eta[ids] * norm2[pending]
      moved <- .rank1_gap(candidate, old[, pending, drop = FALSE]) > 0
      stalled <- accept & !moved
      status[ids[stalled]] <- "line_search_stalled"
      accept <- accept & moved
      if (any(accept)) {
        accepted <- ids[accept]
        gain <- c(gain, proposed[accept] - values[accepted])
        W[, accepted] <- candidate[, accept, drop = FALSE]
        values[accepted] <- proposed[accept]
        iterations[accepted] <- iterations[accepted] + 1L
        eta[accepted] <- pmin(2 * eta[accepted], 1024)
        accepted_ids <- c(accepted_ids, accepted)
      }
      reject <- !accept & !stalled
      eta[ids[reject]] <- eta[ids[reject]] / 2
      backtracks[ids[reject]] <- backtracks[ids[reject]] + 1L
      pending <- pending[reject]
    }
    if (length(pending)) status[active[pending]] <- "line_search_stalled"
    active <- sort(accepted_ids)
    trace[[iter + 1L]] <- data.frame(iteration = iter, active = length(active),
      objective_mean = mean(values), min_objective_gain = if (length(gain)) min(gain) else 0,
      max_gradient = max(sqrt(norm2)))
    if (iter == 1L || iter %% 100L == 0L || !length(active))
      progress(sprintf("Projected gradient %s: iteration %d, %d active", contrast, iter, length(active)))
    if (!length(active)) break
  }
  c <- rank1_ica_contrast(Xt %*% W, contrast, alpha)
  gradient <- X %*% c$g / ncol(X)
  gradient <- gradient - sweep(W, 2, colSums(W * gradient), "*")
  gradient_norm <- sqrt(colSums(gradient^2))
  status[is.finite(gradient_norm) & gradient_norm < tol] <- "converged"
  diagnostics <- data.frame(start = seq_len(nr), iterations = iterations,
    status = status, initial_objective = initial, objective = colMeans(c$G),
    gradient_norm = gradient_norm, backtracks = backtracks)
  list(W = W, diagnostics = diagnostics,
    trace = do.call(rbind, trace[seq_len(iter + 1L)]))
}

# FastICA with mean curvature on the tangent space of the unit sphere.
tangent_ica_step <- function(U, W, nonlin, obs_norm2 = colSums(U^2)) {
  P <- crossprod(U, W)
  out <- nonlin(P)
  curvature <- (drop(crossprod(obs_norm2, out$gprime)) -
    colSums(out$gprime * P^2)) / (nrow(U) - 1L)
  V <- U %*% out$g - sweep(W, 2, curvature, "*")
  V <- sweep(V, 2, ifelse(colSums(W * V) < 0, -1, 1), "*")
  norms <- sqrt(colSums(V^2))
  good <- is.finite(norms) & norms >= 1e-12
  candidate <- W
  candidate[, good] <- sweep(V[, good, drop = FALSE], 2, norms[good], "/")
  list(W = candidate, good = good, finite = is.finite(norms))
}

fit_rank1_tangent <- function(X, starts, contrast, alpha = 1,
                              maxiter = 5000L, tol = 1e-7) {
  stopifnot(nrow(X) >= 2L)
  W <- .rank1_prepare(X, starts, maxiter, tol)
  nonlin <- function(P) rank1_ica_contrast(P, contrast, alpha)
  obs_norm2 <- colSums(X^2)
  step <- function(W) tangent_ica_step(X, W, nonlin, obs_norm2)
  active <- seq_len(ncol(W))
  lagged <- W
  iterations <- integer(ncol(W))
  trace <- vector("list", maxiter)
  for (iter in seq_len(maxiter)) {
    old <- W[, active, drop = FALSE]
    next_step <- step(old)
    gap <- .rank1_gap(next_step$W, old)
    two_gap <- .rank1_gap(next_step$W, lagged[, active, drop = FALSE])
    lagged[, active] <- old
    W[, active] <- next_step$W
    iterations[active] <- iter
    stop_now <- !next_step$good
    # Give shrinking alternating steps time to reach the one-step tolerance.
    check <- which(next_step$good & (gap < tol |
      (iter > 1L & gap > 1e-5 & two_gap < tol)))
    if (length(check)) {
      candidate <- next_step$W[, check, drop = FALSE]
      forward <- step(candidate)
      next_forward <- step(forward$W)
      one <- .rank1_gap(forward$W, candidate)
      two <- .rank1_gap(next_forward$W, candidate)
      stop_now[check] <- !forward$good | one < tol |
        (next_forward$good & one > 1e-5 & two < tol)
    }
    active <- active[!stop_now]
    trace[[iter]] <- data.frame(iteration = iter, active = length(active),
      max_step = max(gap), median_step = median(gap))
    if (!length(active)) break
  }
  forward <- step(W)
  next_forward <- step(forward$W)
  residual <- .rank1_gap(forward$W, W)
  two_residual <- .rank1_gap(next_forward$W, W)
  status <- ifelse(iterations == maxiter, "iteration_limit", "tolerance_not_met")
  status[forward$good & residual < tol] <- "converged"
  status[forward$good & next_forward$good & residual >= tol &
    two_residual < tol] <- "two_cycle"
  # A safeguarded return of the old direction does not certify convergence.
  status[!forward$good] <- ifelse(forward$finite[!forward$good],
    "zero_update", "numerical_failure")
  L <- crossprod(X, W)
  out <- nonlin(L)
  gradient <- X %*% out$g / ncol(X)
  tangent_gradient <- gradient - sweep(W, 2, colSums(W * gradient), "*")
  diagnostics <- data.frame(start = seq_len(ncol(W)), iterations = iterations,
    status = status, objective = colMeans(out$G),
    gradient_norm = sqrt(colSums(tangent_gradient^2)),
    one_step_residual = residual, two_step_residual = two_residual)
  list(W = W, L = L, diagnostics = diagnostics,
    trace = do.call(rbind, trace[seq_len(iter)]))
}
