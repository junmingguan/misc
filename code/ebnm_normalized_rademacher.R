normalized_rademacher <- function(pi, c = 0) {
  out <- data.frame(pi = pi, c = c)
  class(out) <- c("normalized_rademacher", class(out))
  out
}

bound_normalized_rademacher_pi <- function(pi) {
  min(max(as.numeric(pi), 1e-8), 1 - 1e-8)
}

normalized_rademacher_pi <- function(g_init = NULL, default = 0.5) {
  if (!is.null(g_init) && !is.null(g_init$pi)) {
    default <- g_init$pi
  }
  if (length(default) != 1 || !is.finite(default) ||
      default <= 0 || default >= 1) {
    stop("pi must be a single finite value between 0 and 1.")
  }
  bound_normalized_rademacher_pi(default)
}

normalized_rademacher_c <- function(x, s, g_init = NULL) {
  if (!is.null(g_init) && !is.null(g_init$c)) {
    c_init <- g_init$c
  } else {
    w <- 1 / s^2
    c_init <- sum(w * x) / sum(w)
  }
  if (length(c_init) != 1 || !is.finite(c_init)) {
    stop("c must be a single finite value.")
  }
  as.numeric(c_init)
}

normalized_rademacher_atoms <- function(pi, c) {
  pi <- bound_normalized_rademacher_pi(pi)
  delta_minus <- sqrt((1 - pi) / pi)
  delta_plus <- sqrt(pi / (1 - pi))
  c(theta_minus = c - delta_minus, theta_plus = c + delta_plus)
}

llik_normalized_rademacher <- function(x, s, pi, c) {
  pi <- bound_normalized_rademacher_pi(pi)
  theta <- normalized_rademacher_atoms(pi, c)
  ll_minus <- log(pi)     + stats::dnorm(x, mean = theta["theta_minus"], sd = s, log = TRUE)
  ll_plus  <- log(1 - pi) + stats::dnorm(x, mean = theta["theta_plus"],  sd = s, log = TRUE)
  sum(matrixStats::rowLogSumExps(cbind(ll_minus, ll_plus)))
}

post_prob_plus_normalized_rademacher <- function(x, s, pi, c) {
  pi <- bound_normalized_rademacher_pi(pi)
  theta <- normalized_rademacher_atoms(pi, c)
  ll_minus <- log(pi)     + stats::dnorm(x, mean = theta["theta_minus"], sd = s, log = TRUE)
  ll_plus  <- log(1 - pi) + stats::dnorm(x, mean = theta["theta_plus"],  sd = s, log = TRUE)
  log_norm <- matrixStats::rowLogSumExps(cbind(ll_minus, ll_plus))
  exp(ll_plus - log_norm)
}

opt_normalized_rademacher <- function(x, s, pi_init, c_init, fix_pi = FALSE,
                                     control = NULL) {
  optim_control <- if (is.null(control)) list() else control

  if (fix_pi) {
    opt <- stats::optim(
      par = c_init,
      fn = function(c) -llik_normalized_rademacher(x, s, pi_init, c),
      method = "BFGS",
      control = optim_control
    )

    return(list(
      pi = pi_init,
      c = opt$par,
      loglik = -opt$value,
      convergence = opt$convergence
    ))
  }

  opt <- stats::optim(
    par = c(stats::qlogis(pi_init), c_init),
    fn = function(par) {
      pi <- bound_normalized_rademacher_pi(stats::plogis(par[1]))
      c <- par[2]
      -llik_normalized_rademacher(x, s, pi, c)
    },
    method = "BFGS",
    control = optim_control
  )

  list(
    pi = bound_normalized_rademacher_pi(stats::plogis(opt$par[1])),
    c = opt$par[2],
    loglik = -opt$value,
    convergence = opt$convergence
  )
}

opt_normalized_rademacher_fixed_c <- function(x, s, pi_init, c_fixed,
                                             control = NULL) {
  optim_control <- if (is.null(control)) list() else control

  opt <- stats::optim(
    par = stats::qlogis(pi_init),
    fn = function(eta) {
      pi <- bound_normalized_rademacher_pi(stats::plogis(eta))
      -llik_normalized_rademacher(x, s, pi, c_fixed)
    },
    method = "BFGS",
    control = optim_control
  )

  list(
    pi = bound_normalized_rademacher_pi(stats::plogis(opt$par)),
    c = c_fixed,
    loglik = -opt$value,
    convergence = opt$convergence
  )
}

post_summary_normalized_rademacher <- function(x, s, pi, c) {
  x <- as.vector(x)
  theta <- normalized_rademacher_atoms(pi, c)
  r <- post_prob_plus_normalized_rademacher(x, s, pi, c)

  post_mean <- (1 - r) * theta["theta_minus"] + r * theta["theta_plus"]
  post_second_moment <- (1 - r) * theta["theta_minus"]^2 + r * theta["theta_plus"]^2
  post_var <- post_second_moment - post_mean^2

  data.frame(
    mean = as.vector(post_mean),
    sd = sqrt(pmax(0, post_var)),
    second_moment = as.vector(post_second_moment)
  )
}

normalized_rademacher_fit <- function(x, s, pi_hat, c_hat, loglik) {
  post <- post_summary_normalized_rademacher(x, s, pi_hat, c_hat)

  structure(
    list(
      data = data.frame(x = x, s = s),
      posterior = post,
      fitted_g = normalized_rademacher(pi_hat, c_hat),
      log_likelihood = loglik
    ),
    class = c("list", "ebnm")
  )
}

ebnm_normalized_rademacher <- function(
    x,
    s = 1,
    mode = 0,
    g_init = NULL,
    fix_g = FALSE,
    output = NULL,
    control = NULL
) {
  if (mode != 0) {
    stop("ebnm_normalized_rademacher only supports mode = 0.")
  }

  x <- as.vector(x)
  if (length(s) == 1) s <- rep(s, length(x))
  if (length(s) != length(x) || any(!is.finite(s)) || any(s <= 0)) {
    stop("s must be positive and have length 1 or length(x).")
  }

  pi_init <- normalized_rademacher_pi(g_init)
  c_init <- normalized_rademacher_c(x, s, g_init)

  if (!fix_g) {
    opt <- opt_normalized_rademacher(
      x = x,
      s = s,
      pi_init = pi_init,
      c_init = c_init,
      fix_pi = FALSE,
      control = control
    )
    pi_hat <- opt$pi
    c_hat <- opt$c
    loglik <- opt$loglik
  } else {
    pi_hat <- pi_init
    c_hat <- c_init
    loglik <- llik_normalized_rademacher(x, s, pi_hat, c_hat)
  }

  normalized_rademacher_fit(x, s, pi_hat, c_hat, loglik)
}

ebnm_normalized_rademacher_fixed_pi <- function(pi_fixed) {
  if (length(pi_fixed) != 1 || !is.finite(pi_fixed) ||
      pi_fixed <= 0 || pi_fixed >= 1) {
    stop("pi_fixed must be a single finite value between 0 and 1.")
  }

  function(x, s = 1, mode = 0, g_init = NULL, fix_g = FALSE,
           output = NULL, control = NULL) {
    if (mode != 0) {
      stop("ebnm_normalized_rademacher_fixed_pi only supports mode = 0.")
    }

    x <- as.vector(x)
    if (length(s) == 1) s <- rep(s, length(x))
    if (length(s) != length(x) || any(!is.finite(s)) || any(s <= 0)) {
      stop("s must be positive and have length 1 or length(x).")
    }

    c_init <- normalized_rademacher_c(x, s, g_init)

    if (!fix_g) {
      opt <- opt_normalized_rademacher(
        x = x,
        s = s,
        pi_init = pi_fixed,
        c_init = c_init,
        fix_pi = TRUE,
        control = control
      )
      c_hat <- opt$c
      loglik <- opt$loglik
    } else {
      c_hat <- c_init
      loglik <- llik_normalized_rademacher(x, s, pi_fixed, c_hat)
    }

    normalized_rademacher_fit(x, s, pi_fixed, c_hat, loglik)
  }
}

ebnm_normalized_rademacher_fixed_c <- function(c_fixed) {
  if (length(c_fixed) != 1 || !is.finite(c_fixed)) {
    stop("c_fixed must be a single finite value.")
  }

  function(x, s = 1, mode = 0, g_init = NULL, fix_g = FALSE,
           output = NULL, control = NULL) {
    if (mode != 0) {
      stop("ebnm_normalized_rademacher_fixed_c only supports mode = 0.")
    }

    x <- as.vector(x)
    if (length(s) == 1) s <- rep(s, length(x))
    if (length(s) != length(x) || any(!is.finite(s)) || any(s <= 0)) {
      stop("s must be positive and have length 1 or length(x).")
    }

    pi_init <- normalized_rademacher_pi(g_init)

    if (!fix_g) {
      opt <- opt_normalized_rademacher_fixed_c(
        x = x,
        s = s,
        pi_init = pi_init,
        c_fixed = c_fixed,
        control = control
      )
      pi_hat <- opt$pi
      loglik <- opt$loglik
    } else {
      pi_hat <- pi_init
      loglik <- llik_normalized_rademacher(x, s, pi_hat, c_fixed)
    }

    normalized_rademacher_fit(x, s, pi_hat, c_fixed, loglik)
  }
}
