shifted_rademacher <- function(pi, c = 0) {
  out <- data.frame(pi = pi, c = c)
  class(out) <- c("shifted_rademacher", class(out))
  out
}

shifted_rademacher_pi <- function(g_init = NULL, default = 0.5) {
  if (!is.null(g_init) && !is.null(g_init$pi)) {
    default <- g_init$pi
  }
  if (length(default) != 1 || !is.finite(default) || default <= 0 || default >= 1) {
    stop("pi must be a single finite value between 0 and 1.")
  }
  min(max(as.numeric(default), 1e-8), 1 - 1e-8)
}

shifted_rademacher_c <- function(x, s, pi, g_init = NULL) {
  if (!is.null(g_init) && !is.null(g_init$c)) {
    c_init <- g_init$c
  } else {
    w <- 1 / s^2
    c_init <- sum(w * x) / sum(w) - (2 * pi - 1)
  }
  if (length(c_init) != 1 || !is.finite(c_init)) {
    stop("c must be a single finite value.")
  }
  as.numeric(c_init)
}

llik_shifted_rademacher <- function(x, s, pi, c) {
  ll_minus <- log(1 - pi) + stats::dnorm(x, mean = c - 1, sd = s, log = TRUE)
  ll_plus  <- log(pi)     + stats::dnorm(x, mean = c + 1, sd = s, log = TRUE)
  sum(matrixStats::rowLogSumExps(cbind(ll_minus, ll_plus)))
}

post_prob_plus_shifted_rademacher <- function(x, s, pi, c) {
  ll_minus <- log(1 - pi) + stats::dnorm(x, mean = c - 1, sd = s, log = TRUE)
  ll_plus  <- log(pi)     + stats::dnorm(x, mean = c + 1, sd = s, log = TRUE)
  log_norm <- matrixStats::rowLogSumExps(cbind(ll_minus, ll_plus))
  exp(ll_plus - log_norm)
}

opt_shifted_rademacher <- function(x, s, pi_init, c_init, fix_pi = FALSE,
                                   control = NULL) {
  optim_control <- if (is.null(control)) list() else control

  if (fix_pi) {
    opt <- stats::optim(
      par = c_init,
      fn = function(c) -llik_shifted_rademacher(x, s, pi_init, c),
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
      pi <- stats::plogis(par[1])
      c <- par[2]
      -llik_shifted_rademacher(x, s, pi, c)
    },
    method = "BFGS",
    control = optim_control
  )

  list(
    pi = stats::plogis(opt$par[1]),
    c = opt$par[2],
    loglik = -opt$value,
    convergence = opt$convergence
  )
}

post_summary_shifted_rademacher <- function(x, s, pi, c) {
  x <- as.vector(x)
  r <- post_prob_plus_shifted_rademacher(x, s, pi, c)

  post_mean <- c + 2 * r - 1
  post_second_moment <- r * (c + 1)^2 + (1 - r) * (c - 1)^2
  post_var <- post_second_moment - post_mean^2

  data.frame(
    mean = as.vector(post_mean),
    sd = sqrt(pmax(0, post_var)),
    second_moment = as.vector(post_second_moment)
  )
}

shifted_rademacher_fit <- function(x, s, pi_hat, c_hat, loglik) {
  post <- post_summary_shifted_rademacher(x, s, pi_hat, c_hat)

  structure(
    list(
      data = data.frame(x = x, s = s),
      posterior = post,
      fitted_g = shifted_rademacher(pi_hat, c_hat),
      log_likelihood = loglik
    ),
    class = c("list", "ebnm")
  )
}

ebnm_generalized_shifted_rademacher <- function(
    x,
    s = 1,
    mode = 0,
    g_init = NULL,
    fix_g = FALSE,
    output = NULL,
    control = NULL
) {
  if (mode != 0) {
    stop("ebnm_generalized_shifted_rademacher only supports mode = 0.")
  }

  x <- as.vector(x)
  if (length(s) == 1) s <- rep(s, length(x))
  if (length(s) != length(x) || any(!is.finite(s)) || any(s <= 0)) {
    stop("s must be positive and have length 1 or length(x).")
  }

  pi_init <- shifted_rademacher_pi(g_init)
  c_init <- shifted_rademacher_c(x, s, pi_init, g_init)

  if (!fix_g) {
    opt <- opt_shifted_rademacher(
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
    loglik <- llik_shifted_rademacher(x, s, pi_hat, c_hat)
  }

  shifted_rademacher_fit(x, s, pi_hat, c_hat, loglik)
}

ebnm_shifted_rademacher <- function(
    x,
    s = 1,
    mode = 0,
    g_init = NULL,
    fix_g = FALSE,
    output = NULL,
    control = NULL
) {
  if (mode != 0) {
    stop("ebnm_shifted_rademacher only supports mode = 0.")
  }

  x <- as.vector(x)
  if (length(s) == 1) s <- rep(s, length(x))
  if (length(s) != length(x) || any(!is.finite(s)) || any(s <= 0)) {
    stop("s must be positive and have length 1 or length(x).")
  }

  pi_hat <- 0.5
  c_init <- shifted_rademacher_c(x, s, pi_hat, g_init)

  if (!fix_g) {
    opt <- opt_shifted_rademacher(
      x = x,
      s = s,
      pi_init = pi_hat,
      c_init = c_init,
      fix_pi = TRUE,
      control = control
    )
    c_hat <- opt$c
    loglik <- opt$loglik
  } else {
    c_hat <- c_init
    loglik <- llik_shifted_rademacher(x, s, pi_hat, c_hat)
  }

  shifted_rademacher_fit(x, s, pi_hat, c_hat, loglik)
}

ebnm_shifted_rademacher_fixed_pi <- function(pi_fixed) {
  if (length(pi_fixed) != 1 || !is.finite(pi_fixed) ||
      pi_fixed <= 0 || pi_fixed >= 1) {
    stop("pi_fixed must be a single finite value between 0 and 1.")
  }

  function(x, s = 1, mode = 0, g_init = NULL, fix_g = FALSE,
           output = NULL, control = NULL) {
    if (mode != 0) {
      stop("ebnm_shifted_rademacher_fixed_pi only supports mode = 0.")
    }

    x <- as.vector(x)
    if (length(s) == 1) s <- rep(s, length(x))
    if (length(s) != length(x) || any(!is.finite(s)) || any(s <= 0)) {
      stop("s must be positive and have length 1 or length(x).")
    }

    c_init <- shifted_rademacher_c(x, s, pi_fixed, g_init)

    if (!fix_g) {
      opt <- opt_shifted_rademacher(
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
      loglik <- llik_shifted_rademacher(x, s, pi_fixed, c_hat)
    }

    shifted_rademacher_fit(x, s, pi_fixed, c_hat, loglik)
  }
}
