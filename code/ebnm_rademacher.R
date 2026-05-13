rademacher <- function(pi, m = 1) {
  out <- data.frame(pi = pi, m = m)
  class(out) <- c("rademacher", class(out))
  out
}

llik_rademacher <- function(x, s, pi, m = 1) {
  ll_plus  <- log(pi)     + dnorm(x, mean = +m, sd = s, log = TRUE)
  ll_minus <- log(1 - pi) + dnorm(x, mean = -m, sd = s, log = TRUE)
  sum(matrixStats::rowLogSumExps(cbind(ll_plus, ll_minus)))
}

opt_rademacher <- function(x, s, pi_init, m = 1) {
  stats::optim(
    par    = pi_init,
    fn     = function(pi) -llik_rademacher(x, s, pi, m),
    method = "L-BFGS-B",
    lower  = 1e-8,
    upper  = 1 - 1e-8
  )
}


post_summary_rademacher <- function(x, s, pi, m = 1) {
  x <- as.vector(x)

  log_odds <- log(pi) - log(1 - pi)

  post_mean <- m * tanh( (m * x / s^2) + 0.5 * log_odds )

  post_var <- m^2 - post_mean^2
  post_sd  <- sqrt(pmax(0, post_var))

  data.frame(
    mean = as.vector(post_mean),
    sd   = as.vector(post_sd),
    second_moment = rep(m^2, length(x))
  )
}

ebnm_generalized_rademacher <- function(
    x,
    s = 1,
    mode = 0,
    g_init = NULL,
    fix_g = FALSE,
    output = NULL,
    control = NULL,
    m = 1
) {
  if (mode != 0) {
    stop("ebnm_generalized_rademacher only supports mode = 0.")
  }

  if (missing(m) && !is.null(g_init) && !is.null(g_init$m)) {
    m <- g_init$m
  }
  if (length(m) != 1 || !is.finite(m) || m < 0) {
    stop("m must be a single non-negative finite value.")
  }

  if (length(s) == 1) s <- rep(s, length(x))

  if (!is.null(g_init)) {
    pi_init <- g_init$pi
  } else {
    pi_init <- 0.5
  }

  if (!fix_g) {
    opt     <- opt_rademacher(x, s, pi_init, m)
    pi_hat  <- opt$par
    loglik  <- -opt$value
  } else {
    pi_hat  <- pi_init
    loglik  <- llik_rademacher(x, s, pi_hat, m)
  }

  post <- post_summary_rademacher(x, s, pi_hat, m)

  structure(
    list(
      data           = data.frame(x = x, s = s, m = m),
      posterior      = post,
      fitted_g       = rademacher(pi_hat, m),
      log_likelihood = loglik
    ),
    class = c("list", "ebnm")
  )
}


ebnm_rademacher <- function(
    x,
    s = 1,
    mode = 0,
    g_init = NULL,
    fix_g = FALSE,
    output = NULL,
    control = NULL
) {
  ebnm_generalized_rademacher(
    x = x,
    s = s,
    mode = mode,
    m = 1,
    g_init = rademacher(pi = 0.5, m = 1),
    fix_g = TRUE
  )
}
