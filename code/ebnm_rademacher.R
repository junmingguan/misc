rademacher <- function(w) {
  out <- data.frame(w = w, a = 1)
  class(out) <- c("rademacher", class(out))
  out
}

llik_rademacher <- function(x, s, w) {
  ll_plus  <- log(w)     + dnorm(x, mean = +1, sd = s, log = TRUE)
  ll_minus <- log(1 - w) + dnorm(x, mean = -1, sd = s, log = TRUE)
  sum(matrixStats::rowLogSumExps(cbind(ll_plus, ll_minus)))
}

opt_rademacher <- function(x, s, w_init) {
  stats::optim(
    par    = w_init,
    fn     = function(w) -llik_rademacher(x, s, w),
    method = "L-BFGS-B",
    lower  = 1e-8,
    upper  = 1 - 1e-8
  )
}


post_summary_rademacher <- function(x, s, w) {
  x <- as.vector(x)

  log_odds <- log(w) - log(1 - w)

  post_mean <- tanh( (x / s^2) + 0.5 * log_odds )

  post_var <- 1 - post_mean^2
  post_sd  <- sqrt(pmax(0, post_var))

  data.frame(
    mean = as.vector(post_mean),
    sd   = as.vector(post_sd),
    second_moment = rep(1, length(x))
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
  if (mode != 0) {
    stop("ebnm_rademacher only supports mode = 0.")
  }

  if (length(s) == 1) s <- rep(s, length(x))

  if (!is.null(g_init)) {
    w_init <- g_init$w
  } else {
    w_init <- 0.5
  }

  if (!fix_g) {
    opt    <- opt_rademacher(x, s, w_init)
    w_hat  <- opt$par
    loglik <- -opt$value
  } else {
    w_hat  <- w_init
    loglik <- llik_rademacher(x, s, w_hat)
  }

  post <- post_summary_rademacher(x, s, w_hat)

  structure(
    list(
      data           = data.frame(x = x, s = s),
      posterior      = post,
      fitted_g       = rademacher(w_hat),
      log_likelihood = loglik
    ),
    class = c("list", "ebnm")
  )
}


ebnm_symm_rademacher <- function(x, s, g_init, fix_g, output) {
  ebnm_rademacher(
    x = x,
    s = s,
    g_init = rademacher(w = 0.5),
    fix_g = TRUE
  )
}
