point_rademacher <- function(pi_minus, pi0, pi_plus, m = 1) {
  out <- data.frame(
    pi_minus = pi_minus,
    pi0 = pi0,
    pi_plus = pi_plus,
    m = m
  )
  class(out) <- c("point_rademacher", class(out))
  out
}

point_rademacher_probs <- function(g_init = NULL) {
  if (is.null(g_init)) {
    return(c(pi_minus = 1 / 3, pi0 = 1 / 3, pi_plus = 1 / 3))
  }

  probs <- c(
    pi_minus = g_init$pi_minus,
    pi0 = g_init$pi0,
    pi_plus = g_init$pi_plus
  )
  if (any(!is.finite(probs)) || any(probs <= 0) || abs(sum(probs) - 1) > 1e-6) {
    stop("g_init must contain positive pi_minus, pi0, and pi_plus values summing to 1.")
  }
  probs / sum(probs)
}

point_rademacher_m <- function(m, g_init = NULL) {
  if (missing(m) && !is.null(g_init) && !is.null(g_init$m)) {
    m <- g_init$m
  }
  if (length(m) != 1 || !is.finite(m) || m < 0) {
    stop("m must be a single non-negative finite value.")
  }
  m
}

softmax_point_rademacher <- function(eta) {
  logits <- c(eta[1], 0, eta[2])
  exp_logits <- exp(logits - max(logits))
  probs <- exp_logits / sum(exp_logits)
  names(probs) <- c("pi_minus", "pi0", "pi_plus")
  probs
}

constrained_point_rademacher_probs <- function(pi0) {
  pi0 <- as.numeric(pi0)
  c(pi_minus = (1 - pi0) / 2, pi0 = pi0, pi_plus = (1 - pi0) / 2)
}

llik_point_rademacher <- function(x, s, probs, m = 1) {
  ll_minus <- log(probs["pi_minus"]) + dnorm(x, mean = -m, sd = s, log = TRUE)
  ll_zero  <- log(probs["pi0"])      + dnorm(x, mean = 0,  sd = s, log = TRUE)
  ll_plus  <- log(probs["pi_plus"])  + dnorm(x, mean = +m, sd = s, log = TRUE)
  sum(matrixStats::rowLogSumExps(cbind(ll_minus, ll_zero, ll_plus)))
}

opt_point_generalized_rademacher <- function(x, s, probs_init, m = 1) {
  eta_init <- log(c(probs_init["pi_minus"], probs_init["pi_plus"]) / probs_init["pi0"])
  stats::optim(
    par = eta_init,
    fn = function(eta) {
      probs <- softmax_point_rademacher(eta)
      -llik_point_rademacher(x, s, probs, m)
    },
    method = "BFGS"
  )
}

opt_point_rademacher <- function(x, s, pi0_init, m = 1) {
  stats::optim(
    par = pi0_init,
    fn = function(pi0) {
      probs <- constrained_point_rademacher_probs(pi0)
      -llik_point_rademacher(x, s, probs, m)
    },
    method = "L-BFGS-B",
    lower = 1e-8,
    upper = 1 - 1e-8
  )
}

post_summary_point_rademacher <- function(x, s, probs, m = 1) {
  x <- as.vector(x)

  ll_minus <- log(probs["pi_minus"]) + dnorm(x, mean = -m, sd = s, log = TRUE)
  ll_zero  <- log(probs["pi0"])      + dnorm(x, mean = 0,  sd = s, log = TRUE)
  ll_plus  <- log(probs["pi_plus"])  + dnorm(x, mean = +m, sd = s, log = TRUE)
  log_norm <- matrixStats::rowLogSumExps(cbind(ll_minus, ll_zero, ll_plus))

  post_minus <- exp(ll_minus - log_norm)
  post_plus <- exp(ll_plus - log_norm)

  post_mean <- m * (post_plus - post_minus)
  post_second_moment <- m^2 * (post_minus + post_plus)
  post_var <- post_second_moment - post_mean^2

  data.frame(
    mean = as.vector(post_mean),
    sd = sqrt(pmax(0, post_var)),
    second_moment = as.vector(post_second_moment)
  )
}

point_rademacher_fit <- function(x, s, probs_hat, m, loglik) {
  post <- post_summary_point_rademacher(x, s, probs_hat, m)

  structure(
    list(
      data = data.frame(x = x, s = s, m = m),
      posterior = post,
      fitted_g = point_rademacher(
        pi_minus = probs_hat["pi_minus"],
        pi0 = probs_hat["pi0"],
        pi_plus = probs_hat["pi_plus"],
        m = m
      ),
      log_likelihood = loglik
    ),
    class = c("list", "ebnm")
  )
}

ebnm_point_generalized_rademacher <- function(
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
    stop("ebnm_point_generalized_rademacher only supports mode = 0.")
  }

  m <- point_rademacher_m(m, g_init)
  if (length(s) == 1) s <- rep(s, length(x))

  probs_init <- point_rademacher_probs(g_init)

  if (!fix_g) {
    opt <- opt_point_generalized_rademacher(x, s, probs_init, m)
    probs_hat <- softmax_point_rademacher(opt$par)
    loglik <- -opt$value
  } else {
    probs_hat <- probs_init
    loglik <- llik_point_rademacher(x, s, probs_hat, m)
  }

  point_rademacher_fit(x, s, probs_hat, m, loglik)
}

ebnm_point_rademacher <- function(
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
    stop("ebnm_point_rademacher only supports mode = 0.")
  }

  m <- point_rademacher_m(m, g_init)
  if (length(s) == 1) s <- rep(s, length(x))

  probs_init <- point_rademacher_probs(g_init)
  pi0_init <- probs_init["pi0"]

  if (!fix_g) {
    opt <- opt_point_rademacher(x, s, pi0_init, m)
    pi0_hat <- opt$par
    probs_hat <- constrained_point_rademacher_probs(pi0_hat)
    loglik <- -opt$value
  } else {
    pi0_hat <- probs_init["pi0"]
    probs_hat <- constrained_point_rademacher_probs(pi0_hat)
    loglik <- llik_point_rademacher(x, s, probs_hat, m)
  }

  point_rademacher_fit(x, s, probs_hat, m, loglik)
}
