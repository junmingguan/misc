ebnm_normal_fixed_sd <- function(x,
                               s = 1,
                               mode = 0,
                               g_init = NULL,
                               fix_g = FALSE,
                               output = NULL,
                               control = NULL,
                               prior_sd = 1,
                               prior_mean = 0) {
  if (length(prior_sd) != 1 || !is.finite(prior_sd) || prior_sd < 0) {
    stop("prior_sd must be a single non-negative finite value.")
  }
  if (length(prior_mean) != 1 || !is.finite(prior_mean)) {
    stop("prior_mean must be a single finite value.")
  }

  g_fixed <- ashr::normalmix(pi = 1, mean = prior_mean, sd = prior_sd)

  ebnm::ebnm_normal(
    x = x,
    s = s,
    mode = mode,
    g_init = g_fixed,
    fix_g = TRUE,
    output = output,
    control = control
  )
}

make_ebnm_normal_fixed_sd <- function(prior_sd, prior_mean = 0) {
  force(prior_sd)
  force(prior_mean)

  function(x,
           s = 1,
           mode = 0,
           g_init = NULL,
           fix_g = FALSE,
           output = NULL,
           control = NULL) {
    ebnm_normal_fixed_sd(
      x = x,
      s = s,
      mode = mode,
      output = output,
      control = control,
      prior_sd = prior_sd,
      prior_mean = prior_mean
    )
  }
}
