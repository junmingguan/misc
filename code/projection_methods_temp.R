# Projection methods copied from analysis/pancreas_endocrinogenesis.Rmd.
# Includes rank-r extensions for Log-Phi, cubic, skew, and integrated log-cosh ICA.
# ICA expects whitened data (PCs x cells); EBproj expects an orthonormal
# matrix of left singular vectors (cells x PCs).

fastica_update = function(U, W) {
  P  <- t(U) %*% W          # n x n_starts: source estimates
  G  <- tanh(P)
  G2 <- 1 - G^2
  W  <- U %*% G - sweep(W, 2, colSums(G2), "*")
  # Add epsilon (1e-15) to prevent 0/0
  sweep(W, 2, sqrt(colSums(W^2)) + 1e-15, "/")
}

gradient_minica_update = function(U, W, lr = 0.1) {
  # n_samples = ncol(U); n_features = nrow(U); n_starts = ncol(W)

  P <- t(U) %*% W          # n_samples x n_starts: source estimates
  G <- tanh(P)             # First derivative of log-cosh

  # Calculate the gradient
  # We divide by ncol(U) to average over samples, keeping the learning rate
  # stable regardless of your dataset size.
  grad <- (U %*% G) / ncol(U)

  # Gradient descent step (minimize log-cosh)
  # Note: To MAXIMIZE log-cosh (standard for super-Gaussian sources), use + instead of -
  W <- W - lr * grad

  # Normalize to project back to the unit sphere (norm = 1)
  # Add epsilon (1e-15) to prevent 0/0
  sweep(W, 2, sqrt(colSums(W^2)) + 1e-15, "/")
}

polar <- function(W) {
  eig <- eigen(t(W) %*% W, symmetric=TRUE)
  Ainvhalf <- eig$vectors %*% diag(1/sqrt(pmax(eig$values, 1e-14)), nrow = ncol(W)) %*% t(eig$vectors)
  W %*% Ainvhalf
}

fastica_update_rankr <- function(U, W) {
  P  <- t(U) %*% W           # n x r
  G  <- tanh(P)              # n x r
  G2 <- 1 - G^2             # n x r
  W  <- U %*% G - sweep(W, 2, colSums(G2), "*")
  polar(W)
}

gradient_minica_update_rankr <- function(U, W, lr = 0.1) {
  P    <- t(U) %*% W          # n x r
  G    <- tanh(P)             # n x r
  grad <- (U %*% G) / ncol(U)
  W    <- W - lr * grad
  polar(W)
}

fastica_update_tlc <- function(U, W, lambda = 1) {
  P  <- t(U) %*% W
  G  <- tanh(P) + 2 * lambda * abs(P)
  G2 <- 1 - tanh(P)^2 + 2 * lambda * sign(P)
  W  <- U %*% G - sweep(W, 2, colSums(G2), "*")
  sweep(W, 2, sqrt(colSums(W^2)) + 1e-15, "/")
}

fastica_update_skew <- function(U, W) {
  P  <- t(U) %*% W
  G  <- 2 * abs(P)
  G2 <- 2 * sign(P)
  W  <- U %*% G - sweep(W, 2, colSums(G2), "*")
  sweep(W, 2, sqrt(colSums(W^2)) + 1e-15, "/")
}

# Contrast G(x) = x*abs(x)/2; use sign(0) = 0 in the second derivative.
fastica_update_skew_rankr <- function(U, W) {
  P  <- t(U) %*% W
  G  <- abs(P)
  G2 <- sign(P)
  W  <- U %*% G - sweep(W, 2, colSums(G2), "*")
  polar(W)
}

fastica_update_skew_rankr_gs <- function(U, W) {
  W <- fastica_update_skew(U, W)
  P <- t(U) %*% W
  objective <- colMeans(P * abs(P) / 2)
  gs_reorth(W, objective)
}

# Stable log(cosh(x)), used as g = G' for integrated log-cosh ICA.
logcosh <- function(x) abs(x) + log1p(exp(-2 * abs(x))) - log(2)

# G(x) = integral_0^x log(cosh(t)) dt. Trapezoidal grid for objective values.
integrated_logcosh <- function(x) {
  grid <- seq(0, max(abs(x)) + 0.001, by = 0.001)
  y <- logcosh(grid)
  area <- c(0, cumsum((head(y, -1) + tail(y, -1)) * 0.001 / 2))
  value <- sign(x) * approx(grid, area, xout = abs(x))$y
  dim(value) <- dim(x)
  value
}

fastica_update_ilogcosh <- function(U, W) {
  P <- t(U) %*% W
  W <- U %*% logcosh(P) - sweep(W, 2, colSums(tanh(P)), "*")
  sweep(W, 2, sqrt(colSums(W^2)) + 1e-15, "/")
}

fastica_update_ilogcosh_rankr <- function(U, W) {
  P <- t(U) %*% W
  W <- U %*% logcosh(P) - sweep(W, 2, colSums(tanh(P)), "*")
  polar(W)
}

fastica_update_ilogcosh_rankr_gs <- function(U, W) {
  W <- fastica_update_ilogcosh(U, W)
  objective <- colMeans(integrated_logcosh(t(U) %*% W))
  gs_reorth(W, objective)
}

fastica_update_tlc_rankr <- function(U, W, lambda = 1) {
  P  <- t(U) %*% W
  G  <- tanh(P) + 2 * lambda * abs(P)
  G2 <- 1 - tanh(P)^2 + 2 * lambda * sign(P)
  W  <- U %*% G - sweep(W, 2, colSums(G2), "*")
  polar(W)
}

fastica_update_logphi <- function(U, W, alpha = 2) {
  P  <- t(U) %*% W
  u  <- alpha * P
  h  <- exp(dnorm(u, log = TRUE) - pnorm(u, log.p = TRUE))
  G  <- alpha * h
  G2 <- alpha^2 * (-u * h - h^2)
  W  <- U %*% G - sweep(W, 2, colSums(G2), "*")
  sweep(W, 2, sqrt(colSums(W^2)) + 1e-15, "/")
}

fastica_update_logphi_rankr <- function(U, W, alpha = 2) {
  P  <- t(U) %*% W
  u  <- alpha * P
  h  <- exp(dnorm(u, log = TRUE) - pnorm(u, log.p = TRUE))
  G  <- alpha * h
  G2 <- alpha^2 * (-u * h - h^2)
  W  <- U %*% G - sweep(W, 2, colSums(G2), "*")
  polar(W)
}

# Matthew's objective-ordered Gram-Schmidt, returning the original column order.
# As in his code, residuals with norm <= 1e-10 are left unnormalized.
gs_reorth <- function(W, objective) {
  ord <- order(objective, decreasing = TRUE)
  V <- W[, ord, drop = FALSE]
  for (j in seq_len(ncol(V))) {
    if (j > 1)
      for (l in seq_len(j - 1))
        V[, j] <- V[, j] - sum(V[, j] * V[, l]) * V[, l]
    nrm <- sqrt(sum(V[, j]^2))
    if (nrm > 1e-10) V[, j] <- V[, j] / nrm
  }
  V[, order(ord), drop = FALSE]
}

fastica_update_logphi_rankr_gs <- function(U, W, alpha = 2) {
  W <- fastica_update_logphi(U, W, alpha)
  objective <- colMeans(pnorm(alpha * (t(U) %*% W), log.p = TRUE))
  gs_reorth(W, objective)
}

# Contrast G(x) = x^3, with derivatives 3*x^2 and 6*x.
fastica_update_cubic <- function(U, W) {
  P  <- t(U) %*% W
  G  <- 3 * P^2
  G2 <- 6 * P
  W  <- U %*% G - sweep(W, 2, colSums(G2), "*")
  sweep(W, 2, sqrt(colSums(W^2)) + 1e-15, "/")
}

fastica_update_cubic_rankr <- function(U, W) {
  P  <- t(U) %*% W
  G  <- 3 * P^2
  G2 <- 6 * P
  W  <- U %*% G - sweep(W, 2, colSums(G2), "*")
  polar(W)
}

fastica_update_cubic_rankr_gs <- function(U, W) {
  W <- fastica_update_cubic(U, W)
  objective <- colMeans((t(U) %*% W)^3)
  gs_reorth(W, objective)
}

# Matthew's goldenplus schedule with objective-ordered GS and our unit priors.
# Requires eb_proj_newton() and the normalized Bernoulli/centered Exp solvers.
eb_proj_goldenplus_rankr_gs <- function(
    U, W, prior = c("bernoulli", "exponential"), n_warmup = 50, n_full = 100) {
  prior <- match.arg(prior)
  n <- nrow(U)
  r <- ncol(W)
  stopifnot(nrow(W) == ncol(U), r <= ncol(U), n > 2 * ncol(U))
  c_golden <- 2 / (1 + sqrt(5))
  prior_scale <- c_golden
  fixed_prior <- prior == "exponential"

  # Keep fitted_g on the unit-prior scale when moving between phases.
  scaled_prior <- function(x, s, g_init = NULL, fix_g = FALSE, output = NULL) {
    fit <- if (fixed_prior) {
      ebnm_centered_exponential(x / prior_scale, s / prior_scale)
    } else {
      ebnm_normalized_bernoulli(x / prior_scale, s / prior_scale,
        g_init = g_init, fix_g = fix_g, output = output)
    }
    fit$posterior$mean <- prior_scale * fit$posterior$mean
    fit$posterior$sd <- prior_scale * fit$posterior$sd
    fit$posterior$second_moment <- prior_scale^2 * fit$posterior$second_moment
    fit$log_likelihood <- as.numeric(fit$log_likelihood) - length(x) * log(prior_scale)
    fit
  }

  fits <- lapply(seq_len(r), function(j)
    eb_proj_newton(U, ebnm_fn = scaled_prior, w_init = W[, j],
      g_init = if (fixed_prior) NULL else 0.5, sigma2 = c_golden,
      fix_w = TRUE, fix_g = TRUE, fix_tau = TRUE, maxiter = 1))

  for (phase in 1:3) {
    prior_scale <- if (phase < 3) c_golden else 1
    n_iter <- if (phase < 3) n_warmup else n_full
    for (iter in seq_len(n_iter)) {
      for (j in seq_len(r)) {
        fit <- fits[[j]]
        fits[[j]] <- eb_proj_newton(U, ebnm_fn = scaled_prior,
          w_init = fit$w, g_init = fit$fitted_g, tau = fit$tau,
          hessian = if (phase == 1) "none" else "ica",
          fix_g = TRUE, fix_tau = TRUE, maxiter = 1)
      }
      W <- gs_reorth(sqrt(n) * sapply(fits, `[[`, "w"),
                     sapply(fits, `[[`, "objective"))
      # Refresh prior and precision at the directions AFTER orthogonalization.
      for (j in seq_len(r)) {
        fit <- fits[[j]]
        fits[[j]] <- eb_proj_newton(U, ebnm_fn = scaled_prior,
          w_init = W[, j], g_init = fit$fitted_g, tau = fit$tau,
          fix_w = TRUE, fix_g = fixed_prior, fix_tau = phase < 3, maxiter = 1)
      }
    }
  }
  list(loadings = sapply(fits, `[[`, "x"),
       objective = sapply(fits, `[[`, "objective"), fits = fits)
}

log_Z_binary = function(x, y_0, y_1, p, sigma2) {
  # Use -(x-y_j)^2/(2*sigma2) <= 0 to avoid overflow; add back x^2/(2*sigma2) at the end
  lp0 = log(1 - p) - (x - y_0)^2 / (2*sigma2)
  lp1 = log(p)     - (x - y_1)^2 / (2*sigma2)
  m   = pmax(lp0, lp1)
  log(exp(lp0 - m) + exp(lp1 - m)) + m + x^2 / (2*sigma2)
}

post_prob1 = function(x, y_0, y_1, p, sigma2) {
  # x^2/(2*sigma2) cancels in the ratio so omit it
  lp0 = log(1 - p) - (x - y_0)^2 / (2*sigma2)
  lp1 = log(p)     - (x - y_1)^2 / (2*sigma2)
  m   = pmax(lp0, lp1)
  e0  = exp(lp0 - m);  e1 = exp(lp1 - m)
  e1 / (e0 + e1)
}

post_mean_binary = function(x, y_0, y_1, p, sigma2) {
  pi1 = post_prob1(x, y_0, y_1, p, sigma2)
  (1 - pi1)*y_0 + pi1*y_1
}

post_mean_deriv_binary = function(x, y_0, y_1, p, sigma2) {
  pi1 = post_prob1(x, y_0, y_1, p, sigma2)
  pi1 * (1 - pi1) * (y_1 - y_0)^2 / sigma2
}

ebproj_init = function(U, b=0, c=1, tau=NULL, sigma2=NULL, w=NULL, d.init=NULL) {
  n = nrow(U); k = ncol(U); nu = n - k - 2*b
  if (!is.null(tau) && !is.null(sigma2)) stop("specify at most one of tau and sigma2")
  if (is.null(tau) && is.null(sigma2)) tau = 1
  if (is.null(sigma2)) sigma2 = n / (nu * tau)
  if (is.null(tau))    tau    = n / (nu * sigma2)
  if (is.null(d.init)) d.init = rep(n, k)
  d = rep(n, k)
  if (is.null(w)) w = sqrt(d.init) * rnorm(k)
  w = w / sqrt(sum(d * w^2))
  list(U=U, d=d, b=b, n=n, k=k, nu=nu, sigma2=sigma2, w=w, c=c,
       y_0=-1, y_1=1, p=0.5, objective=NULL, tau=tau)
}

ebproj_calc_x = function(fit) as.vector(fit$U %*% (fit$d * fit$w))

ebproj_objective = function(fit) {
  ey0 = fit$c * fit$y_0;  ey1 = fit$c * fit$y_1
  x   = ebproj_calc_x(fit)
  J   = sum(log_Z_binary(x, ey0, ey1, fit$p, fit$sigma2))
  (fit$nu/2) * (log(fit$tau) - fit$tau + 1) + J
}

ebproj_update_w = function(fit, hessian="ica", eps=1e-6) {
  ey0 = fit$c * fit$y_0;  ey1 = fit$c * fit$y_1
  x        = ebproj_calc_x(fit)
  Mx       = post_mean_binary(x, ey0, ey1, fit$p, fit$sigma2)
  Mpx      = post_mean_deriv_binary(x, ey0, ey1, fit$p, fit$sigma2)
  w_target = as.vector(t(fit$U) %*% Mx)
  S_diag   = as.vector((fit$U^2) %*% fit$d)
  H_diag   = switch(hessian,
    diag  = { c_j = as.vector(t(fit$U^2) %*% Mpx); c_j * fit$d },
    trace = { c_t = sum(Mpx * S_diag) / sum(fit$d); c_t * fit$d },
    iso   = { c_i = sum(Mpx * S_diag) / fit$k; rep(c_i, fit$k) },
    ica   = { c_c = sum(Mpx); rep(c_c, fit$k) },
    none  = rep(0, fit$k),
    stop("hessian must be 'diag', 'trace', 'iso', 'ica', or 'none'")
  )
  lambda  = max(sum(fit$d * fit$w * w_target), max(H_diag) + eps)
  w_new   = (w_target - H_diag * fit$w) / (lambda - H_diag)
  scale   = sqrt(sum(fit$d * w_new^2))
  if (scale < 1e-10) w_new = w_target   # fall back to gradient step if Newton collapses
  fit$w   = w_new / sqrt(sum(fit$d * w_new^2))
  fit
}

ebproj_update_g = function(fit) {
  ey0 = fit$c * fit$y_0;  ey1 = fit$c * fit$y_1
  x   = ebproj_calc_x(fit)
  pi1 = post_prob1(x, ey0, ey1, fit$p, fit$sigma2)
  pi0 = 1 - pi1
  fit$p  = mean(pi1)
  s0 = sum(pi0); s1 = sum(pi1)
  y0_new = if (s0 > 1e-10) sum(pi0 * x) / s0 else fit$c * fit$y_0
  y1_new = if (s1 > 1e-10) sum(pi1 * x) / s1 else fit$c * fit$y_1
  var_g   = (1 - fit$p) * y0_new^2 + fit$p * y1_new^2 -
            ((1 - fit$p) * y0_new + fit$p * y1_new)^2
  if (var_g < 1e-10) return(fit)   # degenerate posterior; keep current g
  fit$c   = sqrt(var_g)
  fit$y_0 = y0_new / fit$c
  fit$y_1 = y1_new / fit$c
  fit
}

ebproj_update_g0 = function(fit) {
  x = ebproj_calc_x(fit)
  obj_p = function(logit_p) {
    pp = 1 / (1 + exp(-logit_p))
    y1 =  sqrt((1 - pp) / pp);  y0 = -sqrt(pp / (1 - pp))
    sum(log_Z_binary(x, fit$c * y0, fit$c * y1, pp, fit$sigma2))
  }
  res     = optimize(obj_p, c(-5, 5), maximum=TRUE)
  pp      = 1 / (1 + exp(-res$maximum))
  fit$p   = pp
  fit$y_1 =  sqrt((1 - pp) / pp)
  fit$y_0 = -sqrt(pp / (1 - pp))
  fit
}

ebproj_update_tau = function(fit) {
  ey0 = fit$c * fit$y_0;  ey1 = fit$c * fit$y_1
  x    = ebproj_calc_x(fit)
  pi1  = post_prob1(x, ey0, ey1, fit$p, fit$sigma2)
  vbar = post_mean_binary(x, ey0, ey1, fit$p, fit$sigma2)
  Uv         = as.vector(t(fit$U) %*% vbar)
  Uq         = sum((1 - pi1) * ey0^2 + pi1 * ey1^2)
  rho        = sum(fit$d * Uv^2) / (fit$n * Uq)
  rho        = min(rho, 1 - 1/fit$n)   # cap tau at n to avoid numerical blow-up
  fit$tau    = 1 / (1 - rho)
  fit$sigma2 = fit$n / (fit$nu * fit$tau)
  fit
}

ebproj_fit = function(fit, hessian="ica", fix_w=FALSE, g_update="g",
                      fix_tau=FALSE, max_iter=200, tol=1e-6, verbose=FALSE) {
  fit$objective = ebproj_objective(fit)
  for (iter in seq_len(max_iter)) {
    if (!fix_tau) fit = ebproj_update_tau(fit)
    fit = switch(g_update,
      g    = ebproj_update_g(fit),
      g0   = ebproj_update_g0(fit),
      none = fit,
      stop("g_update must be 'g', 'g0', or 'none'")
    )
    if (!fix_w) fit = ebproj_update_w(fit, hessian=hessian)
    obj_new = ebproj_objective(fit)
    if (verbose) cat(sprintf("iter %3d  obj=%.4f  tau=%.2f  p=%.3f  c=%.4f\n",
                             iter, obj_new, fit$tau, fit$p, fit$c))
    if (!is.null(fit$objective) && abs(obj_new - fit$objective) < tol) {
      fit$objective = obj_new; break
    }
    fit$objective = obj_new
  }
  fit$iter = iter
  if (iter == max_iter) warning("ebproj_fit reached max_iter")
  fit
}
