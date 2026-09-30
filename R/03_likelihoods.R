.nll_mcar <- function(theta_par, data, p, g, covariance_type, entropy_floor = 1e-300) {
  theta <- .unpack_theta(theta_par, p, g, covariance_type)
  lc <- .model_components(data$x, theta, covariance_type)
  lj <- .known_label_logjoint(data$z, lc$log_joint)
  obs <- !data$m
  val <- sum(lj[obs]) + sum(lc$log_marginal[data$m])
  if (!is.finite(val)) return(1e100)
  -val
}

.nll_mar <- function(par, data, p, g, covariance_type, entropy_floor = 1e-300) {
  dt <- .theta_length(p, g, covariance_type)
  theta <- .unpack_theta(par[seq_len(dt)], p, g, covariance_type)
  xi0 <- par[dt + 1L]
  xi1 <- exp(par[dt + 2L])
  lc <- .model_components(data$x, theta, covariance_type)
  lj <- .known_label_logjoint(data$z, lc$log_joint)
  en <- .posterior_entropy(data$x, theta, covariance_type, entropy_floor)
  eta <- xi0 + xi1 * log(en)
  lq <- .log_sigmoid(eta)
  l1q <- .log1m_sigmoid(eta)
  obs <- !data$m; mis <- data$m
  val <- sum(lj[obs] + l1q[obs]) + sum(lc$log_marginal[mis] + lq[mis])
  if (!is.finite(val)) return(1e100)
  -val
}

.nll_mixed_observed <- function(par, data, p, g, covariance_type, alpha,
                                entropy_floor = 1e-300) {
  dt <- .theta_length(p, g, covariance_type)
  theta <- .unpack_theta(par[seq_len(dt)], p, g, covariance_type)
  xi0 <- par[dt + 1L]
  xi1 <- exp(par[dt + 2L])
  lc <- .model_components(data$x, theta, covariance_type)
  lj <- .known_label_logjoint(data$z, lc$log_joint)
  en <- .posterior_entropy(data$x, theta, covariance_type, entropy_floor)
  eta <- xi0 + xi1 * log(en)
  lq <- .log_sigmoid(eta)
  l1q <- .log1m_sigmoid(eta)

  iobs <- data$source == "observed"
  imcar <- data$source == "mcar"
  imar <- data$source == "mar"
  if (any(iobs & data$m) || any((imcar | imar) & !data$m)) return(1e100)

  val <- 0
  if (any(iobs)) val <- val + sum(lj[iobs] + log1p(-alpha) + l1q[iobs])
  if (any(imcar)) {
    if (alpha <= 0) return(1e100)
    val <- val + sum(lc$log_marginal[imcar] + log(alpha))
  }
  if (any(imar)) val <- val + sum(lc$log_marginal[imar] + log1p(-alpha) + lq[imar])
  if (!is.finite(val)) return(1e100)
  -val
}

.nll_mixed_latent <- function(par, data, p, g, covariance_type,
                              entropy_floor = 1e-300) {
  dt <- .theta_length(p, g, covariance_type)
  theta <- .unpack_theta(par[seq_len(dt)], p, g, covariance_type)
  xi0 <- par[dt + 1L]
  xi1 <- exp(par[dt + 2L])
  alpha <- par[dt + 3L]
  if (!is.finite(alpha) || alpha < 0 || alpha >= 1) return(1e100)

  lc <- .model_components(data$x, theta, covariance_type)
  lj <- .known_label_logjoint(data$z, lc$log_joint)
  en <- .posterior_entropy(data$x, theta, covariance_type, entropy_floor)
  eta <- xi0 + xi1 * log(en)
  lq <- .log_sigmoid(eta)
  l1q <- .log1m_sigmoid(eta)
  obs <- !data$m; mis <- data$m
  log1ma <- log1p(-alpha)
  lmiss <- if (alpha == 0) lq else .logsumexp2(rep(log(alpha), length(lq)), log1ma + lq)
  val <- sum(lj[obs] + log1ma + l1q[obs]) +
    sum(lc$log_marginal[mis] + lmiss[mis])
  if (!is.finite(val)) return(1e100)
  -val
}
