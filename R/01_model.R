.n_chol_par <- function(p) p * (p + 1L) / 2L

.chol_positions <- function(p) {
  out <- matrix(0L, .n_chol_par(p), 2L)
  h <- 1L
  for (i in seq_len(p)) {
    for (j in seq_len(i)) {
      out[h, ] <- c(i, j)
      h <- h + 1L
    }
  }
  out
}

.pack_chol <- function(S) {
  p <- nrow(S)
  L <- t(chol(S))
  pos <- .chol_positions(p)
  out <- numeric(nrow(pos))
  for (h in seq_len(nrow(pos))) {
    i <- pos[h, 1L]; j <- pos[h, 2L]
    out[h] <- if (i == j) log(L[i, j]) else L[i, j]
  }
  out
}

.unpack_chol <- function(par, p) {
  pos <- .chol_positions(p)
  if (length(par) != nrow(pos)) stop("Internal covariance parameter length mismatch.", call. = FALSE)
  L <- matrix(0, p, p)
  for (h in seq_len(nrow(pos))) {
    i <- pos[h, 1L]; j <- pos[h, 2L]
    L[i, j] <- if (i == j) exp(par[h]) else par[h]
  }
  L %*% t(L)
}

.softmax_baseline <- function(eta) {
  v <- c(eta, 0)
  m <- max(v)
  ex <- exp(v - m)
  ex / sum(ex)
}

.theta_length <- function(p, g, covariance_type) {
  mix <- max(0L, g - 1L)
  covn <- if (covariance_type == "equal") .n_chol_par(p) else g * .n_chol_par(p)
  mix + g * p + covn
}

.pack_theta <- function(theta, covariance_type) {
  pi <- as.numeric(theta$pi)
  g <- length(pi)
  mu <- as.matrix(theta$mu)
  p <- ncol(mu)
  if (nrow(mu) != g) stop("Internal 'mu' dimension mismatch.", call. = FALSE)
  eta <- if (g == 1L) numeric(0) else log(pi[seq_len(g - 1L)] / pi[g])
  out <- c(eta, as.numeric(t(mu)))
  if (covariance_type == "equal") {
    out <- c(out, .pack_chol(theta$sigma))
  } else {
    for (k in seq_len(g)) {
      Sk <- matrix(theta$sigma[, , k], nrow = p, ncol = p)
      out <- c(out, .pack_chol(Sk))
    }
  }
  out
}

.unpack_theta <- function(par, p, g, covariance_type) {
  need <- .theta_length(p, g, covariance_type)
  if (length(par) != need) stop("Internal model parameter length mismatch.", call. = FALSE)
  h <- 1L
  if (g == 1L) {
    pi <- 1
  } else {
    eta <- par[h:(h + g - 2L)]
    h <- h + g - 1L
    pi <- .softmax_baseline(eta)
  }
  mu_vec <- par[h:(h + g * p - 1L)]
  h <- h + g * p
  mu <- matrix(mu_vec, nrow = g, ncol = p, byrow = TRUE)
  nc <- .n_chol_par(p)
  if (covariance_type == "equal") {
    sigma <- .unpack_chol(par[h:(h + nc - 1L)], p)
  } else {
    sigma <- array(NA_real_, dim = c(p, p, g))
    for (k in seq_len(g)) {
      sigma[, , k] <- .unpack_chol(par[h:(h + nc - 1L)], p)
      h <- h + nc
    }
  }
  list(pi = pi, mu = mu, sigma = sigma)
}

.theta_bounds <- function(p, g, covariance_type) {
  d <- .theta_length(p, g, covariance_type)
  lo <- rep(-Inf, d); hi <- rep(Inf, d)
  mixn <- max(0L, g - 1L)
  if (mixn > 0L) {
    lo[seq_len(mixn)] <- -8
    hi[seq_len(mixn)] <- 8
  }
  cov_start <- mixn + g * p + 1L
  pos <- .chol_positions(p)
  blocks <- if (covariance_type == "equal") 1L else g
  nc <- nrow(pos)
  for (b in seq_len(blocks)) {
    for (j in seq_len(nc)) {
      idx <- cov_start + (b - 1L) * nc + j - 1L
      if (pos[j, 1L] == pos[j, 2L]) {
        lo[idx] <- -5; hi[idx] <- 5
      } else {
        lo[idx] <- -20; hi[idx] <- 20
      }
    }
  }
  list(lower = lo, upper = hi)
}

.dmvn_log <- function(y, mu, sigma) {
  y <- as.matrix(y)
  p <- ncol(y)
  n <- nrow(y)
  if (!length(mu) || length(mu) != p || any(!is.finite(mu)) ||
      any(!is.finite(sigma))) return(rep(-Inf, n))

  R <- tryCatch(chol(sigma), error = function(e) NULL)
  if (is.null(R) || any(!is.finite(R)) || any(diag(R) <= 0)) {
    return(rep(-Inf, n))
  }

  z <- sweep(y, 2L, mu, "-")
  if (any(!is.finite(z))) return(rep(-Inf, n))
  zz <- tryCatch(
    backsolve(R, t(z), transpose = TRUE),
    error = function(e) NULL
  )
  if (is.null(zz) || any(!is.finite(zz))) return(rep(-Inf, n))

  quad <- colSums(zz^2)
  logdet <- 2 * sum(log(diag(R)))
  ans <- -0.5 * (p * log(2 * pi) + logdet + quad)
  ans[!is.finite(ans)] <- -Inf
  ans
}

.model_components <- function(x, theta, covariance_type) {
  x <- as.matrix(x)
  g <- length(theta$pi)
  n <- nrow(x)
  lj <- matrix(NA_real_, n, g)
  for (k in seq_len(g)) {
    S <- if (covariance_type == "equal") theta$sigma else
      matrix(theta$sigma[, , k], nrow = ncol(x), ncol = ncol(x))
    lj[, k] <- log(theta$pi[k]) + .dmvn_log(x, theta$mu[k, ], S)
  }
  lm <- .row_logsumexp(lj)
  post <- matrix(1 / g, nrow = n, ncol = g)
  good <- is.finite(lm)
  if (any(good)) {
    post[good, ] <- exp(lj[good, , drop = FALSE] - lm[good])
    rs <- rowSums(post[good, , drop = FALSE])
    okrs <- is.finite(rs) & rs > 0
    if (any(okrs)) {
      post[which(good)[okrs], ] <- post[which(good)[okrs], , drop = FALSE] / rs[okrs]
    }
  }
  list(log_joint = lj, log_marginal = lm, posterior = post)
}

.posterior_entropy <- function(x, theta, covariance_type, entropy_floor = 1e-300) {
  post <- .model_components(x, theta, covariance_type)$posterior
  val <- -rowSums(ifelse(post > 0, post * log(post), 0))
  pmax(val, entropy_floor)
}

.known_label_logjoint <- function(z, log_joint) {
  out <- rep(NA_real_, length(z))
  obs <- !is.na(z)
  if (any(obs)) out[obs] <- log_joint[cbind(which(obs), z[obs])]
  out
}
