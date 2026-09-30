#' Initialize an SSLfmm model
#'
#' Constructs stable starting values for SSLfmm model fitting. Class-specific labelled
#' means are used when available; otherwise k-means centers provide fallbacks.
#'
#' @param x Numeric feature matrix or data.frame.
#' @param y Class labels in 1:g (or factor/character), with NA for missing labels.
#' @param g Number of mixture components. Inferred from observed labels when possible.
#' @param covariance_type Either "equal" or "unequal".
#' @param ridge Positive diagonal regularization added to covariance starts.
#' @param seed Optional random seed used only for k-means fallback initialization.
#' @return A list with pi, mu, sigma, g, covariance_type and encoded label levels.
#' @export
initialize_sslfmm <- function(x, y, g = NULL,
                              covariance_type = c("equal", "unequal"),
                              ridge = 1e-4, seed = NULL) {
  covariance_type <- match.arg(covariance_type)
  if (length(ridge) != 1L || !is.numeric(ridge) || !is.finite(ridge) || ridge <= 0)
    stop("'ridge' must be a single positive finite number.", call. = FALSE)
  x <- .check_x(x)
  lab <- .encode_labels(y, g)
  z <- lab$z; g <- lab$g
  if (length(z) != nrow(x)) stop("'y' must have one entry per row of 'x'.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(x); p <- ncol(x)
  obs <- !is.na(z)
  global_mu <- colMeans(x)
  global_cov <- .safe_cov(x, ridge = ridge)

  # k-means is only a fallback; observed class labels anchor class identities.
  km_centers <- matrix(rep(global_mu, times = g), nrow = g, byrow = TRUE)
  if (g > 1L && n >= g) {
    km <- try(stats::kmeans(x, centers = g, nstart = 10L), silent = TRUE)
    if (!inherits(km, "try-error")) km_centers <- km$centers
  }

  mu <- matrix(NA_real_, g, p)
  used_center <- rep(FALSE, g)
  for (k in seq_len(g)) {
    ik <- obs & z == k
    if (any(ik)) {
      mu[k, ] <- colMeans(x[ik, , drop = FALSE])
      if (g > 1L) {
        d <- rowSums((sweep(km_centers, 2L, mu[k, ], "-"))^2)
        d[used_center] <- Inf
        j <- which.min(d)
        if (length(j) && is.finite(d[j])) used_center[j] <- TRUE
      }
    }
  }
  for (k in seq_len(g)) {
    if (anyNA(mu[k, ])) {
      cand <- which(!used_center)
      j <- if (length(cand)) cand[1L] else k
      mu[k, ] <- km_centers[j, ]
      if (j <= length(used_center)) used_center[j] <- TRUE
    }
  }

  counts <- tabulate(z[obs], nbins = g)
  pi_hat <- (counts + 0.5) / (sum(counts) + 0.5 * g)
  if (!any(obs)) pi_hat <- rep(1 / g, g)

  if (covariance_type == "equal") {
    resid <- matrix(NA_real_, sum(obs), p)
    if (any(obs)) {
      xo <- x[obs, , drop = FALSE]
      zo <- z[obs]
      for (i in seq_len(nrow(xo))) resid[i, ] <- xo[i, ] - mu[zo[i], ]
      if (nrow(resid) > 1L) {
        S <- crossprod(resid) / nrow(resid)
        S <- .ensure_spd(S, ridge = ridge)
      } else S <- global_cov
    } else S <- global_cov
    sigma <- S
  } else {
    sigma <- array(NA_real_, dim = c(p, p, g))
    for (k in seq_len(g)) {
      ik <- obs & z == k
      sigma[, , k] <- if (sum(ik) >= 2L) .safe_cov(x[ik, , drop = FALSE], ridge) else global_cov
    }
  }

  feature_names <- colnames(x)
  if (is.null(feature_names)) feature_names <- paste0("x", seq_len(p))
  names(pi_hat) <- lab$levels[seq_len(g)]
  rownames(mu) <- lab$levels[seq_len(g)]
  colnames(mu) <- feature_names

  structure(list(pi = pi_hat, mu = mu, sigma = sigma,
                 g = g, p = p, covariance_type = covariance_type,
                 label_levels = lab$levels), class = "SSLfmm_init")
}

.initial_xi <- function(x, theta, covariance_type, response,
                        eligible = rep(TRUE, nrow(x)), entropy_floor = 1e-300) {
  en <- .posterior_entropy(x, theta, covariance_type, entropy_floor)
  yy <- as.integer(response[eligible])
  xx <- log(en[eligible])
  if (length(yy) < 10L || length(unique(yy)) < 2L) {
    pr <- if (length(yy)) mean(yy) else 0.2
    pr <- min(max(pr, 1e-3), 1 - 1e-3)
    return(c(xi0 = stats::qlogis(pr), eta_xi = 0))
  }
  fit <- try(stats::glm(yy ~ xx, family = stats::binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) {
    pr <- min(max(mean(yy), 1e-3), 1 - 1e-3)
    return(c(xi0 = stats::qlogis(pr), eta_xi = 0))
  }
  cf <- stats::coef(fit)
  if (length(cf) < 2L || any(!is.finite(cf))) {
    pr <- min(max(mean(yy), 1e-3), 1 - 1e-3)
    return(c(xi0 = stats::qlogis(pr), eta_xi = 0))
  }
  slope <- unname(cf[2L])
  if (!is.finite(slope) || slope <= 0) slope <- 1
  c(xi0 = unname(cf[1L]), eta_xi = log(max(slope, 1e-3)))
}
