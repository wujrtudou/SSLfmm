.validate_mixture <- function(pi, mu, sigma, covariance_type = NULL, ridge = 0) {
  pi <- as.numeric(pi)
  if (!length(pi) || any(!is.finite(pi)) || any(pi <= 0)) {
    stop("'pi' must contain positive finite mixing proportions.", call. = FALSE)
  }
  if (abs(sum(pi) - 1) > 1e-8) stop("'pi' must sum to 1.", call. = FALSE)
  g <- length(pi)

  if (is.data.frame(mu)) {
    if (!all(vapply(mu, is.numeric, logical(1)))) {
      stop("'mu' must be a numeric p x g matrix, with one column per mixture component.", call. = FALSE)
    }
  } else if (!is.matrix(mu) || !is.numeric(mu)) {
    stop("'mu' must be a numeric p x g matrix, with one column per mixture component.", call. = FALSE)
  }
  mu <- as.matrix(mu)
  storage.mode(mu) <- "double"
  if (length(dim(mu)) != 2L || nrow(mu) < 1L || ncol(mu) != g) {
    stop(sprintf("'mu' must be a p x g matrix with %d column%s (one per mixture component).",
                 g, if (g == 1L) "" else "s"), call. = FALSE)
  }
  if (any(!is.finite(mu))) stop("'mu' must contain only finite values.", call. = FALSE)
  p <- nrow(mu)

  if (!is.null(covariance_type)) {
    covariance_type <- match.arg(covariance_type, c("equal", "unequal"))
  }
  ns <- .normalize_sigma(sigma, p, g, covariance_type = covariance_type, ridge = ridge)
  list(pi = pi, mu = mu, sigma = ns$sigma,
       covariance_type = ns$covariance_type, p = p, g = g)
}

.rmvnorm_chol <- function(n, mu, sigma) {
  p <- length(mu)
  if (n <= 0L) return(matrix(numeric(0), 0L, p))
  R <- chol(sigma)
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  sweep(z %*% R, 2L, mu, "+")
}

#' Simulate from a Gaussian finite mixture
#'
#' @param n Number of observations.
#' @param pi Positive mixing proportions of length g that sum to one.
#' @param mu A numeric p x g matrix of component means.
#' @param sigma Either a shared p x p covariance matrix or a p x p x g array.
#' @param seed_number Optional random seed.
#' @return A data.frame with columns x1,...,xp and truth.
#' @export
rmix <- function(n = 200L, pi, mu, sigma, seed_number = NULL) {
  if (length(n) != 1L || !is.finite(n) || n < 1 || n != as.integer(n)) {
    stop("'n' must be a positive integer.", call. = FALSE)
  }
  setup <- .validate_mixture(pi, mu, sigma)
  if (!is.null(seed_number)) set.seed(seed_number)
  n <- as.integer(n); g <- setup$g; p <- setup$p
  truth <- sample.int(g, size = n, replace = TRUE, prob = setup$pi)
  x <- matrix(NA_real_, n, p)
  for (k in seq_len(g)) {
    ik <- which(truth == k)
    if (!length(ik)) next
    S <- if (setup$covariance_type == "equal") setup$sigma else
      matrix(setup$sigma[, , k], nrow = p, ncol = p)
    x[ik, ] <- .rmvnorm_chol(length(ik), setup$mu[, k], S)
  }
  out <- as.data.frame(x)
  names(out) <- paste0("x", seq_len(p))
  out$truth <- truth
  out
}

.calibrate_xi0_empirical <- function(entropy, xi1, target) {
  if (length(target) != 1L || !is.finite(target) || target <= 0 || target >= 1) {
    stop("'mar_rate' must lie strictly between 0 and 1.", call. = FALSE)
  }
  entropy <- as.numeric(entropy)
  if (!length(entropy) || any(!is.finite(entropy)) || any(entropy <= 0)) {
    stop("Entropy values used to calibrate 'xi0' must be positive and finite.", call. = FALSE)
  }
  if (length(xi1) != 1L || !is.finite(xi1) || xi1 <= 0) {
    stop("'xi1' must be a single positive number.", call. = FALSE)
  }

  # Write eta_i = xi0 + offset_i. Because plogis() is monotone,
  # qlogis(target) - max(offset) and qlogis(target) - min(offset)
  # bracket the empirical calibration equation. Constant entropy (including
  # g = 1) therefore has an immediate analytic solution.
  offset <- xi1 * log(entropy)
  qtarget <- stats::qlogis(target)
  lo <- qtarget - max(offset)
  hi <- qtarget - min(offset)

  if (!is.finite(lo) || !is.finite(hi)) {
    stop("Could not calibrate 'xi0': the entropy offset is not finite.", call. = FALSE)
  }

  if (abs(hi - lo) <= 100 * .Machine$double.eps * max(1, abs(lo), abs(hi))) {
    return((lo + hi) / 2)
  }

  f <- function(xi0) mean(stats::plogis(xi0 + offset)) - target
  flo <- f(lo)
  fhi <- f(hi)
  tol <- 1e-12
  if (abs(flo) <= tol) return(lo)
  if (abs(fhi) <= tol) return(hi)
  if (!is.finite(flo) || !is.finite(fhi) || flo > 0 || fhi < 0) {
    stop("Could not calibrate 'xi0' to the requested 'mar_rate'.", call. = FALSE)
  }
  stats::uniroot(f, interval = c(lo, hi), tol = 1e-10)$root
}

#' Simulate semi-supervised finite-mixture data with label missingness
#'
#' @description
#' Simulates semi-supervised Gaussian finite-mixture data with class-label
#' missingness under complete-case, MCAR, entropy-dependent MAR, or mixed
#' MCAR/MAR mechanisms.
#'
#' Generates Gaussian finite-mixture features and hides class labels under complete
#' observation (cc), MCAR, entropy-dependent MAR, or mixed MCAR/MAR missingness.
#' For mixed missingness, the MCAR trigger has probability alpha and the potential
#' MAR trigger has observation-specific probability q_i. MCAR takes precedence when
#' both triggers occur; this gives the same observed probabilities as the sequential
#' construction alpha + (1 - alpha) q_i while retaining the latent strata needed for
#' reproducible simulation diagnostics.
#'
#' @param n Number of observations.
#' @param pi Positive mixing proportions summing to one.
#' @param mu A numeric p x g matrix of component means.
#' @param sigma Shared p x p covariance matrix or p x p x g covariance array.
#' @param mechanism One of "cc", "mcar", "mar", or "mixed".
#' @param alpha MCAR probability for the mixed mechanism. For mechanism="mcar", this
#'   is the missing-label rate when missing_rate is NULL.
#' @param missing_rate Optional MCAR missing-label rate for mechanism="mcar" only.
#' @param xi0 MAR intercept. If NULL, it is calibrated to mar_rate on the generated data.
#' @param xi1 Positive entropy slope for MAR missingness.
#' @param mar_rate Target mean MAR probability used when xi0 is NULL.
#' @param seed Optional random seed.
#' @param entropy_floor Positive entropy floor.
#' @return A list with data, true_setup, groups, probs and raw. See the help page for
#'   the exact, stable field names and column order.
#' @export
simulate_sslfmm <- function(n = 500L, pi, mu, sigma,
                            mechanism = c("mixed", "mar", "mcar", "cc"),
                            alpha = 0.10, missing_rate = NULL,
                            xi0 = NULL, xi1 = 3, mar_rate = 0.30,
                            seed = NULL, entropy_floor = 1e-300) {
  mechanism <- match.arg(mechanism)
  setup <- .validate_mixture(pi, mu, sigma)
  if (!is.null(seed)) set.seed(seed)

  if (length(entropy_floor) != 1L || !is.numeric(entropy_floor) ||
      !is.finite(entropy_floor) || entropy_floor <= 0) {
    stop("'entropy_floor' must be a single positive finite number.", call. = FALSE)
  }

  if (!is.null(missing_rate) && mechanism != "mcar") {
    stop("'missing_rate' is only used when mechanism = 'mcar'.", call. = FALSE)
  }

  if (mechanism == "mixed" || (mechanism == "mcar" && is.null(missing_rate))) {
    if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) || alpha < 0 || alpha > 1) {
      stop("'alpha' must lie in [0, 1] when used as an MCAR probability.", call. = FALSE)
    }
  }

  if (mechanism %in% c("mar", "mixed")) {
    if (length(xi1) != 1L || !is.numeric(xi1) || !is.finite(xi1) || xi1 <= 0) {
      stop("'xi1' must be a single positive finite number for MAR or mixed simulation.", call. = FALSE)
    }
    if (!is.null(xi0) && (length(xi0) != 1L || !is.numeric(xi0) || !is.finite(xi0))) {
      stop("'xi0' must be NULL or a single finite number.", call. = FALSE)
    }
  }

  rate <- NULL
  if (mechanism == "mcar") {
    rate <- missing_rate %||% alpha
    if (length(rate) != 1L || !is.numeric(rate) || !is.finite(rate) || rate < 0 || rate > 1) {
      stop("The MCAR missing-label rate must lie in [0, 1].", call. = FALSE)
    }
  }

  raw <- rmix(n = n, pi = setup$pi, mu = setup$mu, sigma = setup$sigma)
  x <- as.matrix(raw[paste0("x", seq_len(setup$p))])
  theta <- list(pi = setup$pi, mu = t(setup$mu), sigma = setup$sigma)
  entropy <- .posterior_entropy(x, theta, setup$covariance_type, entropy_floor)

  if (is.null(xi0) && mechanism %in% c("mar", "mixed")) {
    xi0 <- .calibrate_xi0_empirical(entropy, xi1, mar_rate)
  }
  if (is.null(xi0)) xi0 <- NA_real_
  q <- if (mechanism %in% c("mar", "mixed")) {
    stats::plogis(xi0 + xi1 * log(entropy))
  } else {
    rep(0, nrow(x))
  }

  # potential_mar records the latent MAR Bernoulli draw before MCAR precedence.
  # This lets the returned groups reproduce the documented mar_group / obs_group /
  # mcar_in_mar / mcar_in_obs decomposition without changing the observed mechanism.
  potential_mar <- rep(FALSE, nrow(x))
  if (mechanism %in% c("mar", "mixed")) {
    potential_mar <- stats::runif(nrow(x)) < q
  }

  latent_missing <- rep(FALSE, nrow(x)) # true MCAR-channel trigger
  if (mechanism == "mcar") {
    latent_missing <- stats::runif(nrow(x)) < rate
  } else if (mechanism == "mixed") {
    latent_missing <- stats::runif(nrow(x)) < alpha
  }

  mar_missing <- potential_mar & !latent_missing
  observed_missing <- switch(
    mechanism,
    cc = rep(FALSE, nrow(x)),
    mcar = latent_missing,
    mar = potential_mar,
    mixed = latent_missing | potential_mar
  )

  source <- rep("observed", nrow(x))
  source[mar_missing] <- "mar"
  source[latent_missing] <- "mcar"

  label <- raw$truth
  label[observed_missing] <- NA_integer_

  # Stable documented column order. The first block retains the
  # x1..xp, en, missing, label, truth interface. The remaining columns expose
  # explicit missingness indicators, source labels, MAR probabilities, and entropy.
  dat <- raw[paste0("x", seq_len(setup$p))]
  dat$en <- entropy
  dat$missing <- observed_missing
  dat$label <- label
  dat$truth <- raw$truth
  dat$observed_missing <- observed_missing
  dat$latent_missing <- latent_missing
  dat$missing_source <- source
  dat$prob_mar <- q
  dat$entropy <- entropy

  # raw starts from the exact rmix() output and is augmented rather than replaced.
  raw_aug <- raw
  raw_aug$en <- entropy
  raw_aug$missing <- observed_missing
  raw_aug$label <- label
  raw_aug$observed_missing <- observed_missing
  raw_aug$latent_missing <- latent_missing
  raw_aug$missing_source <- source
  raw_aug$potential_mar <- potential_mar
  raw_aug$prob_mar <- q
  raw_aug$entropy <- entropy

  groups <- list(
    mar_group = which(potential_mar),
    obs_group = which(!potential_mar),
    mcar_in_mar = which(latent_missing & potential_mar),
    mcar_in_obs = which(latent_missing & !potential_mar),
    observed = which(!observed_missing),
    mcar = which(latent_missing),
    mar = which(mar_missing),
    missing = which(observed_missing)
  )

  list(
    data = dat,
    true_setup = list(
      pi = setup$pi,
      mu = setup$mu,
      sigma = setup$sigma,
      covariance_type = setup$covariance_type,
      mechanism = mechanism,
      alpha = if (mechanism == "mixed") alpha else 0,
      missing_rate = if (mechanism == "mcar") rate else NULL,
      xi = if (mechanism %in% c("mar", "mixed")) c(xi0 = xi0, xi1 = xi1) else NULL
    ),
    groups = groups,
    probs = q,
    raw = raw_aug
  )
}

#' Simulate data under the mixed missingness mechanism
#'
#' Convenience wrapper around simulate_sslfmm(mechanism = "mixed").
#'
#' @inheritParams simulate_sslfmm
#' @return The same documented list structure as simulate_sslfmm().
#' @export
simulate_mixed_missingness <- function(n = 500L, pi, mu, sigma,
                                       alpha = 0.10, xi0 = NULL, xi1 = 3,
                                       mar_rate = 0.30, seed = NULL,
                                       entropy_floor = 1e-300) {
  simulate_sslfmm(n = n, pi = pi, mu = mu, sigma = sigma,
                  mechanism = "mixed", alpha = alpha, xi0 = xi0,
                  xi1 = xi1, mar_rate = mar_rate, seed = seed,
                  entropy_floor = entropy_floor)
}
