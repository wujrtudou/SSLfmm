#' Simulate a Gaussian mixture dataset with mixed missingness (MAR + MCAR)
#'
#' @param n Integer; sample size.
#' @param pi Numeric vector; mixing proportions (sum to 1).
#' @param mu Matrix (p x K); component means, columns = components.
#' @param sigma Array (p x p x K); component covariance matrices.
#' @param xi0 Numeric; MAR logit intercept.
#' @param xi1 Numeric; MAR logit slope on entropy.
#' @param alpha Numeric in \code{[0,1]}; MCAR rate applied within both MAR and observed groups.
#' @param seed_id Integer; seed passed to rmix() (your generator).
#'
#' @details
#' Requires user-provided functions:
#'   - rmix(n, pi, mu, sigma, seed_number)
#'   - get_entropy(dat, n, p, g, paralist)
#'
#' Missingness mechanism codes:
#'   - 0 = fully observed
#'   - 1 = MCAR
#'   - 2 = MAR (entropy-based)
#'
#' @return A list with:
#'   - data: data.frame with columns x1..xp, en, missing, label, truth
#'   - true_setup: list(pi, mu, sigma)
#'   - groups: list(mar_group, obs_group, mcar_in_mar, mcar_in_obs)
#'   - probs: vector prob_mar
#'   - raw: original rmix output `dat` augmented with en and labels
#'
#' @export
simulate_mixed_missingness <- function(
    n = 500,
    pi,
    mu,
    sigma,
    xi0 = 2,
    xi1 = 3,
    alpha = 0.1,
    seed_id = 123
) {
  # --- Basic checks ---
  K <- length(pi)
  p <- nrow(mu)
  stopifnot(abs(sum(pi) - 1) < 1e-8, p >= 1, dim(sigma)[1] == p, dim(sigma)[2] == p, dim(sigma)[3] == K)
  if (alpha < 0 || alpha > 1) stop("alpha must be in [0,1].")
  
  true_setup <- list(pi = pi, mu = mu, sigma = sigma)
  
  # --- Generate mixture and entropy ---
  dat <- rmix(n = n, pi = pi, mu = mu, sigma = sigma, seed_number = seed_id)
  dat$en <- get_entropy(dat = dat$Y, n = n, p = p, g = K, paralist = true_setup)
  
  # --- Entropy-based MAR probability and assignment ---
  logit_p <- xi0 + xi1 * log(dat$en)
  prob_mar <- exp(logit_p) / (1 + exp(logit_p))
  
  is_mar_flag <- rbinom(n, 1, prob_mar) == 1
  mar_group <- which(is_mar_flag)
  obs_group <- which(!is_mar_flag)
  
  # --- Inject MCAR within both MAR and observed groups (rate = alpha) ---

  mcar_in_mar <- if (length(mar_group)) sample(mar_group, size = floor(alpha * length(mar_group))) else integer(0)

  mcar_in_obs <- if (length(obs_group)) sample(obs_group, size = floor(alpha * length(obs_group))) else integer(0)
  
  # --- Labels and mechanism codes ---
  final_label <- dat$clust
  final_label[mar_group] <- NA
  final_label[mcar_in_mar] <- NA
  final_label[mcar_in_obs] <- NA
  
  missing_mechanism <- integer(n) # 0 fully observed
  missing_mechanism[c(mcar_in_mar, mcar_in_obs)] <- 1L
  missing_mechanism[setdiff(mar_group, mcar_in_mar)] <- 2L
  
  # --- Assemble data.frame ---
  x_names <- paste0("y", 1:p)
  data <- cbind(dat$Y, en = dat$en, missing = missing_mechanism,
                label = final_label, truth = dat$clust)
  colnames(data) <- c(x_names, "en", "missing", "label", "truth")
  data <-as.data.frame(data)
  
  return(data = data)
}
