#' Negative Log-Likelihood for Semi-supervised FMM with a Mixed-Missingness Mechanism
#'
#' Computes the negative log-likelihood for a semi-supervised Gaussian mixture
#' model under a mixed missingness mechanism (MCAR + entropy-based MAR).
#' Assumes a ** covariance matrix** \eqn{\Sigma} across all mixture
#' components.
#'
#' @param theta Numeric vector of packed model parameters to be unpacked by `unpack_fn`.
#' @param Y Numeric matrix of observations (n x p).
#' @param m_j Integer or logical vector of length n indicating missingness:
#'   0 for observed (labeled block), 1 for unlabeled/missingness block.
#' @param Z Integer vector of length n with class labels for labeled samples
#'   (1..g); use `NA` for unlabeled rows.
#' @param d2_yj Numeric vector of length n with the entropy-like score used in the
#'   MAR mechanism (e.g., posterior entropy or any scalar proxy).
#' @param xi Numeric length-2 vector \code{c(xi0, xi1)} for the logistic MAR model
#'   \eqn{q_j = logistic(xi0 + xi1 * d2_yj)}.
#' @param alpha_k Numeric scalar in (0,1), the MCAR mixing proportion in the
#'   missingness mechanism.
#' @param unpack_fn Function that takes `theta` and returns a list with elements:
#'   \describe{
#'     \item{\code{pi}}{Numeric vector of length g with mixture weights.}
#'     \item{\code{mu}}{List of length g; each element is a numeric mean vector (length p).}
#'     \item{\code{sigma}}{Shared covariance matrix (p x p).}
#'   }
#'
#' @return A single numeric value: the negative log-likelihood.
#'
#' @details
#' The total log-likelihood is composed of three parts:
#' \enumerate{
#'   \item Labeled samples (\eqn{m_j=0}) with observed class labels \eqn{Z_j}.
#'   \item Unlabeled samples attributed to MCAR with probability mass \eqn{m_{1j}}.
#'   \item Unlabeled samples attributed to MAR with probability mass \eqn{m_{2j}}.
#' }
#' The MAR probability for each sample is \eqn{q_j = \mathrm{logistic}(xi0 + xi1\, d2\_yj)}.
#' Internally, the function uses a numerically stable \code{logSumExp}.
#'
#' @note
#' This implementation is for the **equal covariance** case (shared \eqn{\Sigma}).
#'
#' @examples
#' \dontrun{
#' # Minimal example (illustrative only):
#' library(mvtnorm)
#' set.seed(1)
#' n <- 20; p <- 2; g <- 2
#' Y <- matrix(rnorm(n*p), n, p)
#' Z <- sample(c(1:g, rep(NA, n - g)), n, replace = TRUE)
#' m_j <- ifelse(is.na(Z), 1L, 0L)
#' d2_yj <- runif(n)
#' xi <- c(-1, 2)
#' alpha_k <- 0.4
#' unpack_fn <- function(theta) {
#'   list(pi = c(0.6, 0.4),
#'        mu = list(c(0,0), c(1,1)),
#'        sigma = diag(p))
#' }
#' theta <- numeric(1) # not used in this toy unpack_fn
#' neg_loglik(theta, Y, m_j, Z, d2_yj, xi, alpha_k, unpack_fn)
#' }
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats plogis
#' @export
neg_loglik <- function(theta, Y, m_j, Z, d2_yj, xi, alpha_k, unpack_fn) {
  eps <- 1e-10

  # --- helpers ---
  logSumExp <- function(x) {
    m <- max(x)
    if (!is.finite(m)) return(m)
    m + log(sum(exp(x - m)))
  }
  is_sym_pd <- function(S) {
    if (!is.matrix(S)) return(FALSE)
    if (nrow(S) != ncol(S)) return(FALSE)
    if (!isTRUE(all.equal(S, t(S), tolerance = 1e-8))) return(FALSE)
    ok <- TRUE
    tryCatch({ chol(S) }, error = function(e) { ok <<- FALSE })
    ok
  }

  # --- unpack parameters ---
  param <- unpack_fn(theta)
  pi_k  <- pmax(as.numeric(param$pi), eps)
  mu_k  <- param$mu                 # list of g mean vectors
  Sigma_raw <- param$sigma          # either a matrix (shared) or a list of length g
  alpha <- as.numeric(alpha_k)
  xi0   <- xi[1]
  xi1   <- xi[2]

  g <- length(pi_k)
  n <- nrow(Y)
  p <- ncol(Y)

  # --- sanity checks (lightweight) ---
  if (length(m_j) != n) stop("m_j must have length nrow(Y).")
  if (length(d2_yj) != n) stop("d2_yj must have length nrow(Y).")
  if (length(xi) != 2L) stop("xi must be length-2: c(xi0, xi1).")
  if (!is.list(mu_k) || length(mu_k) != g) stop("param$mu must be a list of length g.")
  if (any(pi_k <= 0)) stop("All mixture weights must be positive.")
  if (abs(sum(pi_k) - 1) > 1e-6) warning("Mixture weights do not sum to 1 (tolerance 1e-6).")
  if (!(alpha > 0 && alpha < 1)) stop("alpha_k must be in (0,1).")
  if (ncol(Y) != p) stop("Y must have p columns.")
  if (!all(vapply(mu_k, length, 1L) == p)) stop("Each mean vector must have length p.")

  # --- detect covariance mode & normalize to list-of-g matrices ---
  if (is.matrix(Sigma_raw)) {
    # shared covariance
    if (!is_sym_pd(Sigma_raw)) stop("Shared covariance matrix must be symmetric positive definite.")
    Sigma_list <- rep(list(Sigma_raw), g)
    shared_cov <- TRUE
  } else if (is.list(Sigma_raw) && length(Sigma_raw) == g) {
    # heteroscedastic
    if (!all(vapply(Sigma_raw, function(S) nrow(S) == p && ncol(S) == p, logical(1))))
      stop("Each component covariance must be a p x p matrix.")
    if (!all(vapply(Sigma_raw, is_sym_pd, logical(1))))
      stop("All component covariances must be symmetric positive definite.")
    Sigma_list <- Sigma_raw
    shared_cov <- FALSE
  } else {
    stop("param$sigma must be either a p times  p matrix (shared) or a list of length g of p  times  p matrices.")
  }

  # --- MAR probabilities ---
  q_j   <- stats::plogis(xi0 + xi1 * d2_yj)
  denom <- alpha + (1 - alpha) * q_j
  denom <- pmax(denom, eps)
  m1j <- alpha / denom
  m2j <- (1 - alpha) * q_j / denom

  # Only applicable where m_j == 1
  m1j_k <- ifelse(m_j == 1, m1j, 0)
  m2j_k <- ifelse(m_j == 1, m2j, 0)

  # --- 1) Labeled part (m_j == 0 and Z observed) ---
  idx_obs <- which(m_j == 0 & !is.na(Z))
  L_labelled <- 0
  if (length(idx_obs) > 0) {
    L_labelled <- sum(vapply(idx_obs, function(j) {
      k <- Z[j]
      if (is.na(k) || k < 1 || k > g) return(0)
      mv_log <- mvtnorm::dmvnorm(Y[j, ], mean = mu_k[[k]], sigma = Sigma_list[[k]], log = TRUE)
      # observed block prob: 1 - alpha - (1 - alpha) * q_j
      log(pi_k[k] + eps) + mv_log + log(1 - alpha - (1 - alpha) * q_j[j] + eps)
    }, numeric(1)))
  }

  # --- Common mixture log-density for a row ---
  mix_logdens <- function(y_row) {
    logSumExp(vapply(seq_len(g), function(k) {
      mv_log <- mvtnorm::dmvnorm(y_row, mean = mu_k[[k]], sigma = Sigma_list[[k]], log = TRUE)
      log(pi_k[k] + eps) + mv_log
    }, numeric(1)))
  }

  # --- 2) MCAR part (m_j == 1) & 3) MAR part (m_j == 1) ---
  idx_miss <- which(m_j == 1)
  L_mcar <- 0
  L_mar  <- 0
  if (length(idx_miss) > 0) {
    L_mcar <- sum(vapply(idx_miss, function(j) {
      log_mix <- mix_logdens(Y[j, ])
      m1j_k[j] * (log_mix + log(alpha + eps))
    }, numeric(1)))
    L_mar <- sum(vapply(idx_miss, function(j) {
      log_mix <- mix_logdens(Y[j, ])
      m2j_k[j] * (log_mix + log((1 - alpha) * q_j[j] + eps))
    }, numeric(1)))
  }

  L_total <- L_labelled + L_mcar + L_mar
  -L_total
}
