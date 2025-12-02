#' Compute Theoretical Bayes' Error for a Binary Gaussian Mixture
#'
#' @description
#' Computes the Bayes classification error rate for a two-component Gaussian
#' mixture given \code{mu_hat}, \code{Sigma_hat}, and \code{pi_hat}. If
#' \code{mu_hat} is supplied as a list of length 2, it is converted to a
#' \eqn{p \times 2} matrix internally.
#'
#' @param mu_hat Either a numeric matrix of size \eqn{p \times 2} whose columns
#'   are component means, or a list of two numeric vectors.
#' @param Sigma_hat Numeric \eqn{p \times p} covariance matrix shared
#'   across components.
#' @param pi_hat Numeric vector of length 2 with mixing proportions
#'   \eqn{(\pi_1, \pi_2)} that are non-negative and sum to 1.
#'
#' @return
#' A numeric scalar giving the theoretical Bayes classification error rate.
#'
#' @details
#' The linear discriminant is
#' \deqn{
#'   \beta_1 = \Sigma^{-1}(\mu_1 - \mu_2), \qquad
#'   \beta_0 = -\frac{1}{2}(\mu_1 + \mu_2)^\top \Sigma^{-1}(\mu_1 - \mu_2)
#'            + \log(\pi_1 / \pi_2)
#' }
#' and the Bayes error is
#' \deqn{
#'   \mathrm{Err}
#'   = \sum_{k = 1}^2 \pi_k \,
#'     \Phi\left(
#'       \frac{(-1)^k \{\beta_0 + \beta_1^\top \mu_k\}}
#'            {\|\beta_1\|}
#'     \right),
#' }
#' where \eqn{\Phi} is the standard normal cdf.
#'
#' @examples
#' mu_hat <- matrix(c(1, 0, -1, 0), nrow = 2)  # columns are mu1, mu2
#' Sigma_hat <- diag(2)
#' pi_hat <- c(0.5, 0.5)
#' error_beta_classification(mu_hat, Sigma_hat, pi_hat)
#'
#' @export
error_beta_classification <- function(mu_hat, Sigma_hat, pi_hat) {
  # Accept list input for mu_hat and convert to matrix p x 2
  if (is.list(mu_hat)) {
    if (length(mu_hat) != 2L) stop("mu_hat list must have length 2.")
    mu_hat <- do.call(cbind, lapply(mu_hat, function(v) as.numeric(v)))
  }
  if (!is.matrix(mu_hat)) stop("mu_hat must be a matrix (p x 2) or a list of two vectors.")
  if (ncol(mu_hat) != 2L) stop("mu_hat must have 2 columns (two components).")
  p <- nrow(mu_hat)

  # Basic checks for Sigma_hat and pi_hat
  if (!is.matrix(Sigma_hat) || any(dim(Sigma_hat) != c(p, p))) {
    stop("Sigma_hat must be a p x p matrix matching nrow(mu_hat).")
  }
  if (length(pi_hat) != 2L) stop("pi_hat must be length 2.")
  if (any(pi_hat < 0)) stop("pi_hat must be non-negative.")
  s <- sum(pi_hat)
  if (!is.finite(s) || abs(s - 1) > 1e-8) {
    pi_hat <- pi_hat / s
    if (abs(sum(pi_hat) - 1) > 1e-8) stop("pi_hat must sum to 1.")
  }

  # Ensure Sigma_hat is symmetric and positive definite
  Sigma_hat <- 0.5 * (Sigma_hat + t(Sigma_hat))
  R <- tryCatch(chol(Sigma_hat), error = function(e) NULL)
  if (is.null(R)) stop("Sigma_hat must be symmetric positive definite.")

  # Efficient inverse via Cholesky
  Sigma_inv <- chol2inv(R)

  delta_mu <- mu_hat[, 1] - mu_hat[, 2]            # mu1 - mu2
  beta1 <- as.numeric(Sigma_inv %*% delta_mu)      # p-vector
  beta0 <- as.numeric(
    -0.5 * crossprod(mu_hat[, 1] + mu_hat[, 2], Sigma_inv %*% delta_mu) +
      log(pi_hat[1] / pi_hat[2])
  )

  norm_beta1 <- sqrt(sum(beta1 * beta1))
  if (norm_beta1 == 0) {
    # Means identical under shared Sigma -> error is min(pi)
    return(min(pi_hat))
  }

  # Compute error = sum_k pi_k * Phi( (-1)^k * (beta0 + beta1^T mu_k)/||beta1|| )
  v1 <- (beta0 + sum(beta1 * mu_hat[, 1])) / norm_beta1
  v2 <- (-beta0 - sum(beta1 * mu_hat[, 2])) / norm_beta1
  err <- pi_hat[1] * stats::pnorm(v1) + pi_hat[2] * stats::pnorm(v2)

  as.numeric(err)
}
