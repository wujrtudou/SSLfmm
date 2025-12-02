#' Complete-data warm-up initialization for semi-supervised FMM with mixed-missingness mechanisms
#'
#' Uses both labeled and unlabeled subsets of the data to obtain quick initial
#' estimates for mixture parameters and missingness mechanism parameters
#' (\code{alpha}, \code{xi}) via a warm-up EM procedure.
#'
#' @param data A data frame containing:
#'   \itemize{
#'     \item The first \code{p} columns: numeric variables used in the FMM.
#'     \item A column \code{missing}: indicator (0 = labeled, 1 = unlabeled/missing).
#'     \item A column \code{obers}: class labels for labeled rows (\code{1:g}); \code{NA} for unlabeled.
#'   }
#' @param g Integer, number of mixture components (default \code{2}).
#' @param ncov Integer, covariance structure: \code{1} = shared (equal), \code{2} = class-specific (unequal).
#' @param alpha_init Numeric in (0,1), initial MCAR proportion (default \code{0.01}).
#' @param warm_up_iter Integer, number of warm-up EM iterations (default \code{200}).
#' @param tol Convergence tolerance on \code{alpha} (default \code{1e-6}).
#'
#' @return A list with initial values from \code{\link{EM_FMM_SemiSupervised_Initial}}:
#' \itemize{
#'   \item \code{pi} - mixture weights.
#'   \item \code{mu} - list of component mean vectors.
#'   \item \code{Sigma} - covariance matrix/matrices.
#'   \item \code{alpha} - MCAR proportion.
#'   \item \code{xi} - logistic regression coefficients for MAR mechanism.
#' }
#'
#' @details
#' - This function first calls \code{\link{initialestimate}} to get initial \eqn{\pi}, \eqn{\mu}, \eqn{\Sigma}.
#' - Then it calls \code{\link{EM_FMM_SemiSupervised_Initial}} with these values for a short warm-up run.
#' - Covariance structure (\code{equal} vs. \code{unequal}) is determined by \code{ncov}.
#'
#' @examples
#' \dontrun{
#' # Suppose df contains p features, a 'missing' column, and 'obers' labels
#' res_init <- EM_FMM_SemiSupervised_Complete_Initial(df, g = 3, ncov = 2)
#' }
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom matrixStats rowLogSumExps
#' @export
EM_FMM_SemiSupervised_Complete_Initial <- function(data, g = 2, ncov = 1,
                                                   alpha_init = 0.01,
                                                   warm_up_iter = 200,
                                                   tol = 1e-6) {
  
  p <- ncol(data) - 2L  # assumes last 2 cols are 'missing' and 'obers'
  
  Y_labelled <- data[data$missing == 0, c(1:p)]
  Y_unlabelled <- data[data$missing != 0, c(1:p)]
  Z_labelled <- data[data$missing == 0, "z"]
  
  # Initial values from labeled subset
  inits <- initialestimate(dat = data[, 1:p], zm = data$z, g = g, ncov = ncov)
  pi_init <- inits$pi
  mu_init <- lapply(1:ncol(inits$mu), function(i) inits$mu[, i])
  Sigma_init <- inits$sigma
  
  # Warm-up EM to refine alpha and xi
  init_res <- EM_FMM_SemiSupervised_Initial(
    Y_labelled, Z_labelled, Y_unlabelled,
    g = g, pi_init = pi_init, mu_init = mu_init, Sigma_init = Sigma_init,
    alpha_init = alpha_init, warm_up_iter = warm_up_iter, tol = tol
  )
  
  return(list(
    pi = init_res$pi,
    mu = init_res$mu,
    Sigma = init_res$Sigma,
    alpha = init_res$alpha,
    xi = init_res$xi
  ))
}
