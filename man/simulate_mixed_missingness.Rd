\name{simulate_mixed_missingness}
\alias{simulate_mixed_missingness}
\title{Simulate Semi-Supervised Data Under Mixed Label Missingness}
\usage{
simulate_mixed_missingness(
  n = 500L, pi, mu, sigma,
  alpha = 0.10, xi0 = NULL, xi1 = 3,
  mar_rate = 0.30, seed = NULL,
  entropy_floor = 1e-300
)
}
\arguments{
\item{n}{Positive integer number of observations.}
\item{pi}{Positive finite mixing proportions of length \eqn{g} that sum to one.}
\item{mu}{Numeric \eqn{p \times g} matrix of component means.}
\item{sigma}{Shared \eqn{p \times p} covariance matrix or \eqn{p \times p \times g} component-specific covariance array.}
\item{alpha}{MCAR-channel probability in \eqn{[0,1]}.}
\item{xi0}{MAR intercept; when \code{NULL}, calibrated to \code{mar_rate}.}
\item{xi1}{Positive finite entropy slope.}
\item{mar_rate}{Target mean MAR probability when \code{xi0 = NULL}.}
\item{seed}{Optional random seed.}
\item{entropy_floor}{Positive finite entropy floor.}
}
\description{
Convenience wrapper for \code{simulate_sslfmm(..., mechanism = "mixed")}. It uses the same mixed missingness mechanism and the same stable return-value contract.
}
\value{
The same five-component list documented for \code{simulate_sslfmm()}: \code{data}, \code{true_setup}, \code{groups}, \code{probs}, and \code{raw}.
}
\examples{
mu <- matrix(c(-1, 1), nrow = 1, ncol = 2)
sim <- simulate_mixed_missingness(
  n = 40,
  pi = c(0.5, 0.5),
  mu = mu,
  sigma = matrix(1, 1, 1),
  alpha = 0.10,
  mar_rate = 0.25,
  seed = 1
)
stopifnot(
  is.list(sim),
  is.data.frame(sim$data),
  identical(sim$data$missing, sim$data$observed_missing)
)
}
\seealso{\code{simulate_sslfmm()}, \code{rmix()}}
