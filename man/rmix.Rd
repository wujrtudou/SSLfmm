\name{rmix}
\alias{rmix}
\title{Simulate from a Gaussian Finite Mixture}
\usage{
rmix(n = 200L, pi, mu, sigma, seed_number = NULL)
}
\arguments{
\item{n}{Positive integer number of observations.}
\item{pi}{Positive finite mixing proportions of length \eqn{g} that sum to one.}
\item{mu}{Numeric \eqn{p \times g} matrix of component means; each column is one component.}
\item{sigma}{Either one shared numeric \eqn{p \times p} symmetric positive-definite covariance matrix or a numeric \eqn{p \times p \times g} array containing one symmetric positive-definite covariance matrix per component. For \eqn{p=1}, a length-one scalar is also accepted as a shared variance.}
\item{seed_number}{Optional random seed.}
}
\description{
Generates observations from a Gaussian finite mixture without introducing label missingness.
}
\value{
A data.frame with feature columns \code{x1}, ..., \code{xp} followed by integer component label \code{truth}.
}
\examples{
mu <- matrix(c(-1, 1), nrow = 1, ncol = 2)
dat <- rmix(
  n = 20,
  pi = c(0.5, 0.5),
  mu = mu,
  sigma = matrix(1, 1, 1),
  seed_number = 1
)
head(dat)
}
