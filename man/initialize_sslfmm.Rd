\name{initialize_sslfmm}
\alias{initialize_sslfmm}
\title{Initialize SSLfmm Model Parameters}
\usage{
initialize_sslfmm(
  x, y, g = NULL,
  covariance_type = c("equal", "unequal"),
  ridge = 1e-4, seed = NULL
)
}
\arguments{
\item{x}{Numeric feature matrix or data.frame with observations in rows and features in columns.}
\item{y}{Class labels with one entry per row of \code{x}; use \code{NA} for missing labels.}
\item{g}{Positive integer number of mixture components. Inferred from observed labels when possible.}
\item{covariance_type}{\code{"equal"} for a shared covariance matrix or \code{"unequal"} for component-specific covariance matrices.}
\item{ridge}{Single positive finite diagonal regularization value for covariance starts.}
\item{seed}{Optional seed used only for k-means fallback initialization.}
}
\description{
Constructs stable starting values for SSLfmm model fitting. Observed class labels anchor component identities; k-means is used only as a fallback for components without adequate labelled observations.
}
\value{
An object of class \code{"SSLfmm_init"}, implemented as a list containing \code{pi}, \code{mu}, \code{sigma}, \code{g}, \code{p}, \code{covariance_type}, and \code{label_levels}. The mean matrix \code{mu} is \eqn{g \times p}, with one component per row. For equal covariance, \code{sigma} is \eqn{p \times p}; for unequal covariance, it is \eqn{p \times p \times g}.
}
\examples{
x <- rbind(matrix(rnorm(20, -1), 10, 2),
           matrix(rnorm(20, 1), 10, 2))
y <- c(rep(1, 7), rep(NA, 3), rep(2, 7), rep(NA, 3))
init <- initialize_sslfmm(x, y, g = 2,
                           covariance_type = "unequal", seed = 1)
init$pi
init$mu
dim(init$sigma)
}
