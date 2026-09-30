\name{predict.SSLfmm}
\alias{predict.SSLfmm}
\alias{print.SSLfmm}
\alias{summary.SSLfmm}
\title{Methods for SSLfmm Fits}
\usage{
\method{predict}{SSLfmm}(object, newdata,
  type = c("class", "posterior", "entropy", "all"), ...)
\method{print}{SSLfmm}(x, ...)
\method{summary}{SSLfmm}(object, ...)
}
\arguments{
\item{object,x}{An object of class \code{"SSLfmm"}.}
\item{newdata}{Numeric feature matrix or data.frame with the same number of feature columns as the fitted model.}
\item{type}{Prediction output: \code{"class"}, \code{"posterior"}, \code{"entropy"}, or \code{"all"}.}
\item{...}{Unused.}
}
\description{
Standard prediction, printing, and summary methods for fitted SSLfmm objects.
}
\value{
For \code{predict()}, \code{"class"} returns predicted labels in the same broad label type used for fitting; \code{"posterior"} returns an \eqn{n \times g} numeric matrix with class-labelled columns; \code{"entropy"} returns a numeric vector; and \code{"all"} returns a list with components \code{class}, \code{posterior}, and \code{entropy}. \code{print()} returns the fit invisibly. \code{summary()} returns a compact fitted-parameter list.
}
\examples{
mu <- matrix(c(-2, 2), nrow = 1, ncol = 2)
s <- simulate_sslfmm(50, c(0.5, 0.5), mu, matrix(1, 1, 1),
                      mechanism = "mcar", missing_rate = 0.2, seed = 1)
x <- as.matrix(s$data["x1"])
f <- fit_sslfmm(x, s$data$label, g = 2, method = "cc")
predict(f, x[1:5, , drop = FALSE], type = "posterior")
predict(f, x[1:5, , drop = FALSE], type = "all")
}
