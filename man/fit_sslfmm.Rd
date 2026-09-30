\name{fit_sslfmm}
\alias{fit_sslfmm}
\title{Fit a Semi-Supervised Gaussian Finite Mixture Model}
\usage{
fit_sslfmm(
  x, y, g = NULL,
  method = c("mixed", "mar", "mcar", "cc"),
  covariance_type = c("equal", "unequal"),
  indicator = c("unknown", "known"),
  missing_source = NULL, init = NULL, n_starts = 5L,
  control = list(), entropy_floor = 1e-300,
  alpha_starts = c(0, 0.05, 0.15, 0.30),
  alpha_upper = 1 - 1e-8, ridge = 1e-4, seed = NULL
)
}
\arguments{
\item{x}{Numeric feature matrix or data.frame with observations in rows and features in columns. Feature values must be fully observed.}
\item{y}{Class labels with one entry per row of \code{x}; use \code{NA} for missing labels. Numeric labels must be positive integers in \code{1:g}; factor and character labels are also accepted.}
\item{g}{Positive integer number of mixture components. It is inferred from observed labels when possible.}
\item{method}{One of \code{"cc"}, \code{"mcar"}, \code{"mar"}, or \code{"mixed"}.}
\item{covariance_type}{\code{"equal"} for one covariance matrix shared by all components or \code{"unequal"} for one covariance matrix per component.}
\item{indicator}{For \code{method = "mixed"}, whether the MCAR/MAR source indicator is \code{"unknown"} or \code{"known"}.}
\item{missing_source}{Required only for \code{method = "mixed", indicator = "known"}. Supply \code{"mcar"}/\code{"mar"} on missing-label rows, or a logical/0-1 MCAR-source indicator. Values on labelled rows are ignored.}
\item{init}{Optional output from \code{initialize_sslfmm()} or a list containing \code{pi}, \code{mu}, and \code{sigma}. The fitting orientation for \code{mu} is \eqn{g \times p} (one component per row). A \eqn{p \times g} matrix, matching the simulation input orientation, is also accepted when \eqn{p \ne g} and therefore unambiguous.}
\item{n_starts}{Positive integer number of fitting starts.}
\item{control}{A list controlling numerical fitting.}
\item{entropy_floor}{Single positive finite lower bound used before taking the logarithm of posterior entropy.}
\item{alpha_starts}{Finite starting values for the latent MCAR-channel probability in the latent mixed model. Every value must lie between zero and \code{alpha_upper}.}
\item{alpha_upper}{Single finite upper bound for alpha, strictly between zero and one.}
\item{ridge}{Single positive finite covariance regularization value used for initialization and complete-case fitting.}
\item{seed}{Optional random seed for initialization and perturbed starts.}
}
\description{
Fits semi-supervised Gaussian finite-mixture classifiers under complete-case, MCAR, entropy-dependent MAR, or mixed MCAR/MAR label-missingness formulations. For mixed models, the MCAR/MAR source may be known or unknown.
}
\details{
For \code{method = "mcar"}, unlabelled observations contribute the Gaussian-mixture marginal density. For \code{method = "mar"}, label missingness follows a logistic model in log posterior entropy. For \code{method = "mixed"}, the missing-label probability combines an MCAR channel and the entropy-dependent MAR channel. With a known source indicator, alpha is estimated from the recorded MCAR-source assignments; with an unknown source indicator, alpha is estimated jointly with the other model parameters.

Entropy-dependent \code{"mar"} and \code{"mixed"} fitting requires \code{g >= 2}. With one component posterior entropy is constant and the entropy slope is not identifiable. Simulation with \code{g = 1} remains supported as a boundary/API check.
}
\value{
An object of class \code{"SSLfmm"}. Important components include:
\itemize{
\item \code{pi}: fitted mixing proportions, named by class;
\item \code{mu}: fitted \eqn{g \times p} mean matrix, one component per row;
\item \code{sigma}: one \eqn{p \times p} matrix for equal covariance or a \eqn{p \times p \times g} array for unequal covariance;
\item \code{xi} and \code{alpha}: fitted missingness parameters when applicable;
\item \code{posterior}: fitted class-posterior matrix with class-labelled columns;
\item \code{en} and \code{entropy}: identical fitted posterior-entropy vectors;
\item \code{missing} and \code{observed_missing}: identical observed label-missingness indicators;
\item \code{latent_missing}: recorded MCAR-source indicator for known-source mixed fitting, or \code{NA} when that source is unknown;
\item \code{latent_missing_probability}: posterior probability of the MCAR source for each row in the unknown-source mixed model;
\item \code{loglik}, \code{convergence}, \code{optimizer}, \code{start_summary}, and \code{start_diagnostics}: fit and start diagnostics.
}
}
\examples{
set.seed(1)
mu <- matrix(c(-2, 2), nrow = 1, ncol = 2)
sim <- simulate_sslfmm(
  n = 60,
  pi = c(0.5, 0.5),
  mu = mu,
  sigma = matrix(1, 1, 1),
  mechanism = "mcar",
  missing_rate = 0.20,
  seed = 1
)
x <- as.matrix(sim$data["x1"])

# A complete-case fit is deterministic and demonstrates the main interface.
fit <- fit_sslfmm(x, sim$data$label, g = 2,
                  method = "cc", covariance_type = "equal")
fit
predict(fit, x[1:5, , drop = FALSE], type = "posterior")
}
