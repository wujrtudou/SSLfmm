\name{simulate_sslfmm}
\alias{simulate_sslfmm}
\title{Simulate Semi-Supervised Finite-Mixture Data}
\description{
Simulates semi-supervised Gaussian finite-mixture data with class-label
missingness under complete-case, MCAR, entropy-dependent MAR, or mixed
MCAR/MAR mechanisms.
}
\usage{
simulate_sslfmm(
  n = 500L, pi, mu, sigma,
  mechanism = c("mixed", "mar", "mcar", "cc"),
  alpha = 0.10, missing_rate = NULL,
  xi0 = NULL, xi1 = 3, mar_rate = 0.30,
  seed = NULL, entropy_floor = 1e-300
)
}
\arguments{
\item{n}{Positive integer number of observations.}
\item{pi}{Positive finite mixing proportions of length \eqn{g} that sum to one.}
\item{mu}{Numeric \eqn{p \times g} matrix of component means; each column is one component.}
\item{sigma}{Either one shared numeric \eqn{p \times p} symmetric positive-definite covariance matrix or a numeric \eqn{p \times p \times g} array with one covariance matrix per component. For \eqn{p=1}, a length-one scalar is also accepted as a shared variance.}
\item{mechanism}{One of \code{"cc"}, \code{"mcar"}, \code{"mar"}, or \code{"mixed"}.}
\item{alpha}{MCAR-channel probability for \code{mechanism = "mixed"}. For \code{mechanism = "mcar"}, this is the missing-label rate when \code{missing_rate} is \code{NULL}.}
\item{missing_rate}{Optional MCAR missing-label rate used only for \code{mechanism = "mcar"}.}
\item{xi0}{MAR intercept. If \code{NULL}, it is calibrated on the generated feature sample so that the mean MAR probability equals \code{mar_rate}.}
\item{xi1}{Positive finite slope multiplying log posterior entropy for MAR or mixed simulation.}
\item{mar_rate}{Target average MAR probability when \code{xi0 = NULL}; must lie strictly between zero and one.}
\item{seed}{Optional random seed.}
\item{entropy_floor}{Positive finite lower bound used before taking the logarithm of posterior entropy.}
}
\details{
For mixed missingness, let \eqn{q_i} denote the entropy-dependent MAR probability and let the MCAR trigger have probability \eqn{\alpha}. The observed missing-label probability is
\deqn{\alpha + (1-\alpha)q_i.}
The simulator records both the true MCAR trigger and a latent potential-MAR trigger. MCAR takes precedence when both occur. This is distributionally equivalent to generating the MAR trigger only among observations not selected by MCAR, while also making the documented group decomposition reproducible.

When \eqn{g=1}, posterior entropy is constant. Simulation remains supported for API/boundary testing, and automatic \code{xi0} calibration uses the resulting constant-probability analytic solution. Entropy-dependent MAR slope estimation with \eqn{g=1} is not identifiable and is therefore not supported by \code{fit_sslfmm()}.
}
\value{
Always a list with exactly the top-level components \code{data}, \code{true_setup}, \code{groups}, \code{probs}, and \code{raw}.

\itemize{
\item \code{data}: a data.frame whose leading columns are, in order, \code{x1}, ..., \code{xp}, \code{en}, \code{missing}, \code{label}, and \code{truth}. It then contains \code{observed_missing}, \code{latent_missing}, \code{missing_source}, \code{prob_mar}, and \code{entropy}. \code{en} and \code{entropy} are identical; \code{missing} and \code{observed_missing} are identical. \code{latent_missing} is the true MCAR-channel trigger in the simulation.
\item \code{true_setup}: a list beginning with \code{pi}, \code{mu}, and \code{sigma}, followed by covariance/missingness metadata and the true missingness parameters when applicable.
\item \code{groups}: a list beginning with \code{mar_group}, \code{obs_group}, \code{mcar_in_mar}, and \code{mcar_in_obs}. \code{mar_group} is the latent potential-MAR stratum, \code{obs_group} its complement, and the two \code{mcar_in_*} entries are intersections with the MCAR trigger. Additional entries \code{observed}, \code{mcar}, \code{mar}, and \code{missing} give directly useful row indices.
\item \code{probs}: numeric vector of MAR probabilities \eqn{q_i}; it is identical to \code{data$prob_mar}.
\item \code{raw}: the original \code{rmix()} data.frame augmented with \code{en}, \code{missing}, \code{label}, the explicit missingness indicators/source, and MAR probabilities.
}
}
\examples{
mu <- matrix(c(-1.5, 1.5), nrow = 1, ncol = 2)
Sigma <- matrix(1, 1, 1)

s_cc <- simulate_sslfmm(30, c(0.5, 0.5), mu, Sigma,
                         mechanism = "cc", seed = 1)
s_mcar <- simulate_sslfmm(30, c(0.5, 0.5), mu, Sigma,
                           mechanism = "mcar", missing_rate = 0.25, seed = 2)
s_mar <- simulate_sslfmm(30, c(0.5, 0.5), mu, Sigma,
                          mechanism = "mar", mar_rate = 0.25, seed = 3)
s_mixed <- simulate_sslfmm(30, c(0.5, 0.5), mu, Sigma,
                            mechanism = "mixed", alpha = 0.10,
                            mar_rate = 0.25, seed = 4)

names(s_mixed)
head(s_mixed$data)
}
