\name{SSLfmm-package}
\alias{SSLfmm-package}
\alias{SSLfmm}
\title{SSLfmm: Semi-Supervised Finite Mixture Models with Mixed Label Missingness}
\description{
Semi-supervised Gaussian finite mixture models for partially observed class labels under complete-case, MCAR, entropy-dependent MAR, and mixed MCAR/MAR formulations. The mixed formulation allows the missingness source to be observed or latent. The package supports equal or component-specific covariance matrices, simulation, prediction, classification performance assessment, and entropy-based diagnostics.
}
\details{
The intended user-facing functions are \code{fit_sslfmm()}, \code{initialize_sslfmm()}, \code{rmix()}, \code{simulate_sslfmm()}, \code{simulate_mixed_missingness()}, \code{classification_performance()}, together with standard S3 \code{predict()}, \code{print()}, and \code{summary()} methods for fitted objects. The package also includes the semi-synthetic \code{blood_transfusion} data used in the software-paper application. Low-level likelihood, parameter-packing, Cholesky, and entropy helper functions are internal and are not exported.
}
