\name{blood_transfusion}
\alias{blood_transfusion}
\docType{data}
\title{Semi-Synthetic Blood Transfusion Data}
\description{
A semi-synthetic partially observed version of the Blood Transfusion
classification data used to illustrate the missing-label models in
\pkg{SSLfmm}.
}
\usage{
data("blood_transfusion")
}
\format{
A data frame with 748 observations and 7 variables:
\describe{
  \item{\code{id}}{Observation identifier.}
  \item{\code{truth}}{Complete reference class label, coded 1 or 2.}
  \item{\code{observed}}{Class label supplied during fitting; \code{NA} denotes a missing label.}
  \item{\code{missing_indicator}}{Label status/source: \code{observed}, \code{mar}, or \code{mcar}.}
  \item{\code{Recency}}{Standardized recency measurement.}
  \item{\code{Frequency}}{Standardized donation-frequency measurement.}
  \item{\code{Time}}{Standardized time measurement.}
}
}
\details{
The complete labels in \code{truth} are retained only for evaluation in the
semi-synthetic case study. The response supplied to \code{fit_sslfmm()} is
\code{observed}.
}
\examples{
data("blood_transfusion")
table(blood_transfusion$missing_indicator)
}
\keyword{datasets}
