\name{classification_performance}
\alias{classification_performance}
\title{Classification Performance Metrics}
\usage{
classification_performance(truth, predicted, posterior = NULL)
}
\arguments{
\item{truth}{True class labels.}
\item{predicted}{Predicted class labels; must have the same length as \code{truth}.}
\item{posterior}{Optional numeric matrix with one row per observation and nonnegative finite class probabilities. If columns are named, names must be unique and include every true class used in the calculation.}
}
\description{
Computes overall and per-class classification performance after removing observations for which \code{truth} or \code{predicted} is missing. Macro averages are taken over classes represented in the retained truth labels. If a represented true class is never predicted, its precision and F1 score are defined as zero and remain in the macro average.
}
\value{
A list with:
\itemize{
\item \code{metrics}: named numeric vector containing accuracy, error rate, balanced accuracy, macro precision, macro recall, and macro F1; when \code{posterior} is supplied, log loss is added;
\item \code{per_class}: data.frame with class support, precision, recall, and F1;
\item \code{confusion_matrix}: truth-by-predicted contingency table.
}
Posterior rows are normalized internally to unit sum after validation; negative or non-finite values are rejected rather than silently corrected.
}
\examples{
truth <- c(1, 1, 2, 2)
pred <- c(1, 2, 2, 2)
post <- matrix(c(.8, .2,
                 .4, .6,
                 .1, .9,
                 .3, .7), ncol = 2, byrow = TRUE)
colnames(post) <- c("1", "2")
classification_performance(truth, pred, post)
}
