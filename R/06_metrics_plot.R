.safe_div <- function(a, b, zero = NA_real_) {
  ifelse(b > 0, a / b, zero)
}

#' Classification performance metrics
#'
#' @param truth True class labels.
#' @param predicted Predicted class labels.
#' @param posterior Optional matrix of class posterior probabilities, one row per observation.
#' @return A list with overall metrics, per-class metrics and a confusion matrix.
#' @details Macro averages are taken over classes represented in `truth` after
#'   removing incomplete truth/prediction pairs. If a represented true class is
#'   never predicted, its precision and F1 score are defined as zero rather than
#'   omitted from the macro average.
#' @export
classification_performance <- function(truth, predicted, posterior = NULL) {
  if (length(truth) != length(predicted)) stop("'truth' and 'predicted' must have the same length.", call. = FALSE)
  keep <- !is.na(truth) & !is.na(predicted)
  if (!any(keep)) stop("No complete truth/prediction pairs are available.", call. = FALSE)
  tt <- as.character(truth[keep])
  pp <- as.character(predicted[keep])
  lev <- sort(unique(c(tt, pp)))
  cm <- table(factor(tt, levels = lev), factor(pp, levels = lev),
              dnn = c("truth", "predicted"))
  tp <- diag(cm)
  support <- rowSums(cm)
  predicted_n <- colSums(cm)

  # Macro metrics are defined over classes that occur in the truth labels.
  # A true class that receives no predictions must contribute zero precision
  # and zero F1, rather than being silently removed with na.rm = TRUE.
  truth_class <- support > 0
  recall <- .safe_div(tp, support)
  precision <- .safe_div(tp, predicted_n, zero = 0)
  f1 <- .safe_div(2 * precision * recall, precision + recall, zero = 0)

  accuracy <- sum(tp) / sum(cm)
  metrics <- c(accuracy = accuracy,
               error_rate = 1 - accuracy,
               balanced_accuracy = mean(recall[truth_class]),
               macro_precision = mean(precision[truth_class]),
               macro_recall = mean(recall[truth_class]),
               macro_f1 = mean(f1[truth_class]))

  if (!is.null(posterior)) {
    posterior <- as.matrix(posterior)
    if (nrow(posterior) != length(truth)) stop("'posterior' must have one row per observation.", call. = FALSE)
    post <- posterior[keep, , drop = FALSE]
    if (!is.numeric(post) || any(!is.finite(post))) {
      stop("'posterior' must contain only finite numeric probabilities.", call. = FALSE)
    }
    if (any(post < 0)) {
      stop("'posterior' probabilities must be nonnegative.", call. = FALSE)
    }
    if (is.null(colnames(post))) {
      if (ncol(post) != length(lev)) {
        stop("Without column names, 'posterior' must have one column per class represented by 'truth'/'predicted'.", call. = FALSE)
      }
      colnames(post) <- lev
    } else if (anyDuplicated(colnames(post))) {
      stop("'posterior' column names must be unique class labels.", call. = FALSE)
    }
    idx <- match(tt, colnames(post))
    if (anyNA(idx)) stop("'posterior' must contain a column for every true class.", call. = FALSE)
    rs <- rowSums(post)
    if (any(!is.finite(rs)) || any(rs <= 0)) stop("Each 'posterior' row must contain positive finite probability mass.", call. = FALSE)
    post <- post / rs
    ptrue <- post[cbind(seq_along(idx), idx)]
    log_loss <- -mean(log(pmax(ptrue, 1e-15)))
    metrics <- c(metrics, log_loss = log_loss)
  }

  per_class <- data.frame(class = lev, support = as.integer(support),
                          precision = as.numeric(precision), recall = as.numeric(recall),
                          f1 = as.numeric(f1), row.names = NULL)
  list(metrics = metrics, per_class = per_class, confusion_matrix = cm)
}
