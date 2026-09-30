.prepare_fit_data <- function(x, y, g = NULL, source = NULL) {
  x <- .check_x(x)
  lab <- .encode_labels(y, g)
  if (length(lab$z) != nrow(x)) stop("'y' must have one entry per row of 'x'.", call. = FALSE)
  z <- lab$z
  m <- is.na(z)
  out <- list(x = x, z = z, m = m, g = lab$g, label_levels = lab$levels)
  if (!is.null(source)) out$source <- source
  out
}

.normalize_source <- function(source, m) {
  n <- length(m)
  if (is.null(source)) {
    stop("For mixed missingness with indicator = 'known', supply 'missing_source'.", call. = FALSE)
  }
  if (length(source) != n) stop("'missing_source' must have one entry per observation.", call. = FALSE)

  # A logical or 0/1 vector is interpreted as a known MCAR-source indicator:
  # TRUE/1 = MCAR source; FALSE/0 = MAR source on rows whose labels are missing.
  if (is.logical(source) || (is.numeric(source) && all(is.na(source) | source %in% c(0, 1)))) {
    v <- as.logical(source)
    s <- rep("observed", n)
    if (any(m & is.na(v))) stop("The known source indicator cannot be NA on rows with missing labels.", call. = FALSE)
    s[m & v] <- "mcar"
    s[m & !v] <- "mar"
    return(s)
  }

  s <- tolower(as.character(source))
  s[!m & (is.na(s) | s == "" | s == "none" | s == "0")] <- "observed"
  s[!m] <- "observed"
  aliases <- c("mcar" = "mcar", "mar" = "mar", "observed" = "observed",
               "cc" = "observed", "labelled" = "observed", "labeled" = "observed")
  bad <- is.na(s) | !(s %in% names(aliases))
  if (any(bad)) {
    stop("'missing_source' must identify each row as 'observed', 'mcar', or 'mar', or be a logical 0/1 MCAR-source indicator.", call. = FALSE)
  }
  s <- unname(aliases[s])
  if (any(m & s == "observed")) {
    stop("Rows with missing labels must have source 'mcar' or 'mar'.", call. = FALSE)
  }
  s
}

.closed_form_cc <- function(data, covariance_type, ridge = 1e-6) {
  x <- data$x; z <- data$z; g <- data$g; p <- ncol(x)
  obs <- !is.na(z)
  if (!any(obs)) stop("Complete-case fitting requires at least one observed class label.", call. = FALSE)
  counts <- tabulate(z[obs], nbins = g)
  if (any(counts == 0L)) {
    stop("Complete-case fitting requires at least one labelled observation from every component.", call. = FALSE)
  }
  mu <- matrix(NA_real_, g, p)
  for (k in seq_len(g)) mu[k, ] <- colMeans(x[obs & z == k, , drop = FALSE])
  pi_hat <- counts / sum(counts)

  if (covariance_type == "equal") {
    xo <- x[obs, , drop = FALSE]
    zo <- z[obs]
    resid <- xo
    for (i in seq_len(nrow(xo))) resid[i, ] <- xo[i, ] - mu[zo[i], ]
    S <- crossprod(resid) / nrow(resid)
    sigma <- .ensure_spd(S, ridge = ridge)
  } else {
    sigma <- array(NA_real_, c(p, p, g))
    fallback <- .safe_cov(x[obs, , drop = FALSE], ridge)
    for (k in seq_len(g)) {
      ik <- obs & z == k
      if (sum(ik) >= 2L) {
        centered <- sweep(x[ik, , drop = FALSE], 2L, mu[k, ], "-")
        S <- crossprod(centered) / sum(ik)
        sigma[, , k] <- .ensure_spd(S, ridge = ridge)
      } else {
        sigma[, , k] <- fallback
      }
    }
  }
  theta <- list(pi = as.numeric(pi_hat), mu = mu, sigma = sigma)
  lc <- .model_components(x[obs, , drop = FALSE], theta, covariance_type)
  zz <- z[obs]
  ll <- sum(lc$log_joint[cbind(seq_along(zz), zz)])
  ss <- data.frame(
    start = 1L,
    initial_objective = -ll,
    objective = -ll,
    convergence = 0L,
    iterations = 0L,
    optimizer = "closed-form",
    message = "",
    stringsAsFactors = FALSE
  )
  diag <- data.frame(
    start = 1L, initial_objective = -ll, optimizer = "closed-form",
    status = "closed_form", message = "No numerical optimization was required.",
    stringsAsFactors = FALSE
  )
  list(theta = theta, loglik = ll, convergence = 0L, optimizer = "closed-form",
       start_summary = ss, start_diagnostics = diag)
}

.fit_mcar_internal <- function(data, init, covariance_type, n_starts, control,
                               entropy_floor, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- ncol(data$x); g <- data$g
  p0 <- .pack_theta(init, covariance_type)
  bd <- .theta_bounds(p, g, covariance_type)
  starts <- c(list(p0), replicate(max(0L, as.integer(n_starts) - 1L),
                                  .perturb_start(p0), simplify = FALSE))
  ms <- .run_multistart_nlminb(starts, .nll_mcar, lower = bd$lower, upper = bd$upper,
                               control = control, data = data, p = p, g = g,
                               covariance_type = covariance_type,
                               entropy_floor = entropy_floor)
  best <- ms$best
  list(par = best$par,
       theta = .unpack_theta(best$par, p, g, covariance_type),
       loglik = -best$objective, convergence = best$convergence,
       optimizer = best$optimizer %||% "nlminb",
       start_summary = ms$summary,
       start_diagnostics = ms$diagnostics)
}

.xi_bounds <- function(theta_bounds, include_alpha = FALSE, alpha_upper = 1 - 1e-8) {
  lo <- c(theta_bounds$lower, xi0 = -40, eta_xi = -7)
  hi <- c(theta_bounds$upper, xi0 = 40, eta_xi = 7)
  if (include_alpha) {
    lo <- c(lo, alpha = 0)
    hi <- c(hi, alpha = alpha_upper)
  }
  list(lower = lo, upper = hi)
}

#' Fit a semi-supervised finite mixture model
#'
#' Fits semi-supervised Gaussian finite-mixture classifiers under complete-case, MCAR,
#' entropy-dependent MAR, or mixed MCAR/MAR label-missingness formulations. For mixed
#' models, the MCAR/MAR source may be known or unknown.
#'
#' @param x Numeric feature matrix or data.frame; feature values must be fully observed.
#' @param y Class labels in 1:g (or factor/character), with NA for unlabelled observations.
#' @param g Number of mixture components. Inferred from observed labels when possible.
#' @param method One of "cc", "mcar", "mar", or "mixed".
#' @param covariance_type "equal" for one covariance matrix shared by all components or
#'   "unequal" for one covariance matrix per component.
#' @param indicator For method="mixed", whether the MCAR/MAR source indicator is
#'   "unknown" or "known".
#' @param missing_source For method="mixed" and indicator="known", a vector containing
#'   "mcar" or "mar" for rows with missing labels. Values on labelled rows are ignored.
#' @param init Optional initialization returned by initialize_sslfmm() or a list with pi,
#'   mu and sigma.
#' @param n_starts Number of fitting starts.
#' @param control Control list for numerical fitting.
#' @param entropy_floor Small positive lower bound for posterior entropy.
#' @param alpha_starts Starting values for the latent MCAR mixing probability in the mixed model.
#' @param alpha_upper Upper bound for alpha, strictly below 1.
#' @param ridge Covariance regularization used for initialization and complete-case fitting.
#' @param seed Optional random seed for initialization and perturbed starts.
#' @return An object of class "SSLfmm".
#' @export
fit_sslfmm <- function(x, y, g = NULL,
                        method = c("mixed", "mar", "mcar", "cc"),
                        covariance_type = c("equal", "unequal"),
                        indicator = c("unknown", "known"),
                        missing_source = NULL,
                        init = NULL,
                        n_starts = 5L,
                        control = list(),
                        entropy_floor = 1e-300,
                        alpha_starts = c(0, 0.05, 0.15, 0.30),
                        alpha_upper = 1 - 1e-8,
                        ridge = 1e-4,
                        seed = NULL) {
  call <- match.call()
  method <- match.arg(method)
  covariance_type <- match.arg(covariance_type)
  indicator <- match.arg(indicator)
  if (!is.numeric(entropy_floor) || length(entropy_floor) != 1L ||
      !is.finite(entropy_floor) || entropy_floor <= 0) {
    stop("'entropy_floor' must be a single positive finite number.", call. = FALSE)
  }
  if (!is.numeric(n_starts) || length(n_starts) != 1L || !is.finite(n_starts) ||
      n_starts < 1L || abs(n_starts - round(n_starts)) > sqrt(.Machine$double.eps)) {
    stop("'n_starts' must be a positive integer.", call. = FALSE)
  }
  if (length(ridge) != 1L || !is.numeric(ridge) || !is.finite(ridge) || ridge <= 0) {
    stop("'ridge' must be a single positive finite number.", call. = FALSE)
  }
  if (!is.list(control)) {
    stop("'control' must be a list suitable for stats::nlminb().", call. = FALSE)
  }

  label_type <- if (is.factor(y)) "factor" else if (is.character(y)) "character" else "integer"
  data <- .prepare_fit_data(x, y, g)
  if (!any(!data$m)) stop("At least one observed class label is required.", call. = FALSE)
  g <- data$g; p <- ncol(data$x)

  if (g == 1L && method %in% c("mar", "mixed")) {
    stop(
      sprintf(
        "Method '%s' requires g >= 2 for fitting because posterior entropy is identically zero when g = 1; use g = 1 only for simulation boundary checks or fit method = 'cc'/'mcar'.",
        method
      ),
      call. = FALSE
    )
  }

  if (method %in% c("mar", "mixed") && (all(data$m) || !any(data$m))) {
    stop(sprintf("Method '%s' requires both labelled and unlabelled observations.", method), call. = FALSE)
  }

  if (is.null(init)) {
    init <- initialize_sslfmm(data$x, data$z, g = g,
                              covariance_type = covariance_type,
                              ridge = ridge, seed = seed)
  } else {
    if (!is.list(init) || !all(c("pi", "mu", "sigma") %in% names(init))) {
      stop("'init' must be a list containing 'pi', 'mu', and 'sigma'.", call. = FALSE)
    }
    init_pi <- as.numeric(init$pi)
    if (length(init_pi) != g || any(!is.finite(init_pi)) || any(init_pi <= 0) ||
        abs(sum(init_pi) - 1) > 1e-8) {
      stop(sprintf("'init$pi' must contain %d positive finite probabilities summing to 1.", g), call. = FALSE)
    }
    if (!is.matrix(init$mu) || !is.numeric(init$mu)) {
      stop("'init$mu' must be a numeric matrix.", call. = FALSE)
    }
    init_mu <- init$mu
    if (identical(dim(init_mu), c(g, p))) {
      # Fitting/initialization orientation: one component per row.
    } else if (g != p && identical(dim(init_mu), c(p, g))) {
      # Also accept the simulation orientation when it is unambiguous.
      init_mu <- t(init_mu)
    } else {
      stop(sprintf("'init$mu' must have dimension %d x %d (g x p)%s.",
                   g, p, if (g != p) sprintf(" or %d x %d (p x g)", p, g) else ""),
           call. = FALSE)
    }
    if (any(!is.finite(init_mu))) stop("'init$mu' must contain only finite values.", call. = FALSE)
    ns <- .normalize_sigma(init$sigma, p, g, covariance_type, ridge = ridge)
    init <- list(pi = init_pi, mu = init_mu, sigma = ns$sigma)
  }

  xi <- NULL; alpha <- NULL
  internal_par <- NULL

  if (method == "cc") {
    fit <- .closed_form_cc(data, covariance_type, ridge = ridge)
  } else {
    mcar_start_count <- if (method == "mcar") as.integer(n_starts) else
      max(2L, min(as.integer(n_starts), 3L))
    mcar_fit <- .fit_mcar_internal(data, init, covariance_type,
                                   n_starts = mcar_start_count,
                                   control = control, entropy_floor = entropy_floor,
                                   seed = seed)
    if (method == "mcar") {
      fit <- mcar_fit
      internal_par <- fit$par
    } else if (method == "mar") {
      p0xi <- .initial_xi(data$x, mcar_fit$theta, covariance_type,
                          response = data$m, entropy_floor = entropy_floor)
      base <- c(mcar_fit$par, p0xi)
      bd <- .xi_bounds(.theta_bounds(p, g, covariance_type), include_alpha = FALSE)
      starts <- c(list(base), replicate(max(0L, as.integer(n_starts) - 1L),
                                        .perturb_start(base), simplify = FALSE))
      ms <- .run_multistart_nlminb(starts, .nll_mar, lower = bd$lower, upper = bd$upper,
                                   control = control, data = data, p = p, g = g,
                                   covariance_type = covariance_type,
                                   entropy_floor = entropy_floor)
      best <- ms$best
      dt <- .theta_length(p, g, covariance_type)
      fit <- list(par = best$par,
                  theta = .unpack_theta(best$par[seq_len(dt)], p, g, covariance_type),
                  loglik = -best$objective, convergence = best$convergence,
                  optimizer = best$optimizer %||% "nlminb",
                  start_summary = ms$summary,
                  start_diagnostics = ms$diagnostics)
      xi <- c(xi0 = best$par[dt + 1L], xi1 = exp(best$par[dt + 2L]))
      alpha <- NULL
      internal_par <- best$par
    } else if (indicator == "known") {
      data$source <- .normalize_source(missing_source, data$m)
      alpha <- mean(data$source == "mcar")
      if (alpha >= 1) stop("Known mixed-source data cannot have every row assigned to the MCAR source.", call. = FALSE)
      eligible <- data$source != "mcar"
      response <- data$source == "mar"
      p0xi <- .initial_xi(data$x, mcar_fit$theta, covariance_type,
                          response = response, eligible = eligible,
                          entropy_floor = entropy_floor)
      base <- c(mcar_fit$par, p0xi)
      bd <- .xi_bounds(.theta_bounds(p, g, covariance_type), include_alpha = FALSE)
      starts <- c(list(base), replicate(max(0L, as.integer(n_starts) - 1L),
                                        .perturb_start(base), simplify = FALSE))
      ms <- .run_multistart_nlminb(starts, .nll_mixed_observed,
                                   lower = bd$lower, upper = bd$upper, control = control,
                                   data = data, p = p, g = g,
                                   covariance_type = covariance_type, alpha = alpha,
                                   entropy_floor = entropy_floor)
      best <- ms$best
      dt <- .theta_length(p, g, covariance_type)
      fit <- list(par = best$par,
                  theta = .unpack_theta(best$par[seq_len(dt)], p, g, covariance_type),
                  loglik = -best$objective, convergence = best$convergence,
                  optimizer = best$optimizer %||% "nlminb",
                  start_summary = ms$summary,
                  start_diagnostics = ms$diagnostics)
      xi <- c(xi0 = best$par[dt + 1L], xi1 = exp(best$par[dt + 2L]))
      internal_par <- best$par
    } else {
      p0xi <- .initial_xi(data$x, mcar_fit$theta, covariance_type,
                          response = data$m, entropy_floor = entropy_floor)
      if (!is.numeric(alpha_upper) || length(alpha_upper) != 1L || !is.finite(alpha_upper) ||
          alpha_upper <= 0 || alpha_upper >= 1) {
        stop("'alpha_upper' must be a single finite number strictly between 0 and 1.", call. = FALSE)
      }
      aa <- as.numeric(alpha_starts)
      if (!length(aa) || any(!is.finite(aa))) {
        stop("'alpha_starts' must contain finite numeric starting values.", call. = FALSE)
      }
      if (any(aa < 0 | aa > alpha_upper)) {
        stop("Every value in 'alpha_starts' must lie between 0 and 'alpha_upper'.", call. = FALSE)
      }
      aa <- unique(aa)
      aa <- aa[seq_len(min(length(aa), as.integer(n_starts)))]
      base_theta_xi <- c(mcar_fit$par, p0xi)
      starts <- lapply(aa, function(a) c(base_theta_xi, alpha = a))
      while (length(starts) < as.integer(n_starts)) {
        starts[[length(starts) + 1L]] <- .perturb_start(starts[[1L]], 0.03)
      }
      bd <- .xi_bounds(.theta_bounds(p, g, covariance_type), include_alpha = TRUE,
                       alpha_upper = alpha_upper)
      ms <- .run_multistart_nlminb(starts, .nll_mixed_latent,
                                   lower = bd$lower, upper = bd$upper, control = control,
                                   data = data, p = p, g = g,
                                   covariance_type = covariance_type,
                                   entropy_floor = entropy_floor)
      best <- ms$best
      dt <- .theta_length(p, g, covariance_type)
      fit <- list(par = best$par,
                  theta = .unpack_theta(best$par[seq_len(dt)], p, g, covariance_type),
                  loglik = -best$objective, convergence = best$convergence,
                  optimizer = best$optimizer %||% "nlminb",
                  start_summary = ms$summary,
                  start_diagnostics = ms$diagnostics)
      xi <- c(xi0 = best$par[dt + 1L], xi1 = exp(best$par[dt + 2L]))
      alpha <- best$par[dt + 3L]
      fit$start_summary$alpha_solution <- vapply(ms$all, function(z) z$par[dt + 3L], numeric(1))
      internal_par <- best$par
    }
  }

  theta <- fit$theta
  feature_names <- colnames(data$x)
  if (is.null(feature_names)) feature_names <- paste0("x", seq_len(p))
  names(theta$pi) <- data$label_levels[seq_len(g)]
  rownames(theta$mu) <- data$label_levels[seq_len(g)]
  colnames(theta$mu) <- feature_names
  components <- .model_components(data$x, theta, covariance_type)
  colnames(components$posterior) <- data$label_levels[seq_len(g)]
  entropy <- -rowSums(ifelse(components$posterior > 0,
                             components$posterior * log(components$posterior), 0))
  mar_probability <- NULL
  missing_probability <- NULL
  latent_missing <- NULL
  latent_missing_probability <- NULL
  if (!is.null(xi)) {
    mar_probability <- stats::plogis(xi["xi0"] + xi["xi1"] * log(pmax(entropy, entropy_floor)))
  }
  if (method == "mar") {
    missing_probability <- mar_probability
  } else if (method == "mixed") {
    missing_probability <- alpha + (1 - alpha) * mar_probability
    if (indicator == "known") {
      latent_missing <- data$source == "mcar"
      latent_missing_probability <- as.numeric(latent_missing)
    } else {
      latent_missing <- rep(NA, length(data$m))
      latent_missing_probability <- rep(0, length(data$m))
      denom <- pmax(missing_probability[data$m], .Machine$double.xmin)
      latent_missing_probability[data$m] <- alpha / denom
    }
  }
  out <- list(call = call,
              method = method,
              covariance_type = covariance_type,
              indicator = if (method == "mixed") indicator else NA_character_,
              g = g, p = p,
              pi = theta$pi, mu = theta$mu, sigma = theta$sigma,
              theta = theta, xi = xi, alpha = alpha,
              loglik = fit$loglik, convergence = fit$convergence,
              optimizer = fit$optimizer %||% (if (method == "cc") "closed-form" else NA_character_),
              start_summary = fit$start_summary,
              start_diagnostics = fit$start_diagnostics %||% NULL,
              posterior = components$posterior,
              en = entropy,
              entropy = entropy,
              missing = data$m,
              observed_missing = data$m,
              latent_missing = latent_missing,
              latent_missing_probability = latent_missing_probability,
              mar_probability = mar_probability,
              missing_probability = missing_probability,
              missing_source = if (!is.null(data$source)) data$source else NULL,
              label_levels = data$label_levels,
              label_type = label_type,
              internal_par = internal_par)
  class(out) <- "SSLfmm"
  out
}

#' @export
print.SSLfmm <- function(x, ...) {
  cat("SSLfmm fit\n")
  cat("  method:             ", x$method, "\n", sep = "")
  cat("  covariance:         ", x$covariance_type, "\n", sep = "")
  if (x$method == "mixed") cat("  source indicator:    ", x$indicator, "\n", sep = "")
  cat("  components:          ", x$g, "\n", sep = "")
  cat("  log-likelihood:      ", format(x$loglik, digits = 6), "\n", sep = "")
  cat("  convergence code:    ", x$convergence, "\n", sep = "")
  if (!is.null(x$optimizer)) cat("  optimizer:           ", x$optimizer, "\n", sep = "")
  if (x$method == "mixed" && !is.null(x$alpha)) cat("  alpha:               ", format(x$alpha, digits = 5), "\n", sep = "")
  if (!is.null(x$xi)) cat("  xi:                  ", paste(format(x$xi, digits = 5), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' @export
summary.SSLfmm <- function(object, ...) {
  ans <- list(method = object$method,
              covariance_type = object$covariance_type,
              indicator = object$indicator,
              g = object$g,
              pi = object$pi,
              mu = object$mu,
              sigma = object$sigma,
              xi = object$xi,
              alpha = object$alpha,
              loglik = object$loglik,
              convergence = object$convergence,
              optimizer = object$optimizer,
              missing_rate = mean(object$observed_missing))
  class(ans) <- "summary.SSLfmm"
  ans
}

#' Predict classes, posterior probabilities, or entropy
#'
#' @param object An SSLfmm fit.
#' @param newdata Numeric feature matrix or data.frame.
#' @param type One of "class", "posterior", "entropy", or "all".
#' @param ... Unused.
#' @return Predicted class labels, posterior probabilities, entropy, or a list containing all three.
#' @export
predict.SSLfmm <- function(object, newdata, type = c("class", "posterior", "entropy", "all"), ...) {
  type <- match.arg(type)
  x <- .check_x(newdata)
  if (ncol(x) != object$p) stop("'newdata' has the wrong number of columns.", call. = FALSE)
  comp <- .model_components(x, object$theta, object$covariance_type)
  post <- comp$posterior
  colnames(post) <- object$label_levels[seq_len(object$g)]
  cls_idx <- max.col(post, ties.method = "first")
  if (object$label_type == "integer") {
    cls <- as.integer(cls_idx)
  } else if (object$label_type == "factor") {
    cls <- factor(object$label_levels[cls_idx], levels = object$label_levels)
  } else {
    cls <- object$label_levels[cls_idx]
  }
  en <- -rowSums(ifelse(post > 0, post * log(post), 0))
  if (type == "class") return(cls)
  if (type == "posterior") return(post)
  if (type == "entropy") return(en)
  list(class = cls, posterior = post, entropy = en)
}
