`%||%` <- function(x, y) if (is.null(x)) y else x

.softplus <- function(x) {
  out <- numeric(length(x))
  hi <- x > 30
  lo <- x < -30
  mid <- !(hi | lo)
  out[hi] <- x[hi] + log1p(exp(-x[hi]))
  out[lo] <- exp(x[lo])
  out[mid] <- log1p(exp(x[mid]))
  out
}

.log_sigmoid <- function(x) -.softplus(-x)
.log1m_sigmoid <- function(x) -.softplus(x)

.row_logsumexp <- function(M) {
  M <- as.matrix(M)
  n <- nrow(M)
  mx <- apply(M, 1L, max)
  out <- rep(-Inf, n)
  good <- is.finite(mx)
  if (any(good)) {
    Mg <- M[good, , drop = FALSE]
    mg <- mx[good]
    out[good] <- mg + log(rowSums(exp(sweep(Mg, 1L, mg, "-"))))
  }
  out
}

.logsumexp2 <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  out[!is.finite(m)] <- m[!is.finite(m)]
  out
}

.check_x <- function(x) {
  if (is.data.frame(x)) {
    numeric_cols <- vapply(x, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      stop("'x' must be a numeric matrix or a data.frame containing only numeric columns.", call. = FALSE)
    }
  } else if (!is.matrix(x) && !is.numeric(x)) {
    stop("'x' must be a numeric matrix or data.frame.", call. = FALSE)
  } else if (is.matrix(x) && !is.numeric(x)) {
    stop("'x' must be a numeric matrix or data.frame.", call. = FALSE)
  }

  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (length(dim(x)) != 2L || nrow(x) < 1L || ncol(x) < 1L) {
    stop("'x' must be a non-empty numeric matrix or data.frame.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("'x' must contain only finite numeric values; feature missingness is not supported.", call. = FALSE)
  }
  x
}

.check_g <- function(g) {
  if (length(g) != 1L || !is.finite(g) || g < 1 || g != as.integer(g)) {
    stop("'g' must be a positive integer.", call. = FALSE)
  }
  as.integer(g)
}

.encode_labels <- function(y, g = NULL) {
  if (is.null(y)) stop("'y' is required and may contain NA for missing class labels.", call. = FALSE)
  if (is.factor(y) || is.character(y)) {
    yy <- as.character(y)
    lev <- sort(unique(yy[!is.na(yy)]))
    if (!length(lev)) {
      if (is.null(g)) stop("'g' must be supplied when all labels are missing.", call. = FALSE)
      lev <- as.character(seq_len(.check_g(g)))
    }
    if (!is.null(g) && length(lev) > g) stop("Observed labels contain more than 'g' classes.", call. = FALSE)
    map <- stats::setNames(seq_along(lev), lev)
    z <- rep(NA_integer_, length(yy))
    ok <- !is.na(yy)
    z[ok] <- unname(map[yy[ok]])
    if (anyNA(z[ok])) stop("Could not encode one or more observed labels.", call. = FALSE)
    gg <- g %||% length(lev)
    gg <- .check_g(gg)
    if (length(lev) < gg) lev <- c(lev, paste0("class", seq.int(length(lev) + 1L, gg)))
    return(list(z = z, levels = lev, g = gg))
  }

  if (!is.numeric(y) && !is.integer(y)) {
    stop("'y' must be integer-like class labels, a factor, or character vector.", call. = FALSE)
  }
  yy_num <- as.numeric(y)
  ok_num <- !is.na(yy_num)
  if (any(ok_num & abs(yy_num - round(yy_num)) > sqrt(.Machine$double.eps))) {
    stop("Numeric labels must be integer-like values.", call. = FALSE)
  }
  z <- as.integer(yy_num)
  ok <- !is.na(z)
  if (any(ok & (z < 1L))) stop("Numeric labels must be positive integers starting at 1.", call. = FALSE)
  gg <- g %||% if (any(ok)) max(z[ok]) else NA_integer_
  if (!is.finite(gg)) stop("'g' must be supplied when all labels are missing.", call. = FALSE)
  gg <- .check_g(gg)
  if (any(ok & z > gg)) stop("Observed numeric labels must lie in 1:'g'.", call. = FALSE)
  list(z = z, levels = as.character(seq_len(gg)), g = gg)
}

.safe_cov <- function(x, ridge = 1e-6) {
  x <- as.matrix(x)
  p <- ncol(x)
  if (nrow(x) <= 1L) {
    v <- rep(1, p)
    S <- diag(v, p)
  } else if (p == 1L) {
    v <- stats::var(as.numeric(x))
    if (!is.finite(v) || v <= 0) v <- 1
    S <- matrix(v, 1L, 1L)
  } else {
    S <- stats::cov(x)
    if (any(!is.finite(S))) S <- diag(1, p)
  }
  S <- (as.matrix(S) + t(as.matrix(S))) / 2
  S + diag(ridge, p)
}

.ensure_spd <- function(S, ridge = 1e-8, name = "sigma") {
  S <- as.matrix(S)
  if (!is.numeric(S)) stop(sprintf("'%s' must be numeric.", name), call. = FALSE)
  storage.mode(S) <- "double"
  if (nrow(S) != ncol(S) || nrow(S) < 1L) stop(sprintf("'%s' must be square.", name), call. = FALSE)
  if (any(!is.finite(S))) stop(sprintf("'%s' must contain only finite values.", name), call. = FALSE)
  asym <- max(abs(S - t(S)))
  tol <- 100 * .Machine$double.eps * max(1, max(abs(S)))
  if (!is.finite(asym) || asym > tol) {
    stop(sprintf("'%s' must be symmetric.", name), call. = FALSE)
  }
  S <- (S + t(S)) / 2
  p <- nrow(S)
  if (!is.finite(ridge) || ridge < 0) stop("'ridge' must be nonnegative.", call. = FALSE)
  if (ridge == 0) {
    if (inherits(try(chol(S), silent = TRUE), "try-error"))
      stop(sprintf("'%s' must be positive definite.", name), call. = FALSE)
    return(S)
  }
  add <- ridge
  for (i in 0:8) {
    candidate <- S + diag(add, p)
    if (!inherits(try(chol(candidate), silent = TRUE), "try-error")) return(candidate)
    add <- add * 10
  }
  stop(sprintf("'%s' could not be made positive definite using the requested ridge regularization.", name), call. = FALSE)
}

.normalize_sigma <- function(sigma, p, g, covariance_type = NULL, ridge = 0) {
  if (is.null(dim(sigma))) {
    if (p == 1L && length(sigma) == 1L) sigma <- matrix(as.numeric(sigma), 1L, 1L)
    else stop("'sigma' must be a p x p matrix or a p x p x g array; for p = 1, a length-one scalar is also accepted as a shared variance.", call. = FALSE)
  }

  ds <- dim(sigma)
  if (length(ds) == 2L) {
    if (!identical(as.integer(ds), c(p, p))) {
      stop(sprintf("A matrix 'sigma' must have dimension %d x %d.", p, p), call. = FALSE)
    }
    S <- .ensure_spd(sigma, ridge = ridge, name = "sigma")
    if (is.null(covariance_type)) covariance_type <- "equal"
    if (!identical(covariance_type, "equal")) {
      arr <- array(NA_real_, dim = c(p, p, g))
      for (k in seq_len(g)) arr[, , k] <- S
      return(list(sigma = arr, covariance_type = "unequal"))
    }
    return(list(sigma = S, covariance_type = "equal"))
  }

  if (length(ds) == 3L) {
    if (!identical(as.integer(ds), c(p, p, g))) {
      stop(sprintf("An array 'sigma' must have dimension %d x %d x %d.", p, p, g), call. = FALSE)
    }
    arr <- array(NA_real_, dim = c(p, p, g))
    for (k in seq_len(g)) {
      Sk <- matrix(sigma[, , k], nrow = p, ncol = p)
      arr[, , k] <- .ensure_spd(Sk, ridge = ridge,
                                name = sprintf("sigma[,,%d]", k))
    }
    if (is.null(covariance_type)) covariance_type <- "unequal"
    if (identical(covariance_type, "equal")) {
      ref <- arr[, , 1L]
      same <- all(vapply(seq_len(g), function(k) isTRUE(all.equal(ref, arr[, , k], tolerance = 1e-10)), logical(1)))
      if (!same) stop("'covariance_type = \"equal\"' requires identical covariance matrices across components.", call. = FALSE)
      return(list(sigma = ref, covariance_type = "equal"))
    }
    return(list(sigma = arr, covariance_type = "unequal"))
  }

  stop("'sigma' must be a p x p matrix or a p x p x g array; for p = 1, a length-one scalar is also accepted as a shared variance.", call. = FALSE)
}

.run_multistart_nlminb <- function(starts, objective, lower, upper, control = list(), ...) {
  if (!length(starts)) stop("At least one optimization start is required.", call. = FALSE)
  dots <- list(...)
  d <- length(starts[[1L]])
  if (length(lower) != d || length(upper) != d) {
    stop("Internal optimization bounds do not match the parameter-vector length.", call. = FALSE)
  }
  if (any(lower > upper)) stop("Internal optimization lower bounds exceed upper bounds.", call. = FALSE)

  penalty <- 1e100
  eval_objective <- function(par) {
    val <- tryCatch(
      do.call(objective, c(list(as.numeric(par)), dots)),
      error = function(e) penalty
    )
    if (length(val) != 1L || !is.numeric(val) || !is.finite(val) || val >= penalty) {
      return(penalty)
    }
    as.numeric(val)
  }

  # nlminb is the primary optimizer, matching the original research code.  Some R/BLAS
  # combinations can nevertheless abort a start during internal finite differencing even
  # when the objective is finite.  In that case use bounded L-BFGS-B on the *same*
  # observed-data objective.  This is an optimizer fallback, not a change to the model.
  optim_control <- list(maxit = as.integer(control$iter.max %||% control$eval.max %||% 1000L))
  if (!is.null(control$trace)) optim_control$trace <- as.integer(control$trace)

  fits <- vector("list", length(starts))
  diagnostics <- vector("list", length(starts))

  standardize_nlminb <- function(z) {
    list(
      par = as.numeric(z$par),
      objective = as.numeric(z$objective),
      convergence = as.integer(z$convergence),
      iterations = z$iterations %||% NA_integer_,
      message = z$message %||% "",
      optimizer = "nlminb"
    )
  }

  standardize_optim <- function(z) {
    list(
      par = as.numeric(z$par),
      objective = as.numeric(z$value),
      convergence = as.integer(z$convergence),
      iterations = unname(z$counts[["function"]] %||% NA_integer_),
      message = z$message %||% "",
      optimizer = "L-BFGS-B"
    )
  }

  for (i in seq_along(starts)) {
    if (length(starts[[i]]) != d || any(!is.finite(starts[[i]]))) {
      diagnostics[[i]] <- data.frame(
        start = i, initial_objective = NA_real_, optimizer = NA_character_,
        status = "invalid_start",
        message = "Starting vector has the wrong length or contains non-finite values.",
        stringsAsFactors = FALSE
      )
      next
    }

    st <- pmin(pmax(as.numeric(starts[[i]]), lower), upper)
    init_raw <- tryCatch(
      do.call(objective, c(list(st), dots)),
      error = function(e) e
    )
    if (inherits(init_raw, "error")) {
      diagnostics[[i]] <- data.frame(
        start = i, initial_objective = NA_real_, optimizer = NA_character_,
        status = "initial_objective_error",
        message = conditionMessage(init_raw), stringsAsFactors = FALSE
      )
      next
    }
    if (length(init_raw) != 1L || !is.numeric(init_raw) || !is.finite(init_raw) || init_raw >= penalty) {
      diagnostics[[i]] <- data.frame(
        start = i,
        initial_objective = if (length(init_raw) == 1L && is.numeric(init_raw)) as.numeric(init_raw) else NA_real_,
        optimizer = NA_character_, status = "invalid_initial_objective",
        message = "Initial objective is non-finite or numerically invalid.",
        stringsAsFactors = FALSE
      )
      next
    }
    init_obj <- as.numeric(init_raw)

    nl_fit <- tryCatch(
      stats::nlminb(
        start = st,
        objective = eval_objective,
        lower = lower,
        upper = upper,
        control = control
      ),
      error = function(e) e
    )

    nl_ok <- !inherits(nl_fit, "error") && length(nl_fit$objective) == 1L &&
      is.finite(nl_fit$objective) && nl_fit$objective < penalty

    if (nl_ok) {
      fits[[i]] <- standardize_nlminb(nl_fit)
      diagnostics[[i]] <- data.frame(
        start = i, initial_objective = init_obj, optimizer = "nlminb",
        status = "returned", message = nl_fit$message %||% "",
        stringsAsFactors = FALSE
      )
      next
    }

    nl_msg <- if (inherits(nl_fit, "error")) conditionMessage(nl_fit) else
      sprintf("nlminb returned invalid objective: %s", format(nl_fit$objective, digits = 8))

    op_fit <- tryCatch(
      stats::optim(
        par = st,
        fn = eval_objective,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = optim_control
      ),
      error = function(e) e
    )
    op_ok <- !inherits(op_fit, "error") && length(op_fit$value) == 1L &&
      is.finite(op_fit$value) && op_fit$value < penalty

    if (op_ok) {
      fits[[i]] <- standardize_optim(op_fit)
      diagnostics[[i]] <- data.frame(
        start = i, initial_objective = init_obj, optimizer = "L-BFGS-B",
        status = "fallback_returned",
        message = paste0("nlminb failed: ", nl_msg,
                         if (nzchar(op_fit$message %||% "")) paste0("; L-BFGS-B: ", op_fit$message) else ""),
        stringsAsFactors = FALSE
      )
    } else {
      op_msg <- if (inherits(op_fit, "error")) conditionMessage(op_fit) else
        sprintf("L-BFGS-B returned invalid objective: %s", format(op_fit$value, digits = 8))
      diagnostics[[i]] <- data.frame(
        start = i, initial_objective = init_obj, optimizer = "nlminb -> L-BFGS-B",
        status = "both_failed",
        message = paste0("nlminb: ", nl_msg, "; L-BFGS-B: ", op_msg),
        stringsAsFactors = FALSE
      )
    }
  }

  diagnostics <- do.call(rbind, diagnostics)
  ok <- vapply(fits, function(z) {
    is.list(z) && length(z$objective) == 1L && is.finite(z$objective) && z$objective < penalty
  }, logical(1))

  if (!any(ok)) {
    detail <- paste(
      sprintf(
        "start %d: initial objective=%s; optimizer=%s; %s; %s",
        diagnostics$start,
        ifelse(is.finite(diagnostics$initial_objective),
               format(diagnostics$initial_objective, digits = 10), "NA"),
        ifelse(is.na(diagnostics$optimizer), "none", diagnostics$optimizer),
        diagnostics$status,
        diagnostics$message
      ),
      collapse = "\n"
    )
    stop(
      paste0(
        "All direct-likelihood optimization starts failed. Diagnostics:\n",
        detail
      ),
      call. = FALSE
    )
  }

  cand <- fits[ok]
  conv <- vapply(cand, function(z) isTRUE(z$convergence == 0L), logical(1))
  pool <- if (any(conv)) cand[conv] else cand
  best <- pool[[which.min(vapply(pool, `[[`, numeric(1), "objective"))]]
  summary <- do.call(rbind, lapply(which(ok), function(i) {
    z <- fits[[i]]
    data.frame(
      start = i,
      initial_objective = diagnostics$initial_objective[i],
      objective = z$objective,
      convergence = z$convergence,
      iterations = z$iterations %||% NA_integer_,
      optimizer = z$optimizer %||% NA_character_,
      message = z$message %||% "",
      stringsAsFactors = FALSE
    )
  }))
  list(best = best, all = cand, summary = summary, diagnostics = diagnostics)
}

.perturb_start <- function(x, scale = 0.05) x + stats::rnorm(length(x), sd = scale)
