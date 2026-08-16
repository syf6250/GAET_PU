# Simple DETM-only wrapper for sourced PUEM authors' functions.
#
# In your main R script, source/load these BEFORE this wrapper:
#   library(glmnet)
#   source("PUEM/R/DETMfull.R")
#   source("PUEM/R/DETMnull.R")
#   source("puem_detm_wrapper_simple.R")
#
# Main function:
#   fit <- puem_detm_analyze(posdata, negdata, pi_init = 0.3,
#                            pihalf = TRUE, pi_input = NULL)
#
# posdata = numeric matrix/data.frame for labeled/source positives.
# negdata = numeric matrix/data.frame for target unlabeled observations.
#
# This file does NOT install or source the authors' code; it assumes
# DETMfull() and DETMnull() already exist in the global environment.

as_numeric_matrix <- function(x, name = "x") {
  if (is.data.frame(x)) {
    bad <- !vapply(x, is.numeric, logical(1))
    if (any(bad)) {
      stop(name, " must contain only numeric columns. Non-numeric columns: ",
           paste(names(x)[bad], collapse = ", "))
    }
    x <- as.matrix(x)
  } else {
    x <- as.matrix(x)
  }
  storage.mode(x) <- "double"
  if (anyNA(x)) stop(name, " contains NA. Impute or remove missing values first.")
  x
}

stable_detm_prob <- function(x, alpha1, beta1, alpha2, beta2, pi) {
  x <- as_numeric_matrix(x, "x")
  eta1 <- as.numeric(alpha1 + x %*% beta1)
  eta2 <- as.numeric(alpha2 + x %*% beta2)
  loga <- log(pi) + eta1
  logb <- log1p(-pi) + eta2
  1 / (1 + exp(logb - loga))
}

swap_detm_components <- function(fit) {
  fit$pival <- 1 - fit$pival
  tmp_alp <- fit$alp1;  fit$alp1 <- fit$alp2;  fit$alp2 <- tmp_alp
  tmp_beta <- fit$beta1; fit$beta1 <- fit$beta2; fit$beta2 <- tmp_beta
  fit$relabeled <- TRUE
  fit
}

relabel_detm_by_pihalf <- function(fit, pihalf = TRUE) {
  fit$relabeled <- FALSE
  if (!is.logical(pihalf) || length(pihalf) != 1 || is.na(pihalf)) {
    stop("pihalf must be TRUE or FALSE.")
  }
  if (is.null(fit$pival) || is.na(fit$pival)) return(fit)

  # pihalf = TRUE: the oriented/true pi is assumed to be < 0.5.
  # pihalf = FALSE: the oriented/true pi is assumed to be > 0.5.
  need_swap <- (pihalf && fit$pival > 0.5) || (!pihalf && fit$pival < 0.5)
  if (need_swap) fit <- swap_detm_components(fit)
  fit
}

transform_detm_to_original_scale <- function(detm, center, scalev) {
  b1 <- detm$beta1 / scalev
  b2 <- detm$beta2 / scalev
  a1 <- detm$alp1 - sum(detm$beta1 * center / scalev)
  a2 <- detm$alp2 - sum(detm$beta2 * center / scalev)
  list(pival = detm$pival, alp1 = a1, beta1 = b1, alp2 = a2, beta2 = b2,
       loglik = detm$loglik, penal_loglik = detm$penal_loglik, itenum = detm$itenum,
       relabeled = isTRUE(detm$relabeled))
}

make_detm_parameter_table <- function(detm, covariate_names) {
  data.frame(
    parameter = c("pival", "alp1", paste0("beta1_", covariate_names),
                  "alp2", paste0("beta2_", covariate_names)),
    estimate = c(detm$pival, detm$alp1, detm$beta1, detm$alp2, detm$beta2),
    row.names = NULL
  )
}

fit_best_detmfull <- function(pos, unl, pi_init, alp1_init, beta1_init, alp2_init, beta2_init,
                              n_starts = 10, seed = NULL, jitter_alpha = 0.5,
                              jitter_beta = 0.5, jitter_pi = 0.15,
                              maxite = 1000, eps = 1e-4) {
  if (!is.null(seed)) set.seed(seed)
  p <- ncol(pos)
  starts <- vector("list", n_starts)
  starts[[1]] <- list(pi = pi_init, alp1 = alp1_init, beta1 = beta1_init,
                      alp2 = alp2_init, beta2 = beta2_init)
  if (n_starts >= 2) {
    for (s in 2:n_starts) {
      pi_s <- pi_init + stats::runif(1, -jitter_pi, jitter_pi)
      pi_s <- min(max(pi_s, 1 / max(nrow(unl), 2)), 1 - 1 / max(nrow(unl), 2))
      starts[[s]] <- list(
        pi = pi_s,
        alp1 = alp1_init + stats::runif(1, -jitter_alpha, jitter_alpha),
        beta1 = beta1_init + stats::runif(p, -jitter_beta, jitter_beta),
        alp2 = alp2_init + stats::runif(1, -jitter_alpha, jitter_alpha),
        beta2 = beta2_init + stats::runif(p, -jitter_beta, jitter_beta)
      )
    }
  }

  best <- list(penal_loglik = -Inf)
  all_fits <- vector("list", n_starts)
  successful <- rep(FALSE, n_starts)
  converged <- rep(FALSE, n_starts)

  for (s in seq_len(n_starts)) {
    st <- starts[[s]]
    fit_s <- tryCatch(
      DETMfull(pos, unl, st$pi, st$alp1, st$beta1, st$alp2, st$beta2,
               maxite = maxite, eps = eps),
      error = function(e) list(penal_loglik = NA_real_, loglik = NA_real_,
                               itenum = NA_integer_, error = conditionMessage(e))
    )

    finite_fit <- !is.null(fit_s$penal_loglik) &&
      length(fit_s$penal_loglik) == 1 &&
      is.finite(fit_s$penal_loglik)

    # Convergence flag used by the wrapper:
    # finite penal_loglik and, if the authors' DETMfull() reports itenum,
    # the algorithm stops before reaching maxite.
    it_ok <- is.null(fit_s$itenum) ||
      (length(fit_s$itenum) == 1 && !is.na(fit_s$itenum) && fit_s$itenum < maxite)

    successful[s] <- finite_fit
    converged[s] <- finite_fit && it_ok
    fit_s$successful <- successful[s]
    fit_s$converged <- converged[s]
    fit_s$start_id <- s
    all_fits[[s]] <- fit_s

    # Choose the best only among converged starts.
    if (converged[s] && fit_s$penal_loglik > best$penal_loglik) {
      best <- fit_s
      best$best_start <- s
    }
  }

  if (!any(converged)) {
    return(list(best = NULL, starts = starts, all_fits = all_fits,
                convergence = FALSE,
                n_successful_starts = sum(successful),
                n_converged_starts = 0,
                successful_by_start = successful,
                convergence_by_start = converged))
  }

  list(best = best, starts = starts, all_fits = all_fits,
       convergence = TRUE,
       n_successful_starts = sum(successful),
       n_converged_starts = sum(converged),
       successful_by_start = successful,
       convergence_by_start = converged)
}

compute_detm_fixed_pi_elr <- function(pos, unl, full_fit, oriented_fit, pi0,
                                      maxite = 1000, eps = 1e-4) {
  if (is.null(pi0) || length(pi0) != 1 || is.na(pi0) || pi0 <= 0 || pi0 >= 1) {
    stop("pi0 must be one number strictly between 0 and 1.")
  }

  null_fit <- tryCatch(
    DETMnull(pos, unl, pi0,
                   oriented_fit$alp1, oriented_fit$beta1,
                   oriented_fit$alp2, oriented_fit$beta2,
                   maxite = maxite, eps = eps),
    error = function(e) list(loglik = NA_real_, penal_loglik = NA_real_, error = conditionMessage(e))
  )

  elr <- if (is.na(null_fit$loglik)) NA_real_ else 2 * (full_fit$loglik - null_fit$loglik)
  list(pi0 = pi0, elr = elr, null_fit = null_fit)
}

compute_detm_pi_profile <- function(pos, unl, full_fit, oriented_fit, ci_grid, ci_level,
                                    maxite = 1000, eps = 1e-4) {
  cutoff <- stats::qchisq(ci_level, df = 1)
  vals <- vapply(ci_grid, function(pi0) {
    tmp <- compute_detm_fixed_pi_elr(pos = pos, unl = unl, full_fit = full_fit,
                                     oriented_fit = oriented_fit, pi0 = pi0,
                                     maxite = maxite, eps = eps)
    tmp$elr
  }, numeric(1))

  profile <- data.frame(pi = ci_grid, elr = vals)
  keep <- ci_grid[!is.na(vals) & vals <= cutoff]
  ci <- if (length(keep) > 0) range(keep) else c(NA_real_, NA_real_)
  list(confidence_interval = ci, profile = profile, cutoff = cutoff, level = ci_level)
}

puem_detm_analyze <- function(posdata,
                              negdata,
                              pi_init = 0.3,
                              pihalf = TRUE,
                              pi_input = NULL,
                              standardize = TRUE,
                              ci_level = 0.95,
                              ci_grid = NULL,
                              ci_grid_size = 201,
                              classification_threshold = 0.5,
                              n_starts = 10,
                              seed = 2024,
                              maxite = 1000,
                              eps = 1e-4) {
  if (!exists("DETMfull", mode = "function") || !exists("DETMnull", mode = "function")) {
    stop("DETMfull() and DETMnull() are not found. In your main script, run: ",
         "source('PUEM/R/DETMfull.R'); source('PUEM/R/DETMnull.R') before sourcing/calling this wrapper.")
  }
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package glmnet is not installed. Install glmnet first.")
  }
  if (!is.numeric(pi_init) || length(pi_init) != 1 || pi_init <= 0 || pi_init >= 1) {
    stop("pi_init must be one number strictly between 0 and 1.")
  }
  if (!is.logical(pihalf) || length(pihalf) != 1 || is.na(pihalf)) {
    stop("pihalf must be TRUE or FALSE.")
  }
  if (!is.null(pi_input) && (!is.numeric(pi_input) || length(pi_input) != 1 ||
                             pi_input <= 0 || pi_input >= 1)) {
    stop("pi_input must be NULL or one number strictly between 0 and 1.")
  }
  if (!is.numeric(n_starts) || n_starts < 1) stop("n_starts must be at least 1.")
  n_starts <- as.integer(n_starts)
  if (!is.null(seed)) set.seed(seed)

  pos0 <- as_numeric_matrix(posdata, "posdata")
  unl0 <- as_numeric_matrix(negdata, "negdata")
  if (ncol(pos0) != ncol(unl0)) stop("posdata and negdata must have the same number of columns.")
  p <- ncol(pos0)
  covariate_names <- colnames(pos0)
  if (is.null(covariate_names)) covariate_names <- paste0("X", seq_len(p))

  pos <- pos0
  unl <- unl0
  center <- rep(0, p)
  scalev <- rep(1, p)
  names(center) <- names(scalev) <- covariate_names

  if (standardize) {
    allx <- rbind(pos, unl)
    center <- colMeans(allx)
    scalev <- apply(allx, 2, stats::sd)
    scalev[is.na(scalev) | scalev == 0] <- 1
    pos <- scale(pos, center = center, scale = scalev)
    unl <- scale(unl, center = center, scale = scalev)
  }
  colnames(pos) <- colnames(unl) <- covariate_names

  # Practical initialization: alpha1=0, beta1=0; alpha2,beta2 from a logistic
  # regression separating source-positive rows from target-unlabeled rows.
  init_dat <- data.frame(y = c(rep(0, nrow(pos)), rep(1, nrow(unl))), rbind(pos, unl))
  init_fit <- stats::glm(y ~ ., data = init_dat, family = stats::binomial())
  cf <- stats::coef(init_fit)
  cf[is.na(cf)] <- 0
  alp1_init <- 0
  beta1_init <- rep(0, p)
  alp2_init <- unname(cf[1])
  beta2_init <- unname(cf[-1])

  detm_full <- fit_best_detmfull(
    pos = pos, unl = unl, pi_init = pi_init,
    alp1_init = alp1_init, beta1_init = beta1_init,
    alp2_init = alp2_init, beta2_init = beta2_init,
    n_starts = n_starts, seed = seed, maxite = maxite, eps = eps
  )

  if (!isTRUE(detm_full$convergence)) {
    out <- list(
      call = match.call(),
      convergence = FALSE,
      input = list(n_source_positive = nrow(pos),
                   m_target_unlabeled = nrow(unl),
                   p = p,
                   covariate_names = covariate_names,
                   standardized = standardize,
                   center = center,
                   scale = scalev,
                   classification_threshold = classification_threshold,
                   n_starts = n_starts,
                   pihalf = pihalf,
                   pi_input = pi_input),
      init = list(pi_init = pi_init,
                  alp1 = alp1_init,
                  beta1 = beta1_init,
                  alp2 = alp2_init,
                  beta2 = beta2_init),
      convergence_detail = list(
        n_successful_starts = detm_full$n_successful_starts,
        n_converged_starts = detm_full$n_converged_starts,
        successful_by_start = detm_full$successful_by_start,
        convergence_by_start = detm_full$convergence_by_start,
        rule = "finite penal_loglik and itenum < maxite"
      ),
      pi_estimation = NA_real_,
      pi_confidence_interval = list(
        ci = c(NA_real_, NA_real_),
        level = ci_level,
        cutoff = stats::qchisq(ci_level, df = 1),
        method = "not computed because DETMfull did not converge"
      ),
      pi_profile = NULL,
      pi_input_check = if (is.null(pi_input)) NULL else list(
        pi_input = pi_input,
        covered_by_CI = NA,
        elr = NA_real_,
        cutoff = stats::qchisq(ci_level, df = 1),
        level = ci_level,
        method = "not computed because DETMfull did not converge"
      ),
      parameter_estimates = NULL,
      predictions_unlabeled = NULL,
      detm = NULL,
      detm_raw_before_optional_relabeling = NULL,
      all_detm_starts = detm_full$all_fits,
      start_values = detm_full$starts,
      predict = function(newdata, threshold = classification_threshold) {
        stop("Prediction is unavailable because DETMfull did not converge.")
      }
    )
    class(out) <- "puem_detm_analysis"
    return(out)
  }

  detm_raw <- detm_full$best
  detm <- relabel_detm_by_pihalf(detm_raw, pihalf = pihalf)

  # Construct default CI grid on the oriented side determined by pihalf.
  if (is.null(ci_grid)) {
    eps_pi <- 1 / max(nrow(unl), 2)
    if (pihalf) {
      lo <- max(eps_pi, detm$pival - 0.30)
      hi <- min(0.5 - eps_pi, detm$pival + 0.30)
    } else {
      lo <- max(0.5 + eps_pi, detm$pival - 0.30)
      hi <- min(1 - eps_pi, detm$pival + 0.30)
    }
    if (lo >= hi) {
      lo <- if (pihalf) eps_pi else 0.5 + eps_pi
      hi <- if (pihalf) 0.5 - eps_pi else 1 - eps_pi
    }
    ci_grid <- seq(lo, hi, length.out = ci_grid_size)
  }
  ci_grid <- sort(unique(as.numeric(ci_grid)))
  ci_grid <- ci_grid[ci_grid > 0 & ci_grid < 1]
  if (pihalf) ci_grid <- ci_grid[ci_grid < 0.5]
  if (!pihalf) ci_grid <- ci_grid[ci_grid > 0.5]
  if (length(ci_grid) == 0) stop("ci_grid has no values on the side implied by pihalf.")

  pi_prof <- compute_detm_pi_profile(
    pos = pos, unl = unl, full_fit = detm_raw,
    oriented_fit = detm, ci_grid = ci_grid, ci_level = ci_level,
    maxite = maxite, eps = eps
  )

  pi_input_result <- NULL
  if (!is.null(pi_input)) {
    side_ok <- if (pihalf) pi_input < 0.5 else pi_input > 0.5
    fixed_pi <- compute_detm_fixed_pi_elr(
      pos = pos, unl = unl, full_fit = detm_raw,
      oriented_fit = detm, pi0 = pi_input,
      maxite = maxite, eps = eps
    )
    cutoff <- stats::qchisq(ci_level, df = 1)
    pi_input_result <- list(
      pi_input = pi_input,
      side_assumption = if (pihalf) "pi < 0.5" else "pi > 0.5",
      side_ok = side_ok,
      elr = fixed_pi$elr,
      cutoff = cutoff,
      level = ci_level,
      covered_by_CI = if (is.na(fixed_pi$elr)) NA else as.logical(fixed_pi$elr < cutoff),
      method = "Authors' Figure 2 coverage check: 2 * (DETMfull loglik - DETMnull fixed-pi loglik) < qchisq(level, 1)",
      null_fit = fixed_pi$null_fit
    )
  }

  predict_fun <- function(newdata, threshold = classification_threshold) {
    x <- as_numeric_matrix(newdata, "newdata")
    if (ncol(x) != p) stop("newdata must have the same number of columns as the training data.")
    if (standardize) x <- scale(x, center = center, scale = scalev)
    colnames(x) <- covariate_names
    pr <- stable_detm_prob(x, detm$alp1, detm$beta1, detm$alp2, detm$beta2, detm$pival)
    data.frame(row_id = seq_len(nrow(x)),
               prob_target_positive = pr,
               predicted_label = as.integer(pr >= threshold),
               threshold = threshold,
               model = "DETM")
  }

  predictions_unlabeled <- predict_fun(unl0, threshold = classification_threshold)

  detm_original <- if (standardize) transform_detm_to_original_scale(detm, center, scalev) else detm

  parameter_estimates <- list(
    standardized = make_detm_parameter_table(detm, covariate_names),
    original_scale = make_detm_parameter_table(detm_original, covariate_names)
  )

  out <- list(
    call = match.call(),
    convergence = TRUE,
    convergence_detail = list(
      n_successful_starts = detm_full$n_successful_starts,
      n_converged_starts = detm_full$n_converged_starts,
      best_start = detm_raw$best_start,
      successful_by_start = detm_full$successful_by_start,
      convergence_by_start = detm_full$convergence_by_start,
      rule = "finite penal_loglik and itenum < maxite"
    ),
    input = list(n_source_positive = nrow(pos),
                 m_target_unlabeled = nrow(unl),
                 p = p,
                 covariate_names = covariate_names,
                 standardized = standardize,
                 center = center,
                 scale = scalev,
                 classification_threshold = classification_threshold,
                 n_starts = n_starts,
                 pihalf = pihalf,
                 pi_input = pi_input),
    init = list(pi_init = pi_init,
                alp1 = alp1_init,
                beta1 = beta1_init,
                alp2 = alp2_init,
                beta2 = beta2_init),
    label_switching = list(
      pihalf = pihalf,
      side_assumption = if (pihalf) "pi < 0.5" else "pi > 0.5",
      raw_pi = detm_raw$pival,
      oriented_pi = detm$pival,
      switched = isTRUE(detm$relabeled)
    ),
    pi_estimation = detm$pival,
    pi_confidence_interval = list(
      ci = pi_prof$confidence_interval,
      level = ci_level,
      cutoff = pi_prof$cutoff,
      method = "DETM profile empirical likelihood ratio using DETMnull over ci_grid on the side implied by pihalf"
    ),
    pi_profile = pi_prof$profile,
    pi_input_check = pi_input_result,
    parameter_estimates = parameter_estimates,
    predictions_unlabeled = predictions_unlabeled,
    detm = detm,
    detm_raw_before_optional_relabeling = detm_raw,
    all_detm_starts = detm_full$all_fits,
    start_values = detm_full$starts,
    predict = predict_fun
  )
  class(out) <- "puem_detm_analysis"
  out
}

print.puem_detm_analysis <- function(x, ...) {
  cat("PUEM DETM-only analysis\n")
  cat("  n source positive: ", x$input$n_source_positive, "\n", sep = "")
  cat("  m target unlabeled: ", x$input$m_target_unlabeled, "\n", sep = "")
  cat("  p covariates: ", x$input$p, "\n", sep = "")
  cat("  n DETM starts: ", x$input$n_starts, "\n", sep = "")
  cat("  convergence: ", x$convergence, "\n", sep = "")
  if (!isTRUE(x$convergence)) {
    cat("  n converged starts: ", x$convergence_detail$n_converged_starts, "\n", sep = "")
    return(invisible(x))
  }
  cat("  pihalf assumption: ", if (x$input$pihalf) "pi < 0.5" else "pi > 0.5", "\n", sep = "")
  cat("  raw pi_hat: ", signif(x$label_switching$raw_pi, 5), "\n", sep = "")
  cat("  oriented DETM pi_hat: ", signif(x$pi_estimation, 5), "\n", sep = "")
  cat("  label switched: ", x$label_switching$switched, "\n", sep = "")
  cat("  DETM ", 100 * x$pi_confidence_interval$level, "% CI for pi: ",
      paste(signif(x$pi_confidence_interval$ci, 5), collapse = " to "), "\n", sep = "")
  if (!is.null(x$pi_input_check)) {
    cat("  pi_input: ", signif(x$pi_input_check$pi_input, 5),
        "; ELR: ", signif(x$pi_input_check$elr, 5),
        "; covered: ", x$pi_input_check$covered_by_CI, "\n", sep = "")
  }
  invisible(x)
}
