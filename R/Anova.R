
# Mapping between mod$survreg and the fitting function
.survstan_fit_fun <- function(survreg) {
  switch(survreg,
    aft = survstan::aftreg,
    ph  = survstan::phreg,
    po  = survstan::poreg,
    ah  = survstan::ahreg,
    yp  = survstan::ypreg,
    eh  = survstan::ehreg,
    stop("Unknown survstan regression type: ", survreg)
  )
}

# Helper: extract term names from model formula (without intercept)
.survstan_term_names <- function(mod) {
  tt <- terms(mod$formula)
  nms <- attr(tt, "term.labels")
  nms
}

# Helper: returns the relative (lower-order) terms of a given term
.relatives <- function(term, names, factors) {
  is.relative <- function(term1, term2) {
    all(!(factors[, term1] & (!factors[, term2])))
  }
  if (length(names) == 1) return(NULL)
  which.term <- which(term == names)
  (seq_along(names))[-which.term][
    sapply(names[-which.term], function(term2) is.relative(term, term2))
  ]
}

#' Anova method for survstan models (Type II and III Likelihood Ratio Tests)
#'
#' @aliases Anova.survstan
#' @description Computes Type II or Type III likelihood ratio tests (LRT) for
#'   parametric survival regression models fitted with the \pkg{survstan} package.
#'   Requires the \pkg{car} package.
#' @importFrom stats terms model.matrix pchisq model.response
#' @importFrom car Anova
#' @export
#' @param mod a fitted model of class \code{survstan}.
#' @param type type of test: \code{"II"} (default) or \code{"III"} (also accepts \code{2} or \code{3}).
#' @param test.statistic test statistic to use: \code{"LR"} (likelihood ratio, default) or \code{"Wald"}.
#' @param ... further arguments passed to or from other methods.
#' @return an object of class \code{c("anova", "data.frame")}.
#' @examples
#' \donttest{
#' library(survstan)
#' data(ovarian, package = "survival")
#' ovarian$rx <- as.factor(ovarian$rx)
#' fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
#' car::Anova(fit, type = "II")
#' car::Anova(fit, type = "III")
#' }
Anova.survstan <- function(mod, type = c("II", "III", 2, 3),
                            test.statistic = c("LR", "Wald"), ...) {
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)

  switch(type,
    II  = switch(test.statistic,
                 LR   = .Anova_II_LR_survstan(mod, ...),
                 Wald = car::Anova.default(mod, type = "II",
                                           test.statistic = "Chisq",
                                           vcov. = vcov(mod), ...)),
    III = switch(test.statistic,
                 LR   = .Anova_III_LR_survstan(mod, ...),
                 Wald = car::Anova.default(mod, type = "III",
                                           test.statistic = "Chisq",
                                           vcov. = vcov(mod), ...)),
    "2" = switch(test.statistic,
                 LR   = .Anova_II_LR_survstan(mod, ...),
                 Wald = car::Anova.default(mod, type = "II",
                                           test.statistic = "Chisq",
                                           vcov. = vcov(mod), ...)),
    "3" = switch(test.statistic,
                 LR   = .Anova_III_LR_survstan(mod, ...),
                 Wald = car::Anova.default(mod, type = "III",
                                           test.statistic = "Chisq",
                                           vcov. = vcov(mod), ...))
  )
}

# ── Type II LR ────────────────────────────────────────────────────────────────

.Anova_II_LR_survstan <- function(mod, ...) {
  fit_fun  <- .survstan_fit_fun(mod$survreg)
  baseline <- mod$baseline
  data     <- mod$mf                          # model frame stored in the fitted object
  y        <- stats::model.response(data)     # Surv object

  fac   <- attr(terms(mod$formula), "factors")
  names <- .survstan_term_names(mod)
  n.terms <- length(names)

  X_full <- stats::model.matrix(mod$formula, data = data)
  asgn_full <- attr(X_full, "assign")        # 0 for intercept, 1..k for terms

  # Strip intercept column so indices align with asgn values
  X    <- X_full[, asgn_full != 0, drop = FALSE]
  asgn <- asgn_full[asgn_full != 0]

  which.nms <- function(name) which(asgn == which(names == name))

  df <- df_terms_survstan(mod)
  LR <- p <- rep(0, n.terms)

  loglik_full <- mod$loglik

  for (term in seq_len(n.terms)) {
    rels      <- names[.relatives(names[term], names, fac)]
    exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))

    # Model 1: without the term AND its relatives (numerator)
    X1   <- X[, -exclude.1, drop = FALSE]
    if (ncol(X1) == 0) {
      mod1 <- fit_fun(y ~ 1, data = data.frame(y = y), baseline = baseline, init = 0)
    } else {
      mod1 <- fit_fun(y ~ X1 - 1, data = data.frame(y = y, X1),
                      baseline = baseline, init = 0)
    }
    ll1  <- mod1$loglik

    # Model 2: without only the relatives (denominator); if no relatives, use full model
    if (length(rels) == 0) {
      ll2 <- loglik_full
    } else {
      exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
      X2   <- X[, -exclude.2, drop = FALSE]
      if (ncol(X2) == 0) {
        mod2 <- fit_fun(y ~ 1, data = data.frame(y = y), baseline = baseline, init = 0)
      } else {
        mod2 <- fit_fun(y ~ X2 - 1, data = data.frame(y = y, X2),
                        baseline = baseline, init = 0)
      }
      ll2 <- mod2$loglik
    }

    LR[term] <- -2 * (ll1 - ll2)
    p[term]  <- stats::pchisq(LR[term], df[term], lower.tail = FALSE)
  }

  result <- data.frame(LR, df, p)
  row.names(result) <- names
  names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c(
    paste0("Analysis of Deviance Table (Type II tests)\n"),
    paste0("Model: ", mod$baseline, "(", mod$survreg, ")"),
    paste0("Response: ", deparse(mod$formula[[2]]))
  )
  result
}

# ── Type III LR ───────────────────────────────────────────────────────────────

.Anova_III_LR_survstan <- function(mod, ...) {
  fit_fun  <- .survstan_fit_fun(mod$survreg)
  baseline <- mod$baseline
  data     <- mod$mf
  y        <- stats::model.response(data)

  names   <- .survstan_term_names(mod)
  n.terms <- length(names)

  X_full <- stats::model.matrix(mod$formula, data = data)
  asgn_full <- attr(X_full, "assign")

  # Strip intercept column so indices align with asgn values
  X    <- X_full[, asgn_full != 0, drop = FALSE]
  asgn <- asgn_full[asgn_full != 0]

  which.nms <- function(name) which(asgn == which(names == name))

  df <- df_terms_survstan(mod)
  LR <- p <- rep(0, n.terms)

  loglik_full <- mod$loglik

  for (term in seq_len(n.terms)) {
    exclude <- which.nms(names[term])
    X0      <- X[, -exclude, drop = FALSE]
    if (ncol(X0) == 0) {
      mod0 <- fit_fun(y ~ 1, data = data.frame(y = y), baseline = baseline, init = 0)
    } else {
      mod0 <- fit_fun(y ~ X0 - 1, data = data.frame(y = y, X0),
                      baseline = baseline, init = 0)
    }
    ll0     <- mod0$loglik

    LR[term] <- -2 * (ll0 - loglik_full)
    p[term]  <- stats::pchisq(LR[term], df[term], lower.tail = FALSE)
  }

  result <- data.frame(LR, df, p)
  row.names(result) <- names
  names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c(
    paste0("Analysis of Deviance Table (Type III tests)\n"),
    paste0("Model: ", mod$baseline, "(", mod$survreg, ")"),
    paste0("Response: ", deparse(mod$formula[[2]]))
  )
  result
}

# ── Helper: degrees of freedom per term ───────────────────────────────────────

#' Degrees of freedom per term for a survstan model
#'
#' @param mod a fitted model of class \code{survstan}.
#' @return a named integer vector of degrees of freedom for each term.
#' @keywords internal
df_terms_survstan <- function(mod) {
  X    <- stats::model.matrix(mod$formula, data = mod$mf)
  asgn <- attr(X, "assign")
  nms  <- .survstan_term_names(mod)
  # count columns assigned to each term (excluding intercept, assign == 0)
  asgn_no_int  <- asgn[asgn != 0]
  term_indices <- seq_along(nms)  # term 1, 2, ...
  df <- sapply(term_indices, function(i) sum(asgn_no_int == i))
  names(df) <- nms
  df
}
