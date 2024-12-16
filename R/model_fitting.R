#' Fit Models and Extract P-values
#' @keywords internal
#' @param dat_rr The sampled dataset
#' @return Numeric p-value or NA
fit_models_and_get_pvalue <- function(dat_rr) {
  # Fit full and null models
  fit_full <- try(
    lme4::glmer(y ~ x1 * minority + (1 | trial_id),
                family = binomial, data = dat_rr),
    silent = TRUE
  )
  fit_null <- try(
    lme4::glmer(y ~ x1 + minority + (1 | trial_id),
                family = binomial, data = dat_rr),
    silent = TRUE
  )

  if (inherits(fit_full, "try-error") ||
      inherits(fit_null, "try-error") ||
      lme4::isSingular(fit_full) ||
      AIC(fit_full) == Inf) {
    return(NA_real_)
  } else {
    anova_res <- anova(fit_null, fit_full)
    p_val <- anova_res$`Pr(>Chisq)`[2]
    return(p_val)
  }
}
