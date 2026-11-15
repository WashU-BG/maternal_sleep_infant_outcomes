library(dplyr)
library(lmerTest)
library(car)

#NOTE: this is done manually so you will need to update covariates, independent variables, and wheter you are running mulitple sensitivty tests
#sensitivity covars:

library(dplyr)
library(lmerTest)
library(car)

PheWAS_baseline <- read.csv("~/Downloads/PheWAS_baseline_Z.csv")

PheWAS_baseline[PheWAS_baseline == "NA"] <- NA
PheWAS_baseline[PheWAS_baseline == ""] <- NA

PheWAS_baseline$child_sex <- as.factor(PheWAS_baseline$child_sex)
PheWAS_baseline$family_id <- as.factor(PheWAS_baseline$family_id)

PheWAS_baseline$study_type <- factor(
  PheWAS_baseline$study_type,
  levels = c("Cuddel", "eLABE")
)

PheWAS_baseline$mom_cannabis_use <- as.factor(PheWAS_baseline$mom_cannabis_use)

PheWAS_baseline$child_birthlength_cm <- as.numeric(PheWAS_baseline$child_birthlength_cm)

PheWAS_baseline <- PheWAS_baseline %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.))))

model_sleepqual3 <- lmer(
  itsea_competence_tscore ~
    sleepqual_t1 + sleepqual_t2 + sleepqual_t3 +
    mom_cannabis_use +
    maternal_age_birth + income_needs_demo +   
    (1 | study_type),
  data = PheWAS_baseline,
  REML = FALSE
)



get_coef_table <- function(mod, terms_keep) {
  sm <- summary(mod)$coefficients
  coef_tab <- as.data.frame(sm)
  coef_tab$term <- rownames(coef_tab)
  
  coef_tab <- coef_tab %>%
    filter(term %in% terms_keep) %>%
    transmute(
      row_type   = "coefficient",
      label      = term,
      estimate   = Estimate,
      std_error  = `Std. Error`,
      conf_low   = Estimate - 1.96 * `Std. Error`,
      conf_high  = Estimate + 1.96 * `Std. Error`,
      statistic  = `t value`,
      df         = df,
      p_value    = `Pr(>|t|)`,
      test       = "t"
    )
  
  coef_tab
}

get_contrast_row <- function(mod, cvec, label, hyp_string, test = "Chisq") {
  beta <- lme4::fixef(mod)
  V <- as.matrix(vcov(mod))
  
  est <- as.numeric(sum(cvec * beta))
  se  <- as.numeric(sqrt(t(cvec) %*% V %*% cvec))
  ci_low  <- est - 1.96 * se
  ci_high <- est + 1.96 * se
  
  lh <- car::linearHypothesis(mod, hyp_string, test = test)
  tab <- as.data.frame(lh)
  out_row <- tab[2, , drop = FALSE]
  
  data.frame(
    row_type   = "contrast",
    label      = label,
    estimate   = est,
    std_error  = se,
    conf_low   = ci_low,
    conf_high  = ci_high,
    statistic  = if ("Chisq" %in% names(out_row)) out_row$Chisq else out_row$F,
    df         = if ("Df" %in% names(out_row)) out_row$Df else out_row$`num Df`,
    p_value    = if ("Pr(>Chisq)" %in% names(out_row)) out_row$`Pr(>Chisq)` else out_row$`Pr(>F)`,
    test       = test,
    stringsAsFactors = FALSE
  )
}

get_joint_row <- function(mod, hyp, label, test = "Chisq") {
  lh <- car::linearHypothesis(mod, hyp, test = test)
  tab <- as.data.frame(lh)
  out_row <- tab[2, , drop = FALSE]
  
  data.frame(
    row_type   = "joint_test",
    label      = label,
    estimate   = NA_real_,
    std_error  = NA_real_,
    conf_low   = NA_real_,
    conf_high  = NA_real_,
    statistic  = if ("Chisq" %in% names(out_row)) out_row$Chisq else out_row$F,
    df         = if ("Df" %in% names(out_row)) out_row$Df else out_row$`num Df`,
    p_value    = if ("Pr(>Chisq)" %in% names(out_row)) out_row$`Pr(>Chisq)` else out_row$`Pr(>F)`,
    test       = test,
    stringsAsFactors = FALSE
  )
}

coef_table <- get_coef_table(
  model_sleepqual3,
  terms_keep = c("sleepqual_t1", "sleepqual_t2", "sleepqual_t3")
)

beta_names <- names(lme4::fixef(model_sleepqual3))

make_cvec <- function(term1, term2, beta_names) {
  v <- rep(0, length(beta_names))
  names(v) <- beta_names
  v[term1] <- 1
  v[term2] <- -1
  v
}

c_t1_t2 <- make_cvec("sleepqual_t1", "sleepqual_t2", beta_names)
c_t1_t3 <- make_cvec("sleepqual_t1", "sleepqual_t3", beta_names)
c_t2_t3 <- make_cvec("sleepqual_t2", "sleepqual_t3", beta_names)

contrast_table <- bind_rows(
  get_contrast_row(
    model_sleepqual3,
    cvec = c_t1_t2,
    label = "sleepqual_t1 - sleepqual_t2",
    hyp_string = "sleepqual_t1 = sleepqual_t2"
  ),
  get_contrast_row(
    model_sleepqual3,
    cvec = c_t1_t3,
    label = "sleepqual_t1 - sleepqual_t3",
    hyp_string = "sleepqual_t1 = sleepqual_t3"
  ),
  get_contrast_row(
    model_sleepqual3,
    cvec = c_t2_t3,
    label = "sleepqual_t2 - sleepqual_t3",
    hyp_string = "sleepqual_t2 = sleepqual_t3"
  )
)

joint_table <- get_joint_row(
  model_sleepqual3,
  hyp = c("sleepqual_t1 = sleepqual_t2", "sleepqual_t1 = sleepqual_t3"),
  label = "sleepqual_t1 = sleepqual_t2 = sleepqual_t3"
)

results_table <- bind_rows(coef_table, contrast_table, joint_table) %>%
  mutate(
    p_adj_bh = p.adjust(p_value, method = "BH")
  )

print(results_table)

