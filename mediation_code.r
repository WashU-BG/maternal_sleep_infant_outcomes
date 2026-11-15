library(lavaan)
library(dplyr)

set.seed(68)


# EDIT ONLY THIS BLOCK EACH TIME- currently for itsea but can be easily modified to switch to bayley (m)vcars will change

x_var <- "whiirs_t3"
y_var <- "itsea_internal_tscore"
run_label <- "itsea_t1_agechecl00"

m_covars <- c( #this is for brain
  "whiirs_t1",
  "whiirs_t2",
  "mrinotes_ga",
  "mom_cannabis_use",
  "maternal_age_birth",
  "income_needs_demo",
  "child_sex",
  "study_type"
)

y_covars <- c( 
  "whiirs_t1",
  "whiirs_t2",
  "mom_cannabis_use",
  "maternal_age_birth",
  "income_needs_demo",
  "child_sex",
  "study_type"
)
#loop for all imaging based on modality- remember that mri_notes_ga is only for all non vol modalities 
mediators_to_run <- unique(c(
  "intracranial_volume_volume", "subcortical_gray_matter_volume", "supratentorial_not_ventricle_volume", "total_brain_volume_volume", "total_cerebral_white_matter_volume", 
  "total_cortical_gray_matter_volume"
  
  
))

PheWAS_baseline <- read.csv(
  "~/Downloads/PheWAS_baseline_with_ALL_totals_v12.csv",
  stringsAsFactors = FALSE
)

PheWAS_baseline[PheWAS_baseline == "NA"] <- NA
PheWAS_baseline[PheWAS_baseline == ""] <- NA

#recode if needed
PheWAS_baseline$study_type <- dplyr::case_when(
  PheWAS_baseline$study_type == "Cuddel" ~ 0,
  PheWAS_baseline$study_type == "eLABE" ~ 1,
  TRUE ~ NA_real_
)

#improve run time if not using HTCF
all_needed <- unique(c(x_var, y_var, mediators_to_run, m_covars, y_covars))
missing_vars <- setdiff(all_needed, names(PheWAS_baseline))

if (length(missing_vars) > 0) {
  warning("Missing variables: ", paste(missing_vars, collapse = ", "))
}

all_needed_present <- intersect(all_needed, names(PheWAS_baseline))
dat <- PheWAS_baseline[, all_needed_present, drop = FALSE]


dat <- dat %>%
  mutate(across(
    everything(),
    ~ {
      if (is.character(.x)) {
        suppressWarnings(type.convert(.x, as.is = TRUE))
      } else {
        .x
      }
    }
  ))


vars_to_scale <- unique(c(
  x_var, y_var, mediators_to_run, #double check these!! sleep vars would be removed for cannabis anlayses
  "whiirs_t1", "whiirs_t2",
  "mrinotes_ga",
  "maternal_age_birth", "income_needs_demo"
))

vars_to_scale <- intersect(vars_to_scale, names(dat))
vars_to_scale <- vars_to_scale[sapply(dat[vars_to_scale], is.numeric)]

dat <- dat %>%
  mutate(across(all_of(vars_to_scale), ~ as.numeric(scale(.))))


run_one_mediation <- function(data, x, m, y, m_covars, y_covars) {
  
  needed <- unique(c(x, m, y, m_covars, y_covars))
  needed <- needed[needed %in% names(data)]
  
  d <- data[, needed, drop = FALSE]
  #sanity check 
  if (!(x %in% names(d))) stop(paste("Predictor missing:", x))
  if (!(m %in% names(d))) stop(paste("Mediator missing:", m))
  if (!(y %in% names(d))) stop(paste("Outcome missing:", y))
  
  if (!is.numeric(d[[x]]) || !is.numeric(d[[m]]) || !is.numeric(d[[y]])) {
    stop("x, mediator, and y must be numeric for this ML/FIML model.")
  }
  
  m_rhs <- paste(c(paste0("a*", x), m_covars[m_covars %in% names(d)]), collapse = " + ")
  y_rhs <- paste(c(paste0("b*", m), paste0("cprime*", x), y_covars[y_covars %in% names(d)]), collapse = " + ")
  
  med_model <- paste0(
    m, " ~ ", m_rhs, "\n",
    y, " ~ ", y_rhs, "\n",
    "indirect := a*b\n",
    "total := cprime + indirect\n"
  )
  
  fit <- tryCatch(
    sem(
      med_model,
      data = d,
      se = "bootstrap",
      bootstrap = 5000,
      missing = "fiml",
      estimator = "ML",
      fixed.x = FALSE,
      em.h1.iter.max = 10000
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(data.frame(
      predictor = x,
      mediator = m,
      outcome = y,
      param = NA,
      lhs = NA,
      op = NA,
      rhs = NA,
      est = NA,
      se = NA,
      z = NA,
      pvalue = NA,
      ci.lower = NA,
      ci.upper = NA,
      std.all = NA,
      n_obs = NA,
      converged = FALSE,
      error_msg = fit$message
    ))
  }
  
  pe <- parameterEstimates(
    fit,
    ci = TRUE,
    boot.ci.type = "bca.simple",
    standardized = TRUE
  )
  
  keep <- (pe$label %in% c("a", "b", "cprime")) |
    (pe$lhs %in% c("indirect", "total") & pe$op == ":=")
  
  out <- pe[keep, ] %>%
    mutate(
      param = ifelse(lhs %in% c("indirect", "total"), lhs, label),
      predictor = x,
      mediator = m,
      outcome = y,
      n_obs = lavInspect(fit, "nobs"),
      converged = lavInspect(fit, "converged"),
      error_msg = NA_character_
    ) %>%
    select(
      predictor, mediator, outcome,
      param, lhs, op, rhs,
      est, se, z, pvalue,
      ci.lower, ci.upper, std.all,
      n_obs, converged, error_msg
    )
  
  out
}


all_results <- vector("list", length(mediators_to_run))

for (i in seq_along(mediators_to_run)) {
  m_var <- mediators_to_run[i]
  
  cat("Running", i, "of", length(mediators_to_run), ":", m_var, "\n")
  
  all_results[[i]] <- tryCatch(
    run_one_mediation(
      data = dat,
      x = x_var,
      m = m_var,
      y = y_var,
      m_covars = m_covars,
      y_covars = y_covars
    ),
    error = function(e) {
      data.frame(
        predictor = x_var,
        mediator = m_var,
        outcome = y_var,
        param = NA,
        lhs = NA,
        op = NA,
        rhs = NA,
        est = NA,
        se = NA,
        z = NA,
        pvalue = NA,
        ci.lower = NA,
        ci.upper = NA,
        std.all = NA,
        n_obs = NA,
        converged = FALSE,
        error_msg = e$message
      )
    }
  )
}

final_results <- bind_rows(all_results)


timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S") # this will time stamp ver number as i sometimes forget to rename files

out_file <- file.path(
  "~/Downloads",
  paste0("mediation_results_", run_label, "_", timestamp, ".csv")
)

write.csv(final_results, out_file, row.names = FALSE)

