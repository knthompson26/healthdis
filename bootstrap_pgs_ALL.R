library(purrr)
library(tidyr)
library(boot)
library(tidyverse)
library(dplyr)
library(rlang)
library(magrittr)

set.seed(123)
remove(list = ls())

# load data
load("/home/thom1336/healthdis/data/dat_pheno_pgs_W2.RData")
load("/home/thom1336/healthdis/data/dat_pheno_pgs_W4.RData")

# list of mediators
mediators <- c("family", "school", "friends", "neigh", "religion", "social")

# list of waves and their corresponding data and outcome names
waves <- list(
  W2 = list(data = dat_pheno_pgs_W2, outcome = "CESD_W2"),
  W4 = list(data = dat_pheno_pgs_W4, outcome = "CESD_W4")
)

# list of pgs variables
pgs_vars <- c("pgs_reg2", "pgs_reg3") 

# function to run counterfactual calculations and Monte Carlo simulations within bootstrapping 
bootstrap_stat <- function(data, indices, wave_name, mediators, pgs_var) {
  
  #--------------------------------------------
  # Bootstrap sample (original sample)
  #--------------------------------------------
  dat <- data[indices, ]
  dat$original <- 1
  dat$UM <- rnorm(nrow(dat))  # mediator error (original)
  dat$UY <- rnorm(nrow(dat))  # outcome error (original)
  
  outcome_var <- waves[[wave_name]]$outcome
  
  #--------------------------------------------
  # Monte Carlo expansion (inner loop)
  #--------------------------------------------
  dat_expanded <- dat[rep(1:nrow(dat), each = 1000), ]  
  dat_expanded$original <- rep(c(1, rep(0, 999)), times = nrow(dat))
  dat_expanded$UM <- rnorm(nrow(dat_expanded))
  dat_expanded$UY <- rnorm(nrow(dat_expanded))
  
  #--------------------------------------------
  # Counterfactual loop
  #--------------------------------------------
  
  # create a list to store model results
  model_list_reg <- list()
  
  #----------------------------------------------------
  # Step 1: Variable for social (M) when pgs_var = 0
  #----------------------------------------------------
  
  for (m in mediators) {
    # fit model
    modelM <- lm(as.formula(paste(m, "~", pgs_var)), data = dat)  # fit on original bootstrapped sample
    
    # store model in list
    model_list_reg[[paste0("model_", m, "_", pgs_var)]] <- modelM
    
    # extract intercept and residual SD
    intercept <- coef(modelM)[["(Intercept)"]]
    resid_sd <- summary(modelM)$sigma
    
    # simulate mediator under pgs_var = 0 (on expanded)
    new_var <- paste0(m, "_", pgs_var, "0")
    dat_expanded[[new_var]] <- intercept + resid_sd * dat_expanded$UM # applied to all expanded simulations

  }
  
  #------------------------------------
  # Step 2: Model for Dep (Y)
  #------------------------------------
  
  walk(mediators, function(m) {
    # construct formula for the model
    formula_str <- paste0(
      outcome_var, " ~ ",
      pgs_var, " + ", 
      m, " + ", 
      m, ":", pgs_var, " + ", 
      m, "2 + ",
      "H1GI1Y +", pgs_var, ":H1GI1Y"
    )
    
    formula_obj <- as.formula(formula_str)
    modelY <- lm(formula_obj, data = dat)  # fit outcome model on original sample
    model_nameY <- paste0("model_", outcome_var, "_", pgs_var, "_", m)
    assign(model_nameY, modelY, envir = .GlobalEnv)
    
  })
  
  #-----------------------------------------------------------------------
  # Step 2A: Variable for Y, when pgs_var = 0 and M=m_0, using model from Step 2 
  #----------------------------------------------------------------------
  
  walk(mediators, function(m) {
    model_obj <- get(paste0("model_", outcome_var, "_", pgs_var, "_", m))
    coefs <- coef(model_obj)
    rmse <- summary(model_obj)$sigma
    pgs0_var <- sym(paste0(m, "_", pgs_var, "0"))
    output_var_00 <- sym(paste0(outcome_var, "_00_", m))
    
    dat_expanded <<- dat_expanded %>%  # apply prediction values to all expanded data sets
      mutate(
        !!output_var_00 :=
          coefs[["(Intercept)"]] +                                        # intercept
          (coefs[[pgs_var]] %||% 0)           * 0 +                      # main effect of pgs_var = 0
          coefs[[m]]                        * !!pgs0_var +               # main effect of mediator
          (coefs[[paste0(pgs_var, ":", m)]] %||% 0) * 0 +                # interaction of pgs_var × mediator when pgs_var = 0
          coefs[[paste0(m, "2")]]           * (!!pgs0_var)^2 +           # quadratic effect of mediator
          (coefs[["H1GI1Y"]] %||% 0)        * H1GI1Y +                    # H1GI1Y main effect
          (coefs[[paste0(pgs_var, ":H1GI1Y")]] %||% 0)    * H1GI1Y * 0 + # interaction: pgs_var × H1GI1Y when pgs_var = 0
          rmse * UY                                                       # residual
      )
    
    cat("Created:", paste0(outcome_var, "_00_", m), "\n")
  })
  
  #----------------------------------------------------------------------
  # Step 2B: Variable for Y, when pgs_var = 1 and M=m_0, using model from Step 2 
  #----------------------------------------------------------------------
  
  walk(mediators, function(m) {
    model_obj <- get(paste0("model_", outcome_var, "_", pgs_var, "_", m))
    coefs <- coef(model_obj)
    rmse <- summary(model_obj)$sigma
    
    pgs0_var   <- sym(paste0(m, "_", pgs_var, "0"))
    output_var_10 <- sym(paste0(outcome_var, "_10_", m))
    
    dat_expanded <<- dat_expanded %>%  # apply prediction values to all expanded data sets
      mutate(
        !!output_var_10 :=
          coefs[["(Intercept)"]] +                                    # intercept
          coefs[[pgs_var]] +                                          # main effect of pgs_var = 1
          coefs[[m]]             * !!pgs0_var +                       # main effect of mediator
          (coefs[[paste0(pgs_var, ":", m)]] %||% 0) * !!pgs0_var +    # interaction: pgs_var × mediator
          coefs[[paste0(m, "2")]] * (!!pgs0_var)^2 +                  # quadratic effect of mediator
          coefs[["H1GI1Y"]]      * H1GI1Y +                           # H1GI1Y main effect
          coefs[[paste0(pgs_var, ":H1GI1Y")]]  * H1GI1Y +             # interaction: pgs_var × H1GI1Y
          rmse * UY                                                   # residual
      )
    
    cat("Created:", paste0(outcome_var, "_10_", m), "\n")
  })
  
  #----------------------------------------------------------------------  
  # Step 3: Model for the Total Causal Effect (TCE), Y ~ X
  #----------------------------------------------------------------------
  
  model <- lm(as.formula(paste0(outcome_var, " ~ ", pgs_var, " + H1GI1Y +", pgs_var, ":H1GI1Y")), 
              data = dat) # computed on original data
  
  coefs <- coef(model)
  
  #----------------------------------------------------------------------
  # Step 3A: Variables for Dep (Y), when pgs_var = 0 (no mediator involved) 
  #----------------------------------------------------------------------
  
  output_var_0 <- sym(paste0(outcome_var, "_0")) # when pgs_var = 0
  dat_expanded <- dat_expanded %>% # applied to expanded data
    mutate(
      !!output_var_0 :=
        coefs["(Intercept)"] + 
        (coefs[[pgs_var]] %||% 0) * 0 + 
        coefs["H1GI1Y"] * H1GI1Y +
        (coefs[[paste0(pgs_var, ":H1GI1Y")]]) * H1GI1Y * 0
    )
  
  #----------------------------------------------------------------------
  # Step 3B: Variables for Dep (Y), when pgs_var = 1 (no mediator involved) 
  #----------------------------------------------------------------------
  
  output_var_1 <- sym(paste0(outcome_var, "_1")) # when pgs_var = 1
  dat_expanded <- dat_expanded %>% # applied to expanded data
    mutate(
      !!output_var_1 :=
        coefs["(Intercept)"] + 
        coefs[[pgs_var]] + 
        coefs["H1GI1Y"] * H1GI1Y +
        (coefs[[paste0(pgs_var, ":H1GI1Y")]])  * H1GI1Y 
    )
  
  cat("Created:", output_var_0, "\n")
  cat("Created:", output_var_1, "\n")
  
  #--------------------------------------------
  # Step 4: Compute differences: de, tce, dif
  #--------------------------------------------
  
  walk(mediators, function(m) {
    de_var <- sym(paste0("de_", pgs_var, "_", m))
    tce_var <- sym(paste0("tce_", pgs_var, "_", m))
    dif_var <- sym(paste0("dif_", pgs_var, "_", m))
    dep_00 <- sym(paste0(outcome_var, "_00_", m))
    dep_10 <- sym(paste0(outcome_var, "_10_", m))
    dep_0  <- sym(paste0(outcome_var, "_0"))
    dep_1  <- sym(paste0(outcome_var, "_1"))
    
    dat_expanded <<- dat_expanded %>% # calculated for each expanded dataset
      mutate(
        !!de_var := !!dep_10 - !!dep_00,
        !!tce_var := !!dep_1 - !!dep_0,
        !!dif_var := !!tce_var - !!de_var
      )
  })
  
  #-----------------------------------------------------------------
  # Calulate mean de, tce, and dif across all expanded data sets
  #-----------------------------------------------------------------
  
  effects <- map(mediators, function(m) {
    c(
      CDM = mean(dat_expanded[[paste0("de_", pgs_var, "_", m)]], na.rm = TRUE),
      TCE = mean(dat_expanded[[paste0("tce_", pgs_var, "_", m)]], na.rm = TRUE),
      DIF = mean(dat_expanded[[paste0("dif_", pgs_var, "_", m)]], na.rm = TRUE)
    )
  })
  
  names_list <- unlist(map2(mediators, seq_along(mediators),
                            ~ c(paste0("CDM", .y), paste0("TCE", .y), paste0("DIF", .y))))
  out <- unlist(effects)
  names(out) <- names_list
  return(out)
}

#----------------------------------
# Run Bootstrapping (outer loop)
#----------------------------------

for (wave in names(waves)) {
  for (pgs_var in pgs_vars) {  # NEW LINE: Nested loop over each pgs variable
    cat("\nBootstrapping wave:", wave, "for pgs variable:", pgs_var, "\n")
    
    # This re-samples individuals 1000 times to capture sampling variability and calculate confidence intervals 
    boot_results <- boot(
      data = waves[[wave]]$data,
      statistic = function(d, i) bootstrap_stat(d, i, wave_name = wave, mediators = mediators, pgs_var = pgs_var), # loops through pgs_var
      R = 1000
    )
    
    estimates <- boot_results$t
    col_names <- names(boot_results$t0)
    
    # Calculate mean, SE, and percentile confidence intervals
    results <- tibble(
      Effect  = col_names,
      Mean    = boot_results$t0,
      SE      = apply(estimates, 2, sd, na.rm = TRUE),
      CI_Low  = apply(estimates, 2, quantile, probs = 0.025, na.rm = TRUE),
      CI_High = apply(estimates, 2, quantile, probs = 0.975, na.rm = TRUE),
      Wave    = wave,
      pgs    = pgs_var  # Add the pgs variable to the results
    ) %>%
      # Extract mediator number and effect type (CDM/TCE/DIF)
      mutate(
        Mediator_Num = as.numeric(gsub("[^0-9]", "", Effect)),
        Effect_Type  = gsub("[0-9]", "", Effect),
        Mediator     = mediators[Mediator_Num],
        Wave         = wave
      ) %>%
      dplyr::select(Wave, pgs, Mediator, Effect_Type, Mean, SE, CI_Low, CI_High) # final table
    
    saveRDS(results, file = paste0("/home/thom1336/healthdis/data/pgs/bootstrap_results_ALL_", wave, "_", pgs_var, ".rds"))
  #  saveRDS(boot_results, file = paste0("/home/thom1336/healthdis/data/pgs/boot_object_ALL_", wave, "_", pgs_var, ".rds"))
    
    cat("Saved results for wave:", wave, "and pgs variable:", pgs_var, "\n")
  }
}
