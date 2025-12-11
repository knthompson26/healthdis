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
load("/home/thom1336/healthdis/data/dat_pheno_W2.RData")
load("/home/thom1336/healthdis/data/dat_pheno_W4.RData")

# list of mediators
mediators <- c("family", "school", "friends", "neigh", "religion", "social")

# list of waves and their corresponding data and outcome names
waves <- list(
  W2 = list(data = dat_pheno_W2, outcome = "CESD_W2"),
  W4 = list(data = dat_pheno_W4, outcome = "CESD_W4")
)

# list of racial variables
race_vars <- c("black", "asian", "natam") 

# function to run counterfactual calculations and Monte Carlo simulations within bootstrapping 
bootstrap_stat <- function(data, indices, wave_name, mediators, race_var) {
  
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
  # Step 1: Variable for social (M) when race_var = 0
  #----------------------------------------------------
  
  for (m in mediators) {
    # fit model
    modelM <- lm(as.formula(paste(m, "~", race_var)), data = dat)  # fit on original bootstrapped sample
    
    # store model in list
    model_list_reg[[paste0("model_", m, "_", race_var)]] <- modelM
    
    # extract intercept and residual SD
    intercept <- coef(modelM)[["(Intercept)"]]
    resid_sd <- summary(modelM)$sigma
    
    # simulate mediator under race_var = 0 (on expanded)
    new_var <- paste0(m, "_", race_var, "0")
    dat_expanded[[new_var]] <- intercept + resid_sd * dat_expanded$UM # applied to all expanded simulations
  }
  
  #------------------------------------
  # Step 2: Model for Dep (Y)
  #------------------------------------
  
  walk(mediators, function(m) {
    # construct formula for the model
    formula_str <- paste0(
      outcome_var, " ~ ",                   # y
      race_var, " + ",                      # x
      m, " + ",                             # m
      m, ":", race_var, " + ",              # m*x
      m, "2 + ",                            # m2
      "H1GI1Y + ",                          # c1
      "SES + ",                             # c2
      "sex + ",                             # c3
      "happyP + drinkP_bin + health_bin + violence + ",   # new confounders
      
      # race x covariate
      race_var, ":H1GI1Y + ",               # c1*x
      race_var, ":SES + ",             
      race_var, ":sex + ",                  
      race_var, ":happyP + ",                 
      race_var, ":drinkP_bin + ",               
      race_var, ":health_bin + ",                 
      race_var, ":violence + ",                 
      
      # covariate x covariate 
      "H1GI1Y:SES + ",                     
      "H1GI1Y:sex + ",                     
      "H1GI1Y:happyP + ",                     
      "H1GI1Y:drinkP_bin + ",                     
      "H1GI1Y:health_bin + ",                     
      "H1GI1Y:violence + ",                     
      "SES:sex + ",                          
      "SES:happyP + ",                          
      "SES:drinkP_bin + ",                          
      "SES:health_bin + ",                          
      "SES:violence + ",   
      "sex:happyP + ",
      "sex:drinkP_bin + ",
      "sex:health_bin + ",
      "sex:violence + ",
      "happyP:drinkP_bin + ",
      "happyP:health_bin + ",
      "happyP:violence + ",
      "drinkP_bin:health_bin + ",
      "drinkP_bin:violence + ",
      "health_bin:violence"
    )
    
    formula_obj <- as.formula(formula_str)
    modelY <- lm(formula_obj, data = dat)  # fit outcome model on original sample
    model_nameY <- paste0("model_", outcome_var, "_", race_var, "_", m)
    assign(model_nameY, modelY, envir = .GlobalEnv)
  })
  
  #-----------------------------------------------------------------------
  # Step 2A: Variable for Y, when race_var = 0 and M=m_0, using model from Step 2 
  #----------------------------------------------------------------------
  
  walk(mediators, function(m) {
    model_obj <- get(paste0("model_", outcome_var, "_", race_var, "_", m))
    coefs <- coef(model_obj)
    rmse <- summary(model_obj)$sigma
    race0_var <- sym(paste0(m, "_", race_var, "0"))
    output_var_00 <- sym(paste0(outcome_var, "_00_", m))
    
    dat_expanded <<- dat_expanded %>%  # apply prediction values to all expanded data sets
      mutate(
        !!output_var_00 :=
          # intercept
          (coefs[["(Intercept)"]] %||% 0) +
          
          # main effect of race_var = 0
          (coefs[[race_var]] %||% 0) * 0 +
          
          # mediator terms
          (coefs[[m]] %||% 0)                        * !!race0_var +
          (coefs[[paste0(race_var, ":", m)]] %||% 0) * 0 * !!race0_var +
          (coefs[[paste0(m, "2")]] %||% 0)           * (!!race0_var)^2 +
          
          # covariate main effects
          (coefs[["H1GI1Y"]]    %||% 0) * H1GI1Y +
          (coefs[["SES"]]       %||% 0) * SES +
          (coefs[["sex"]]       %||% 0) * sex +
          (coefs[["happyP"]]    %||% 0) * happyP +
          (coefs[["drinkP_bin"]]%||% 0) * drinkP_bin +
          (coefs[["health_bin"]]%||% 0) * health_bin +
          (coefs[["violence"]]  %||% 0) * violence +
          
          # interactions with race_var (X = 0 → all vanish)
          (coefs[[paste0(race_var, ":H1GI1Y")]]   %||% 0) * H1GI1Y   * 0 +
          (coefs[[paste0(race_var, ":SES")]]      %||% 0) * SES      * 0 +
          (coefs[[paste0(race_var, ":sex")]]      %||% 0) * sex      * 0 +
          (coefs[[paste0(race_var, ":happyP")]]   %||% 0) * happyP   * 0 +
          (coefs[[paste0(race_var, ":drinkP_bin")]] %||% 0) * drinkP_bin * 0 +
          (coefs[[paste0(race_var, ":health_bin")]] %||% 0) * health_bin * 0 +
          (coefs[[paste0(race_var, ":violence")]] %||% 0) * violence * 0 +
          
          # covariate–covariate interactions
          (coefs[["H1GI1Y:SES"]]        %||% 0) * H1GI1Y * SES +
          (coefs[["H1GI1Y:sex"]]        %||% 0) * H1GI1Y * sex +
          (coefs[["H1GI1Y:happyP"]]     %||% 0) * H1GI1Y * happyP +
          (coefs[["H1GI1Y:drinkP_bin"]] %||% 0) * H1GI1Y * drinkP_bin +
          (coefs[["H1GI1Y:health_bin"]] %||% 0) * H1GI1Y * health_bin +
          (coefs[["H1GI1Y:violence"]]   %||% 0) * H1GI1Y * violence +
          
          (coefs[["SES:sex"]]           %||% 0) * SES * sex +
          (coefs[["SES:happyP"]]        %||% 0) * SES * happyP +
          (coefs[["SES:drinkP_bin"]]    %||% 0) * SES * drinkP_bin +
          (coefs[["SES:health_bin"]]    %||% 0) * SES * health_bin +
          (coefs[["SES:violence"]]      %||% 0) * SES * violence +
          
          (coefs[["sex:happyP"]]        %||% 0) * sex * happyP +
          (coefs[["sex:drinkP_bin"]]    %||% 0) * sex * drinkP_bin +
          (coefs[["sex:health_bin"]]    %||% 0) * sex * health_bin +
          (coefs[["sex:violence"]]      %||% 0) * sex * violence +
          
          (coefs[["happyP:drinkP_bin"]] %||% 0) * happyP * drinkP_bin +
          (coefs[["happyP:health_bin"]] %||% 0) * happyP * health_bin +
          (coefs[["happyP:violence"]]   %||% 0) * happyP * violence +
          
          (coefs[["drinkP_bin:health_bin"]] %||% 0) * drinkP_bin * health_bin +
          (coefs[["drinkP_bin:violence"]]   %||% 0) * drinkP_bin * violence +
          
          (coefs[["health_bin:violence"]]   %||% 0) * health_bin * violence +
          
          # residual term 
          rmse * UY
      )
  })
  
  #----------------------------------------------------------------------
  # Step 2B: Variable for Y, when race_var = 1 and M=m_0, using model from Step 2 
  #----------------------------------------------------------------------
  
  walk(mediators, function(m) {
    model_obj <- get(paste0("model_", outcome_var, "_", race_var, "_", m))
    coefs <- coef(model_obj)
    rmse <- summary(model_obj)$sigma
    
    race0_var   <- sym(paste0(m, "_", race_var, "0"))
    output_var_10 <- sym(paste0(outcome_var, "_10_", m))
    
    dat_expanded <<- dat_expanded %>%
      mutate(
        !!output_var_10 :=
          # intercept
          (coefs[["(Intercept)"]] %||% 0) +
          
          # main effect of race_var = 1
          (coefs[[race_var]] %||% 0) * 1 +
          
          # mediator terms
          (coefs[[m]] %||% 0)                        * !!race0_var +
          (coefs[[paste0(race_var, ":", m)]] %||% 0) * 1 * !!race0_var +
          (coefs[[paste0(m, "2")]] %||% 0)           * (!!race0_var)^2 +
          
          # covariate main effects
          (coefs[["H1GI1Y"]]    %||% 0) * H1GI1Y +
          (coefs[["SES"]]       %||% 0) * SES +
          (coefs[["sex"]]       %||% 0) * sex +
          (coefs[["happyP"]]    %||% 0) * happyP +
          (coefs[["drinkP_bin"]]%||% 0) * drinkP_bin +
          (coefs[["health_bin"]]%||% 0) * health_bin +
          (coefs[["violence"]]  %||% 0) * violence +
          
          # interactions with race_var (X = 1)
          (coefs[[paste0(race_var, ":H1GI1Y")]]   %||% 0) * H1GI1Y +
          (coefs[[paste0(race_var, ":SES")]]      %||% 0) * SES +
          (coefs[[paste0(race_var, ":sex")]]      %||% 0) * sex +
          (coefs[[paste0(race_var, ":happyP")]]   %||% 0) * happyP +
          (coefs[[paste0(race_var, ":drinkP_bin")]] %||% 0) * drinkP_bin +
          (coefs[[paste0(race_var, ":health_bin")]] %||% 0) * health_bin +
          (coefs[[paste0(race_var, ":violence")]] %||% 0) * violence +
          
          # covariate–covariate interactions
          (coefs[["H1GI1Y:SES"]]        %||% 0) * H1GI1Y * SES +
          (coefs[["H1GI1Y:sex"]]        %||% 0) * H1GI1Y * sex +
          (coefs[["H1GI1Y:happyP"]]     %||% 0) * H1GI1Y * happyP +
          (coefs[["H1GI1Y:drinkP_bin"]] %||% 0) * H1GI1Y * drinkP_bin +
          (coefs[["H1GI1Y:health_bin"]] %||% 0) * H1GI1Y * health_bin +
          (coefs[["H1GI1Y:violence"]]   %||% 0) * H1GI1Y * violence +
          
          (coefs[["SES:sex"]]           %||% 0) * SES * sex +
          (coefs[["SES:happyP"]]        %||% 0) * SES * happyP +
          (coefs[["SES:drinkP_bin"]]    %||% 0) * SES * drinkP_bin +
          (coefs[["SES:health_bin"]]    %||% 0) * SES * health_bin +
          (coefs[["SES:violence"]]      %||% 0) * SES * violence +
          
          (coefs[["sex:happyP"]]        %||% 0) * sex * happyP +
          (coefs[["sex:drinkP_bin"]]    %||% 0) * sex * drinkP_bin +
          (coefs[["sex:health_bin"]]    %||% 0) * sex * health_bin +
          (coefs[["sex:violence"]]      %||% 0) * sex * violence +
          
          (coefs[["happyP:drinkP_bin"]] %||% 0) * happyP * drinkP_bin +
          (coefs[["happyP:health_bin"]] %||% 0) * happyP * health_bin +
          (coefs[["happyP:violence"]]   %||% 0) * happyP * violence +
          
          (coefs[["drinkP_bin:health_bin"]] %||% 0) * drinkP_bin * health_bin +
          (coefs[["drinkP_bin:violence"]]   %||% 0) * drinkP_bin * violence +
          
          (coefs[["health_bin:violence"]]   %||% 0) * health_bin * violence +
          
          # residual
          rmse * UY
      )
  })
  
  #----------------------------------------------------------------------  
  # Step 3: Model for the Total Causal Effect (TCE), Y ~ X
  #----------------------------------------------------------------------
  
  model <- lm(
    as.formula(
      paste0(
        outcome_var, " ~ ",
        race_var, " + ",
        "H1GI1Y + SES + sex + ",
        race_var, ":H1GI1Y + ",
        race_var, ":SES + ",
        race_var, ":sex + ",
        "H1GI1Y:SES + H1GI1Y:sex + SES:sex"
      )
    ),
    data = dat  # computed on original data
  )
  
  coefs <- coef(model)
  
  #----------------------------------------------------------------------
  # Step 3A: Variables for Dep (Y), when race_var = 0 (no mediator involved) 
  #----------------------------------------------------------------------
  
  output_var_0 <- sym(paste0(outcome_var, "_0")) # when race_var = 0
  dat_expanded <- dat_expanded %>% # applied to expanded data
    mutate(
      !!output_var_0 :=
        (coefs[["(Intercept)"]] %||% 0) +
        (coefs[[race_var]] %||% 0) * 0 +
        (coefs[["H1GI1Y"]] %||% 0) * H1GI1Y +
        (coefs[["SES"]]    %||% 0) * SES +
        (coefs[["sex"]]    %||% 0) * sex +
        (coefs[[paste0(race_var, ":H1GI1Y")]] %||% 0) * H1GI1Y * 0 +
        (coefs[[paste0(race_var, ":SES")]]    %||% 0) * SES     * 0 +
        (coefs[[paste0(race_var, ":sex")]]    %||% 0) * sex     * 0 +
        (coefs[["H1GI1Y:SES"]] %||% 0)  * H1GI1Y * SES +
        (coefs[["H1GI1Y:sex"]] %||% 0)  * H1GI1Y * sex +
        (coefs[["SES:sex"]]    %||% 0)  * SES     * sex
    )
  
  #----------------------------------------------------------------------
  # Step 3B: Variables for Dep (Y), when race_var = 1 (no mediator involved) 
  #----------------------------------------------------------------------
  
  output_var_1 <- sym(paste0(outcome_var, "_1")) # when race_var = 1
  dat_expanded <- dat_expanded %>% # applied to expanded data
    mutate(
      !!output_var_1 :=
        (coefs[["(Intercept)"]] %||% 0) +
        (coefs[[race_var]] %||% 0) +
        (coefs[["H1GI1Y"]] %||% 0) * H1GI1Y +
        (coefs[["SES"]]    %||% 0) * SES +
        (coefs[["sex"]]    %||% 0) * sex +
        (coefs[[paste0(race_var, ":H1GI1Y")]] %||% 0) * H1GI1Y +
        (coefs[[paste0(race_var, ":SES")]]    %||% 0) * SES +
        (coefs[[paste0(race_var, ":sex")]]    %||% 0) * sex +
        (coefs[["H1GI1Y:SES"]] %||% 0)  * H1GI1Y * SES +
        (coefs[["H1GI1Y:sex"]] %||% 0)  * H1GI1Y * sex +
        (coefs[["SES:sex"]]    %||% 0)  * SES     * sex
    )
  
  #--------------------------------------------
  # Step 4: Compute differences: de, tce, dif
  #--------------------------------------------
  
  walk(mediators, function(m) {
    de_var <- sym(paste0("de_", race_var, "_", m))
    tce_var <- sym(paste0("tce_", race_var, "_", m))
    dif_var <- sym(paste0("dif_", race_var, "_", m))
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
      CDM = mean(dat_expanded[[paste0("de_", race_var, "_", m)]], na.rm = TRUE),
      TCE = mean(dat_expanded[[paste0("tce_", race_var, "_", m)]], na.rm = TRUE),
      DIF = mean(dat_expanded[[paste0("dif_", race_var, "_", m)]], na.rm = TRUE)
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
  for (race_var in race_vars) {  # NEW LINE: Nested loop over each race variable
    cat("\nBootstrapping wave:", wave, "for race variable:", race_var, "\n")
    
    # This re-samples individuals 1000 times to capture sampling variability and calculate confidence intervals 
    boot_results <- boot(
      data = waves[[wave]]$data,
      statistic = function(d, i) bootstrap_stat(d, i, wave_name = wave, mediators = mediators, race_var = race_var), # loops through race_var
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
      Race    = race_var  # Add the race variable to the results
    ) %>%
      # Extract mediator number and effect type (CDM/TCE/DIF)
      mutate(
        Mediator_Num = as.numeric(gsub("[^0-9]", "", Effect)),
        Effect_Type  = gsub("[0-9]", "", Effect),
        Mediator     = mediators[Mediator_Num],
        Wave         = wave
      ) %>%
      dplyr::select(Wave, Race, Mediator, Effect_Type, Mean, SE, CI_Low, CI_High) # final table
    
    saveRDS(results, file = paste0("/home/thom1336/healthdis/data/confound/all/race/bootstrap_results_", wave, "_", race_var, ".rds"))
  #  saveRDS(boot_results, file = paste0("/home/thom1336/healthdis/data/confound/all/race/boot_object_", wave, "_", race_var, ".rds"))
    
    cat("Saved results for wave:", wave, "and race variable:", race_var, "\n")
  }
}
