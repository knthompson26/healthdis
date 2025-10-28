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
load("/home/thom1336/healthdis/data/dat_pheno_W2_bin.RData")
load("/home/thom1336/healthdis/data/dat_pheno_W4_bin.RData")

# list of mediators - all social items 
sc_items_all <- c(
  # school
  "H1ED19",  # You feel close to the people at your school
  "H1ED20",  # You feel like you are part of your school
  "H1ED21",  # Students at your school are prejudiced
  "H1ED22",  # You are happy to be at your school
  "H1ED23",  # The teachers at your school treat students fairly
  "H1ED24",  # You feel safe in your school
  "H1PR2",   # How much do you feel that your teachers care about you
  
  # family
  "H1WP9",   # How close do you feel to your mother
  "H1WP10",  # How much do you think your mother cares about you
  "H1WP13",  # How close do you feel to your father
  "H1WP14",  # How much do you think your father cares about you
  "H1PF1",   # Most of the time, your mother is warm and loving towards you
  "H1PF4",   # You are satisfied with the way your mother and you communicate with each other
  "H1PF5",   # Overall, you are satisfied with your relationship with your mother
  "H1PR3",   # How much do you feel that your parents care about you
  "H1PR5",   # How much do you feel that people in your family understand you
  
  # friends
  "H16A",    # Did you go to X house during the past seven days - male/female friend 1
  "H17A",    # Did you meet X after school to hang out or go somewhere during the past seven days - male/female friend 1
  "H18A",    # Did you spend time with X during the past weekend - male/female friend 1
  "H19A",    # Did you talk to X about a problem during the past seven days - male/female friend 1
  "H110A",   # Did you talk to X on the telephone during the past seven days - male/female friend 1
  "H1RR1",   # Have you had a special romantic relationship with any one
  "H1DA7",   # During the past week, how many times did you just hang out with friends
  "H1PR4",   # How much do you feel that your friends care about you
  
  # neighbourhood
  "H1NB1",   # You know most of the people in your neighbourhood
  "H1NB2",   # In the past month, you have stopped on the street to talk with someone who lives in your neighbourhood
  "H1NB3",   # People in this neighbourhood look out for each other
  "H1RE7",   # Many places of worship have activities for teenagers, in the past 12 months, how often did you attend such youth activities
  
  # general social support and inclusion
  "H1PF35",  # You feel socially accepted
  "H1PF36",  # You feel loved and wanted
  "H1PR1"    # How much do you feel that adults care about you
)

mediators <- sc_items_all

# list of waves and their corresponding data and outcome names
waves <- list(
  W2 = list(data = dat_pheno_W2_bin, outcome = "CESD_W2"),
  W4 = list(data = dat_pheno_W4_bin, outcome = "CESD_W4")
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
 # dat$UM <- rnorm(nrow(dat))  # mediator error (original) - removed
  dat$UY <- rnorm(nrow(dat))  # outcome error (original)
  
  outcome_var <- waves[[wave_name]]$outcome
  
  #--------------------------------------------
  # Monte Carlo expansion (inner loop)
  #--------------------------------------------
  dat_expanded <- dat[rep(1:nrow(dat), each = 1000), ]  
  dat_expanded$original <- rep(c(1, rep(0, 999)), times = nrow(dat))
 # dat_expanded$UM <- rnorm(nrow(dat_expanded)) - removed for binary 
  dat_expanded$UY <- rnorm(nrow(dat_expanded))
  
  #--------------------------------------------
  # Counterfactual loop
  #--------------------------------------------
  
  # create a list to store model results
  model_list_reg <- list()
  
  #----------------------------------------------------
  # Step 1: Variable for BINARY social (M) when race_var = 0
  #----------------------------------------------------
  
  for (m in mediators) {
    
    # model fit - BINARY so logistic model 
    fitM <- try(                                          # tries to run without crashing whole script
      glm(as.formula(paste(m, "~", race_var)),            # build glm formula "mediator ~ race"
          data = dat,                                     # fit using the bootstrapped dataset
          family = binomial(),                            # specify logistic regression (binomial with logit link)
          control = glm.control(maxit = 50)),             # give more room to try more iterations              
      silent = TRUE)                                      # suppress warnings so code won’t crash mid-bootstrap
    
    # check if glm() gave an error or did not converge
    if (inherits(fitM, "try-error")) {
      message(paste(Sys.time(), "Logistic model failed for mediator:", m))
    } else if (!fitM$converged) {
      message(paste(Sys.time(), "Logistic model did NOT converge for mediator:", m))
    }
    
    # predict P(M = 1 | sex = 0) on expanded data
    new0 <- dat_expanded                                 # copy the expanded Monte-Carlo dataset
    new0[[race_var]] <- 0                                 # set exposure = 0 for everyone (counterfactual condition) - make sure calling var [[]]
    p0 <- as.numeric(                                    # convert predicted probabilities to numeric vector
      predict(fitM, newdata = new0, type = "response")   # logistic inverse link -> P(M=1|sex=0)
    )
    
    # simulate mediator draws under sex=0
    new_var <- paste0(m, "_", race_var, "0")             # name of the new simulated mediator variable
    dat_expanded[[new_var]] <- rbinom(                   # simulate binary outcomes from Bernoulli(p0)
      n = length(p0),                                    # number of draws = rows in expanded data
      size = 1,                                          # each trial produces a single 0/1 draw
      prob = p0)                                         # individual-specific success probability from logistic model
  }  
  
  #------------------------------------
  # Step 2: Model for Dep (Y)
  #------------------------------------
  
  walk(mediators, function(m) {
    # construct formula for the model
    formula_str <- paste0(
      outcome_var, " ~ ",
      race_var, " + ", 
      m, " + ", 
      m, ":", race_var, " + ", 
    #  m, "2 + ",
      "H1GI1Y +", race_var, ":H1GI1Y"
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
          coefs[["(Intercept)"]] +                                        # intercept
          (coefs[[race_var]] %||% 0)           * 0 +                      # main effect of race_var = 0
          coefs[[m]]                        * !!race0_var +               # main effect of mediator
          (coefs[[paste0(race_var, ":", m)]] %||% 0) * 0 +                # interaction of race_var × mediator when race_var = 0
       #   coefs[[paste0(m, "2")]]           * (!!race0_var)^2 +           # quadratic effect of mediator
          (coefs[["H1GI1Y"]] %||% 0)        * H1GI1Y +                    # H1GI1Y main effect
          (coefs[[paste0(race_var, ":H1GI1Y")]] %||% 0)    * H1GI1Y * 0 + # interaction: race_var × H1GI1Y when race_var = 0
          rmse * UY                                                       # residual
      )
    
    cat("Created:", paste0(outcome_var, "_00_", m), "\n")
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
    
    dat_expanded <<- dat_expanded %>%  # apply prediction values to all expanded data sets
      mutate(
        !!output_var_10 :=
          coefs[["(Intercept)"]] +                                  # intercept
          coefs[[race_var]] +                                       # main effect of race_var = 1
          coefs[[m]]             * !!race0_var +                    # main effect of mediator
          (coefs[[paste0(race_var, ":", m)]] %||% 0) * !!race0_var +         # interaction: race_var × mediator
       #   coefs[[paste0(m, "2")]] * (!!race0_var)^2 +               # quadratic effect of mediator
          coefs[["H1GI1Y"]]      * H1GI1Y +                         # H1GI1Y main effect
          coefs[[paste0(race_var, ":H1GI1Y")]]  * H1GI1Y +          # interaction: race_var × H1GI1Y
          rmse * UY                                                 # residual
      )
    
    cat("Created:", paste0(outcome_var, "_10_", m), "\n")
  })
  
  #----------------------------------------------------------------------  
  # Step 3: Model for the Total Causal Effect (TCE), Y ~ X
  #----------------------------------------------------------------------
  
  model <- lm(as.formula(paste0(outcome_var, " ~ ", race_var, " + H1GI1Y +", race_var, ":H1GI1Y")), 
              data = dat) # computed on original data
  
  coefs <- coef(model)
  
  #----------------------------------------------------------------------
  # Step 3A: Variables for Dep (Y), when race_var = 0 (no mediator involved) 
  #----------------------------------------------------------------------
  
  output_var_0 <- sym(paste0(outcome_var, "_0")) # when race_var = 0
  dat_expanded <- dat_expanded %>% # applied to expanded data
    mutate(
      !!output_var_0 :=
        coefs["(Intercept)"] + 
        (coefs[[race_var]] %||% 0) * 0 + 
        coefs["H1GI1Y"] * H1GI1Y +
        (coefs[[paste0(race_var, ":H1GI1Y")]])  * H1GI1Y * 0
    )
  
  #----------------------------------------------------------------------
  # Step 3B: Variables for Dep (Y), when race_var = 1 (no mediator involved) 
  #----------------------------------------------------------------------
  
  output_var_1 <- sym(paste0(outcome_var, "_1")) # when race_var = 1
  dat_expanded <- dat_expanded %>% # applied to expanded data
    mutate(
      !!output_var_1 :=
        coefs["(Intercept)"] + 
        coefs[[race_var]] + 
        coefs["H1GI1Y"] * H1GI1Y +
        (coefs[[paste0(race_var, ":H1GI1Y")]]) * H1GI1Y 
    )
  
  cat("Created:", output_var_0, "\n")
  cat("Created:", output_var_1, "\n")
  
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
  for (race_var in race_vars) {  # loop over each race variable
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
    
    saveRDS(results, file = paste0("/home/thom1336/healthdis/data/sensitivity/item/race/bootstrap_results_", wave, "_", race_var, ".rds"))
  #  saveRDS(boot_results, file = paste0("/home/thom1336/healthdis/data/sensitivity/item/race/boot_object_", wave, "_", race_var, ".rds"))
    
    cat("Saved results for wave:", wave, "and race variable:", race_var, "\n")
  }
}
