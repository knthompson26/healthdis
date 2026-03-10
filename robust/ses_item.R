library(purrr)
library(tidyr)
library(boot)
library(tidyverse)
library(dplyr)
library(rlang)
library(magrittr)

set.seed(123)
remove(list = ls())

# load data - bin for binary 
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

# list of ses variables
ses_vars <- c("SESlow", "SESmid") 

# function to run counterfactual calculations and Monte Carlo simulations within bootstrapping 
bootstrap_stat <- function(data, indices, wave_name, mediators, ses_var) {
  
  #--------------------------------------------
  # Bootstrap sample (original sample)
  #--------------------------------------------
  dat <- data[indices, ]
  dat$original <- 1
 # dat$UM <- rnorm(nrow(dat))  # mediator error (original) - we dont have a cont mediator anymore so we dont need this
  dat$UY <- rnorm(nrow(dat))  # outcome error (original)
  
  outcome_var <- waves[[wave_name]]$outcome
  
  #--------------------------------------------
  # Monte Carlo expansion (inner loop)
  #--------------------------------------------
  dat_expanded <- dat[rep(1:nrow(dat), each = 1000), ]  
  dat_expanded$original <- rep(c(1, rep(0, 999)), times = nrow(dat))
 # dat_expanded$UM <- rnorm(nrow(dat_expanded)) - we dont have a cont mediator anymore so we dont need this
  dat_expanded$UY <- rnorm(nrow(dat_expanded))
  
  #--------------------------------------------
  # Counterfactual loop
  #--------------------------------------------
  
  # create a list to store model results
  model_list_reg <- list()
  
  #----------------------------------------------------
  # Step 1: Variable for BINARY social (M) when ses_var = 0
  #----------------------------------------------------
  
  for (m in mediators) {
    
    # model fit - BINARY so logistic model 
    fitM <- try(                                          # tries to run without crashing whole script
      glm(as.formula(paste(m, "~", ses_var)),             # build glm formula "mediator ~ ses"
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
    new0[[ses_var]] <- 0                                 # set exposure = 0 for everyone (counterfactual condition) - make sure calling var [[]]
    p0 <- as.numeric(                                    # convert predicted probabilities to numeric vector
      predict(fitM, newdata = new0, type = "response")   # logistic inverse link -> P(M=1|sex=0)
    )
    
    # simulate mediator draws under sex=0
    new_var <- paste0(m, "_", ses_var, "0")              # name of the new simulated mediator variable
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
      outcome_var, " ~ ",                                           # y
      ses_var, " + ",                                               # x
      m, " + ",                                                     # m
      m, ":", ses_var, " + ",                                       # m*x

      # covariates (main effects)
      "H1GI1Y + ",                                                  # c1
      "black_cov + natam_cov + asian_cov + ",                       # c2 (dummies)
      "sex + ",                                                     # c3
      "happyP + drinkP_bin + health_bin + violence + ",             # new confounders
      
      # ses_var × covariates
      ses_var, ":H1GI1Y + ",
      ses_var, ":black_cov + ",
      ses_var, ":natam_cov + ",
      ses_var, ":asian_cov + ",
      ses_var, ":sex + ",
      ses_var, ":happyP + ",
      ses_var, ":drinkP_bin + ",
      ses_var, ":health_bin + ",
      ses_var, ":violence + ",
      
      # covariate–covariate interactions
      "H1GI1Y:black_cov + H1GI1Y:natam_cov + H1GI1Y:asian_cov + ",
      "H1GI1Y:sex + ",
      "H1GI1Y:happyP + ",
      "H1GI1Y:drinkP_bin + ",
      "H1GI1Y:health_bin + ",
      "H1GI1Y:violence + ",
      #   "sex:black_cov + sex:natam_cov + sex:asian_cov + ",
      "sex:happyP + ",
      "sex:drinkP_bin + ",
      "sex:health_bin + ",
      "sex:violence + ",
      #  "happyP:black_cov + happyP:natam_cov + happyP:asian_cov + ",
      "happyP:drinkP_bin + ",
      "happyP:health_bin + ",
      "happyP:violence + ",
      #    "drinkP_bin:black_cov + drinkP_bin:natam_cov + drinkP_bin:asian_cov + ",
      "drinkP_bin:health_bin + ",
      "drinkP_bin:violence + ",
      #   "health_bin:black_cov + health_bin:natam_cov + health_bin:asian_cov + ",
      "health_bin:violence"
      #    "violence:black_cov + violence:natam_cov + violence:asian_cov"
    )
    
    formula_obj <- as.formula(formula_str)
    modelY <- lm(formula_obj, data = dat)  # fit outcome model on original sample
    model_nameY <- paste0("model_", outcome_var, "_", ses_var, "_", m)
    assign(model_nameY, modelY, envir = .GlobalEnv)
    
  })
  
  #-----------------------------------------------------------------------
  # Step 2A: Variable for Y, when ses_var = 0 and M=m_0, using model from Step 2 
  #----------------------------------------------------------------------
  
  walk(mediators, function(m) {
    model_obj <- get(paste0("model_", outcome_var, "_", ses_var, "_", m))
    coefs <- coef(model_obj)
    rmse <- summary(model_obj)$sigma
    ses0_var <- sym(paste0(m, "_", ses_var, "0"))
    output_var_00 <- sym(paste0(outcome_var, "_00_", m))
    
    dat_expanded <<- dat_expanded %>%
      mutate(
        !!output_var_00 :=
          # intercept
          (coefs[["(Intercept)"]] %||% 0) +
          
          # X = ses_var = 0
          (coefs[[ses_var]] %||% 0) * 0 +
          
          # mediator terms
          (coefs[[m]] %||% 0)                          * !!ses0_var +
          (coefs[[paste0(ses_var, ":", m)]] %||% 0)    * 0 * !!ses0_var +
          
          # covariate main effects
          (coefs[["H1GI1Y"]]     %||% 0) * H1GI1Y +
          (coefs[["black_cov"]]  %||% 0) * black_cov +
          (coefs[["natam_cov"]]  %||% 0) * natam_cov +
          (coefs[["asian_cov"]]  %||% 0) * asian_cov +
          (coefs[["sex"]]        %||% 0) * sex +
          (coefs[["happyP"]]     %||% 0) * happyP +
          (coefs[["drinkP_bin"]] %||% 0) * drinkP_bin +
          (coefs[["health_bin"]] %||% 0) * health_bin +
          (coefs[["violence"]]   %||% 0) * violence +
          
          # ses_var × covariate interactions (X = 0 → vanish)
          (coefs[[paste0(ses_var, ":H1GI1Y")]]     %||% 0) * 0 * H1GI1Y +
          (coefs[[paste0(ses_var, ":black_cov")]]  %||% 0) * 0 * black_cov +
          (coefs[[paste0(ses_var, ":natam_cov")]]  %||% 0) * 0 * natam_cov +
          (coefs[[paste0(ses_var, ":asian_cov")]]  %||% 0) * 0 * asian_cov +
          (coefs[[paste0(ses_var, ":sex")]]        %||% 0) * 0 * sex +
          (coefs[[paste0(ses_var, ":happyP")]]     %||% 0) * 0 * happyP +
          (coefs[[paste0(ses_var, ":drinkP_bin")]] %||% 0) * 0 * drinkP_bin +
          (coefs[[paste0(ses_var, ":health_bin")]] %||% 0) * 0 * health_bin +
          (coefs[[paste0(ses_var, ":violence")]]   %||% 0) * 0 * violence +
          
          # covariate–covariate interactions (expanded)
          (coefs[["H1GI1Y:black_cov"]]    %||% 0) * H1GI1Y * black_cov +
          (coefs[["H1GI1Y:natam_cov"]]    %||% 0) * H1GI1Y * natam_cov +
          (coefs[["H1GI1Y:asian_cov"]]    %||% 0) * H1GI1Y * asian_cov +
          (coefs[["H1GI1Y:sex"]]          %||% 0) * H1GI1Y * sex +
          (coefs[["H1GI1Y:happyP"]]       %||% 0) * H1GI1Y * happyP +
          (coefs[["H1GI1Y:drinkP_bin"]]   %||% 0) * H1GI1Y * drinkP_bin +
          (coefs[["H1GI1Y:health_bin"]]   %||% 0) * H1GI1Y * health_bin +
          (coefs[["H1GI1Y:violence"]]     %||% 0) * H1GI1Y * violence +
          
          #    (coefs[["sex:black_cov"]]       %||% 0) * sex * black_cov +
          #    (coefs[["sex:natam_cov"]]       %||% 0) * sex * natam_cov +
          #    (coefs[["sex:asian_cov"]]       %||% 0) * sex * asian_cov +
          (coefs[["sex:happyP"]]          %||% 0) * sex * happyP +
          (coefs[["sex:drinkP_bin"]]      %||% 0) * sex * drinkP_bin +
          (coefs[["sex:health_bin"]]      %||% 0) * sex * health_bin +
          (coefs[["sex:violence"]]        %||% 0) * sex * violence +
          
          #  (coefs[["happyP:black_cov"]]    %||% 0) * happyP * black_cov +
          #  (coefs[["happyP:natam_cov"]]    %||% 0) * happyP * natam_cov +
          #  (coefs[["happyP:asian_cov"]]    %||% 0) * happyP * asian_cov +
          (coefs[["happyP:drinkP_bin"]]   %||% 0) * happyP * drinkP_bin +
          (coefs[["happyP:health_bin"]]   %||% 0) * happyP * health_bin +
          (coefs[["happyP:violence"]]     %||% 0) * happyP * violence +
          
          #   (coefs[["drinkP_bin:black_cov"]]  %||% 0) * drinkP_bin * black_cov +
          #   (coefs[["drinkP_bin:natam_cov"]]  %||% 0) * drinkP_bin * natam_cov +
          #   (coefs[["drinkP_bin:asian_cov"]]  %||% 0) * drinkP_bin * asian_cov +
          (coefs[["drinkP_bin:health_bin"]] %||% 0) * drinkP_bin * health_bin +
          (coefs[["drinkP_bin:violence"]]   %||% 0) * drinkP_bin * violence +
          
          #   (coefs[["health_bin:black_cov"]]  %||% 0) * health_bin * black_cov +
          #   (coefs[["health_bin:natam_cov"]]  %||% 0) * health_bin * natam_cov +
          #   (coefs[["health_bin:asian_cov"]]  %||% 0) * health_bin * asian_cov +
          (coefs[["health_bin:violence"]]   %||% 0) * health_bin * violence +
          
          #   (coefs[["violence:black_cov"]]    %||% 0) * violence * black_cov +
          #   (coefs[["violence:natam_cov"]]    %||% 0) * violence * natam_cov +
          #   (coefs[["violence:asian_cov"]]    %||% 0) * violence * asian_cov +
          
          # residual
          rmse * UY
      )
  })
  
  #----------------------------------------------------------------------
  # Step 2B: Variable for Y, when ses_var = 1 and M=m_0, using model from Step 2 
  #----------------------------------------------------------------------
  
  walk(mediators, function(m) {
    model_obj <- get(paste0("model_", outcome_var, "_", ses_var, "_", m))
    coefs <- coef(model_obj)
    rmse <- summary(model_obj)$sigma
    
    ses0_var   <- sym(paste0(m, "_", ses_var, "0"))
    output_var_10 <- sym(paste0(outcome_var, "_10_", m))
    
    dat_expanded <<- dat_expanded %>%
      mutate(
        !!output_var_10 :=
          # intercept
          (coefs[["(Intercept)"]] %||% 0) +
          
          # X = ses_var = 1
          (coefs[[ses_var]] %||% 0) * 1 +
          
          # mediator terms
          (coefs[[m]] %||% 0)                          * !!ses0_var +
          (coefs[[paste0(ses_var, ":", m)]] %||% 0)    * 1 * !!ses0_var +
          
          # covariate main effects
          (coefs[["H1GI1Y"]]     %||% 0) * H1GI1Y +
          (coefs[["black_cov"]]  %||% 0) * black_cov +
          (coefs[["natam_cov"]]  %||% 0) * natam_cov +
          (coefs[["asian_cov"]]  %||% 0) * asian_cov +
          (coefs[["sex"]]        %||% 0) * sex +
          (coefs[["happyP"]]     %||% 0) * happyP +
          (coefs[["drinkP_bin"]] %||% 0) * drinkP_bin +
          (coefs[["health_bin"]] %||% 0) * health_bin +
          (coefs[["violence"]]   %||% 0) * violence +
          
          # ses_var × covariate interactions (X = 1)
          (coefs[[paste0(ses_var, ":H1GI1Y")]]     %||% 0) * 1 * H1GI1Y +
          (coefs[[paste0(ses_var, ":black_cov")]]  %||% 0) * 1 * black_cov +
          (coefs[[paste0(ses_var, ":natam_cov")]]  %||% 0) * 1 * natam_cov +
          (coefs[[paste0(ses_var, ":asian_cov")]]  %||% 0) * 1 * asian_cov +
          (coefs[[paste0(ses_var, ":sex")]]        %||% 0) * 1 * sex +
          (coefs[[paste0(ses_var, ":happyP")]]     %||% 0) * 1 * happyP +
          (coefs[[paste0(ses_var, ":drinkP_bin")]] %||% 0) * 1 * drinkP_bin +
          (coefs[[paste0(ses_var, ":health_bin")]] %||% 0) * 1 * health_bin +
          (coefs[[paste0(ses_var, ":violence")]]   %||% 0) * 1 * violence +
          
          # covariate–covariate interactions (expanded)
          (coefs[["H1GI1Y:black_cov"]]    %||% 0) * H1GI1Y * black_cov +
          (coefs[["H1GI1Y:natam_cov"]]    %||% 0) * H1GI1Y * natam_cov +
          (coefs[["H1GI1Y:asian_cov"]]    %||% 0) * H1GI1Y * asian_cov +
          (coefs[["H1GI1Y:sex"]]          %||% 0) * H1GI1Y * sex +
          (coefs[["H1GI1Y:happyP"]]       %||% 0) * H1GI1Y * happyP +
          (coefs[["H1GI1Y:drinkP_bin"]]   %||% 0) * H1GI1Y * drinkP_bin +
          (coefs[["H1GI1Y:health_bin"]]   %||% 0) * H1GI1Y * health_bin +
          (coefs[["H1GI1Y:violence"]]     %||% 0) * H1GI1Y * violence +
          
          #    (coefs[["sex:black_cov"]]       %||% 0) * sex * black_cov +
          #    (coefs[["sex:natam_cov"]]       %||% 0) * sex * natam_cov +
          #    (coefs[["sex:asian_cov"]]       %||% 0) * sex * asian_cov +
          (coefs[["sex:happyP"]]          %||% 0) * sex * happyP +
          (coefs[["sex:drinkP_bin"]]      %||% 0) * sex * drinkP_bin +
          (coefs[["sex:health_bin"]]      %||% 0) * sex * health_bin +
          (coefs[["sex:violence"]]        %||% 0) * sex * violence +
          
          #    (coefs[["happyP:black_cov"]]    %||% 0) * happyP * black_cov +
          #    (coefs[["happyP:natam_cov"]]    %||% 0) * happyP * natam_cov +
          #    (coefs[["happyP:asian_cov"]]    %||% 0) * happyP * asian_cov +
          (coefs[["happyP:drinkP_bin"]]   %||% 0) * happyP * drinkP_bin +
          (coefs[["happyP:health_bin"]]   %||% 0) * happyP * health_bin +
          (coefs[["happyP:violence"]]     %||% 0) * happyP * violence +
          
          #   (coefs[["drinkP_bin:black_cov"]]  %||% 0) * drinkP_bin * black_cov +
          #   (coefs[["drinkP_bin:natam_cov"]]  %||% 0) * drinkP_bin * natam_cov +
          #   (coefs[["drinkP_bin:asian_cov"]]  %||% 0) * drinkP_bin * asian_cov +
          (coefs[["drinkP_bin:health_bin"]] %||% 0) * drinkP_bin * health_bin +
          (coefs[["drinkP_bin:violence"]]   %||% 0) * drinkP_bin * violence +
          
          #    (coefs[["health_bin:black_cov"]]  %||% 0) * health_bin * black_cov +
          #    (coefs[["health_bin:natam_cov"]]  %||% 0) * health_bin * natam_cov +
          #    (coefs[["health_bin:asian_cov"]]  %||% 0) * health_bin * asian_cov +
          (coefs[["health_bin:violence"]]   %||% 0) * health_bin * violence +
          
          #    (coefs[["violence:black_cov"]]    %||% 0) * violence * black_cov +
          #    (coefs[["violence:natam_cov"]]    %||% 0) * violence * natam_cov +
          #    (coefs[["violence:asian_cov"]]    %||% 0) * violence * asian_cov +
          
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
        ses_var, " + ",
        "H1GI1Y + black_cov + natam_cov + asian_cov + sex + ",
        ses_var, ":H1GI1Y + ",
        ses_var, ":black_cov + ",
        ses_var, ":natam_cov + ",
        ses_var, ":asian_cov + ",
        ses_var, ":sex + ",
        "H1GI1Y:black_cov + H1GI1Y:natam_cov + H1GI1Y:asian_cov + ",
        "H1GI1Y:sex + ",
        "black_cov:sex + natam_cov:sex + asian_cov:sex"
      )
    ),
    data = dat
  ) # computed on original data
  
  coefs <- coef(model)
  
  #----------------------------------------------------------------------
  # Step 3A: Variables for Dep (Y), when ses_var = 0 (no mediator involved) 
  #----------------------------------------------------------------------
  
  output_var_0 <- sym(paste0(outcome_var, "_0")) # when ses_var = 0
  dat_expanded <- dat_expanded %>%
    mutate(
      !!output_var_0 :=
        (coefs[["(Intercept)"]] %||% 0) +
        (coefs[[ses_var]] %||% 0) * 0 +
        (coefs[["H1GI1Y"]]     %||% 0) * H1GI1Y +
        (coefs[["black_cov"]]  %||% 0) * black_cov +
        (coefs[["natam_cov"]]  %||% 0) * natam_cov +
        (coefs[["asian_cov"]]  %||% 0) * asian_cov +
        (coefs[["sex"]]        %||% 0) * sex +
        (coefs[[paste0(ses_var, ":H1GI1Y")]]  %||% 0) * 0 * H1GI1Y +
        (coefs[[paste0(ses_var, ":black_cov")]] %||% 0) * 0 * black_cov +
        (coefs[[paste0(ses_var, ":natam_cov")]] %||% 0) * 0 * natam_cov +
        (coefs[[paste0(ses_var, ":asian_cov")]] %||% 0) * 0 * asian_cov +
        (coefs[[paste0(ses_var, ":sex")]]      %||% 0) * 0 * sex +
        (coefs[["H1GI1Y:black_cov"]] %||% 0) * H1GI1Y * black_cov +
        (coefs[["H1GI1Y:natam_cov"]] %||% 0) * H1GI1Y * natam_cov +
        (coefs[["H1GI1Y:asian_cov"]] %||% 0) * H1GI1Y * asian_cov +
        (coefs[["H1GI1Y:sex"]]       %||% 0) * H1GI1Y * sex +
        (coefs[["black_cov:sex"]]    %||% 0) * black_cov * sex +
        (coefs[["natam_cov:sex"]]    %||% 0) * natam_cov * sex +
        (coefs[["asian_cov:sex"]]    %||% 0) * asian_cov * sex
    )
  
  #----------------------------------------------------------------------
  # Step 3B: Variables for Dep (Y), when ses_var = 1 (no mediator involved) 
  #----------------------------------------------------------------------
  
  output_var_1 <- sym(paste0(outcome_var, "_1")) # when ses_var = 1
  dat_expanded <- dat_expanded %>%
    mutate(
      !!output_var_1 :=
        (coefs[["(Intercept)"]] %||% 0) +
        (coefs[[ses_var]] %||% 0) +
        (coefs[["H1GI1Y"]]     %||% 0) * H1GI1Y +
        (coefs[["black_cov"]]  %||% 0) * black_cov +
        (coefs[["natam_cov"]]  %||% 0) * natam_cov +
        (coefs[["asian_cov"]]  %||% 0) * asian_cov +
        (coefs[["sex"]]        %||% 0) * sex +
        (coefs[[paste0(ses_var, ":H1GI1Y")]]  %||% 0) * 1 * H1GI1Y +
        (coefs[[paste0(ses_var, ":black_cov")]] %||% 0) * 1 * black_cov +
        (coefs[[paste0(ses_var, ":natam_cov")]] %||% 0) * 1 * natam_cov +
        (coefs[[paste0(ses_var, ":asian_cov")]] %||% 0) * 1 * asian_cov +
        (coefs[[paste0(ses_var, ":sex")]]      %||% 0) * 1 * sex +
        (coefs[["H1GI1Y:black_cov"]] %||% 0) * H1GI1Y * black_cov +
        (coefs[["H1GI1Y:natam_cov"]] %||% 0) * H1GI1Y * natam_cov +
        (coefs[["H1GI1Y:asian_cov"]] %||% 0) * H1GI1Y * asian_cov +
        (coefs[["H1GI1Y:sex"]]       %||% 0) * H1GI1Y * sex +
        (coefs[["black_cov:sex"]]    %||% 0) * black_cov * sex +
        (coefs[["natam_cov:sex"]]    %||% 0) * natam_cov * sex +
        (coefs[["asian_cov:sex"]]    %||% 0) * asian_cov * sex
    )
  
  #--------------------------------------------
  # Step 4: Compute differences: me, tce, dif
  #--------------------------------------------
  
  walk(mediators, function(m) {
    me_var <- sym(paste0("me_", ses_var, "_", m))
    tce_var <- sym(paste0("tce_", ses_var, "_", m))
    dif_var <- sym(paste0("dif_", ses_var, "_", m))
    dep_00 <- sym(paste0(outcome_var, "_00_", m))
    dep_10 <- sym(paste0(outcome_var, "_10_", m))
    dep_0  <- sym(paste0(outcome_var, "_0"))
    dep_1  <- sym(paste0(outcome_var, "_1"))
    
    dat_expanded <<- dat_expanded %>% # calculated for each expanded dataset
      mutate(
        !!me_var := !!dep_10 - !!dep_00,
        !!tce_var := !!dep_1 - !!dep_0,
        !!dif_var := !!tce_var - !!me_var
      )
  })
  
  #-----------------------------------------------------------------
  # Calulate mean me, tce, and dif across all expanded data sets
  #-----------------------------------------------------------------
  
  effects <- map(mediators, function(m) {
    c(
      ME = mean(dat_expanded[[paste0("me_", ses_var, "_", m)]], na.rm = TRUE),
      TCE = mean(dat_expanded[[paste0("tce_", ses_var, "_", m)]], na.rm = TRUE),
      DIF = mean(dat_expanded[[paste0("dif_", ses_var, "_", m)]], na.rm = TRUE)
    )
  })
  
  names_list <- unlist(map2(mediators, seq_along(mediators),
                            ~ c(paste0("ME", .y), paste0("TCE", .y), paste0("DIF", .y))))
  out <- unlist(effects)
  names(out) <- names_list
  return(out)
}

#----------------------------------
# Run Bootstrapping (outer loop)
#----------------------------------

for (wave in names(waves)) {
  for (ses_var in ses_vars) {  # loop over each ses variable
    cat("\nBootstrapping wave:", wave, "for ses variable:", ses_var, "\n")
    
    # This re-samples individuals 1000 times to capture sampling variability and calculate confidence intervals 
    boot_results <- boot(
      data = waves[[wave]]$data,
      statistic = function(d, i) bootstrap_stat(d, i, wave_name = wave, mediators = mediators, ses_var = ses_var), # loops through ses_var
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
      ses    = ses_var  # add ses variable to the results
    ) %>%
      # Extract mediator number and effect type (ME/TCE/DIF)
      mutate(
        Mediator_Num = as.numeric(gsub("[^0-9]", "", Effect)),
        Effect_Type  = gsub("[0-9]", "", Effect),
        Mediator     = mediators[Mediator_Num],
        Wave         = wave
      ) %>%
      dplyr::select(Wave, ses, Mediator, Effect_Type, Mean, SE, CI_Low, CI_High) # final table
    
    saveRDS(results, file = paste0("/home/thom1336/healthdis/data/sensitivity/item/confound/ses/bootstrap_results_", wave, "_", ses_var, ".rds"))
    cat("Saved results for wave:", wave, "and ses variable:", ses_var, "\n")
  }
}
