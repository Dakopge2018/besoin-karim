# =============================================================================
# UTILITIES FILE (utilities.r)
# Reusable functions between the two main files
# =============================================================================

# === LOADING DEPENDENCIES ===
source("./src/const.r")

# === SECTION 1: DATA PREPARATION ===

load_and_prepare_data <- function() {
  "
  Loads raw data and performs initial cleaning
  "
  # Load data
  data <- readxl::read_excel(CONST$DATA$path_to_data, sheet = CONST$DATA$sheet_name)
  
  # Filter invalid data
  data <- data[!is.na(data$Death_counts) & !is.infinite(data$Death_counts) & data$Death_counts > 0, ]
  data <- na.omit(data)
  data <- data[is.finite(rowSums(data[c('year','No_year', 'Temperature', 'Death_counts', 'Weekly_exposure')])), ]
  
  return(data)
}

generate_trig_covariates <- function(week_no, year_no, period = CONST$HMM$period, degree = CONST$HMM$degree_trans_pol) {
  "
  Generates trigonometric covariates (sin, cos) for seasonality
  "
  trig_covs <- data.frame(time = week_no, trend = year_no)
  for (d in 1:degree) {
    trig_covs[[paste0("cos_", d)]] <- cos(2 * pi * d * week_no / period)
    trig_covs[[paste0("sin_", d)]] <- sin(2 * pi * d * week_no / period)
  }
  return(trig_covs)
}

prepare_hmm_dataframe <- function(data_raw) {
  "
  Prepares the dataframe for HMM fitting with all necessary variables
  "
  # Trigonometric covariates
  trig_covs <- generate_trig_covariates(
    data_raw$No_week, 
    data_raw$No_year,
    CONST$HMM$period, 
    CONST$HMM$degree_trans_pol
  )
  
  # Create final dataframe
  df <- data.frame(
    Code = data_raw$Code,
    death = data_raw$Death_counts,
    exposure = data_raw$Weekly_exposure,
    log_death_rate = log(data_raw$Death_counts / data_raw$Weekly_exposure),
    death_rate = data_raw$Death_counts / data_raw$Weekly_exposure,
    temp_extreme = data_raw$Temp_extrem,
    temp_norm = data_raw$Temperature,
    log_exposure = log(data_raw$Weekly_exposure),
    trend = data_raw$No_year,
    Age_factor = data_raw$Age_group,
    trig_covs
  )
  
  df <- na.omit(df)
  
  # Filter to get only the specified code
  df_sub <- subset(df, Code == CONST$DATA$code_filter)
  
  return(df_sub)
}

# === SECTION 2: INITIAL PARAMETERS ===

calculer_parametres_initiaux_diversifies_simple <- function(data, n_states, iteration = 1, age_ref = NULL, mort_var = "log_death_rate", temp_var = "temp_extreme", mort_dist = "norm") {
  "
  Computes diversified initial parameters for HMM fitting
  Uses 4 different strategies for robustness
  mort_dist: 'norm' or 'pois' to indicate the distribution type
  "
  # If age_ref is specified, filter the data
  if (!is.null(age_ref)) {
    data_ref <- subset(data, Age_factor == age_ref)
    mean_death <- mean(data_ref[[mort_var]], na.rm = TRUE)
    sd_death <- sd(data_ref[[mort_var]], na.rm = TRUE)
    mean_temp <- mean(data_ref[[temp_var]], na.rm = TRUE)
    sd_temp <- sd(data_ref[[temp_var]], na.rm = TRUE)
  } else {
    mean_death <- mean(data[[mort_var]], na.rm = TRUE)
    sd_death <- sd(data[[mort_var]], na.rm = TRUE)
    mean_temp <- mean(data[[temp_var]], na.rm = TRUE)
    sd_temp <- sd(data[[temp_var]], na.rm = TRUE)
  }

  # Temperature split into quantiles
  quantiles_temp <- quantile(
    data[[temp_var]],
    probs = seq(0, 1, length.out = n_states + 1),
    na.rm = TRUE
  )

  # Initialization
  death_means <- numeric(n_states)
  death_sds <- numeric(n_states)
  temp_means <- numeric(n_states)
  temp_sds <- numeric(n_states)

  # Temperature by quantile
  for (i in 1:n_states) {
    if (i == 1) {
      mask <- data[[temp_var]] <= quantiles_temp[i + 1]
    } else if (i == n_states) {
      mask <- data[[temp_var]] > quantiles_temp[i]
    } else {
      mask <- data[[temp_var]] > quantiles_temp[i] & data[[temp_var]] <= quantiles_temp[i + 1]
    }

    temp_means[i] <- mean(data[[temp_var]][mask], na.rm = TRUE)
    temp_sds[i] <- sd(data[[temp_var]][mask], na.rm = TRUE)

    if (is.na(temp_means[i])) temp_means[i] <- mean_temp
    if (is.na(temp_sds[i]) || temp_sds[i] <= 0) temp_sds[i] <- sd_temp
  }

  # Strategies for mortality (4 strategies based on iteration)
  strategy <- iteration %% 4

  if (strategy == 1) {
    # Empirical strategy
    for (i in 1:n_states) {
      if (i == 1) {
        mask <- data[[temp_var]] <= quantiles_temp[i + 1]
      } else if (i == n_states) {
        mask <- data[[temp_var]] > quantiles_temp[i]
      } else {
        mask <- data[[temp_var]] > quantiles_temp[i] & data[[temp_var]] <= quantiles_temp[i + 1]
      }

      observed_mean <- mean(data[[mort_var]][mask], na.rm = TRUE)
      observed_sd <- sd(data[[mort_var]][mask], na.rm = TRUE)
      
      if (is.na(observed_mean) || observed_mean <= 0) observed_mean <- mean_death
      if (is.na(observed_sd) || observed_sd <= 0) observed_sd <- sd_death
      
      death_means[i] <- max(observed_mean * runif(1, 0.9, 1.1), 0.0001)
      death_sds[i] <- max(observed_sd * runif(1, 0.9, 1.1), 0.0001)
    }

  } else if (strategy == 2) {
    # Random strategy
    min_mean <- max(mean_death * 0.5, 0.0001)
    max_mean <- mean_death * 2
    death_means <- runif(n_states, min_mean, max_mean)
    
    min_sd <- max(sd_death * 0.5, 0.0001)
    max_sd <- sd_death * 2
    death_sds <- runif(n_states, min_sd, max_sd)

  } else if (strategy == 3) {
    # Permutation strategy
    empirical_means <- numeric(n_states)
    empirical_sds <- numeric(n_states)
    
    for (i in 1:n_states) {
      if (i == 1) {
        mask <- data[[temp_var]] <= quantiles_temp[i + 1]
      } else if (i == n_states) {
        mask <- data[[temp_var]] > quantiles_temp[i]
      } else {
        mask <- data[[temp_var]] > quantiles_temp[i] & data[[temp_var]] <= quantiles_temp[i + 1]
      }

      empirical_means[i] <- mean(data[[mort_var]][mask], na.rm = TRUE)
      empirical_sds[i] <- sd(data[[mort_var]][mask], na.rm = TRUE)
      
      if (is.na(empirical_means[i]) || empirical_means[i] <= 0) empirical_means[i] <- mean_death
      if (is.na(empirical_sds[i]) || empirical_sds[i] <= 0) empirical_sds[i] <- sd_death
    }

    death_means <- sample(empirical_means) * runif(n_states, 0.8, 1.2)
    death_means <- pmax(death_means, 0.0001)
    
    death_sds <- sample(empirical_sds) * runif(n_states, 0.8, 1.2)
    death_sds <- pmax(death_sds, 0.0001)

  } else {
    # Uniform strategy
    death_means <- mean_death * runif(n_states, 0.7, 1.8)
    death_means <- pmax(death_means, 0.0001)
    
    death_sds <- sd_death * runif(n_states, 0.7, 1.8)
    death_sds <- pmax(death_sds, 0.0001)
  }

  # Cleanup
  death_means <- pmax(death_means, 0.0001)
  death_sds <- pmax(death_sds, 0.0001)
  temp_sds <- pmax(temp_sds, 0.01)

  death_means[!is.finite(death_means)] <- mean_death
  death_sds[!is.finite(death_sds)] <- sd_death
  temp_means[!is.finite(temp_means)] <- mean_temp
  temp_sds[!is.finite(temp_sds)] <- sd_temp
  
  # For Poisson: ensure parameters are strictly positive
  if (mort_dist == "pois") {
    death_means <- pmax(death_means, 0.001)  # Minimum for Poisson lambda
  }
  
  # Create parameter list with dynamic names
  par_list <- list()
  
  # Parameters for mortality variable (adapted to distribution)
  if (mort_dist == "pois") {
    # For Poisson: use 'rate' parameter (must match f_obs formula)
    # hmmTMB expects 'rate' for Poisson with offset
    par_list[[mort_var]] <- list(rate = pmax(death_means, 0.001))
  } else {
    # For Normal: use mean and sd
    par_list[[mort_var]] <- list(mean = as.numeric(death_means), sd = as.numeric(death_sds))
  }
  
  # Parameters for temperature (always Normal)
  par_list[[temp_var]] <- list(mean = as.numeric(temp_means), sd = as.numeric(temp_sds))
  
  return(par_list)
}

# === SECTION 3: VISUALIZATIONS ===

sauvegarder_graphiques <- function(plots, repertoire_sortie = getwd(), nb_etats, code) {
  "
  Saves generated plots to PDF and PNG
  "
  if (!dir.exists(repertoire_sortie)) {
    dir.create(repertoire_sortie, recursive = TRUE)
  }
  
  # Convergence plot
  if (!is.null(plots$convergence)) {
    ggsave(file.path(repertoire_sortie, paste0("nn_convergence_criteres_",nb_etats, code, ".pdf")), 
           plots$convergence, width = 10, height = 6)
  }
  
  # Original data
  if (!is.null(plots$donnees_originales)) {
    pdf(file.path(repertoire_sortie, paste0("nn_donnees_originales.pdf_",nb_etats, code, ".pdf")), 
        width = 8, height = 12)
    do.call(grid.arrange, c(plots$donnees_originales, list(ncol = 1)))
    dev.off()
  }
  
  # Simulated data
  if (!is.null(plots$donnees_simulees)) {
    pdf(file.path(repertoire_sortie, paste0("nn_donnees_simulees.pdf_",nb_etats, code, ".pdf")), 
        width = 8, height = 12)
    do.call(grid.arrange, c(plots$donnees_simulees, list(ncol = 1)))
    dev.off()
  }
  
  cat("✓ Plots saved in:", repertoire_sortie, "\n")
}

create_state_palette <- function(n_states) {
  "
  Creates a color palette adapted to the number of states
  "
  palette <- RColorBrewer::brewer.pal(max(3, n_states), CONST$VIZ$state_palette_name)
  if (n_states > 9) {
    palette <- colorRampPalette(palette)(n_states)
  }
  return(palette)
}

# === SECTION 4: STATISTICAL UTILITIES ===

calculate_aic <- function(loglik, k) {
  "Calculate AIC"
  2 * k - 2 * loglik
}

calculate_bic <- function(loglik, k, n) {
  "Calculate BIC"
  k * log(n) - 2 * loglik
}

extract_model_statistics <- function(hmm_model, data) {
  "
  Extracts main statistics from the HMM model
  "
  loglik <- logLik(hmm_model)[1]
  n_obs <- nrow(data)
  n_params <- length(coef(hmm_model))
  
  return(list(
    loglik = loglik,
    aic = calculate_aic(loglik, n_params),
    bic = calculate_bic(loglik, n_params, n_obs),
    n_params = n_params,
    n_obs = n_obs
  ))
}

# === SECTION 5: SIMULATION ===

simulate_from_model <- function(hmm_model, n_sim, data_template) {
  "
  Simulates data from a fitted HMM model
  "
  tryCatch({
    template_sim <- data_template[1:n_sim, ]
    
    if (nrow(data_template) < n_sim) {
      indices_repeat <- rep(seq_len(nrow(data_template)), length.out = n_sim)
      template_sim <- data_template[indices_repeat, ]
    }
  
    simulations <- hmm_model$simulate(n = n_sim, data = template_sim)
    return(simulations)
  }, error = function(e) {
    cat("⚠️  Error in simulation:", e$message, "\n")
    return(NULL)
  })
}

# === SECTION 6: MODEL COMPARISON ===

compare_models_summary <- function(results_list, n_states_list) {
  "
  Creates a comparative summary of fitted models
  "
  comparison_df <- data.frame(
    n_states = n_states_list,
    loglik = sapply(results_list, function(x) x$resume$meilleur_loglik),
    aic = sapply(results_list, function(x) x$resume$meilleur_aic),
    bic = sapply(results_list, function(x) x$resume$meilleur_bic),
    n_converged = sapply(results_list, function(x) x$resume$n_converges)
  )
  
  return(comparison_df)
}

# === SECTION 7: PRICING - COMMON FUNCTIONS ===

precompute_mortality_coefficients <- function(hmm_model, n_states, age_labels) {
  "
  Pre-computes mortality coefficients to accelerate calculations
  "
  obs_obj <- hmm_model$obs()
  coef <- obs_obj$coeff_fe()
  coef_lookup <- setNames(coef[, 1], rownames(coef))
  
  n_age_groups <- length(age_labels)
  coef_types <- c("(Intercept)", "sin_1", "cos_1", "trend")
  
  mu_matrix <- array(0, dim = c(n_states, length(coef_types), n_age_groups),
                    dimnames = list(
                      state = 1:n_states,
                      coef_type = coef_types,
                      age_group = age_labels
                    ))
  
  for (state in 1:n_states) {
    state_prefix <- paste0("log_death_rate.mean.state", state, ".")
    
    for (age_idx in 1:n_age_groups) {
      age_range <- age_labels[age_idx]
      
      for (coef_idx in 1:length(coef_types)) {
        coef_name <- paste0(state_prefix, coef_types[coef_idx])
        mu_matrix[state, coef_idx, age_idx] <- 
          ifelse(coef_name %in% names(coef_lookup), coef_lookup[coef_name], 0.0)
      }
      
      age_coef_name <- paste0(state_prefix, "Age_factor", age_range)
      if (age_coef_name %in% names(coef_lookup)) {
        mu_matrix[state, 1, age_idx] <- 
          mu_matrix[state, 1, age_idx] + coef_lookup[age_coef_name]
      }
    }
  }
  
  return(mu_matrix)
}

compute_mortality_intensity <- function(age_range, state, week, year, coef_matrix) {
  "
  Computes mortality intensity for given parameters
  "
  sin_val <- sin(2 * pi * week / 52)
  cos_val <- cos(2 * pi * week / 52)
  
  age_idx <- which(CONST$PRICING$age_labels == age_range)
  if (length(age_idx) == 0) return(1e-8)
  
  coefs <- coef_matrix[state, , age_idx[1]]
  
  log_death_rate <- coefs[1] + coefs[2] * sin_val + coefs[3] * cos_val + coefs[4] * year
  exp(log_death_rate)  # Return death_rate (intensity)
}

# === EXPORT AND REPORTING ===

export_results_to_csv <- function(results_df, filename) {
  "Exports results to CSV"
  write.csv(results_df, file = filename, row.names = FALSE)
  cat(sprintf("✓ Results exported: %s\n", filename))
}

export_model_summary <- function(model_results, filename) {
  "Exports model summary"
  summary_text <- sprintf("
MODEL SUMMARY
=============

Number of states: %d
Number of convergences: %d/%d
Best AIC: %.2f
Best BIC: %.2f
Best LogLik: %.2f

Optimal iteration: %d
",
    max(sapply(model_results$tous_resultats, function(x) if(!is.null(x$convergence) && x$convergence) x$iteration else 0)),
    model_results$resume$n_converges,
    model_results$resume$n_ajustements,
    model_results$resume$meilleur_aic,
    model_results$resume$meilleur_bic,
    model_results$resume$meilleur_loglik,
    model_results$resume$iteration_optimale
  )
  
  writeLines(summary_text, con = filename)
  cat(sprintf("✓ Summary exported: %s\n", filename))
}

# === SECTION 8: INTERACTIVE QUESTIONS ===

ask_model_type <- function() {
  "
  Asks user which model type to use
  Returns: 'poisson_normal' or 'normal_normal'
  "
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("HMM MODEL CONFIGURATION\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  cat("Please select the model type:\n")
  cat("  1 - Poisson-Normal (death in Poisson, temp in Normal)\n")
  cat("  2 - Normal-Normal (death in Normal, temp in Normal)\n\n")
  
  choice <- readline("Enter your choice (1 or 2): ")
  
  while (!choice %in% c("1", "2")) {
    cat("⚠️  Invalid choice. Please enter 1 or 2.\n")
    choice <- readline("Enter your choice (1 or 2): ")
  }
  
  if (choice == "1") {
    cat("✓ Model selected: Poisson-Normal\n\n")                  
    return("poisson_normal")
  } else {
    cat("✓ Model selected: Normal-Normal\n\n")
    return("normal_normal")
  }
}

ask_temp_variable <- function() {
  "
  Asks user which temperature variable to use
  Returns: 'temp_norm' or 'temp_extreme'
  "
  cat("Please select the temperature variable:\n")
  cat("  1 - Normal temperature (temp_norm)\n")
  cat("  2 - Extreme temperature (temp_extreme)\n\n")
  
  choice <- readline("Enter your choice (1 or 2): ")
  
  while (!choice %in% c("1", "2")) {
    cat("⚠️  Invalid choice. Please enter 1 or 2.\n")
    choice <- readline("Enter your choice (1 or 2): ")
  }
  
  if (choice == "1") {
    cat("✓ Variable selected: Normal temperature\n\n")
    return("temp_norm")
  } else {
    cat("✓ Variable selected: Extreme temperature\n\n")
    return("temp_extreme")
  }
}

ask_plot_results <- function() {
  "
  Asks user if they want to generate plots
  Returns: TRUE or FALSE
  "
  cat("Do you want to display and save plots?\n")
  cat("  1 - Yes\n")
  cat("  2 - No\n\n")
  
  choice <- readline("Enter your choice (1 or 2): ")
  
  while (!choice %in% c("1", "2")) {
    cat("⚠️  Invalid choice. Please enter 1 or 2.\n")
    choice <- readline("Enter your choice (1 or 2): ")
  }
  
  if (choice == "1") {
    cat("✓ Plots will be displayed and saved.\n\n")
    return(TRUE)
  } else {
    cat("✓ Plots will not be displayed.\n\n")
    return(FALSE)
  }
}

# === SECTION 9: ADVANCED METRICS ===

calculate_metrics_ajustement <- function(donnees_obs, donnees_pred, colonne_obs = "log_death_rate", colonne_pred = "log_death_rate") {
  "Computes model fit quality metrics"
  obs <- donnees_obs[[colonne_obs]]
  pred <- donnees_pred[[colonne_pred]]
  
  min_length <- min(length(obs), length(pred))
  obs <- obs[1:min_length]
  pred <- pred[1:min_length]
  
  mask_valid <- !is.na(obs) & !is.na(pred) & is.finite(obs) & is.finite(pred)
  obs <- obs[mask_valid]
  pred <- pred[mask_valid]
  
  erreurs <- obs - pred
  erreurs_abs <- abs(erreurs)
  
  mae <- mean(erreurs_abs, na.rm = TRUE)
  rmse <- sqrt(mean(erreurs^2, na.rm = TRUE))
  mape <- mean(abs((obs - pred) / (abs(obs) + 0.001)), na.rm = TRUE) * 100
  
  pseudo_r2 <- 1 - (sum(erreurs^2, na.rm = TRUE) / sum((obs - mean(obs))^2, na.rm = TRUE))
  
  return(list(
    mae = mae,
    rmse = rmse,
    mape = mape,
    pseudo_r2 = pseudo_r2,
    n_obs = length(obs)
  ))
}

compare_distributions <- function(donnees_reelles, donnees_simulees, colonne = "log_death_rate") {
  "Tests if real and simulated distributions are similar"
  real <- donnees_reelles[[colonne]]
  sim <- donnees_simulees[[colonne]]
  
  real <- real[!is.na(real) & is.finite(real)]
  sim <- sim[!is.na(sim) & is.finite(sim)]
  
  if (length(sim) > length(real)) {
    sim <- sim[1:length(real)]
  } else if (length(sim) < length(real)) {
    real <- real[1:length(sim)]
  }
  
  stats_reelles <- list(
    mean = mean(real, na.rm = TRUE),
    sd = sd(real, na.rm = TRUE),
    min = min(real, na.rm = TRUE),
    max = max(real, na.rm = TRUE)
  )
  
  stats_simulees <- list(
    mean = mean(sim, na.rm = TRUE),
    sd = sd(sim, na.rm = TRUE),
    min = min(sim, na.rm = TRUE),
    max = max(sim, na.rm = TRUE)
  )
  
  ks_test <- ks.test(real, sim)
  
  return(list(
    stats_reelles = stats_reelles,
    stats_simulees = stats_simulees,
    ks_pvalue = ks_test$p.value,
    ks_statistic = ks_test$statistic
  ))
}

calculate_ljung_box <- function(residus, lags = 10) {
  "Ljung-Box test for residual autocorrelation"
  residus_clean <- residus[!is.na(residus) & is.finite(residus)]
  
  if (length(residus_clean) <= lags) {
    return(list(
      test_statistic = NA,
      p_value = NA,
      autocorr_present = NA
    ))
  }
  
  lb_test <- Box.test(residus_clean, lag = lags, type = "Ljung-Box")
  
  return(list(
    test_statistic = lb_test$statistic,
    p_value = lb_test$p.value,
    autocorr_present = lb_test$p.value < 0.05
  ))
}

cat("✅ Utilities loaded successfully!\n")
