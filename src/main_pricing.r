# =============================================================================
# MAIN FILE 2: PRICING (main_pricing.r)
# Estimates insurance premiums based on different persistence probabilities
# =============================================================================

# === INITIALIZATION ===
source("./src/const.r")
source("./src/utilities.r")

initialize_environment()
theme_set(theme_bw())

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("MAIN PROGRAM 2: LIFE INSURANCE PRICING\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# === INTERACTIVE QUESTIONS ===

ask_model_selection <- function() {
  "
  Asks user which HMM model to use
  Returns: list with model_type and temp_variable
  "
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("HMM MODEL SELECTION\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  cat("Please select the model to use:\n\n")
  cat("MODEL TYPE:\n")
  cat("  1 - Poisson-Normal (death in Poisson)\n")
  cat("  2 - Normal-Normal (death in Normal log)\n\n")
  
  choice_model <- readline("Enter your choice (1 or 2): ")
  
  while (!choice_model %in% c("1", "2")) {
    cat("❌ Invalid choice. Please enter 1 or 2.\n")
    choice_model <- readline("Enter your choice (1 or 2): ")
  }
  
  model_type <- ifelse(choice_model == "1", "poisson_normal", "normal_normal")
  
  cat("\nTEMPERATURE VARIABLE:\n")
  cat("  1 - Normal temperature (temp_norm)\n")
  cat("  2 - Extreme temperature (temp_extreme)\n\n")
  
  choice_temp <- readline("Enter your choice (1 or 2): ")
  
  while (!choice_temp %in% c("1", "2")) {
    cat("❌ Invalid choice. Please enter 1 or 2.\n")
    choice_temp <- readline("Enter your choice (1 or 2): ")
  }
  
  temp_variable <- ifelse(choice_temp == "1", "temp_norm", "temp_extreme")
  
  cat(sprintf("\n✓ Model selected: %s + %s\n\n", model_type, temp_variable))
  
  return(list(
    model_type = model_type,
    temp_variable = temp_variable
  ))
}

# Get model selection
model_selection <- ask_model_selection()

# === HELPER FUNCTION: Create pricing model ===

main_init_pricing_model <- function(hmm_model, n_states) {
  "
  Helper function: Initializes the pricing engine
  Called from main() with fitted HMM models
  "
  
  cat("\n🔧 Initializing pricing model\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  # Create private environment
  env <- new.env()
  env$hmm_model <- hmm_model
  env$n_states <- n_states
  env$risk_free_rate <- CONST$PRICING$risk_free_rate
  env$age_labels <- CONST$PRICING$age_labels
  env$age_breaks <- CONST$PRICING$age_breaks
  env$coef_matrix <- NULL
  env$xi_H <- diag(1, n_states)  # Identity by default
  env$homo <- FALSE
  
  cat(sprintf("✓ Model initialized with %d states\n", n_states))
  cat(sprintf("✓ Risk-free rate: %.2f%%\n", CONST$PRICING$risk_free_rate * 100))
  
  # === PRE-COMPUTING COEFFICIENTS ===
  
  initialize_coefficients <- function() {
    env$coef_matrix <- precompute_mortality_coefficients(
      env$hmm_model, 
      env$n_states, 
      env$age_labels
    )
    cat("✓ Mortality coefficients pre-computed\n")
  }
  
  initialize_coefficients()
  
  # === FUTURE STATES SIMULATION ===
  
  simulate_future_states <- function(n_periods, initial_state = 1) {
    "Simulates future states with stochastic transitions"
    states <- integer(n_periods)
    states[1] <- initial_state
    
    if (env$homo) {
      # Homogeneous mode with xi_H matrix
      trans_matrix <- env$xi_H / rowSums(env$xi_H)
      
      for (t in 2:n_periods) {
        trans_probs <- trans_matrix[states[t-1], ]
        trans_probs <- pmax(trans_probs, 1e-10)
        trans_probs <- trans_probs / sum(trans_probs)
        states[t] <- sample(1:env$n_states, 1, prob = trans_probs)
      }
    } else {
      # Stochastic mode: HMM-based transitions
      for (t in 2:n_periods) {
        tryCatch({
          Q <- as.matrix(env$hmm_model$predict(what = "tpm", t = t)[,,1])
          Q_tilde <- (Q * env$xi_H) / rowSums(Q * env$xi_H)
          trans_probs <- Q_tilde[states[t-1], ]
          trans_probs <- pmax(trans_probs, 1e-10)
          trans_probs <- trans_probs / sum(trans_probs)
          states[t] <- sample(1:env$n_states, 1, prob = trans_probs)
        }, error = function(e) {
          states[t] <<- sample(1:env$n_states, 1, prob = rep(1/env$n_states, env$n_states))
        })
      }
    }
    
    return(states)
  }
  
  # === SURVIVAL CALCULATION ===
  
  compute_survival_probability <- function(age, time_horizon, initial_state = 1L, n_steps = NULL) {
    "Computes integrated survival probability over the period"
    
    if (is.null(n_steps)) {
      n_steps <- max(1, round(time_horizon * 52))
    }
    dt <- time_horizon / n_steps
    
    # Simulate future states
    states_path <- simulate_future_states(n_steps, initial_state)
    
    # Ages and temporal variables
    time_steps <- (0:(n_steps-1)) * dt
    ages_at_time <- age + time_steps
    
    total_weeks <- 0:(n_steps-1)
    years <- 1 + total_weeks %/% 52
    weeks <- 1 + (total_weeks %% 52)
    
    # Determine age groups
    age_group_indices <- findInterval(ages_at_time, env$age_breaks, rightmost.closed = TRUE)
    age_group_indices <- pmax(1, pmin(age_group_indices, length(env$age_labels)))
    age_ranges <- env$age_labels[age_group_indices]
    
    # Compute mortality intensities
    mu_vector <- sapply(1:n_steps, function(i) {
      compute_mortality_intensity(
        age_ranges[i],
        states_path[i],
        weeks[i],
        years[i],
        env$coef_matrix
      )
    })
    
    # Cumulative hazard
    cumulative_hazard <- sum(mu_vector * dt)
    
    # Survival probability
    exp(-cumulative_hazard)
  }
  
  # === CLASSICAL PREMIUM CALCULATION ===
  
  calculate_premium_classical <- function(age_start,
                                         death_benefit,
                                         contract_duration,
                                         initial_state = 1L,
                                         verbose = TRUE) {
    "Computes premium using classical discrete annual method"
    
    n_years <- as.integer(ceiling(contract_duration))
    
    years_vec <- seq_len(n_years)
    ages_vec <- age_start + (years_vec - 1)
    
    # Survival probabilities year by year
    survival_1_year_vec <- sapply(ages_vec, function(age) {
      compute_survival_probability(
        age = age,
        time_horizon = 1.0,
        initial_state = initial_state,
        n_steps = 52
      )
    })
    
    death_probs <- 1.0 - survival_1_year_vec
    survival_probs_to_t <- c(1.0, cumprod(1.0 - death_probs))[1:n_years]
    discount_factors <- (1 + env$risk_free_rate)^(-years_vec)
    contributions <- death_benefit * discount_factors * survival_probs_to_t * death_probs
    premium_total <- sum(contributions)
    
    if (verbose) {
      cat("\n📊 YEAR-BY-YEAR CALCULATION:\n")
      cat(sprintf("%6s %12s %12s %12s %15s\n", "Year", "P(Survival→t)", "P(Death)", "Disc. Factor", "Contribution"))
      cat(paste(rep("-", 65), collapse = ""), "\n")
      
      for (t in seq_len(min(n_years, 10))) {
        cat(sprintf("%6d %12.6f %12.6f %12.6f %15.2f\n", 
                    t, survival_probs_to_t[t], death_probs[t], discount_factors[t], contributions[t]))
      }
      
      if (n_years > 10) cat("   ... (following years hidden)\n")
      cat(paste(rep("-", 65), collapse = ""), "\n")
      cat(sprintf("%6s %12s %12s %12s %15.2f\n", "TOTAL", "", "", "", premium_total))
    }
    
    return(list(
      premium_value = premium_total,
      method = "classical",
      breakdown = data.frame(
        year = years_vec,
        survival_prob = survival_probs_to_t,
        death_prob = death_probs,
        discount_factor = discount_factors,
        contribution = contributions
      ),
      total_death_probability = sum(survival_probs_to_t * death_probs)
    ))
  }
  
  # === PUBLIC INTERFACE ===
  
  return(list(
    compute_survival_probability = compute_survival_probability,
    calculate_premium = function(age, benefit, duration, initial_state = 1, verbose = T) {
      calculate_premium_classical(age, benefit, duration, initial_state, verbose)
    },
    
    set_xi_matrix = function(xi_H, homo = FALSE) {
      if (nrow(xi_H) != env$n_states || ncol(xi_H) != env$n_states) {
        stop(sprintf("xi_H must be a %dx%d matrix", env$n_states, env$n_states))
      }
      env$xi_H <- xi_H
      env$homo <- homo
      cat(sprintf("✓ xi_H matrix updated\n"))
    },
    
    get_info = function() list(
      n_states = env$n_states,
      risk_free_rate = env$risk_free_rate,
      xi_H = env$xi_H,
      homo = env$homo
    )
  ))
}

# === HELPER FUNCTION: Sensitivity analysis ===

main_sensitivity_analysis <- function(pricing_model,
                                      age_start = 65,
                                      death_benefit = 100000,
                                      contract_duration = 40,
                                      initial_state = 1,
                                      n_states = 2,
                                      p22_values = NULL) {
  "
  Helper function: Performs sensitivity analysis
  Computes premiums for different transition matrices
  For 3 states: generates 5 matrices for different p22 values
  "
  
  cat("\n\n🎯 SENSITIVITY ANALYSIS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Probability values to test
  p_values <- c(0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05)
  n_p <- length(p_values)
  
  # For 3 states: test with 5 p22 values
  if (n_states == 3 && is.null(p22_values)) {
    p22_values <- c(0.95, 0.75, 0.5, 0.25, 0.05)
  }
  
  # === 2 STATE CASE: Single matrix ===
  if (n_states == 2) {
    premium_matrix <- matrix(NA, nrow = n_p, ncol = n_p)
    rownames(premium_matrix) <- p_values
    colnames(premium_matrix) <- p_values
    
    cat(sprintf("Computing for %d x %d transition scenarios...\n", n_p, n_p))
    pb <- txtProgressBar(min = 0, max = n_p^2, style = 3)
    counter <- 0
    
    for (i in 1:n_p) {
      for (j in 1:n_p) {
        counter <- counter + 1
        p1 <- p_values[i]
        p2 <- p_values[j]
        
        # 2x2 matrix
        xi_matrix <- matrix(c(p1, 1-p1, 1-p2, p2), nrow = 2, byrow = TRUE)
        
        pricing_model$set_xi_matrix(xi_matrix, homo = TRUE)
        
        result <- pricing_model$calculate_premium(
          age = age_start,
          benefit = death_benefit,
          duration = contract_duration,
          initial_state = initial_state,
          verbose = FALSE
        )
        
        premium_matrix[i, j] <- result$premium_value / death_benefit * 1000  # Taux (pour 1000€)
        setTxtProgressBar(pb, counter)
      }
    }
    close(pb)
    
    cat("\n\n✓ Sensitivity analysis completed\n\n")
    return(premium_matrix)
  }
  
  # === 3 STATE CASE: 5 matrices for 5 p22 values ===
  if (n_states == 3) {
    matrices_list <- list()
    
    total_calcs <- length(p22_values) * n_p^2
    cat(sprintf("Computing for %d p22 values x %d x %d scenarios (Total: %d calculations)...\n", 
                length(p22_values), n_p, n_p, total_calcs))
    pb <- txtProgressBar(min = 0, max = total_calcs, style = 3)
    counter <- 0
    
    for (k in 1:length(p22_values)) {
      p22 <- p22_values[k]
      
      premium_matrix <- matrix(NA, nrow = n_p, ncol = n_p)
      rownames(premium_matrix) <- p_values
      colnames(premium_matrix) <- p_values
      
      for (i in 1:n_p) {
        for (j in 1:n_p) {
          counter <- counter + 1
          p1 <- p_values[i]
          p2 <- p_values[j]
          
          # 3x3 matrix with variable p22
          xi_matrix <- matrix(c(
            p1, (1-p1)/2, (1-p1)/2,
            (1-p22)/2, p22, (1-p22)/2,
            (1-p2)/2, (1-p2)/2, p2
          ), nrow = 3, byrow = TRUE)
          
          pricing_model$set_xi_matrix(xi_matrix, homo = TRUE)
          
          result <- pricing_model$calculate_premium(
            age = age_start,
            benefit = death_benefit,
            duration = contract_duration,
            initial_state = initial_state,
            verbose = FALSE
          )
          
          premium_matrix[i, j] <- result$premium_value / death_benefit * 1000  # Rate (per 1000€)
          setTxtProgressBar(pb, counter)
        }
      }
      
      # Store matrix with descriptive name
      matrix_name <- sprintf("p22_%.2f", p22)
      matrices_list[[matrix_name]] <- premium_matrix
    }
    close(pb)
    
    cat("\n\n✓ Sensitivity analysis completed\n")
    cat(sprintf("✓ 5 matrices generated for p22 = c(%.2f, %.2f, %.2f, %.2f, %.2f)\n\n", 
                p22_values[1], p22_values[2], p22_values[3], p22_values[4], p22_values[5]))
    
    return(matrices_list)
  }
}

# === MAIN FUNCTION: Complete pricing ===

main <- function() {
  "
  Main function:
  1. Loads fitted HMM models
  2. Initializes pricing engines
  3. Computes premiums for different scenarios
  4. Performs sensitivity analyses
  5. Generates reports and visualizations
  "
  
  # === STEP 1: LOAD MODELS ===
  
  cat("\n📊 STEP 1: LOADING HMM MODELS\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  # Create model suffix
  model_suffix <- paste0(
    model_selection$model_type, "_",
    model_selection$temp_variable
  )
  
  file_model_2states <- sprintf("./models/hmm_model_2states_%s.rds", model_suffix)
  file_model_3states <- sprintf("./models/hmm_model_3states_%s.rds", model_suffix)
  
  # Check that models exist
  if (!file.exists(file_model_2states)) {
    cat(sprintf("❌ Error: 2-state model not found\n", file_model_2states))
    cat(sprintf("   Expected file: %s\n", file_model_2states))
    cat("   Please run main_model_comparison.r first with the same configuration\n")
    return(invisible(NULL))
  }
  
  hmm_2states <- readRDS(file_model_2states)
  cat(sprintf("✓ 2-state model loaded: %s\n", file_model_2states))
  
  if (file.exists(file_model_3states)) {
    hmm_3states <- readRDS(file_model_3states)
    cat(sprintf("✓ 3-state model loaded: %s\n", file_model_3states))
  } else {
    cat(sprintf("⚠️  3-state model not found: %s\n", file_model_3states))
    hmm_3states <- NULL
  }
  
  # === STEP 2: INITIALIZE PRICING ENGINES ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("📊 STEP 2: INITIALIZING PRICING ENGINES\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  pricing_2states <- main_init_pricing_model(hmm_2states, 2)
  cat("✓ 2-state pricing engine initialized\n")
  
  if (!is.null(hmm_3states)) {
    pricing_3states <- main_init_pricing_model(hmm_3states, 3)
    cat("✓ 3-state pricing engine initialized\n")
  }
  
  # === STEP 3: SIMPLE PREMIUM CALCULATION ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("📊 STEP 3: PREMIUM CALCULATION (BASE SCENARIOS)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  # Scenarios to test
  scenarios <- data.frame(
    age = c(65, 65, 75, 75),
    duration = c(20, 40, 15, 30),
    benefit = c(100000, 100000, 50000, 50000)
  )
  
  cat("\n🔹 2-STATE MODEL\n")
  results_2 <- list()
  
  for (i in 1:nrow(scenarios)) {
    cat(sprintf("\n  Scenario %d: Age %d, Duration %d years, Benefit %s€\n",
                i, scenarios$age[i], scenarios$duration[i], 
                format(scenarios$benefit[i], big.mark = ",")))
    
    result <- pricing_2states$calculate_premium(
      age = scenarios$age[i],
      benefit = scenarios$benefit[i],
      duration = scenarios$duration[i],
      initial_state = 1,
      verbose = FALSE
    )
    
    cat(sprintf("  → Premium: %.2f€ (Rate: %.2f‰)\n",
                result$premium_value,
                result$premium_value / scenarios$benefit[i] * 1000))
    
    results_2[[i]] <- result
  }
  
  if (!is.null(hmm_3states)) {
    cat("\n\n🔹 3-STATE MODEL\n")
    results_3 <- list()
    
    for (i in 1:nrow(scenarios)) {
      cat(sprintf("\n  Scenario %d: Age %d, Duration %d years, Benefit %s€\n",
                  i, scenarios$age[i], scenarios$duration[i], 
                  format(scenarios$benefit[i], big.mark = ",")))
      
      result <- pricing_3states$calculate_premium(
        age = scenarios$age[i],
        benefit = scenarios$benefit[i],
        duration = scenarios$duration[i],
        initial_state = 2,  # État "moyen"
        verbose = FALSE
      )
      
      cat(sprintf("  → Premium: %.2f€ (Rate: %.2f‰)\n",
                  result$premium_value,
                  result$premium_value / scenarios$benefit[i] * 1000))
      
      results_3[[i]] <- result
    }
  }
  
  # === STEP 4: SENSITIVITY ANALYSES ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("📊 STEP 4: SENSITIVITY ANALYSES\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  # For 2 states
  sensitivity_2 <- main_sensitivity_analysis(
    pricing_2states,
    age_start = 65,
    death_benefit = 100000,
    contract_duration = 40,
    initial_state = 1,
    n_states = 2
  )
  
  cat("\n📋 Sensitivity table (2 states) - Rate per 1000€:\n")
  print(round(sensitivity_2, 4))
  
  # Export matrix with model suffix
  file_sens_2 <- sprintf("./results/sensitivity_2states_%s.csv", model_suffix)
  write.csv(sensitivity_2, file_sens_2)
  cat(sprintf("\n✓ Results exported: %s\n", file_sens_2))
  
  if (!is.null(hmm_3states)) {
    sensitivity_3 <- main_sensitivity_analysis(
      pricing_3states,
      age_start = 65,
      death_benefit = 100000,
      contract_duration = 40,
      initial_state = 2,
      n_states = 3
    )
    
    cat("\n📋 SENSITIVITY TABLES (3 STATES) - Rate per 1000€:\n")
    cat(paste(rep("=", 70), collapse = ""), "\n\n")
    
    for (k in 1:length(sensitivity_3)) {
      matrix_name <- names(sensitivity_3)[k]
      cat(sprintf("\n🔹 %s:\n", matrix_name))
      cat(paste(rep("-", 70), collapse = ""), "\n")
      print(round(sensitivity_3[[k]], 4))
      
      # Export each matrix with model suffix
      filename <- sprintf("./results/sensitivity_3states_%s_%s.csv", matrix_name, model_suffix)
      write.csv(sensitivity_3[[k]], filename)
      cat(sprintf("✓ Exported: %s\n", filename))
    }
  }
  
  # === STEP 5: VISUALIZATIONS ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("📊 STEP 5: GENERATING VISUALIZATIONS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Heatmap sensibilité 2 états
  df_sens2 <- as.data.frame(sensitivity_2)
  df_sens2$p1 <- rownames(sensitivity_2)
  df_sens2_long <- tidyr::pivot_longer(df_sens2, -p1, names_to = "p2", values_to = "taux")
  df_sens2_long$p1 <- as.numeric(df_sens2_long$p1)
  df_sens2_long$p2 <- as.numeric(df_sens2_long$p2)
  
  p_heatmap2 <- ggplot(df_sens2_long, aes(x = p1, y = p2, fill = taux)) +
    geom_tile() +
    scale_fill_gradient(low = "#fee5d9", high = "#a1330d") +
    labs(title = "Premium rate sensitivity (2 states)",
         subtitle = "Age 65, Benefit 100k€, Duration 40 years",
         x = "State 1 persistence probability", 
         y = "State 2 persistence probability",
         fill = "Rate (‰)") +
    theme_minimal()
  
  print(p_heatmap2)
  pdf_name_2 <- sprintf("./results/heatmap_sensitivity_2states_%s.pdf", model_suffix)
  ggsave(pdf_name_2, width = 10, height = 8, dpi = 300)
  cat(sprintf("✓ Plot saved: %s\n", pdf_name_2))
  
  # 3-state sensitivity heatmaps (5 matrices)
  if (!is.null(hmm_3states) && exists("sensitivity_3")) {
    cat("\n\n🔹 3-state heatmaps (5 p22 values):\n")
    
    for (k in 1:length(sensitivity_3)) {
      matrix_name <- names(sensitivity_3)[k]
      sens_matrix <- sensitivity_3[[k]]
      
      df_sens3 <- as.data.frame(sens_matrix)
      df_sens3$p1 <- rownames(sens_matrix)
      df_sens3_long <- tidyr::pivot_longer(df_sens3, -p1, names_to = "p2", values_to = "taux")
      df_sens3_long$p1 <- as.numeric(df_sens3_long$p1)
      df_sens3_long$p2 <- as.numeric(df_sens3_long$p2)
      
      p_heatmap3 <- ggplot(df_sens3_long, aes(x = p1, y = p2, fill = taux)) +
        geom_tile() +
        scale_fill_gradient(low = "#e8f4f8", high = "#08519c") +
        labs(title = sprintf("Premium rate sensitivity (3 states - %s)", matrix_name),
             subtitle = "Age 65, Benefit 100k€, Duration 40 years",
             x = "State 1 persistence probability", 
             y = "State 3 persistence probability",
             fill = "Rate (‰)") +
        theme_minimal()
      
      print(p_heatmap3)
      
      pdf_name <- sprintf("./results/heatmap_sensitivity_3states_%s_%s.pdf", matrix_name, model_suffix)
      ggsave(pdf_name, width = 10, height = 8, dpi = 300)
      cat(sprintf("  ✓ %s\n", pdf_name))
    }
  }
  
  # === FINAL SUMMARY ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("✅ PROGRAM COMPLETED SUCCESSFULLY\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  cat("📁 Generated files:\n")
  cat("  • 2 states: ./results/sensitivity_2states.csv + heatmap PDF\n")
  if (!is.null(hmm_3states)) {
    cat("  • 3 states: 5 matrices (p22=0.95, 0.75, 0.50, 0.25, 0.05)\n")
    cat("             sensitivity_3states_p22_*.csv + heatmap PDF\n")
  }
  cat("\n")
  
  cat("💡 USAGE:\n")
  cat("  • pricing_2states$calculate_premium(age, benefit, duration)\n")
  cat("  • pricing_2states$set_xi_matrix(xi_H, homo=TRUE)\n\n")
  
  return(invisible(list(
    pricing_2states = pricing_2states,
    pricing_3states = if(!is.null(hmm_3states)) pricing_3states else NULL,
    sensitivity_2 = sensitivity_2,
    sensitivity_3 = if(exists("sensitivity_3")) sensitivity_3 else NULL,
    scenarios_results = results_2,
    scenarios_results_3states = if(exists("results_3")) results_3 else NULL
  )))
}

# === PROGRAM EXECUTION ===

results <- main()

cat("💾 Results available in: results\n")
cat("   • results$pricing_2states (pricing object)\n")
cat("   • results$pricing_3states (3-state pricing object)\n")
cat("   • results$scenarios_results (premiums 4 scenarios - 2 states)\n")
cat("   • results$scenarios_results_3states (premiums 4 scenarios - 3 states)\n")
