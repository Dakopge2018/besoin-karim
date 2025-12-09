# =============================================================================
# FICHIER PRINCIPAL 2: TARIFICATION (main_pricing.r)
# Estime les primes d'assurance en fonction des différentes probabilités de persistance
# =============================================================================

# === INITIALISATION ===
source("./src/const.r")
source("./src/utilities.r")

initialize_environment()
theme_set(theme_bw())

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PROGRAMME PRINCIPAL 2: TARIFICATION D'ASSURANCE VIE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# === SOUS-FONCTION: Créer le modèle de tarification ===

main_init_pricing_model <- function(hmm_model, n_states) {
  "
  Fonction auxiliaire: Initialise le moteur de tarification
  Appelée depuis le main() avec les modèles HMM ajustés
  "
  
  cat("\n🔧 Initialisation du modèle de tarification\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  # Créer l'environnement privé
  env <- new.env()
  env$hmm_model <- hmm_model
  env$n_states <- n_states
  env$risk_free_rate <- CONST$PRICING$risk_free_rate
  env$age_labels <- CONST$PRICING$age_labels
  env$age_breaks <- CONST$PRICING$age_breaks
  env$coef_matrix <- NULL
  env$xi_H <- diag(1, n_states)  # Identité par défaut
  env$homo <- FALSE
  
  cat(sprintf("✓ Modèle initialisé avec %d états\n", n_states))
  cat(sprintf("✓ Taux sans risque: %.2f%%\n", CONST$PRICING$risk_free_rate * 100))
  
  # === PRÉ-CALCUL DES COEFFICIENTS ===
  
  initialize_coefficients <- function() {
    env$coef_matrix <- precompute_mortality_coefficients(
      env$hmm_model, 
      env$n_states, 
      env$age_labels
    )
    cat("✓ Coefficients de mortalité pré-calculés\n")
  }
  
  initialize_coefficients()
  
  # === SIMULATION DES ÉTATS FUTURS ===
  
  simulate_future_states <- function(n_periods, initial_state = 1) {
    "Simule les états futurs avec transitions stochastiques"
    states <- integer(n_periods)
    states[1] <- initial_state
    
    if (env$homo) {
      # Mode homogène avec matrice xi_H
      trans_matrix <- env$xi_H / rowSums(env$xi_H)
      
      for (t in 2:n_periods) {
        trans_probs <- trans_matrix[states[t-1], ]
        trans_probs <- pmax(trans_probs, 1e-10)
        trans_probs <- trans_probs / sum(trans_probs)
        states[t] <- sample(1:env$n_states, 1, prob = trans_probs)
      }
    } else {
      # Mode stochastique: transitions basées sur HMM
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
  
  # === CALCUL DE LA SURVIE ===
  
  compute_survival_probability <- function(age, time_horizon, initial_state = 1L, n_steps = NULL) {
    "Calcule la probabilité de survie intégrée sur la période"
    
    if (is.null(n_steps)) {
      n_steps <- max(1, round(time_horizon * 52))
    }
    dt <- time_horizon / n_steps
    
    # Simuler les états futurs
    states_path <- simulate_future_states(n_steps, initial_state)
    
    # Ages et variables temporelles
    time_steps <- (0:(n_steps-1)) * dt
    ages_at_time <- age + time_steps
    
    total_weeks <- 0:(n_steps-1)
    years <- 1 + total_weeks %/% 52
    weeks <- 1 + (total_weeks %% 52)
    
    # Déterminer les groupes d'âge
    age_group_indices <- findInterval(ages_at_time, env$age_breaks, rightmost.closed = TRUE)
    age_group_indices <- pmax(1, pmin(age_group_indices, length(env$age_labels)))
    age_ranges <- env$age_labels[age_group_indices]
    
    # Calculer les intensités de mortalité
    mu_vector <- sapply(1:n_steps, function(i) {
      compute_mortality_intensity(
        age_ranges[i],
        states_path[i],
        weeks[i],
        years[i],
        env$coef_matrix
      )
    })
    
    # Hazard cumulé
    cumulative_hazard <- sum(mu_vector * dt)
    
    # Probabilité de survie
    exp(-cumulative_hazard)
  }
  
  # === CALCUL DE PRIME CLASSIQUE ===
  
  calculate_premium_classical <- function(age_start,
                                         death_benefit,
                                         contract_duration,
                                         initial_state = 1L,
                                         verbose = TRUE) {
    "Calcule la prime par méthode classique discrète annuelle"
    
    n_years <- as.integer(ceiling(contract_duration))
    
    years_vec <- seq_len(n_years)
    ages_vec <- age_start + (years_vec - 1)
    
    # Probabilités de survie année par année
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
      cat("\n📊 CALCUL ANNÉE PAR ANNÉE:\n")
      cat(sprintf("%6s %12s %12s %12s %15s\n", "Année", "P(Survie→t)", "P(Décès)", "Facteur D.", "Contribution"))
      cat(paste(rep("-", 65), collapse = ""), "\n")
      
      for (t in seq_len(min(n_years, 10))) {
        cat(sprintf("%6d %12.6f %12.6f %12.6f %15.2f\n", 
                    t, survival_probs_to_t[t], death_probs[t], discount_factors[t], contributions[t]))
      }
      
      if (n_years > 10) cat("   ... (années suivantes masquées)\n")
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
  
  # === INTERFACE PUBLIQUE ===
  
  return(list(
    compute_survival_probability = compute_survival_probability,
    calculate_premium = function(age, benefit, duration, initial_state = 1, verbose = T) {
      calculate_premium_classical(age, benefit, duration, initial_state, verbose)
    },
    
    set_xi_matrix = function(xi_H, homo = FALSE) {
      if (nrow(xi_H) != env$n_states || ncol(xi_H) != env$n_states) {
        stop(sprintf("xi_H doit être une matrice %dx%d", env$n_states, env$n_states))
      }
      env$xi_H <- xi_H
      env$homo <- homo
      cat(sprintf("✓ Matrice xi_H mise à jour\n"))
    },
    
    get_info = function() list(
      n_states = env$n_states,
      risk_free_rate = env$risk_free_rate,
      xi_H = env$xi_H,
      homo = env$homo
    )
  ))
}

# === SOUS-FONCTION: Analyse de sensibilité ===

main_sensitivity_analysis <- function(pricing_model,
                                      age_start = 65,
                                      death_benefit = 100000,
                                      contract_duration = 40,
                                      initial_state = 1,
                                      n_states = 2) {
  "
  Fonction auxiliaire: Effectue une analyse de sensibilité
  Calcule les primes pour différentes matrices de transition
  "
  
  cat("\n\n🎯 ANALYSE DE SENSIBILITÉ\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Valeurs de probabilités à tester
  p_values <- c(0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05)
  n_p <- length(p_values)
  
  premium_matrix <- matrix(NA, nrow = n_p, ncol = n_p)
  rownames(premium_matrix) <- p_values
  colnames(premium_matrix) <- p_values
  
  cat(sprintf("Calcul pour %d x %d scénarios de transition...\n", n_p, n_p))
  pb <- txtProgressBar(min = 0, max = n_p^2, style = 3)
  counter <- 0
  
  for (i in 1:n_p) {
    for (j in 1:n_p) {
      counter <- counter + 1
      p1 <- p_values[i]
      p2 <- p_values[j]
      
      if (n_states == 2) {
        # Matrice 2x2
        xi_matrix <- matrix(c(p1, 1-p1, 1-p2, p2), nrow = 2, byrow = TRUE)
      } else if (n_states == 3) {
        # Matrice 3x3 (exemple: fixer p22)
        p22 <- 0.5
        xi_matrix <- matrix(c(
          p1, (1-p1)/2, (1-p1)/2,
          (1-p22)/2, p22, (1-p22)/2,
          (1-p2)/2, (1-p2)/2, p2
        ), nrow = 3, byrow = TRUE)
      }
      
      pricing_model$set_xi_matrix(xi_matrix, homo = TRUE)
      
      result <- pricing_model$calculate_premium(
        age_start = age_start,
        death_benefit = death_benefit,
        contract_duration = contract_duration,
        initial_state = initial_state,
        verbose = FALSE
      )
      
      premium_matrix[i, j] <- result$premium_value / death_benefit * 1000  # Taux (pour 1000€)
      setTxtProgressBar(pb, counter)
    }
  }
  close(pb)
  
  cat("\n\n✓ Analyse de sensibilité terminée\n\n")
  
  return(premium_matrix)
}

# === FONCTION PRINCIPALE: Tarification complète ===

main <- function() {
  "
  Fonction principale:
  1. Charge les modèles HMM ajustés
  2. Initialise les moteurs de tarification
  3. Calcule les primes pour différents scénarios
  4. Effectue des analyses de sensibilité
  5. Génère rapports et visualisations
  "
  
  # === ÉTAPE 1: CHARGER LES MODÈLES ===
  
  cat("\n📊 ÉTAPE 1: CHARGEMENT DES MODÈLES HMM\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  # Vérifier que les modèles existent
  if (!file.exists("./models/hmm_model_2states.rds")) {
    cat("❌ Erreur: Modèle 2 états non trouvé\n")
    cat("   Exécutez d'abord main_model_comparison.r\n")
    return(invisible(NULL))
  }
  
  hmm_2states <- readRDS("./models/hmm_model_2states.rds")
  cat("✓ Modèle 2 états chargé\n")
  
  if (file.exists("./models/hmm_model_3states.rds")) {
    hmm_3states <- readRDS("./models/hmm_model_3states.rds")
    cat("✓ Modèle 3 états chargé\n")
  } else {
    hmm_3states <- NULL
  }
  
  # === ÉTAPE 2: INITIALISER LES MODÈLES DE TARIFICATION ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("📊 ÉTAPE 2: INITIALISATION DES MOTEURS DE TARIFICATION\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  pricing_2states <- main_init_pricing_model(hmm_2states, 2)
  cat("✓ Moteur de tarification 2 états initialisé\n")
  
  if (!is.null(hmm_3states)) {
    pricing_3states <- main_init_pricing_model(hmm_3states, 3)
    cat("✓ Moteur de tarification 3 états initialisé\n")
  }
  
  # === ÉTAPE 3: CALCUL DE PRIMES SIMPLES ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("📊 ÉTAPE 3: CALCUL DES PRIMES (SCÉNARIOS DE BASE)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  # Scénarios à tester
  scenarios <- data.frame(
    age = c(65, 65, 75, 75),
    duration = c(20, 40, 15, 30),
    benefit = c(100000, 100000, 50000, 50000)
  )
  
  cat("\n🔹 MODÈLE 2 ÉTATS\n")
  results_2 <- list()
  
  for (i in 1:nrow(scenarios)) {
    cat(sprintf("\n  Scénario %d: Âge %d, Durée %d ans, Capital %s€\n",
                i, scenarios$age[i], scenarios$duration[i], 
                format(scenarios$benefit[i], big.mark = ",")))
    
    result <- pricing_2states$calculate_premium(
      age_start = scenarios$age[i],
      death_benefit = scenarios$benefit[i],
      contract_duration = scenarios$duration[i],
      initial_state = 1,
      verbose = FALSE
    )
    
    cat(sprintf("  → Prime: %.2f€ (Taux: %.2f‰)\n",
                result$premium_value,
                result$premium_value / scenarios$benefit[i] * 1000))
    
    results_2[[i]] <- result
  }
  
  if (!is.null(hmm_3states)) {
    cat("\n\n🔹 MODÈLE 3 ÉTATS\n")
    results_3 <- list()
    
    for (i in 1:nrow(scenarios)) {
      cat(sprintf("\n  Scénario %d: Âge %d, Durée %d ans, Capital %s€\n",
                  i, scenarios$age[i], scenarios$duration[i], 
                  format(scenarios$benefit[i], big.mark = ",")))
      
      result <- pricing_3states$calculate_premium(
        age_start = scenarios$age[i],
        death_benefit = scenarios$benefit[i],
        contract_duration = scenarios$duration[i],
        initial_state = 2,  # État "moyen"
        verbose = FALSE
      )
      
      cat(sprintf("  → Prime: %.2f€ (Taux: %.2f‰)\n",
                  result$premium_value,
                  result$premium_value / scenarios$benefit[i] * 1000))
      
      results_3[[i]] <- result
    }
  }
  
  # === ÉTAPE 4: ANALYSE DE SENSIBILITÉ ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("📊 ÉTAPE 4: ANALYSES DE SENSIBILITÉ\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  # Pour 2 états
  sensitivity_2 <- main_sensitivity_analysis(
    pricing_2states,
    age_start = 65,
    death_benefit = 100000,
    contract_duration = 40,
    initial_state = 1,
    n_states = 2
  )
  
  cat("\n📋 Tableau de sensibilité (2 états) - Taux pour 1000€:\n")
  print(round(sensitivity_2, 4))
  
  # Exporter la matrice
  write.csv(sensitivity_2, "./results/sensitivity_2states.csv")
  cat("\n✓ Résultats exportés: ./results/sensitivity_2states.csv\n")
  
  if (!is.null(hmm_3states)) {
    sensitivity_3 <- main_sensitivity_analysis(
      pricing_3states,
      age_start = 65,
      death_benefit = 100000,
      contract_duration = 40,
      initial_state = 2,
      n_states = 3
    )
    
    cat("\n📋 Tableau de sensibilité (3 états) - Taux pour 1000€:\n")
    print(round(sensitivity_3, 4))
    
    write.csv(sensitivity_3, "./results/sensitivity_3states.csv")
    cat("\n✓ Résultats exportés: ./results/sensitivity_3states.csv\n")
  }
  
  # === ÉTAPE 5: VISUALISATIONS ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("📊 ÉTAPE 5: GÉNÉRATION DES VISUALISATIONS\n")
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
    labs(title = "Sensibilité des taux de prime (2 états)",
         subtitle = "Âge 65, Capital 100k€, Durée 40 ans",
         x = "Probabilité persistance État 1", 
         y = "Probabilité persistance État 2",
         fill = "Taux (‰)") +
    theme_minimal()
  
  print(p_heatmap2)
  ggsave("./results/heatmap_sensitivity_2states.pdf", width = 10, height = 8, dpi = 300)
  cat("✓ Graphique sauvegardé: ./results/heatmap_sensitivity_2states.pdf\n")
  
  # === RÉSUMÉ FINAL ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("✅ PROGRAMME TERMINÉ AVEC SUCCÈS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  cat("📁 Fichiers générés:\n")
  cat("  • Résultats: ./results/sensitivity_*.csv\n")
  cat("  • Graphiques: ./results/heatmap_sensitivity_*.pdf\n\n")
  
  cat("💡 UTILISATION:\n")
  cat("  • pricing_2states$calculate_premium(age, benefit, duration)\n")
  cat("  • pricing_2states$set_xi_matrix(xi_H, homo=TRUE)\n\n")
  
  return(invisible(list(
    pricing_2states = pricing_2states,
    pricing_3states = if(!is.null(hmm_3states)) pricing_3states else NULL,
    sensitivity_2 = sensitivity_2,
    sensitivity_3 = if(exists("sensitivity_3")) sensitivity_3 else NULL,
    scenarios_results = results_2
  )))
}

# === EXÉCUTION DU PROGRAMME ===

results <- main()

cat("💾 Résultats disponibles dans: results\n")
cat("   • results$pricing_2states (objet de tarification)\n")
cat("   • results$sensitivity_2 (matrice de sensibilité)\n\n")
