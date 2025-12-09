# =============================================================================
# FICHIER D'UTILITAIRES (utilities.r)
# Fonctions réutilisables entre les deux fichiers principaux
# =============================================================================

# === CHARGEMENT DES DÉPENDANCES ===
source("./src/const.r")

# === SECTION 1: PRÉPARATION DES DONNÉES ===

load_and_prepare_data <- function() {
  "
  Charge les données brutes et effectue le nettoyage initial
  "
  # Charger les données
  data <- readxl::read_excel(CONST$DATA$path_to_data, sheet = CONST$DATA$sheet_name)
  
  # Filtrer les données invalides
  data <- data[!is.na(data$Death_counts) & !is.infinite(data$Death_counts) & data$Death_counts > 0, ]
  data <- na.omit(data)
  data <- data[is.finite(rowSums(data[c('year','No_year', 'Temperature', 'Death_counts', 'Weekly_exposure')])), ]
  
  return(data)
}

generate_trig_covariates <- function(week_no, year_no, period = CONST$HMM$period, degree = CONST$HMM$degree_trans_pol) {
  "
  Génère les covariables trigonométriques (sin, cos) pour la saisonnalité
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
  Prépare le dataframe pour l'ajustement HMM avec toutes les variables nécessaires
  "
  # Covariables trigonométriques
  trig_covs <- generate_trig_covariates(
    data_raw$No_week, 
    data_raw$No_year,
    CONST$HMM$period, 
    CONST$HMM$degree_trans_pol
  )
  
  # Création du dataframe final
  df <- data.frame(
    Code = 'ITI43',
    death = data_raw$Death_counts,
    Age_factor = as.factor(data_raw$Age_group),
    exposure = data_raw$Weekly_exposure,
    log_death_rate = log(data_raw$Death_counts / data_raw$Weekly_exposure),
    death_rate = data_raw$Death_counts / data_raw$Weekly_exposure,
    temp_extreme = data_raw$Temp_extrem,
    temp_norm = data_raw$Temperature,
    log_exposure = log(data_raw$Weekly_exposure),
    trend = data_raw$No_year,
    trig_covs
  )
  
  df <- na.omit(df)
  return(df)
}

# === SECTION 2: PARAMÈTRES INITIAUX ===

calculer_parametres_initiaux_diversifies_simple <- function(data, n_states, iteration = 1, age_ref = NULL) {
  "
  Calcule des paramètres initiaux diversifiés pour l'ajustement HMM
  Utilise 4 stratégies différentes pour robustesse
  "
  # Si age_ref est spécifié, filtrer les données
  if (!is.null(age_ref)) {
    data_ref <- subset(data, Age_factor == age_ref)
    mean_death <- mean(data_ref$log_death_rate, na.rm = TRUE)
    sd_death <- sd(data_ref$log_death_rate, na.rm = TRUE)
    mean_temp <- mean(data_ref$temp_extreme, na.rm = TRUE)
    sd_temp <- sd(data_ref$temp_extreme, na.rm = TRUE)
  } else {
    mean_death <- mean(data$log_death_rate, na.rm = TRUE)
    sd_death <- sd(data$log_death_rate, na.rm = TRUE)
    mean_temp <- mean(data$temp_extreme, na.rm = TRUE)
    sd_temp <- sd(data$temp_extreme, na.rm = TRUE)
  }

  # Découpage des températures en quantiles
  quantiles_temp <- quantile(
    data$temp_extreme,
    probs = seq(0, 1, length.out = n_states + 1),
    na.rm = TRUE
  )

  # Initialisation
  death_means <- numeric(n_states)
  death_sds <- numeric(n_states)
  temp_means <- numeric(n_states)
  temp_sds <- numeric(n_states)

  # Température par quantile
  for (i in 1:n_states) {
    if (i == 1) {
      mask <- data$temp_extreme <= quantiles_temp[i + 1]
    } else if (i == n_states) {
      mask <- data$temp_extreme > quantiles_temp[i]
    } else {
      mask <- data$temp_extreme > quantiles_temp[i] & data$temp_extreme <= quantiles_temp[i + 1]
    }

    temp_means[i] <- mean(data$temp_extreme[mask], na.rm = TRUE)
    temp_sds[i] <- sd(data$temp_extreme[mask], na.rm = TRUE)

    if (is.na(temp_means[i])) temp_means[i] <- mean_temp
    if (is.na(temp_sds[i]) || temp_sds[i] <= 0) temp_sds[i] <- sd_temp
  }

  # Stratégies pour la mortalité (4 stratégies selon iteration)
  strategy <- iteration %% 4

  if (strategy == 1) {
    # Stratégie empirique
    for (i in 1:n_states) {
      if (i == 1) {
        mask <- data$temp_extreme <= quantiles_temp[i + 1]
      } else if (i == n_states) {
        mask <- data$temp_extreme > quantiles_temp[i]
      } else {
        mask <- data$temp_extreme > quantiles_temp[i] & data$temp_extreme <= quantiles_temp[i + 1]
      }

      observed_mean <- mean(data$log_death_rate[mask], na.rm = TRUE)
      observed_sd <- sd(data$log_death_rate[mask], na.rm = TRUE)
      
      if (is.na(observed_mean) || observed_mean <= 0) observed_mean <- mean_death
      if (is.na(observed_sd) || observed_sd <= 0) observed_sd <- sd_death
      
      death_means[i] <- max(observed_mean * runif(1, 0.9, 1.1), 0.0001)
      death_sds[i] <- max(observed_sd * runif(1, 0.9, 1.1), 0.0001)
    }

  } else if (strategy == 2) {
    # Stratégie aléatoire
    min_mean <- max(mean_death * 0.5, 0.0001)
    max_mean <- mean_death * 2
    death_means <- runif(n_states, min_mean, max_mean)
    
    min_sd <- max(sd_death * 0.5, 0.0001)
    max_sd <- sd_death * 2
    death_sds <- runif(n_states, min_sd, max_sd)

  } else if (strategy == 3) {
    # Stratégie permutation
    empirical_means <- numeric(n_states)
    empirical_sds <- numeric(n_states)
    
    for (i in 1:n_states) {
      if (i == 1) {
        mask <- data$temp_extreme <= quantiles_temp[i + 1]
      } else if (i == n_states) {
        mask <- data$temp_extreme > quantiles_temp[i]
      } else {
        mask <- data$temp_extreme > quantiles_temp[i] & data$temp_extreme <= quantiles_temp[i + 1]
      }

      empirical_means[i] <- mean(data$log_death_rate[mask], na.rm = TRUE)
      empirical_sds[i] <- sd(data$log_death_rate[mask], na.rm = TRUE)
      
      if (is.na(empirical_means[i]) || empirical_means[i] <= 0) empirical_means[i] <- mean_death
      if (is.na(empirical_sds[i]) || empirical_sds[i] <= 0) empirical_sds[i] <- sd_death
    }

    death_means <- sample(empirical_means) * runif(n_states, 0.8, 1.2)
    death_means <- pmax(death_means, 0.0001)
    
    death_sds <- sample(empirical_sds) * runif(n_states, 0.8, 1.2)
    death_sds <- pmax(death_sds, 0.0001)

  } else {
    # Stratégie uniforme
    death_means <- mean_death * runif(n_states, 0.7, 1.8)
    death_means <- pmax(death_means, 0.0001)
    
    death_sds <- sd_death * runif(n_states, 0.7, 1.8)
    death_sds <- pmax(death_sds, 0.0001)
  }

  # Nettoyage
  death_means <- pmax(death_means, 0.0001)
  death_sds <- pmax(death_sds, 0.0001)
  temp_sds <- pmax(temp_sds, 0.01)

  death_means[!is.finite(death_means)] <- mean_death
  death_sds[!is.finite(death_sds)] <- sd_death
  temp_means[!is.finite(temp_means)] <- mean_temp
  temp_sds[!is.finite(temp_sds)] <- sd_temp
  
  return(list(
    log_death_rate = list(mean = death_means, sd = death_sds),
    temp_extreme = list(mean = temp_means, sd = temp_sds)
  ))
}

# === SECTION 3: VISUALISATIONS ===

sauvegarder_graphiques <- function(plots, repertoire_sortie = getwd(), nb_etats, code) {
  "
  Sauvegarde les graphiques générés en PDF et PNG
  "
  if (!dir.exists(repertoire_sortie)) {
    dir.create(repertoire_sortie, recursive = TRUE)
  }
  
  # Graphique de convergence
  if (!is.null(plots$convergence)) {
    ggsave(file.path(repertoire_sortie, paste0("nn_convergence_criteres_",nb_etats, code, ".pdf")), 
           plots$convergence, width = 10, height = 6)
  }
  
  # Données originales
  if (!is.null(plots$donnees_originales)) {
    pdf(file.path(repertoire_sortie, paste0("nn_donnees_originales.pdf_",nb_etats, code, ".pdf")), 
        width = 8, height = 12)
    do.call(grid.arrange, c(plots$donnees_originales, list(ncol = 1)))
    dev.off()
  }
  
  # Données simulées
  if (!is.null(plots$donnees_simulees)) {
    pdf(file.path(repertoire_sortie, paste0("nn_donnees_simulees.pdf_",nb_etats, code, ".pdf")), 
        width = 8, height = 12)
    do.call(grid.arrange, c(plots$donnees_simulees, list(ncol = 1)))
    dev.off()
  }
  
  cat("✓ Graphiques sauvegardés dans:", repertoire_sortie, "\n")
}

create_state_palette <- function(n_states) {
  "
  Crée une palette de couleurs adaptée au nombre d'états
  "
  palette <- RColorBrewer::brewer.pal(max(3, n_states), CONST$VIZ$state_palette_name)
  if (n_states > 9) {
    palette <- colorRampPalette(palette)(n_states)
  }
  return(palette)
}

# === SECTION 4: UTILITAIRES STATISTIQUES ===

calculate_aic <- function(loglik, k) {
  "Calcul de l'AIC"
  2 * k - 2 * loglik
}

calculate_bic <- function(loglik, k, n) {
  "Calcul du BIC"
  k * log(n) - 2 * loglik
}

extract_model_statistics <- function(hmm_model, data) {
  "
  Extrait les statistiques principales du modèle HMM
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
  Simule des données à partir d'un modèle HMM ajusté
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
    cat("⚠️  Erreur dans la simulation:", e$message, "\n")
    return(NULL)
  })
}

# === SECTION 6: COMPARAISON DE MODÈLES ===

compare_models_summary <- function(results_list, n_states_list) {
  "
  Crée un résumé comparatif des modèles ajustés
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

# === SECTION 7: TARIFICATION - FONCTIONS COMMUNES ===

precompute_mortality_coefficients <- function(hmm_model, n_states, age_labels) {
  "
  Pré-calcule les coefficients de mortalité pour accélération des calculs
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
  Calcule l'intensité de mortalité pour des paramètres donnés
  "
  sin_val <- sin(2 * pi * week / 52)
  cos_val <- cos(2 * pi * week / 52)
  
  age_idx <- which(CONST$PRICING$age_labels == age_range)
  if (length(age_idx) == 0) return(1e-8)
  
  coefs <- coef_matrix[state, , age_idx[1]]
  
  log_death_rate <- coefs[1] + coefs[2] * sin_val + coefs[3] * cos_val + coefs[4] * year
  exp(log_death_rate)  # Retourner la death_rate (intensité)
}

# === EXPORT ET REPORTING ===

export_results_to_csv <- function(results_df, filename) {
  "Exporte les résultats en CSV"
  write.csv(results_df, file = filename, row.names = FALSE)
  cat(sprintf("✓ Résultats exportés: %s\n", filename))
}

export_model_summary <- function(model_results, filename) {
  "Exporte un résumé du modèle"
  summary_text <- sprintf("
RÉSUMÉ DU MODÈLE
================

Nombre d'états: %d
Nombre de convergences: %d/%d
Meilleur AIC: %.2f
Meilleur BIC: %.2f
Meilleur LogLik: %.2f

Itération optimale: %d
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
  cat(sprintf("✓ Résumé exporté: %s\n", filename))
}

cat("✅ Utilitaires chargés avec succès!\n")
