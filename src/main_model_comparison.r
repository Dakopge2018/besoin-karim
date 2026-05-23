# =============================================================================
# MAIN FILE 1: MODEL COMPARISON (main_model_comparison.r)
# Estimates and compares models: Poisson-Normal, Normal-Normal 2 and 3 states
# =============================================================================

# === INITIALIZATION ===
source("./src/const.r")
source("./src/utilities.r")

# Initialize environment
initialize_environment()
theme_set(theme_bw())

# === INTERACTIVE QUESTIONS ===
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("MAIN PROGRAM 1: HMM MODEL COMPARISON\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Ask the questions
CONST$INTERACTIVE$model_type <- ask_model_type()
CONST$INTERACTIVE$temp_variable <- ask_temp_variable()
CONST$INTERACTIVE$plot_results <- ask_plot_results()

# Set up distributions according to model type
if (CONST$INTERACTIVE$model_type == "poisson_normal") {
  CONST$HMM$dists <- list(
    death = "pois"
  )
  CONST$HMM$dists[[CONST$INTERACTIVE$temp_variable]] <- "norm"
  mort_dist_name <- "death"
  cat("📊 Distributions: Poisson-Normal (mortality in Poisson)\n\n")
} else {
  CONST$HMM$dists <- list(
    log_death_rate = "norm"
  )
  CONST$HMM$dists[[CONST$INTERACTIVE$temp_variable]] <- "norm"
  mort_dist_name <- "log_death_rate"
  cat("📊 Distributions: Normal-Normal (mortality in log Normal)\n\n")
}

cat(sprintf("🌡️  Temperature variable: %s\n\n", CONST$INTERACTIVE$temp_variable))

# === SUB-FUNCTION: Fit a single HMM model ===

main_fit_hmm_model <- function(data, 
                               nb_etats = 2, 
                               n_simulations = 1000,
                               plot_results = TRUE,
                               n_ajustements = 10) {
  "
  Helper function: Fits an HMM model with given parameters
  Called from main() for each configuration (2 or 3 states)
  "
  
  cat(sprintf("\n🔧 HMM Fitting: %d states (%d attempts)\n", nb_etats, n_ajustements))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  set.seed(CONST$DATA$seed)
  
  # Internal function to fit a single model
  fit_single_model <- function(data, n_states, par_init, iteration) {
    tryCatch({
      # === Create formulas adapted to model type ===
      
      f_transition <- ~ cos_1 + sin_1
      
      # Observation formulas depending on model type
      f_obs <- list()
      
      if (CONST$INTERACTIVE$model_type == "poisson_normal") {
        # Poisson-Normal model: death in Poisson with offset
        f_obs$death <- list(rate = ~ trend + sin_1 + cos_1 + offset(log_exposure))
      } else {
        # Normal-Normal model: log_death_rate in Normal
        f_obs$log_death_rate <- list(mean = ~ trend + sin_1 + cos_1 + Age_factor, 
                                      sd = ~ 1)
      }
      
      # Add temperature formula
      f_obs[[CONST$INTERACTIVE$temp_variable]] <- list(
        mean = ~ sin_1 + cos_1 + trend,
        sd = ~ 1
      )
      
      # Create HMM objects
      hid <- hmmTMB::MarkovChain$new(
        data = data, 
        n_states = n_states, 
        formula = f_transition
      )
      
      obs <- hmmTMB::Observation$new(
        data = data,
        dists = CONST$HMM$dists,
        n_states = n_states,
        par = par_init,
        formulas = f_obs
      )
      
      # Ajuster le modèle
      hmm <- hmmTMB::HMM$new(obs = obs, hid = hid)
      hmm$fit(silent = TRUE, maxit = CONST$HMM$maxit, tol = CONST$HMM$tol)
      
      # Extraire les statistiques
      loglik <- logLik(hmm)[1]
      aic <- AIC(hmm)
      bic <- BIC(hmm)
      
      return(list(
        modele = hmm,
        loglik = loglik,
        aic = aic,
        bic = bic,
        iteration = iteration,
        convergence = TRUE
      ))
      
    }, error = function(e) {
      cat("  ⚠️ Erreur itération", iteration, ":", substr(e$message, 1, 50), "...\n")
      return(list(
        modele = NULL,
        loglik = -Inf,
        aic = Inf,
        bic = Inf,
        iteration = iteration,
        convergence = FALSE
      ))
    })
  }
  
  # === BOUCLE PRINCIPALE D'AJUSTEMENT ===
  
  resultats <- list()
  pb <- txtProgressBar(min = 0, max = n_ajustements, style = 3)
  
  for (i in 1:n_ajustements) {
    # Générer des paramètres initiaux diversifiés
    set.seed(CONST$DATA$seed + i)
    par_init <- calculer_parametres_initiaux_diversifies_simple(
      data, 
      nb_etats, 
      iteration = i,
      mort_var = mort_dist_name,
      temp_var = CONST$INTERACTIVE$temp_variable,
      mort_dist = CONST$HMM$dists[[mort_dist_name]]
    )
    
    # Ajuster le modèle
    resultat <- fit_single_model(data, nb_etats, par_init, i)
    resultats[[i]] <- resultat
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # === SÉLECTION DU MEILLEUR MODÈLE ===
  
  modeles_valides <- resultats[sapply(resultats, function(x) x$convergence)]
  
  if (length(modeles_valides) == 0) {
    stop(sprintf("❌ Aucun modèle à %d états n'a convergé", nb_etats))
  }
  
  cat(sprintf("\n✓ %d modèles ont convergé sur %d tentatives\n", 
              length(modeles_valides), n_ajustements))
  
  # Sélectionner le meilleur par AIC
  aics <- sapply(modeles_valides, function(x) x$aic)
  indice_optimal <- which.min(aics)
  modele_optimal <- modeles_valides[[indice_optimal]]
  
  cat(sprintf("✓ Meilleur modèle (itération %d)\n", modele_optimal$iteration))
  cat(sprintf("  - AIC: %.2f\n", modele_optimal$aic))
  cat(sprintf("  - BIC: %.2f\n", modele_optimal$bic))
  cat(sprintf("  - LogLik: %.2f\n", modele_optimal$loglik))
  
  # === VISUALISATIONS ===
  
  plots <- list()
  
  if (plot_results && CONST$INTERACTIVE$plot_results) {
    cat("  📊 Génération des graphiques...\n")
    
    # Obtenir les états prédits
    etats_predits <- modele_optimal$modele$viterbi()
    data$etat <- etats_predits
    
    # Ajouter la variable temps si elle n'existe pas
    if (!"time" %in% colnames(data)) {
      data$time <- seq_len(nrow(data))
    }
    
    # Palette de couleurs
    etat_palette <- create_state_palette(nb_etats)
    names(etat_palette) <- as.character(1:nb_etats)
    
    # 1. Graphique de convergence des critères
    resultats_ajustement <- data.frame(
      iteration = sapply(modeles_valides, function(x) x$iteration),
      aic = sapply(modeles_valides, function(x) x$aic),
      bic = sapply(modeles_valides, function(x) x$bic),
      loglik = sapply(modeles_valides, function(x) x$loglik)
    )
    
    if (nrow(resultats_ajustement) > 1) {
      p_convergence <- ggplot(resultats_ajustement, aes(x = iteration)) +
        geom_point(aes(y = bic, color = "BIC"), alpha = 0.6, size = 2) +
        geom_point(aes(y = aic, color = "AIC"), alpha = 0.6, size = 2) +
        labs(title = sprintf("Information criteria - %d states", nb_etats),
             x = "Iteration", y = "Criterion value", color = "Criterion") +
        theme_minimal()
      
      plots$convergence <- p_convergence
      print(p_convergence)
    }
    
    # 2. Visualisation des données originales avec états prédits
    p1 <- ggplot(data, aes(x = time, y = factor(etat), color = factor(etat))) +
      geom_point(size = 3) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Predicted hidden states (original data)", 
           x = "Time", y = "State", color = "State") +
      theme_minimal()
    
    # Graphique mortality rate / log_death_rate
    if (CONST$INTERACTIVE$model_type == "poisson_normal") {
      p2 <- ggplot(data, aes(x = time, y = death, color = factor(etat))) +
        geom_point(size = 2) +
        scale_color_manual(values = etat_palette) +
        labs(title = "Death counts by state (original data)", 
             x = "Time", y = "Death counts", color = "State") +
        theme_minimal()
    } else {
      p2 <- ggplot(data, aes(x = time, y = log_death_rate, color = factor(etat))) +
        geom_point(size = 2) +
        scale_color_manual(values = etat_palette) +
        labs(title = "Log mortality rate by state (original data)", 
             x = "Time", y = "Log mortality rate", color = "State") +
        theme_minimal()
    }

    # Graphique température normale
    p3 <- ggplot(data, aes_string(x = "time", y = "temp_norm", color = "factor(etat)")) +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Normal temperature by state (original data)", 
           x = "Time", y = "Temperature (°C)", color = "State") +
      theme_minimal()  
    
    # Graphique température extrême
    p4_extreme <- ggplot(data, aes_string(x = "time", y = "temp_extreme", color = "factor(etat)")) +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Extreme temperature by state (original data)", 
           x = "Time", y = "Temperature (°C)", color = "State") +
      theme_minimal()
    
    plots$donnees_originales <- list(p1, p2, p3, p4_extreme)
    
    # Affichage
    if (length(plots$donnees_originales) > 0) {
      print(do.call(gridExtra::grid.arrange, c(plots$donnees_originales, list(ncol = 1))))
    }
  }
  
  # === RETOUR ===
  
  return(list(
    modele_optimal = modele_optimal$modele,
    tous_resultats = resultats,
    modeles_valides = modeles_valides,
    plots = plots,
    resume = list(
      n_ajustements = n_ajustements,
      n_converges = length(modeles_valides),
      meilleur_aic = modele_optimal$aic,
      meilleur_bic = modele_optimal$bic,
      meilleur_loglik = modele_optimal$loglik,
      iteration_optimale = modele_optimal$iteration,
      nb_etats = nb_etats
    )
  ))
}

# === FONCTION PRINCIPALE: Comparer les modèles ===

main <- function() {
  "
  Fonction principale:
  1. Charge et prépare les données
  2. Estime les modèles pour 2 et 3 états
  3. Compare les résultats
  4. Génère les visualisations et rapports
  "
  
  cat("\n📊 ÉTAPE 1: CHARGEMENT ET PRÉPARATION DES DONNÉES\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  # Charger les données brutes
  data_raw <- load_and_prepare_data()
  cat(sprintf("✓ Données chargées: %d lignes\n", nrow(data_raw)))
  
  # Préparer le dataframe HMM
  df <- prepare_hmm_dataframe(data_raw)
  cat(sprintf("✓ Dataframe HMM créé: %d lignes, %d colonnes\n", nrow(df), ncol(df)))
  
  # Afficher les statistiques descriptives
  cat("\n📈 Statistiques descriptives:\n")
  cat(sprintf("  • Log(taux de mortalité) - Moyenne: %.4f, SD: %.4f\n", 
              mean(df$log_death_rate, na.rm = T), sd(df$log_death_rate, na.rm = T)))
  cat(sprintf("  • Température - Moyenne: %.2f, SD: %.2f\n", 
              mean(df$temp_extreme, na.rm = T), sd(df$temp_extreme, na.rm = T)))
  
  # Sauvegarder le dataframe pour utilisation en tarification
  saveRDS(df, "./data/df_prepared.rds")
  cat("✓ Dataframe sauvegardé: ./data/df_prepared.rds\n")
  
  # === ÉTAPE 2: AJUSTEMENT DES MODÈLES ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("📊 ÉTAPE 2: AJUSTEMENT DES MODÈLES HMM\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")

  dists = CONST$HMM$dists
  
  # Ajuster pour 2 états
  cat("\n\n🔹 MODÈLE À 2 ÉTATS\n")
  results_2states <- main_fit_hmm_model(
    df,
    nb_etats = 2,
    n_simulations = 1000,
    plot_results = TRUE,
    n_ajustements = CONST$HMM$n_ajustements
  )
  
  # Ajuster pour 3 états
  cat("\n\n🔹 MODÈLE À 3 ÉTATS\n")
  results_3states <- main_fit_hmm_model(
    df,
    nb_etats = 3,
    n_simulations = 1000,
    plot_results = TRUE,
    n_ajustements = CONST$HMM$n_ajustements
  )
  
  # === ÉTAPE 3: COMPARAISON DES MODÈLES ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("📊 ÉTAPE 3: COMPARAISON DES MODÈLES\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  comparison_df <- compare_models_summary(
    list(results_2states, results_3states),
    c(2, 3)
  )
  
  cat("📋 TABLEAU COMPARATIF:\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  print(comparison_df)
  cat("\n")
  
  # Analyse
  if (comparison_df$aic[1] < comparison_df$aic[2]) {
    cat("✓ MEILLEUR MODÈLE: 2 états (AIC plus faible)\n")
  } else {
    cat("✓ MEILLEUR MODÈLE: 3 états (AIC plus faible)\n")
  }
  
  # Exporter le tableau de comparaison
  export_results_to_csv(comparison_df, "./results/model_comparison.csv")
  
  # === ÉTAPE 4: VISUALISATIONS ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("📊 ÉTAPE 4: GÉNÉRATION DES VISUALISATIONS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Graphique de convergence pour 2 états
  if (!is.null(results_2states$plots$convergence) && CONST$INTERACTIVE$plot_results) {
    print(results_2states$plots$convergence)
    ggsave("./results/convergence_2states.pdf", width = 10, height = 6, dpi = 300)
  }
  
  # Graphique de convergence pour 3 états
  if (!is.null(results_3states$plots$convergence) && CONST$INTERACTIVE$plot_results) {
    print(results_3states$plots$convergence)
    ggsave("./results/convergence_3states.pdf", width = 10, height = 6, dpi = 300)
  }
  
  # Graphique comparatif AIC/BIC
  p_comparison <- ggplot(comparison_df, aes(x = as.factor(n_states))) +
    geom_bar(aes(y = aic, fill = "AIC"), stat = "identity", position = "dodge", alpha = 0.7) +
    geom_bar(aes(y = bic, fill = "BIC"), stat = "identity", position = "dodge", alpha = 0.7) +
    scale_fill_manual(values = c("AIC" = "#1f77b4", "BIC" = "#ff7f0e")) +
    labs(title = "Comparison of information criteria",
         x = "Number of states", y = "Criterion value") +
    theme_minimal()
  
  if (CONST$INTERACTIVE$plot_results) {
    print(p_comparison)
    ggsave("./results/comparison_criteria.pdf", width = 10, height = 6, dpi = 300)
  }
  
  # === ÉTAPE 5: RAPPORTS ET EXPORT ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("📊 ÉTAPE 5: SAUVEGARDE DES RÉSULTATS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Créer un suffixe descriptif incluant le type de modèle et la variable de température
  model_suffix <- paste0(
    CONST$INTERACTIVE$model_type, "_",
    CONST$INTERACTIVE$temp_variable
  )
  
  # Exporter les modèles optimaux avec noms différenciés
  file_model_2states <- sprintf("./models/hmm_model_2states_%s.rds", model_suffix)
  file_model_3states <- sprintf("./models/hmm_model_3states_%s.rds", model_suffix)
  
  saveRDS(results_2states$modele_optimal, file_model_2states)
  cat(sprintf("✓ Modèle 2 états sauvegardé: %s\n", file_model_2states))
  
  saveRDS(results_3states$modele_optimal, file_model_3states)
  cat(sprintf("✓ Modèle 3 états sauvegardé: %s\n", file_model_3states))
  
  # Exporter les résumés
  export_model_summary(results_2states, "./results/summary_2states.txt")
  export_model_summary(results_3states, "./results/summary_3states.txt")
  
  # Exporter les résultats complets
  saveRDS(list(
    results_2states = results_2states,
    results_3states = results_3states,
    comparison = comparison_df
  ), "./results/all_results.rds")
  
  cat("✓ Tous les résultats sauvegardés dans ./results/\n")
  
  # === RÉSUMÉ FINAL ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("✅ PROGRAMME TERMINÉ AVEC SUCCÈS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  cat("� Configuration utilisée:\n")
  cat(sprintf("  • Modèle: %s\n", CONST$INTERACTIVE$model_type))
  cat(sprintf("  • Variable de température: %s\n", CONST$INTERACTIVE$temp_variable))
  cat(sprintf("  • Graphiques: %s\n", ifelse(CONST$INTERACTIVE$plot_results, "Oui", "Non")))
  
  cat("\n📁 Fichiers générés:\n")
  cat("  • Modèles: ./models/hmm_model_*.rds\n")
  cat("  • Résultats: ./results/\n")
  if (CONST$INTERACTIVE$plot_results) {
    cat("  • Graphiques: ./results/*.pdf\n")
  }
  cat("\n")
  
  cat("📊 Résumé des résultats:\n")
  cat(sprintf("  • Modèle 2 états: AIC = %.2f, BIC = %.2f\n", 
              comparison_df$aic[1], comparison_df$bic[1]))
  cat(sprintf("  • Modèle 3 états: AIC = %.2f, BIC = %.2f\n", 
              comparison_df$aic[2], comparison_df$bic[2]))
  
  # Retour invisibles pour éviter affichage automatique
  invisible(list(
    df = df,
    results_2states = results_2states,
    results_3states = results_3states,
    comparison = comparison_df
  ))
}

# === EXÉCUTION DU PROGRAMME ===

if (!dir.exists("./results")) dir.create("./results", recursive = T)
if (!dir.exists("./models")) dir.create("./models", recursive = T)
if (!dir.exists("./data")) dir.create("./data", recursive = T)

# Exécuter le programme principal
results <- main()

cat("\n💡 POUR UTILISER LES RÉSULTATS:\n")
cat("   results_2states <- results$results_2states\n")
cat("   results_3states <- results$results_3states\n")
cat("   df <- results$df\n")
cat("   comparison <- results$comparison\n\n")

cat("🔧 PARAMÈTRES UTILISÉS:\n")
cat(sprintf("   CONST$INTERACTIVE$model_type = '%s'\n", CONST$INTERACTIVE$model_type))
cat(sprintf("   CONST$INTERACTIVE$temp_variable = '%s'\n", CONST$INTERACTIVE$temp_variable))
cat(sprintf("   CONST$INTERACTIVE$plot_results = %s\n\n", CONST$INTERACTIVE$plot_results))
