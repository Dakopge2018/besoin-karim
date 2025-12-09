# =============================================================================
# FICHIER PRINCIPAL 1: COMPARAISON DE MODÈLES (main_model_comparison.r)
# Estime et compare les modèles: Poisson-Normal, Normal-Normal 2 et 3 états
# =============================================================================

# === INITIALISATION ===
source("./src/const.r")
source("./src/utilities.r")

# Initialiser l'environnement
initialize_environment()
theme_set(theme_bw())

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PROGRAMME PRINCIPAL 1: COMPARAISON DES MODÈLES HMM\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# === SOUS-FONCTION: Ajuster un modèle HMM unique ===

main_fit_hmm_model <- function(data, 
                               nb_etats = 2, 
                               n_simulations = 1000,
                               plot_results = TRUE,
                               n_ajustements = 10) {
  "
  Fonction auxiliaire: Ajuste un modèle HMM avec les paramètres donnés
  Appelée depuis le main() pour chaque configuration (2 ou 3 états)
  "
  
  cat(sprintf("\n🔧 Ajustement HMM: %d états (%d tentatives)\n", nb_etats, n_ajustements))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  set.seed(CONST$DATA$seed)
  
  # Fonction interne pour ajuster un seul modèle
  fit_single_model <- function(data, n_states, par_init, iteration) {
    tryCatch({
      # Créer les objets HMM
      hid <- hmmTMB::MarkovChain$new(
        data = data, 
        n_states = n_states, 
        formula = CONST$FORMULAS$transition
      )
      
      obs <- hmmTMB::Observation$new(
        data = data,
        dists = CONST$HMM$dists,
        n_states = n_states,
        par = par_init,
        formulas = CONST$FORMULAS$observation
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
      iteration = i
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
  
  if (plot_results) {
    # 1. Convergence des critères
    resultats_ajustement <- data.frame(
      iteration = sapply(modeles_valides, function(x) x$iteration),
      aic = sapply(modeles_valides, function(x) x$aic),
      bic = sapply(modeles_valides, function(x) x$bic),
      loglik = sapply(modeles_valides, function(x) x$loglik)
    )
    
    if (nrow(resultats_ajustement) > 1) {
      p_convergence <- ggplot(resultats_ajustement, aes(x = iteration)) +
        geom_line(aes(y = bic, color = "BIC"), size = 1) +
        geom_line(aes(y = aic, color = "AIC"), size = 1) +
        geom_point(aes(y = bic, color = "BIC"), alpha = 0.6, size = 2) +
        geom_point(aes(y = aic, color = "AIC"), alpha = 0.6, size = 2) +
        labs(title = sprintf("Critères d'information - %d états", nb_etats),
             x = "Itération", y = "Valeur du critère", color = "Critère") +
        theme_minimal()
      
      plots$convergence <- p_convergence
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
  if (!is.null(results_2states$plots$convergence)) {
    print(results_2states$plots$convergence)
    ggsave("./results/convergence_2states.pdf", width = 10, height = 6, dpi = 300)
  }
  
  # Graphique de convergence pour 3 états
  if (!is.null(results_3states$plots$convergence)) {
    print(results_3states$plots$convergence)
    ggsave("./results/convergence_3states.pdf", width = 10, height = 6, dpi = 300)
  }
  
  # Graphique comparatif AIC/BIC
  p_comparison <- ggplot(comparison_df, aes(x = as.factor(n_states))) +
    geom_bar(aes(y = aic, fill = "AIC"), stat = "identity", position = "dodge", alpha = 0.7) +
    geom_bar(aes(y = bic, fill = "BIC"), stat = "identity", position = "dodge", alpha = 0.7) +
    scale_fill_manual(values = c("AIC" = "#1f77b4", "BIC" = "#ff7f0e")) +
    labs(title = "Comparaison des critères d'information",
         x = "Nombre d'états", y = "Valeur du critère") +
    theme_minimal()
  
  print(p_comparison)
  ggsave("./results/comparison_criteria.pdf", width = 10, height = 6, dpi = 300)
  
  # === ÉTAPE 5: RAPPORTS ET EXPORT ===
  
  cat("\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("📊 ÉTAPE 5: SAUVEGARDE DES RÉSULTATS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Exporter les modèles optimaux
  saveRDS(results_2states$modele_optimal, "./models/hmm_model_2states.rds")
  cat("✓ Modèle 2 états sauvegardé: ./models/hmm_model_2states.rds\n")
  
  saveRDS(results_3states$modele_optimal, "./models/hmm_model_3states.rds")
  cat("✓ Modèle 3 états sauvegardé: ./models/hmm_model_3states.rds\n")
  
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
  
  cat("📁 Fichiers générés:\n")
  cat("  • Modèles: ./models/hmm_model_*.rds\n")
  cat("  • Résultats: ./results/\n")
  cat("  • Graphiques: ./results/*.pdf\n\n")
  
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
cat("   df <- results$df\n\n")
