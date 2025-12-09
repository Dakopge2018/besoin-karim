#Dataset with temperature (not extreme temperatures)
# Load required libraries
library(readxl)
library(ggplot2)
library(scico)
library(hmmTMB)
library(gridExtra)
library(RColorBrewer)
library(grid)
library(dplyr)
library(pracma) # pour l'intégration numérique
library(expm)   # pour les exponentielles de matrices
pal <- hmmTMB:::hmmTMB_cols


# Set theme for plots
theme_set(theme_bw())
# Importation des données
set.seed(123)
states <- 4
degree_obs_pol <- 1
degree_trans_pol <- 1
period <- 52


#------------------------------------------------
# Importation des données
path_to_data = "C:\\Users\\samue\\OneDrive\\Documents\\cours\\Projet de mémoire\\git\\real_data_code\\data\\main_database.xlsx"
OUTER_IMAGE = "C:\\Users\\samue\\OneDrive\\Documents\\cours\\Projet de mémoire\\git\\real_data_code\\image\\nnnn\\hmm_tmb\\normal"

data <- read_excel(path_to_data, sheet = "database_4")

# Check for and handle missing or infinite values
data <- data[!is.na(data$Death_counts) & !is.infinite(data$Death_counts) & data$Death_counts > 0, ]
data <- na.omit(data)
data <- data[is.finite(rowSums(data[c('year','No_year', 'Temperature', 'Death_counts', 'Weekly_exposure')])), ]

generate_trig_covariates <- function(week_no, year_no, period, degree) {
  trig_covs <- data.frame(time = week_no, trend = year_no)
  for (d in 1:degree) {
    trig_covs[[paste0("cos_", d)]] <- cos(2 * pi * d * week_no / period)
    trig_covs[[paste0("sin_", d)]] <- sin(2 * pi * d * week_no / period)
  }
  return(trig_covs)
}
trig_covs <- generate_trig_covariates(data$No_week, data$No_year,period, degree_trans_pol)

df = data.frame(
    # Code = data$Code,
    Code = 'ITI43',
    death = data$Death_counts,
    Age_factor = as.factor(data$Age_group),
    exposure = data$Weekly_exposure,
    log_death_rate = log(data$Death_counts / data$Weekly_exposure),
    death_rate = data$Death_counts / data$Weekly_exposure,
    temp_extreme = data$Temp_extrem,
    # temp_extreme = data$Temperature,
    temp_norm = data$Temperature,
    log_exposure = log(data$Weekly_exposure),
    trend = data$No_year,
    trig_covs)

df <- na.omit(df)


sauvegarder_graphiques <- function(plots, repertoire_sortie = getwd(), nb_etats, code) {
  if (!dir.exists(repertoire_sortie)) {
    dir.create(repertoire_sortie, recursive = TRUE)
  }
  
  # Sauvegarder le graphique de convergence
  if (!is.null(plots$convergence)) {
    # ggsave(file.path(repertoire_sortie, "convergence_criteres.png"), 
    #        plots$convergence, width = 12, height = 8, dpi = 300)
    ggsave(file.path(repertoire_sortie, paste0("nn_convergence_criteres_",nb_etats, code, ".pdf")), 
           plots$convergence, width = 10, height = 6)
  }
  
  # Sauvegarder les graphiques des données originales
  if (!is.null(plots$donnees_originales)) {
    # png(file.path(repertoire_sortie, "donnees_originales.png"), 
    #     width = 1500, height = 1000)
    # do.call(grid.arrange, c(plots$donnees_originales, list(ncol = 1)))
    # dev.off()
    
    pdf(file.path(repertoire_sortie, paste0("nn_donnees_originales.pdf_",nb_etats, code, ".pdf")), 
        width = 8, height = 12)
    do.call(grid.arrange, c(plots$donnees_originales, list(ncol = 1)))
    dev.off()
  }
  
  # Sauvegarder les graphiques des données simulées
  if (!is.null(plots$donnees_simulees)) {
    # png(file.path(repertoire_sortie, "donnees_simulees.png"), 
    #     width = 1500, height = 1000)
    # do.call(grid.arrange, c(plots$donnees_simulees, list(ncol = 1)))
    # dev.off()
    
    pdf(file.path(repertoire_sortie, paste0("nn_donnees_simulees.pdf_",nb_etats, code, ".pdf")), 
        width = 8, height = 12)
    do.call(grid.arrange, c(plots$donnees_simulees, list(ncol = 1)))
    dev.off()
  }
  
  cat("Graphiques sauvegardés dans:", repertoire_sortie, "\n")
}
# Fonction 1: calculer_parametres_initiaux_diversifies_simple
calculer_parametres_initiaux_diversifies_simple <- function(data, n_states, iteration = 1, age_ref = NULL) {
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
      mask <- data$temp_extreme > quantiles_temp[i] & 
              data$temp_extreme <= quantiles_temp[i + 1]
    }

    temp_means[i] <- mean(data$temp_extreme[mask], na.rm = TRUE)
    temp_sds[i] <- sd(data$temp_extreme[mask], na.rm = TRUE)

    if (is.na(temp_means[i])) temp_means[i] <- mean_temp
    if (is.na(temp_sds[i]) || temp_sds[i] <= 0) temp_sds[i] <- sd_temp
  }

  # Stratégies pour la mortalité
  strategy <- iteration %% 4

  if (strategy == 1) {
    cat("  -> Stratégie empirique\n")
    for (i in 1:n_states) {
      if (i == 1) {
        mask <- data$temp_extreme <= quantiles_temp[i + 1]
      } else if (i == n_states) {
        mask <- data$temp_extreme > quantiles_temp[i]
      } else {
        mask <- data$temp_extreme > quantiles_temp[i] & 
                data$temp_extreme <= quantiles_temp[i + 1]
      }

      observed_mean <- mean(data$log_death_rate[mask], na.rm = TRUE)
      observed_sd <- sd(data$log_death_rate[mask], na.rm = TRUE)
      
      if (is.na(observed_mean) || observed_mean <= 0) {
        observed_mean <- mean_death
      }
      if (is.na(observed_sd) || observed_sd <= 0) {
        observed_sd <- sd_death
      }
      
      death_means[i] <- max(observed_mean * runif(1, 0.9, 1.1), 0.0001)
      death_sds[i] <- max(observed_sd * runif(1, 0.9, 1.1), 0.0001)
    }

  } else if (strategy == 2) {
    cat("  -> Stratégie aléatoire\n")
    min_mean <- max(mean_death * 0.5, 0.0001)
    max_mean <- mean_death * 2
    death_means <- runif(n_states, min_mean, max_mean)
    
    min_sd <- max(sd_death * 0.5, 0.0001)
    max_sd <- sd_death * 2
    death_sds <- runif(n_states, min_sd, max_sd)

  } else if (strategy == 3) {
    cat("  -> Stratégie permutation\n")
    empirical_means <- numeric(n_states)
    empirical_sds <- numeric(n_states)
    
    for (i in 1:n_states) {
      if (i == 1) {
        mask <- data$temp_extreme <= quantiles_temp[i + 1]
      } else if (i == n_states) {
        mask <- data$temp_extreme > quantiles_temp[i]
      } else {
        mask <- data$temp_extreme > quantiles_temp[i] & 
                data$temp_extreme <= quantiles_temp[i + 1]
      }

      empirical_means[i] <- mean(data$log_death_rate[mask], na.rm = TRUE)
      empirical_sds[i] <- sd(data$log_death_rate[mask], na.rm = TRUE)
      
      if (is.na(empirical_means[i]) || empirical_means[i] <= 0) {
        empirical_means[i] <- mean_death
      }
      if (is.na(empirical_sds[i]) || empirical_sds[i] <= 0) {
        empirical_sds[i] <- sd_death
      }
    }

    death_means <- sample(empirical_means) * runif(n_states, 0.8, 1.2)
    death_means <- pmax(death_means, 0.0001)
    
    death_sds <- sample(empirical_sds) * runif(n_states, 0.8, 1.2)
    death_sds <- pmax(death_sds, 0.0001)

  } else {
    cat("  -> Stratégie uniforme\n")
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

# Fonction 2: ajuster_hmm_optimal
ajuster_hmm_optimal <- function(donnees, 
                                 nb_etats = 2, 
                                 n_simulations = 1000,
                                 seed = 123, 
                                 maxit = 1000, 
                                 tol = 1e-6,
                                 plot_results = TRUE, 
                                 n_ajustements = 10,
                                 sauvegarder_plots = FALSE,
                                 repertoire_sortie = getwd(),
                                 age_ref = NULL)
                                  {
  
  # Définir le seed global
  set.seed(seed)
  
  # Fonction pour ajuster un modèle HMM
  ajuster_modele_unique <- function(data, n_states, par_init, iteration) {
    tryCatch({
      # Formules
      f_transition <- ~ cos_1 + sin_1
      f_obs <- list(
        # log_death_rate = list(mean = ~ trend + sin_1 + cos_1, 
        log_death_rate = list(mean = ~ sin_1 + cos_1 + trend + Age_factor, 
                              sd = ~ 1),
        temp_extreme = list(mean = ~ sin_1 + cos_1 + trend,
                           sd = ~ 1)
      )
      
      # Créer les objets MarkovChain et Observation
      hid <- MarkovChain$new(data = data, n_states = n_states, formula = f_transition)
      obs <- Observation$new(data = data,
                            dists = dists,  # Assumé défini globalement
                            n_states = n_states,
                            par = par_init,
                            formulas = f_obs)
      
      # Créer et ajuster le HMM
      hmm <- HMM$new(obs = obs, hid = hid)
      hmm$fit(silent = TRUE, maxit = maxit, tol = tol)
      
      # Calculer les critères de sélection
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
      cat("Erreur dans l'ajustement", iteration, ":", e$message, "\n")
      return(list(
        modele = NULL,
        loglik = -Inf,
        aic = Inf,
        bic = Inf,
        iteration = iteration,
        convergence = FALSE,
        erreur = e$message
      ))
    })
  }
  
  # Fonction de simulation après ajustement
  simuler_donnees <- function(modele_optimal, n_sim, data_template) {
    if (is.null(modele_optimal)) {
      cat("Aucun modèle valide pour la simulation\n")
      return(NULL)
    }
    tryCatch({
      # Créer un template de données pour la simulation avec les bonnes variables
      template_sim <- data_template[1:n_sim, ]
      
      # Si on a moins de lignes que nécessaire, répéter les données
      if (nrow(data_template) < n_sim) {
        indices_repeat <- rep(seq_len(nrow(data_template)), length.out = n_sim)
        template_sim <- data_template[indices_repeat, ]
      }
    
      simulations <- modele_optimal$simulate(n = n_sim, data = template_sim)
      return(simulations)
    }, error = function(e) {
      cat("Erreur dans la simulation:", e$message, "\n")
      return(NULL)
    })
  }
  
  # AJUSTEMENT PRINCIPAL
  cat("Début de l'ajustement avec", n_ajustements, "tentatives...\n")
  
  resultats <- list()
  
  # Barre de progression
  pb <- txtProgressBar(min = 0, max = n_ajustements, style = 3)
  
  for (i in 1:n_ajustements) {
    # Générer des paramètres initiaux avec la fonction diversifiée
    set.seed(seed + i)
    par_init <- calculer_parametres_initiaux_diversifies_simple(
      donnees, 
      nb_etats, 
      iteration = i,
      age_ref = age_ref
    )
    
    # Ajuster le modèle
    resultat <- ajuster_modele_unique(donnees, nb_etats, par_init, i)
    resultats[[i]] <- resultat
    
    # Mettre à jour la barre de progression
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Filtrer les modèles convergés
  modeles_valides <- resultats[sapply(resultats, function(x) x$convergence)]
  
  if (length(modeles_valides) == 0) {
    stop("Aucun modèle n'a convergé. Essayez d'augmenter maxit ou de modifier les paramètres.")
  }
  
  cat("\n", length(modeles_valides), "modèles ont convergé sur", n_ajustements, "tentatives.\n")
  
  # Sélectionner le meilleur modèle (plus faible AIC)
  aics <- sapply(modeles_valides, function(x) x$aic)
  indice_optimal <- which.min(aics)
  modele_optimal <- modeles_valides[[indice_optimal]]
  
  cat("Meilleur modèle: Itération", modele_optimal$iteration, "\n")
  cat("LogLik:", round(modele_optimal$loglik, 2), "\n")
  cat("AIC:", round(modele_optimal$aic, 2), "\n")
  cat("BIC:", round(modele_optimal$bic, 2), "\n")
  
  # Simulation
  cat("\nSimulation de", n_simulations, "observations...\n")
  # Df_2 = age_Factor = 65-69
  df_sub <- subset(df, Age_factor == "65-69")
  df_2 <- df_sub[, c("time", "sin_1", "cos_1", 'Age_factor',"trend", "exposure", "log_exposure")]
  donnees_simulees <- simuler_donnees(modele_optimal$modele, n_simulations, df_2)
  if (!is.null(donnees_simulees)) {
    donnees_simulees$states <- attr(donnees_simulees, "state")
  }
  
  # Création des graphiques avancés
  plots <- list()
  
  if (plot_results && !is.null(donnees_simulees)) {
    # Charger les librairies nécessaires pour les visualisations
    if (!require(RColorBrewer, quietly = TRUE)) {
      install.packages("RColorBrewer")
      library(RColorBrewer)
    }
    if (!require(gridExtra, quietly = TRUE)) {
      install.packages("gridExtra")
      library(gridExtra)
    }
    
    # Obtenir les états prédits pour les données originales
    etats_predits <- modele_optimal$modele$viterbi()
    donnees$etat <- etats_predits
    donnees_2 <- subset(donnees, Age_factor == "+90")


    
    # Créer une variable temporelle si elle n'existe pas
    if (!"date" %in% colnames(donnees_2) && !"time" %in% colnames(donnees_2)) {
      donnees_2$date <- seq_len(nrow(donnees_2))
    }
    if (!"time" %in% colnames(donnees)) {
      donnees_2$time <- donnees_2$date
    }
    
    # Préparer les données simulées
    if (is.data.frame(donnees_simulees)) {
      donnees_simulees$simulation_id <- 1
      donnees_simulees$states <- donnees_simulees$state
      donnees_simulees$time <- seq_len(nrow(donnees_simulees))
    }
    
    # Création d'une palette de couleurs adaptée au nombre d'états
    etat_palette <- brewer.pal(max(3, nb_etats), "Set1")
    if(nb_etats > 9) {
      etat_palette <- colorRampPalette(etat_palette)(nb_etats)
    }
    names(etat_palette) <- as.character(1:nb_etats)
    
    # Graphique de convergence des critères
    resultats_ajustement <- data.frame(
      iteration = sapply(modeles_valides, function(x) x$iteration),
      aic = sapply(modeles_valides, function(x) x$aic),
      bic = sapply(modeles_valides, function(x) x$bic),
      loglik = sapply(modeles_valides, function(x) x$loglik)
    )
    
    if (nrow(resultats_ajustement) > 1) {
      p_convergence <- ggplot(resultats_ajustement, aes(x = iteration)) +
        geom_line(aes(y = bic, color = "BIC")) +
        geom_line(aes(y = aic, color = "AIC")) +
        geom_point(aes(y = bic, color = "BIC"), alpha = 0.6) +
        geom_point(aes(y = aic, color = "AIC"), alpha = 0.6) +
        labs(title = "Évolution des critères d'information lors des ajustements",
             x = "Itération", y = "Valeur du critère", color = "Critère") +
        theme_minimal()
      
      plots$convergence <- p_convergence
      print(p_convergence)
    }
    
    # Visualisation des données originales avec états prédits
    p1 <- ggplot(donnees_2, aes(x = time, y = factor(etat), color = factor(etat))) +
      geom_point(size = 3) +
      scale_color_manual(values = etat_palette) +
      labs(title = "États cachés prédits (données originales)", 
           x = "Temps", y = "État", color = "État") +
      theme_minimal()
    
    p2 <- ggplot(donnees_2, aes(x = time, y = log_death_rate, color = factor(etat))) +
      geom_line() +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Log de Taux de mortalité par état (données originales)", 
           x = "Temps", y = "Log de Taux de mortalité", color = "État") +
      theme_minimal()

    # p3_bis <- ggplot(donnees_2, aes(x = time, y = temp_extreme, color = factor(etat))) +
    #   geom_line() +
    #   geom_point(size = 2) +
    #   scale_color_manual(values = etat_palette) +
    #   labs(title = "Température normale par état (données originales)", 
    #        x = "Temps", y = "Température (°C)", color = "État") +
    #   theme_minimal()  
    
    p3 <- ggplot(donnees_2, aes(x = time, y = temp_extreme, color = factor(etat))) +
      geom_line() +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Température extrême par état (données originales)", 
           x = "Temps", y = "Température (°C)", color = "État") +
      theme_minimal()
    
    # plots$donnees_originales <- list(p1, p2, p3_bis, p3)
    plots$donnees_originales <- list(p1, p2, p3)



    # Visualisation des données simulées
    if (is.data.frame(donnees_simulees)) {
      sim1 <- subset(donnees_simulees, simulation_id == 1)
      
      p4 <- ggplot(sim1, aes(x = time, y = factor(states), color = factor(states))) +
        geom_point(size = 3) +
        scale_color_manual(values = etat_palette) +
        labs(title = "États cachés simulés", 
             x = "Temps", y = "État", color = "État") +
        theme_minimal()
      
      p5 <- ggplot(sim1, aes(x = time, y = log_death_rate, color = factor(states))) +
        geom_line() +
        geom_point(size = 2) +
        scale_color_manual(values = etat_palette) +
        labs(title = "Log de Taux de mortalité simulé", 
             x = "Temps", y = "Log de Taux de mortalité", color = "État") +
        theme_minimal()
      
      p6 <- ggplot(sim1, aes(x = time, y = temp_extreme, color = factor(states))) +
        geom_line() +
        geom_point(size = 2) +
        scale_color_manual(values = etat_palette) +
        labs(title = "Température simulée", 
             x = "Temps", y = "Température (°C)", color = "État") +
        theme_minimal()
      
      plots$donnees_simulees <- list(p4, p5, p6)
      
      # Comparaison entre réel et simulé
      donnees_temp <- donnees
      donnees_temp$type <- "Réel"
      sim1$type <- "Simulé"
      
      colonnes_communes <- intersect(c("time", "log_death_rate", "temp_extreme", "type"), 
                                   intersect(colnames(donnees_temp), colnames(sim1)))
      
      if (length(colonnes_communes) >= 3) {
        donnees_comparaison <- rbind(
          donnees_temp[, colonnes_communes],
          sim1[, colonnes_communes]
        )
        
        p7 <- ggplot(donnees_comparaison, aes(x = time, y = log_death_rate, color = type)) +
          geom_line() +
          geom_point(size = 2) +
          facet_wrap(~type) +
          labs(title = "Comparaison du log de taux de mortalité: réel vs simulé", 
               x = "Temps", y = "Log de Taux de mortalité", color = "Type de données") +
          theme_minimal()
        
        p8 <- ggplot(donnees_comparaison, aes(x = time, y = temp_extreme, color = type)) +
          geom_line() +
          geom_point(size = 2) +
          facet_wrap(~type) +
          labs(title = "Comparaison des températures: réel vs simulé", 
               x = "Temps", y = "Température (°C)", color = "Type de données") +
          theme_minimal()
        
        plots$comparaison <- list(p7, p8)
      }
      
      # Distributions des valeurs par état
      if ("etat" %in% colnames(donnees_2) && "states" %in% colnames(sim1)) {
        y_limits <- range(c(donnees_2$log_death_rate, sim1$log_death_rate), na.rm = TRUE)
        
        p9 <- ggplot(donnees_2, aes(x = factor(etat), y = log_death_rate, fill = factor(etat))) +
          geom_boxplot() +
          scale_fill_manual(values = etat_palette) +
          labs(title = "Distribution du taux de mortalité par état (données réelles)", 
               x = "État", y = "Taux de mortalité", fill = "État") +
          theme_minimal() +
          ylim(y_limits)
        
        p10 <- ggplot(sim1, aes(x = factor(states), y = log_death_rate, fill = factor(states))) +
          geom_boxplot() +
          scale_fill_manual(values = etat_palette) +
          labs(title = "Distribution du Log de taux de mortalité par état (données simulées)", 
               x = "État", y = "Log de Taux de mortalité", fill = "État") +
          theme_minimal() +
          ylim(y_limits)
        
        plots$distributions <- list(p9, p10)
      }
    }
    
    # Affichage des graphiques principaux
    cat("Affichage des visualisations...\n")
    if (length(plots$donnees_originales) > 0) {
      print(do.call(grid.arrange, c(plots$donnees_originales, list(ncol = 1))))
    }
    
    if (length(plots$donnees_simulees) > 0) {
      print(do.call(grid.arrange, c(plots$donnees_simulees, list(ncol = 1))))
    }
    
    if (length(plots$comparaison) > 0) {
      print(do.call(grid.arrange, c(plots$comparaison, list(ncol = 2))))
    }
    
    if (length(plots$distributions) > 0) {
      print(do.call(grid.arrange, c(plots$distributions, list(ncol = 2))))
    }
    
    # Sauvegarder les graphiques si demandé
    if (sauvegarder_plots && length(plots) > 0) {
      sauvegarder_graphiques(plots, repertoire_sortie, nb_etats, donnees_2$Age_factor[1])
    }
  }
  
  # Retourner les résultats
  return(list(
    modele_optimal = modele_optimal$modele,
    tous_resultats = resultats,
    modeles_valides = modeles_valides,
    donnees_simulees = donnees_simulees,
    plots = plots,
    resume = list(
      n_ajustements = n_ajustements,
      n_converges = length(modeles_valides),
      meilleur_aic = modele_optimal$aic,
      meilleur_bic = modele_optimal$bic,
      meilleur_loglik = modele_optimal$loglik,
      iteration_optimale = modele_optimal$iteration
    )
  ))
}

# Exemple d'utilisation:
# Assurez-vous d'avoir défini 'dists' dans votre environnement global
dists <- list(log_death_rate = "norm", temp_extreme = "norm") 


resultats_2 <- ajuster_hmm_optimal(
  donnees = df,
  nb_etats = 2,
  n_simulations = 1000,
  seed = 123,
  maxit = 1000,
  tol = 1e-6,
  plot_results = TRUE,
  n_ajustements = 150,
  sauvegarder_plots = TRUE,  # Nouveau paramètre
  repertoire_sortie = OUTER_IMAGE,
)

# Accéder au meilleur modèle
meilleur_modele <- resultats_3$modele_optimal


hmm_model_3 <- resultats_extreme_3$modele_optimal
hmm_model_2 <- resultats_extreme_2$modele_optimal

# Voir le résumé
print(resultats_extreme_3$resume)

# Coefficients du meilleur modèle
meilleur_modele$coeff_fe()

library(ggplot2)
p_real <- ggplot(donnees_65_69, aes(x = time, y = death)) +
  geom_line(color = "blue") +
  labs(title = "Décès (réel) pour 65-69 ans", x = "Temps", y = "Nombre de décès") +
  theme_minimal()
print(p_real)

# =============================================================================
# TARIFICATION RISK-NEUTRAL POUR MODÈLE LOG_DEATH_RATE
# Compatible avec 2 ou 3 états
# =============================================================================


risk_neutral_pricing_log <- function(hmm_results = NULL) {
  
  # Environnement
  env <- new.env()
  env$hmm_model <- NULL
  env$donnees_originales <- NULL
  env$risk_free_rate <- 0.025
  env$n_states <- 2L
  
  # Paramètres d'ajustement
  env$xi_H <- NULL
  env$homo <- FALSE
  
  # Constantes pour les groupes d'âge
  env$age_breaks <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, Inf)
  env$age_labels <- c(paste0(c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85), 
                            "-", 
                            c(4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79, 84, 89)), 
                      "+90")
  
  # Variables pour l'optimisation
  env$avg_log_exposure <- NULL
  env$coef_matrix <- NULL
  env$coef_lookup <- NULL
  
  # === FONCTIONS D'OPTIMISATION ===
  
  precompute_average_exposures <- function() {
    if (is.null(env$donnees_originales)) {
      env$avg_log_exposure <- c(
        "0-4" = 7.5, "5-9" = 7.6, "10-14" = 7.7, "15-19" = 7.8, "20-24" = 7.9,
        "25-29" = 8.0, "30-34" = 8.1, "35-39" = 8.2, "40-44" = 8.5,
        "45-49" = 8.79, "50-54" = 8.62, "55-59" = 8.48, "60-64" = 8.47,
        "65-69" = 8.31, "70-74" = 8.31, "75-79" = 8.07, "80-84" = 7.76,
        "85-89" = 7.20, "+90" = 6.35
      )
    } else {
      env$avg_log_exposure <- setNames(
        tapply(env$donnees_originales$log_exposure, 
               env$donnees_originales$Age_factor, 
               mean, na.rm = TRUE),
        unique(env$donnees_originales$Age_factor)
      )
      
      default_values <- c(
        "0-4" = 7.5, "5-9" = 7.6, "10-14" = 7.7, "15-19" = 7.8, "20-24" = 7.9,
        "25-29" = 8.0, "30-34" = 8.1, "35-39" = 8.2, "40-44" = 8.5,
        "45-49" = 8.79, "50-54" = 8.62, "55-59" = 8.48, "60-64" = 8.47,
        "65-69" = 8.31, "70-74" = 8.31, "75-79" = 8.07, "80-84" = 7.76,
        "85-89" = 7.20, "+90" = 6.35
      )
      
      missing_ages <- setdiff(names(default_values), names(env$avg_log_exposure))
      env$avg_log_exposure <- c(env$avg_log_exposure, default_values[missing_ages])
    }
    
    cat(sprintf("✓ Expositions moyennes précalculées pour %d groupes d'âge\n", 
                length(env$avg_log_exposure)))
  }
  
  precompute_coefficients <- function() {
    coef <- env$hmm_model$obs()$coeff_fe()
    env$coef_lookup <- setNames(coef[, 1], rownames(coef))
    
    n_age_groups <- length(env$age_labels)
    coef_types <- c("(Intercept)", "sin_1", "cos_1", "trend")
    
    # Matrice des coefficients : [état, type_coef, groupe_âge]
    env$coef_matrix <- array(0, dim = c(env$n_states, length(coef_types), n_age_groups),
                            dimnames = list(
                              state = 1:env$n_states,
                              coef_type = coef_types,
                              age_group = env$age_labels
                            ))
    
    # Remplissage de la matrice
    for (state in 1:env$n_states) {
      state_prefix <- paste0("log_death_rate.mean.state", state, ".")
      
      for (age_idx in 1:n_age_groups) {
        age_range <- env$age_labels[age_idx]
        
        # Coefficients de base
        for (coef_idx in 1:length(coef_types)) {
          coef_name <- paste0(state_prefix, coef_types[coef_idx])
          env$coef_matrix[state, coef_idx, age_idx] <- 
            ifelse(coef_name %in% names(env$coef_lookup), 
                   env$coef_lookup[coef_name], 0.0)
        }
        
        # Coefficient d'âge
        age_coef_name <- paste0(state_prefix, "Age_factor", age_range)
        if (age_coef_name %in% names(env$coef_lookup)) {
          env$coef_matrix[state, 1, age_idx] <- 
            env$coef_matrix[state, 1, age_idx] + env$coef_lookup[age_coef_name]
        }
      }
    }
    
    cat(sprintf("✓ Matrice de coefficients précalculée (%dx%dx%d)\n", 
                dim(env$coef_matrix)[1], dim(env$coef_matrix)[2], dim(env$coef_matrix)[3]))
  }
  
  # === FONCTION VECTORISÉE: Calcul de l'intensité de mortalité ===
  # IMPORTANT: log_death_rate = log(death/exposure) 
  # donc death_rate = exp(log_death_rate) = death/exposure
  # L'intensité de mortalité μ(t) = death_rate
  
  compute_mortality_intensity_vectorized <- function(age_ranges, states, weeks, years) {
    n <- length(age_ranges)
    
    if (length(states) == 1) states <- rep(states, n)
    if (length(weeks) == 1) weeks <- rep(weeks, n)
    if (length(years) == 1) years <- rep(years, n)
    
    # Calcul vectorisé des covariables temporelles
    sin_vals <- sin(2 * pi * weeks / 52)
    cos_vals <- cos(2 * pi * weeks / 52)
    
    intensities <- numeric(n)
    
    for (i in 1:n) {
      age_idx <- which(env$age_labels == age_ranges[i])
      if (length(age_idx) == 0) {
        intensities[i] <- 1e-8
        next
      }
      
      state <- states[i]
      
      # Récupération des coefficients
      coefs <- env$coef_matrix[state, , age_idx[1]]
      
      # Calcul du log_death_rate prédit
      log_death_rate_pred <- coefs[1] +              # (Intercept) + Age_factor
                             coefs[2] * sin_vals[i] + # sin_1
                             coefs[3] * cos_vals[i] + # cos_1
                             coefs[4] * years[i]      # trend
      
      # CONVERSION: log_death_rate -> death_rate = intensité de mortalité
      # death_rate = death/exposure = exp(log_death_rate)
      death_rate <- exp(log_death_rate_pred)
      
      intensities[i] <- death_rate
    }
    
    pmax(intensities, 1e-8)
  }
  
  # === FONCTION: Calcul de la probabilité de survie ===
  
  compute_survival_probability_optimized <- function(age, time_horizon, initial_state = 1L, n_steps = NULL) {
    
    # Vérification de l'état initial
    if (initial_state < 1 || initial_state > env$n_states) {
      warning(sprintf("État initial %d invalide, utilisation de l'état 1", initial_state))
      initial_state <- 1L
    }
    
    if (is.null(n_steps)) {
      n_steps <- max(1, round(time_horizon * 52))
    }
    dt <- time_horizon / n_steps
    
    # Vérification de la matrice xi_H
    if (is.null(env$xi_H)) {
      warning("Matrice xi_H non initialisée, utilisation de la matrice identité")
      env$xi_H <- diag(1, env$n_states)
    }
    
    # Simulation des états futurs
    states_path <- simulate_future_states_optimized(n_steps, initial_state)
    
    # Calcul vectorisé des âges à chaque pas
    time_steps <- (0:(n_steps-1)) * dt
    ages_at_time <- age + time_steps
    
    # Calcul vectorisé des semaines et années
    total_weeks_elapsed <- 0:(n_steps-1)
    current_years <- 1 + total_weeks_elapsed %/% 52
    current_weeks <- 1 + (total_weeks_elapsed %% 52)
    
    # Détermination vectorisée des groupes d'âge
    age_group_indices <- findInterval(ages_at_time, env$age_breaks, rightmost.closed = TRUE)
    age_group_indices <- pmax(1, pmin(age_group_indices, length(env$age_labels)))
    age_ranges <- env$age_labels[age_group_indices]
    
    # Calcul vectorisé des intensités de mortalité
    mu_vector <- compute_mortality_intensity_vectorized(
      age_ranges, states_path, current_weeks, current_years
    )
    
    # Calcul du hazard cumulé
    cumulative_hazard <- sum(mu_vector * dt)
    
    # Probabilité de survie
    survival_prob <- exp(-cumulative_hazard)
    max(survival_prob, 1e-8)
  }
  
  # === Simulation des états futurs ===
  
  simulate_future_states_optimized <- function(n_periods, initial_state = 1L) {
    states <- integer(n_periods)
    states[1] <- initial_state
    
    if (env$homo) {
      # Mode homogène avec matrice xi_H
      trans_matrix <- env$xi_H / rowSums(env$xi_H)
      
      for (t in 2:n_periods) {
        trans_probs <- trans_matrix[states[t-1], ]
        
        # Vérification et correction des probabilités
        if (any(is.na(trans_probs)) || any(trans_probs < 0) || sum(trans_probs) == 0) {
          trans_probs <- rep(1/env$n_states, env$n_states)
        } else {
          trans_probs <- trans_probs / sum(trans_probs)
        }
        
        states[t] <- sample(1:env$n_states, 1, prob = trans_probs)
      }
    } else {
      # Mode non-homogène avec prédictions HMM
      for (t in 2:n_periods) {
        tryCatch({
          Q <- as.matrix(env$hmm_model$predict(what = "tpm", t = t)[,,1])
          Q_tilde <- (Q * env$xi_H) / rowSums(Q * env$xi_H)
          trans_probs <- Q_tilde[states[t-1], ]
          
          # Vérification et correction des probabilités
          if (any(is.na(trans_probs)) || any(trans_probs < 0) || sum(trans_probs) == 0) {
            trans_probs <- rep(1/env$n_states, env$n_states)
          } else {
            trans_probs <- trans_probs / sum(trans_probs)
          }
          
          states[t] <- sample(1:env$n_states, 1, prob = trans_probs)
        }, error = function(e) {
          # En cas d'erreur, utiliser une transition uniforme
          states[t] <<- sample(1:env$n_states, 1, prob = rep(1/env$n_states, env$n_states))
        })
      }
    }
    
    states
  }
  
  # === FONCTIONS PRINCIPALES ===
  
  initialize <- function(hmm_results) {
    env$hmm_model <- hmm_results$modele_optimal
    env$donnees_originales <- df  # Utiliser les données globales
    
    # Détection du nombre d'états
    states <- env$hmm_model$viterbi()
    env$n_states <- max(states)
    
    # Initialisation de la matrice d'ajustement (identité par défaut)
    env$xi_H <- diag(1, env$n_states)
    env$homo <- FALSE
    
    # Pré-calculs
    precompute_average_exposures()
    precompute_coefficients()
    
    cat(sprintf("✓ Modèle initialisé avec %d états (optimisé pour log_death_rate)\n", env$n_states))
  }

  set_transition_adjustment <- function(xi_H, homo = FALSE) {
    # Vérification de la dimension
    if (nrow(xi_H) != env$n_states || ncol(xi_H) != env$n_states) {
      stop(sprintf("xi_H doit être une matrice %dx%d", env$n_states, env$n_states))
    }
    
    # Vérification des valeurs
    if (any(is.na(xi_H)) || any(xi_H < 0)) {
      stop("xi_H ne doit pas contenir de NA ou de valeurs négatives")
    }
    
    # Vérification que les lignes somment à une valeur positive
    row_sums <- rowSums(xi_H)
    if (any(row_sums == 0)) {
      stop("xi_H ne doit pas contenir de lignes qui somment à zéro")
    }
    
    env$xi_H <- xi_H
    env$homo <- homo
    cat(sprintf("✓ Matrice d'ajustement %dx%d mise à jour (homo=%s)\n", 
                env$n_states, env$n_states, homo))
  }

  extract_hmm_parameters <- function() {
    coef <- env$hmm_model$obs()$coeff_fe()
    coef_df <- as.data.frame(coef)
    coef_df$rownames <- rownames(coef)
    
    states <- env$hmm_model$viterbi()
    
    list(
      coefficients = coef_df,
      predicted_states = as.integer(states),
      n_states = env$n_states
    )
  }

  compute_risk_neutral_transitions <- function(t) {
    if (env$homo) {
      Q_tilde <- env$xi_H
      Q_tilde <- Q_tilde / rowSums(Q_tilde)
      return(Q_tilde)
    }
    
    Q <- as.matrix(env$hmm_model$predict(what = "tpm", t = 10)[,,1])
    Q_tilde <- Q * env$xi_H
    Q_tilde <- Q_tilde / rowSums(Q_tilde)
    Q_tilde
  }

  compute_mortality_intensity <- function(age_range, state, week, year) {
    compute_mortality_intensity_vectorized(age_range, state, week, year)[1]
  }
  
  compute_survival_probability <- function(age, time_horizon, initial_state = 1L, n_steps = NULL) {
    compute_survival_probability_optimized(age, time_horizon, initial_state, n_steps)
  }

  # Initialisation automatique
  if (!is.null(hmm_results)) {
    initialize(hmm_results)
  }

  # Interface publique
  list(
    initialize = initialize,
    set_transition_adjustment = set_transition_adjustment,
    compute_survival_probability = compute_survival_probability,
    compute_mortality_intensity = compute_mortality_intensity,
    compute_mortality_intensity_vectorized = compute_mortality_intensity_vectorized,
    compute_survival_probability_optimized = compute_survival_probability_optimized,
    
    get_risk_free_rate = function() env$risk_free_rate,
    set_risk_free_rate = function(rate) {
      env$risk_free_rate <- rate
      cat(sprintf("✓ Taux sans risque: %.4f\n", rate))
    },
    
    get_info = function() {
      list(
        n_states = env$n_states,
        risk_free_rate = env$risk_free_rate,
        xi_H = env$xi_H,
        homo = env$homo,
        avg_exposures_computed = !is.null(env$avg_log_exposure),
        coef_matrix_computed = !is.null(env$coef_matrix),
        n_age_groups = length(env$age_labels)
      )
    }
  )
}

# =============================================================================
# FONCTION DE TARIFICATION PRINCIPALE
# =============================================================================

price_death_benefit_analytical <- function(pricing_model,
                                          age_start,
                                          death_benefit,
                                          contract_duration,
                                          method = "classical",
                                          initial_state = 1L,
                                          verbose = TRUE) {
  
  risk_free_rate <- pricing_model$get_risk_free_rate()
  model_info <- pricing_model$get_info()
  n_states <- model_info$n_states
  
  if (verbose) {
    cat("🧮 Calcul analytique de la prime d'assurance décès\n")
    cat(sprintf("   • Modèle: %d états\n", n_states))
    cat(sprintf("   • Méthode: %s\n", method))
    cat(sprintf("   • Âge: %d ans, Capital: %s, Durée: %.1f ans\n", 
                age_start, format(death_benefit, big.mark = ","), contract_duration))
  }
  
  if (method == "classical") {
    result <- calculate_classical_premium_optimized(pricing_model, age_start, death_benefit, 
                                                   contract_duration, risk_free_rate, 
                                                   initial_state, verbose)
  } else if (method == "continuous") {
    result <- calculate_continuous_premium_optimized(pricing_model, age_start, death_benefit, 
                                                    contract_duration, risk_free_rate, 
                                                    initial_state, verbose)
  } else {
    stop("method doit être 'classical' ou 'continuous'")
  }
  
  return(result)
}

# =============================================================================
# MÉTHODE 1: CALCUL CLASSIQUE DISCRET VECTORISÉ
# =============================================================================

calculate_classical_premium_optimized <- function(pricing_model, age_start, death_benefit, 
                                                 contract_duration, risk_free_rate, 
                                                 initial_state, verbose) {
  
  n_years <- as.integer(ceiling(contract_duration))
  
  years_vec <- seq_len(n_years)
  ages_vec <- age_start + (years_vec - 1)
  
  survival_1_year_vec <- vapply(ages_vec, function(age) {
    pricing_model$compute_survival_probability(
      age = age,
      time_horizon = 1.0,
      initial_state = initial_state
    )
  }, numeric(1))
  
  death_probs <- 1.0 - survival_1_year_vec
  survival_probs_to_t <- c(1.0, cumprod(1.0 - death_probs))[1:n_years]
  discount_factors <- (1 + risk_free_rate)^(-years_vec)
  contributions <- death_benefit * discount_factors * survival_probs_to_t * death_probs
  premium_total <- sum(contributions)
  
  if (verbose) {
    cat("\n📊 Calcul année par année (vectorisé):\n")
    cat(sprintf("%4s %12s %12s %12s %12s\n", "Année", "P(survie→t)", "P(décès t→t+1)", "v^(t+1)", "Contribution"))
    cat(paste(rep("-", 65), collapse = ""), "\n")
    
    for (t in seq_len(min(n_years, 10))) {
      cat(sprintf("%4d %12.6f %12.6f %12.6f %12.2f\n", 
                  t, survival_probs_to_t[t], death_probs[t], discount_factors[t], contributions[t]))
    }
    
    if (n_years > 10) {
      cat("   ... (années suivantes masquées)\n")
    }
    
    cat(paste(rep("-", 65), collapse = ""), "\n")
    cat(sprintf("%4s %12s %12s %12s %12.2f\n", "TOTAL", "", "", "", premium_total))
  }
  
  return(list(
    premium_value = premium_total,
    method = "classical_analytical_optimized",
    breakdown = data.frame(
      year = years_vec,
      survival_prob_to_t = survival_probs_to_t,
      death_prob_t_to_t1 = death_probs,
      discount_factor = discount_factors,
      contribution = contributions
    ),
    total_death_probability = sum(survival_probs_to_t * death_probs),
    parameters = list(
      age_start = age_start,
      death_benefit = death_benefit,
      contract_duration = contract_duration,
      risk_free_rate = risk_free_rate,
      initial_state = initial_state
    )
  ))
}

# =============================================================================
# MÉTHODE 2: CALCUL CONTINU VECTORISÉ
# =============================================================================

calculate_continuous_premium_optimized <- function(pricing_model, age_start, death_benefit, 
                                                  contract_duration, risk_free_rate, 
                                                  initial_state, verbose) {
  
  n_steps <- max(52L, as.integer(contract_duration * 52))
  dt <- contract_duration / n_steps
  
  if (verbose) {
    cat(sprintf("\n🔢 Intégration numérique vectorisée avec %d pas (dt = %.4f ans)\n", n_steps, dt))
  }
  
  times <- seq(0, contract_duration - dt, by = dt)
  n_times <- length(times)
  
  ages_t <- age_start + times
  
  age_breaks <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, Inf)
  age_labels <- c(paste0(c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85), 
                        "-", 
                        c(4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79, 84, 89)), 
                  "+90")
  
  age_group_indices <- findInterval(ages_t, age_breaks, rightmost.closed = TRUE)
  age_group_indices <- pmax(1, pmin(age_group_indices, length(age_labels)))
  age_ranges <- age_labels[age_group_indices]
  
  weeks <- ((seq_along(times) - 1L) %% 52L) + 1L
  years <- ((seq_along(times) - 1L) %/% 52L) + 1L
  states <- rep(initial_state, n_times)
  
  survival_probs <- numeric(n_times)
  for (i in seq_len(n_times)) {
    survival_probs[i] <- pricing_model$compute_survival_probability(
      age = age_start,
      time_horizon = times[i],
      initial_state = initial_state
    )
  }
  
  mu_values <- pricing_model$compute_mortality_intensity_vectorized(
    age_ranges = age_ranges,
    states = states,
    weeks = weeks,
    years = years
  )
  
  discount_factors <- exp(-risk_free_rate * times)
  integrand <- death_benefit * discount_factors * survival_probs * mu_values
  
  if (n_times == 1) {
    premium_value <- integrand[1] * dt
  } else {
    weights <- rep(1.0, n_times)
    weights[1] <- 0.5
    weights[n_times] <- 0.5
    premium_value <- dt * sum(integrand * weights)
  }
  
  total_death_prob <- 1 - pricing_model$compute_survival_probability(
    age = age_start, 
    time_horizon = contract_duration, 
    initial_state = initial_state
  )
  
  if (verbose) {
    cat(sprintf("   • Prime calculée: %.2f€\n", premium_value))
    cat(sprintf("   • Probabilité de décès totale: %.2f%%\n", total_death_prob * 100))
    cat(sprintf("   • Intensité moyenne: %.6f\n", mean(mu_values)))
  }
  
  return(list(
    premium_value = premium_value,
    method = "continuous_analytical_optimized",
    integration_steps = n_times,
    total_death_probability = total_death_prob,
    vectorized_computations = list(
      mu_values = mu_values,
      survival_probs = survival_probs,
      times = times,
      integrand = integrand
    ),
    parameters = list(
      age_start = age_start,
      death_benefit = death_benefit,
      contract_duration = contract_duration,
      risk_free_rate = risk_free_rate,
      initial_state = initial_state
    )
  ))
}

# =============================================================================
# ANALYSE DE SENSIBILITÉ GÉNÉRALISÉE (2 OU 3 ÉTATS)
# =============================================================================

generate_xi_sensitivity_table <- function(pricer, 
                                         age_start = 65, 
                                         death_benefit = 100000, 
                                         contract_duration = 40, 
                                         method = "classical",
                                         initial_state = 1,
                                         p22_val = 0.5) {
  
  model_info <- pricer$get_info()
  n_states <- model_info$n_states
  
  cat(sprintf("🎯 Génération du tableau de sensibilité Xi_H (%d états)\n", n_states))
 cat(paste(rep("=", 50), collapse = ""), "\n")
  
  # Valeurs de probabilités à tester
  p_values <- c(0.95, 0.80, 0.65, 0.50, 0.35, 0.20, 0.05)
  n_p <- length(p_values)
  
  if (n_states == 2) {
    return(generate_xi_table_2_states(pricer, p_values, age_start, death_benefit, 
                                     contract_duration, method, initial_state))
  } else if (n_states == 3) {
    return(generate_xi_table_3_states(pricer, p_values, age_start, death_benefit, 
                                     contract_duration, method, initial_state, p22_val))
  } else {
    stop("Analyse de sensibilité disponible uniquement pour 2 ou 3 états")
  }
}

# --- Cas 2 états ---

generate_xi_table_2_states <- function(pricer, p_values, age_start, death_benefit, 
                                      contract_duration, method, initial_state) {
  
  n <- length(p_values)
  premium_rates <- matrix(NA, nrow = n, ncol = n)
  rownames(premium_rates) <- p_values
  colnames(premium_rates) <- p_values
  
  cat("Calcul en cours (2 états)", flush = TRUE)
  
  for (i in 1:n) {
    for (j in 1:n) {
      p11 <- p_values[i]
      p22 <- p_values[j]
      
      xi_matrix <- matrix(c(p11, 1-p11, 1-p22, p22), nrow = 2, byrow = TRUE)
      
      pricer$set_transition_adjustment(xi_matrix, homo = TRUE)
      
      result <- price_death_benefit_analytical(
        pricing_model = pricer,
        age_start = age_start,
        death_benefit = death_benefit,
        contract_duration = contract_duration,
        method = method,
        initial_state = initial_state,
        verbose = FALSE
      )
      
      premium_rates[i, j] <- result$premium_value / death_benefit
      cat(".", flush = TRUE)
    }
  }
  
  cat(" Terminé!\n")
  return(premium_rates)
}

# --- Cas 3 états ---

generate_xi_table_3_states <- function(pricer, p_values, age_start, death_benefit, 
                                      contract_duration, method, initial_state, p22_val=0.5) {
  
  n <- length(p_values)
  
  # Pour 3 états, on fait varier p11 et p33, avec p22 fixe
  p22_fixed <- p22_val  # Valeur intermédiaire pour l'état 2
  
  premium_rates <- matrix(NA, nrow = n, ncol = n)
  rownames(premium_rates) <- p_values
  colnames(premium_rates) <- p_values
  
  cat(sprintf("Calcul en cours (3 états, p22 fixé à %.2f)", p22_fixed), flush = TRUE)
  
  for (i in 1:n) {
    for (j in 1:n) {
      p11 <- p_values[i]
      p33 <- p_values[j]
      
      # Construction de la matrice 3x3
      # On répartit uniformément les transitions hors-diagonale
      p12 <- (1 - p11) / 2
      p13 <- (1 - p11) / 2
      p21 <- (1 - p22_fixed) / 2
      p23 <- (1 - p22_fixed) / 2
      p31 <- (1 - p33) / 2
      p32 <- (1 - p33) / 2
      
      xi_matrix <- matrix(c(
        p11, p12, p13,
        p21, p22_fixed, p23,
        p31, p32, p33
      ), nrow = 3, byrow = TRUE)
      
      pricer$set_transition_adjustment(xi_matrix, homo = TRUE)
      
      result <- price_death_benefit_analytical(
        pricing_model = pricer,
        age_start = age_start,
        death_benefit = death_benefit,
        contract_duration = contract_duration,
        method = method,
        initial_state = initial_state,
        verbose = FALSE
      )
      
      premium_rates[i, j] <- result$premium_value / death_benefit
      cat(".", flush = TRUE)
    }
  }
  
  cat(" Terminé!\n")
  
  # Ajouter l'information sur p22_fixed dans les attributs
  attr(premium_rates, "p22_fixed") <- p22_fixed
  
  return(premium_rates)
}

# === FONCTION D'AFFICHAGE ===

print_formatted_table <- function(xi_table, n_states = 2) {
  
  p_values <- as.numeric(rownames(xi_table))
  
  if (n_states == 2) {
    cat("\n📊 TABLEAU DES TAUX DE PRIME (2 ÉTATS)\n")
    cat("=====================================\n")
    cat(sprintf("%8s", "p1,1\\p2,2"))
  } else if (n_states == 3) {
    p22_fixed <- attr(xi_table, "p22_fixed")
    cat("\n📊 TABLEAU DES TAUX DE PRIME (3 ÉTATS)\n")
    cat("=====================================\n")
    cat(sprintf("p2,2 fixé à %.2f\n", p22_fixed))
    cat(sprintf("%8s", "p1,1\\p3,3"))
  }
  
  for (p in p_values) {
    cat(sprintf("%8.2f", p))
  }
  cat("\n")
  
  for (i in 1:nrow(xi_table)) {
    cat(sprintf("%8.2f", p_values[i]))
    for (j in 1:ncol(xi_table)) {
      cat(sprintf("%8.4f", xi_table[i, j]))
    }
    cat("\n")
  }
  
  cat("\n")
}

# === FONCTION D'EXPORT ===

create_export_table <- function(xi_table, n_states = 2) {
  
  df <- as.data.frame(xi_table)
  df <- round(df, 4)
  
  df$p_row <- as.numeric(rownames(xi_table))
  df <- df[, c("p_row", colnames(df)[1:(ncol(df)-1)])]
  
  if (n_states == 2) {
    colnames(df)[1] <- "p1.1_p2.2"
  } else if (n_states == 3) {
    colnames(df)[1] <- "p1.1_p3.3"
    p22_fixed <- attr(xi_table, "p22_fixed")
    attr(df, "p22_fixed") <- p22_fixed
  }
  
  return(df)
}

# =============================================================================
# TARIFICATION EN LOT
# =============================================================================

price_death_benefit_batch <- function(pricing_model, 
                                     ages_start,
                                     death_benefits, 
                                     contract_durations,
                                     method = "continuous",
                                     initial_state = 1L,
                                     verbose = TRUE) {
  
  n_policies <- length(ages_start)
  
  if (length(death_benefits) == 1) death_benefits <- rep(death_benefits, n_policies)
  if (length(contract_durations) == 1) contract_durations <- rep(contract_durations, n_policies)
  if (length(initial_state) == 1) initial_state <- rep(initial_state, n_policies)
  
  if (verbose) {
    cat(sprintf("📦 Tarification en lot: %d polices (méthode: %s)\n", n_policies, method))
  }
  
  start_time <- Sys.time()
  
  premiums <- vapply(seq_len(n_policies), function(i) {
    result <- price_death_benefit_analytical(
      pricing_model = pricing_model,
      age_start = ages_start[i],
      death_benefit = death_benefits[i],
      contract_duration = contract_durations[i],
      method = method,
      initial_state = initial_state[i],
      verbose = FALSE
    )
    result$premium_value
  }, numeric(1))
  
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (verbose) {
    cat(sprintf("   ✓ Calculé en %.3f secondes\n", execution_time))
    cat(sprintf("   ✓ Vitesse: %.0f polices/seconde\n", n_policies / execution_time))
    cat(sprintf("   ✓ Prime moyenne: %.2f€\n", mean(premiums)))
    cat(sprintf("   ✓ Prime médiane: %.2f€\n", median(premiums)))
  }
  
  data.frame(
    age_start = ages_start,
    death_benefit = death_benefits,
    contract_duration = contract_durations,
    initial_state = initial_state,
    premium = premiums,
    premium_rate = premiums / death_benefits * 1000
  )
}

# =============================================================================
# BENCHMARK DE PERFORMANCE
# =============================================================================

benchmark_pricing_methods <- function(pricing_model, age_start, death_benefit, 
                                     contract_duration, initial_state = 1L) {
  
  cat("🏃‍♂️ Benchmark des méthodes de tarification\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  
  methods <- c("classical", "continuous")
  results <- list()
  
  for (method in methods) {
    cat(sprintf("\n⏱️  Test méthode: %s\n", method))
    
    start_time <- Sys.time()
    
    result <- price_death_benefit_analytical(
      pricing_model = pricing_model,
      age_start = age_start,
      death_benefit = death_benefit,
      contract_duration = contract_duration,
      method = method,
      initial_state = initial_state,
      verbose = FALSE
    )
    
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    results[[method]] <- list(
      premium = result$premium_value,
      time_seconds = execution_time,
      method_details = result$method
    )
    
    cat(sprintf("   ✓ Prime: %.2f€\n", result$premium_value))
    cat(sprintf("   ✓ Temps: %.4f secondes\n", execution_time))
  }
  
  cat("\n📊 RÉSUMÉ COMPARATIF:\n")
  for (method in methods) {
    r <- results[[method]]
    cat(sprintf("   %s: %.2f€ (%.4fs)\n", 
                toupper(method), r$premium, r$time_seconds))
  }
  
  if (length(methods) == 2) {
    speed_ratio <- results[[methods[2]]]$time_seconds / results[[methods[1]]]$time_seconds
    cat(sprintf("   Ratio de vitesse: %.1fx\n", speed_ratio))
    
    price_diff <- abs(results[[methods[1]]]$premium - results[[methods[2]]]$premium)
    price_diff_pct <- price_diff / mean(c(results[[methods[1]]]$premium, results[[methods[2]]]$premium)) * 100
    cat(sprintf("   Différence de prix: %.2f€ (%.2f%%)\n", price_diff, price_diff_pct))
  }
  
  invisible(results)
}

# =============================================================================
# EXEMPLES D'UTILISATION
# =============================================================================

# --- EXEMPLE 1: Initialisation et tarification simple (2 états) ---

exemple_2_etats <- function() {
  cat(paste(rep("=", 3), collapse = ""), "\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  cat("EXEMPLE 1: TARIFICATION AVEC MODÈLE 2 ÉTATS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  # Initialisation du modèle
  pricing_model_2 <- risk_neutral_pricing_log(resultats_extreme_2)
  pricing_model$set_risk_free_rate(0.01)
  
  # Tarification simple
  cat("\n--- Tarification classique ---\n")
  premium_2states <- price_death_benefit_analytical(
    pricing_model_2, 
    age_start = 65, 
    death_benefit = 100000, 
    contract_duration = 40,
    method = "classical",
    initial_state = 1,
    verbose = TRUE
  )
  
  cat(sprintf("\n💰 Prime calculée: %.2f€\n", premium_2states$premium_value))
  cat(sprintf("📊 Taux de prime: %.2f‰\n", 
              premium_2states$premium_value / 100000 * 1000))
  
  return(pricing_model_2)
}

# --- EXEMPLE 2: Analyse de sensibilité (2 états) ---

exemple_sensibilite_2_etats <- function(pricing_model_2) {
  cat(paste(rep("=", 3), collapse = ""), "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("EXEMPLE 2: ANALYSE DE SENSIBILITÉ (2 ÉTATS)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  xi_table_2 <- generate_xi_sensitivity_table(
    pricer = pricing_model_2,
    age_start = 65,
    death_benefit = 100000,
    contract_duration = 40,
    method = "classical",
    initial_state = 1
  )
  
  print_formatted_table(xi_table_2, n_states = 2)
  
  # Export
  df_export <- create_export_table(xi_table_2, n_states = 2)
  cat("\n✓ Tableau exportable créé (première ligne):\n")
  print(head(df_export, 3))
  
  return(xi_table_2)
}

# --- EXEMPLE 3: Initialisation et tarification (3 états) ---

exemple_3_etats <- function() {
  cat('\n\n\n')
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("EXEMPLE 3: TARIFICATION AVEC MODÈLE 3 ÉTATS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Initialisation du modèle
  pricing_model_3 <- risk_neutral_pricing_log(resultats_extreme_3)
  pricing_model_3$set_risk_free_rate(0.01)
  
  # Tarification simple
  cat("\n--- Tarification classique ---\n")
  premium_3states <- price_death_benefit_analytical(
    pricing_model_3, 
    age_start = 65, 
    death_benefit = 100000, 
    contract_duration = 40,
    method = "classical",
    initial_state = 2,
    verbose = TRUE
  )
  
  cat(sprintf("\n💰 Prime calculée: %.2f€\n", premium_3states$premium_value))
  cat(sprintf("📊 Taux de prime: %.2f‰\n", 
              premium_3states$premium_value / 100000 * 1000))
  
  return(pricing_model_3)
}

# --- EXEMPLE 4: Analyse de sensibilité (3 états) ---

exemple_sensibilite_3_etats <- function(pricing_model_3) {
  cat("\n\n\n")
  
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("EXEMPLE 4: ANALYSE DE SENSIBILITÉ (3 ÉTATS)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  xi_table_3 <- generate_xi_sensitivity_table(
    pricer = pricing_model_3,
    age_start = 65,
    death_benefit = 100000,
    contract_duration = 40,
    method = "classical",
    initial_state = 2
  )
  
  print_formatted_table(xi_table_3, n_states = 3)
  
  # Export
  df_export <- create_export_table(xi_table_3, n_states = 3)
  cat(sprintf("\n✓ Tableau exportable créé (p2,2 = %.2f):\n", 
              attr(df_export, "p22_fixed")))
  print(head(df_export, 3))
  
  return(xi_table_3)
}

# --- EXEMPLE 5: Tarification en lot ---

exemple_batch_pricing <- function(pricing_model) {
  cat("\n\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("EXEMPLE 5: TARIFICATION EN LOT\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Génération d'un portefeuille de polices
  n_policies <- 100
  
  portfolio <- price_death_benefit_batch(
    pricing_model = pricing_model,
    ages_start = sample(60:80, n_policies, replace = TRUE),
    death_benefits = 100000,
    contract_durations = sample(10:25, n_policies, replace = TRUE),
    method = "classical",
    initial_state = 1,
    verbose = TRUE
  )
  
  cat("\n📈 Statistiques du portefeuille:\n")
  cat(sprintf("   • Âge moyen: %.1f ans\n", mean(portfolio$age_start)))
  cat(sprintf("   • Durée moyenne: %.1f ans\n", mean(portfolio$contract_duration)))
  cat(sprintf("   • Prime min/max: %.2f€ / %.2f€\n", 
              min(portfolio$premium), max(portfolio$premium)))
  cat(sprintf("   • Taux de prime moyen: %.2f‰\n", mean(portfolio$premium_rate)))
  
  return(portfolio)
}

# --- EXEMPLE 6: Comparaison des méthodes ---

exemple_benchmark <- function(pricing_model) {

  cat("EXEMPLE 6: BENCHMARK DES MÉTHODES\n")
  
  benchmark_results <- benchmark_pricing_methods(
    pricing_model = pricing_model,
    age_start = 65,
    death_benefit = 100000,
    contract_duration = 40,
    initial_state = 1
  )
  
  return(benchmark_results)
}

# =============================================================================
# SCRIPT PRINCIPAL D'EXÉCUTION
# =============================================================================

run_all_examples <- function() {
  
  cat("\n")
  cat("# DÉMONSTRATION COMPLÈTE DU SYSTÈME DE TARIFICATION\n")
  cat("# Modèle log_death_rate avec 2 et 3 états\n")
  
  # Vérification des modèles HMM
  if (!exists("resultats_extreme_2") || !exists("resultats_extreme_3")) {
    cat("\n⚠️  ATTENTION: Les modèles HMM doivent être ajustés avant la tarification\n")
    cat("   Assurez-vous d'avoir exécuté:\n")
    cat("   - resultats_extreme_2 (modèle 2 états)\n")
    cat("   - resultats_extreme_3 (modèle 3 états)\n")
    return(invisible(NULL))
  }
  
  # Exemples pour 2 états
  pricing_model_2 <- exemple_2_etats()
  xi_table_2 <- exemple_sensibilite_2_etats(pricing_model_2)
  
  # Exemples pour 3 états
  pricing_model_3 <- exemple_3_etats()
  xi_table_3 <- exemple_sensibilite_3_etats(pricing_model_3)
  
  # Exemples additionnels
  portfolio <- exemple_batch_pricing(pricing_model_2)
  benchmark <- exemple_benchmark(pricing_model_2)
  
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("# DÉMONSTRATION TERMINÉE\n")
 cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  return(list(
    pricing_model_2 = pricing_model_2,
    pricing_model_3 = pricing_model_3,
    xi_table_2 = xi_table_2,
    xi_table_3 = xi_table_3,
    portfolio = portfolio,
    benchmark = benchmark
  ))
}

# =============================================================================
# UTILISATION RAPIDE
# =============================================================================

cat("\n📚 GUIDE D'UTILISATION RAPIDE\n")
cat("============================\n\n")
cat("1. Initialiser un modèle de tarification:\n")
pricing_model <- risk_neutral_pricing_log(resultats_2)
pricing_model_3 <- risk_neutral_pricing_log(resultats_3)
cat("2. Calculer une prime:\n")
premium <- price_death_benefit_analytical(pricing_model_3, 65, 100000, 40)
cat("3. Analyse de sensibilité:\n")
xi_table_2 <- generate_xi_sensitivity_table(pricing_model)

xi_table_3 <- generate_xi_sensitivity_table(pricing_model_3, p22_val = 0.05)


cat("4. Tarification en lot:\n")
portfolio <- price_death_benefit_batch(pricing_model, 65, 100000, 40)
cat("5. Exécuter tous les exemples:\n")
results <- run_all_examples()

# =============================================================================
# EXTENSION MONTE CARLO ULTRA-OPTIMISÉE
# Version avec pré-calcul maximal et optimisations agressives
# =============================================================================
# =============================================================================
# EXTENSION MONTE CARLO ULTRA-OPTIMISÉE
# Version avec pré-calcul maximal et optimisations agressives
# =============================================================================

add_monte_carlo_methods <- function(pricing_model, hmm_results) {
  
  hmm_model <- hmm_results$modele_optimal
  donnees_originales <- df
  
  model_info <- pricing_model$get_info()
  n_states <- model_info$n_states
  risk_free_rate <- model_info$risk_free_rate
  
  age_breaks <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, Inf)
  age_labels <- c(paste0(c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85), 
                        "-", 
                        c(4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79, 84, 89)), 
                  "+90")
  
  cat("✓ Modèle Monte Carlo ultra-optimisé\n")
  cat(sprintf("  • Nombre d'états: %d\n", n_states))
  
  # === PRÉ-CALCUL DES PARAMÈTRES ===
  
  extract_normal_parameters <- function() {
    obs_obj <- hmm_model$obs()
    coef <- obs_obj$coeff_fe()
    coef_lookup <- setNames(coef[, 1], rownames(coef))
    
    n_age_groups <- length(age_labels)
    coef_types <- c("(Intercept)", "sin_1", "cos_1", "trend")
    
    mu_matrix <- array(0, dim = c(n_states, length(coef_types), n_age_groups))
    
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
          mu_matrix[state, 1, age_idx] <- mu_matrix[state, 1, age_idx] + coef_lookup[age_coef_name]
        }
      }
    }
    
    sigma_matrix <- array(0, dim = c(n_states, n_age_groups))
    for (state in 1:n_states) {
      state_prefix <- paste0("log_death_rate.sd.state", state, ".")
      sigma_name <- paste0(state_prefix, "(Intercept)")
      
      if (sigma_name %in% names(coef_lookup)) {
        sigma_val <- exp(coef_lookup[sigma_name])
      } else {
        sigma_val <- sd(donnees_originales$log_death_rate, na.rm = TRUE) * 0.5
      }
      sigma_matrix[state, ] <- sigma_val
    }
    
    list(mu_matrix = mu_matrix, sigma_matrix = sigma_matrix)
  }
  
  # === SIMPLIFICATION: Utiliser des transitions HOMOGÈNES ===
  # Éviter les appels coûteux à predict() à chaque pas de temps
  
  get_average_transition_matrix <- function(n_periods = 100) {
    cat("   Calcul de la matrice de transition moyenne...\n")
    
    # Échantillonner quelques matrices de transition
    sample_times <- seq(1, n_periods, length.out = min(20, n_periods))
    Q_sum <- matrix(0, n_states, n_states)
    
    for (t in sample_times) {
      tryCatch({
        Q <- as.matrix(hmm_model$predict(what = "tpm", t = round(t))[,,1])
        Q_sum <- Q_sum + Q
      }, error = function(e) {
        Q_sum <<- Q_sum + matrix(1/n_states, n_states, n_states)
      })
    }
    
    Q_avg <- Q_sum / length(sample_times)
    Q_avg / rowSums(Q_avg)  # Normaliser
  }
  
  # === SIMULATION ULTRA-RAPIDE AVEC TRANSITIONS HOMOGÈNES ===
  
  simulate_states_fast <- function(n_periods, initial_state, n_simulations, Q_matrix) {
    states_matrix <- matrix(initial_state, nrow = n_simulations, ncol = n_periods)
    
    # Créer une matrice de probabilités cumulées pour sample rapide
    Q_cumsum <- t(apply(Q_matrix, 1, cumsum))
    
    for (t in 2:n_periods) {
      U <- matrix(runif(n_simulations), ncol = 1)
      
      for (s in 1:n_states) {
        mask <- states_matrix[, t-1] == s
        if (sum(mask) > 0) {
          U_s <- U[mask, 1]
          states_matrix[mask, t] <- rowSums(sweep(Q_cumsum[s, , drop=FALSE], 1, U_s, ">")) + 1
        }
      }
    }
    
    states_matrix
  }
  
  # === CALCUL ULTRA-VECTORISÉ DES MU ===
  
  compute_mu_vectorized <- function(states_matrix, age_start, time_horizon, params) {
    n_simulations <- nrow(states_matrix)
    n_weeks <- ncol(states_matrix)
    
    # Pré-calcul des variables temporelles
    weeks <- ((0:(n_weeks-1)) %% 52) + 1
    years <- ((0:(n_weeks-1)) %/% 52) + 1
    sin_vals <- sin(2 * pi * weeks / 52)
    cos_vals <- cos(2 * pi * weeks / 52)
    
    ages_at_time <- age_start + (0:(n_weeks-1)) / 52
    age_group_indices <- findInterval(ages_at_time, age_breaks, rightmost.closed = TRUE)
    age_group_indices <- pmax(1, pmin(age_group_indices, length(age_labels)))
    
    # Pré-calculer TOUTES les combinaisons possibles [état, temps]
    # Créer une lookup table 3D: [état, temps, (mu, sigma)]
    mu_lookup <- array(0, dim = c(n_states, n_weeks))
    sigma_lookup <- array(0, dim = c(n_states, n_weeks))
    
    for (state in 1:n_states) {
      for (t in 1:n_weeks) {
        age_idx <- age_group_indices[t]
        coefs <- params$mu_matrix[state, , age_idx]
        
        mu_lookup[state, t] <- coefs[1] + 
                                coefs[2] * sin_vals[t] + 
                                coefs[3] * cos_vals[t] + 
                                coefs[4] * years[t]
        
        sigma_lookup[state, t] <- params$sigma_matrix[state, age_idx]
      }
    }
    
    # Extraction vectorisée optimisée
    # Pour chaque (simulation, temps), extraire mu_lookup[état, temps]
    mu_matrix <- matrix(0, nrow = n_simulations, ncol = n_weeks)
    sigma_matrix <- matrix(0, nrow = n_simulations, ncol = n_weeks)
    
    # Vectorisation par état
    for (state in 1:n_states) {
      state_mask <- states_matrix == state
      mu_matrix[state_mask] <- mu_lookup[state, col(states_matrix)[state_mask]]
      sigma_matrix[state_mask] <- sigma_lookup[state, col(states_matrix)[state_mask]]
    }
    
    list(mu_matrix = mu_matrix, sigma_matrix = sigma_matrix)
  }
  
  # === FONCTION PRINCIPALE ULTRA-OPTIMISÉE ===
  
  price_monte_carlo <- function(age_start,
                               death_benefit,
                               contract_duration,
                               initial_state = 1L,
                               n_simulations = 10000,
                               method = "discrete",  # CHANGÉ: discrete par défaut
                               seed = NULL,
                               verbose = TRUE,
                               xi_H = NULL) {  # NOUVEAU: matrice de transition personnalisée
    
    if (!is.null(seed)) set.seed(seed)
    
    if (verbose) {
      cat("\n🚀 Monte Carlo ULTRA-OPTIMISÉ\n")
      cat(paste(rep("=", 60), collapse = ""), "\n")
      cat(sprintf("   • Simulations: %s\n", format(n_simulations, big.mark = ",")))
      cat(sprintf("   • Durée: %.1f ans\n", contract_duration))
    }
    
    total_start <- Sys.time()
    
    # 1. Extraction paramètres (une seule fois)
    params <- extract_normal_parameters()
    if (verbose) cat("   ✓ Paramètres extraits\n")
    
    # 2. Matrice de transition
    n_weeks <- as.integer(ceiling(contract_duration * 52))
    
    if (!is.null(xi_H)) {
      # Utiliser la matrice fournie
      if (nrow(xi_H) != n_states || ncol(xi_H) != n_states) {
        stop(sprintf("xi_H doit être une matrice %dx%d", n_states, n_states))
      }
      Q_avg <- xi_H / rowSums(xi_H)  # Normaliser
      if (verbose) cat("   ✓ Matrice de transition personnalisée utilisée\n")
    } else {
      # Calculer la matrice moyenne du HMM
      Q_avg <- get_average_transition_matrix(n_weeks)
      if (verbose) cat("   ✓ Matrice de transition moyenne calculée\n")
    }
    
    # 3. Simulation des états (RAPIDE)
    start_sim <- Sys.time()
    states_matrix <- simulate_states_fast(n_weeks, initial_state, n_simulations, Q_avg)
    time_states <- as.numeric(difftime(Sys.time(), start_sim, units = "secs"))
    if (verbose) cat(sprintf("   ✓ États simulés en %.2fs\n", time_states))
    
    # 4. Calcul des mu et sigma (ULTRA-RAPIDE avec lookup)
    start_mu <- Sys.time()
    mu_sigma <- compute_mu_vectorized(states_matrix, age_start, contract_duration, params)
    time_mu <- as.numeric(difftime(Sys.time(), start_mu, units = "secs"))
    if (verbose) cat(sprintf("   ✓ μ et σ calculés en %.2fs\n", time_mu))
    
    # 5. Simulation des log(taux) - VECTORISÉ TOTAL
    start_mort <- Sys.time()
    n_total <- n_simulations * n_weeks
    random_normal <- matrix(rnorm(n_total), nrow = n_simulations, ncol = n_weeks)
    log_mortality_matrix <- mu_sigma$mu_matrix + mu_sigma$sigma_matrix * random_normal
    mortality_matrix <- exp(log_mortality_matrix)
    time_mort <- as.numeric(difftime(Sys.time(), start_mort, units = "secs"))
    if (verbose) cat(sprintf("   ✓ Taux simulés en %.2fs\n", time_mort))
    
    # 6. Calcul des payoffs - VERSION CORRIGÉE
    start_pay <- Sys.time()
    dt <- 1 / 52
    
    if (method == "discrete") {
      # Méthode discrète par ANNÉE (comme dans le code original)
      n_years <- as.integer(ceiling(contract_duration))
      
      # Matrices pour stocker les résultats
      payoffs <- numeric(n_simulations)
      
      for (i in 1:n_simulations) {
        survived_to_t <- 1.0  # Probabilité de survie jusqu'à maintenant
        
        for (year in 1:n_years) {
          week_start <- (year - 1) * 52 + 1
          week_end <- min(year * 52, n_weeks)
          
          # Intensité de mortalité moyenne sur l'année
          # mortality_rate = death/exposure = exp(log_death_rate)
          # Donc on utilise directement mortality_matrix
          avg_intensity <- mean(mortality_matrix[i, week_start:week_end])
          
          # Probabilité de survie sur cette année (continue)
          # P(survie) = exp(-∫ μ(s) ds) ≈ exp(-μ_avg * Δt)
          survival_this_year <- exp(-avg_intensity * 1.0)  # 1 an
          
          # Probabilité de décès cette année
          death_prob_this_year <- 1.0 - survival_this_year
          
          # Contribution au payoff
          # E[Payoff] = death_benefit * P(survie jusqu'à t) * P(décès en t) * v^t
          discount_factor <- (1 + risk_free_rate)^(-year)
          payoffs[i] <- payoffs[i] + death_benefit * survived_to_t * death_prob_this_year * discount_factor
          
          # Mise à jour de la probabilité de survie cumulative
          survived_to_t <- survived_to_t * survival_this_year
        }
      }
      
    } else {
      # Méthode continue (approximation plus fine)
      payoffs <- numeric(n_simulations)
      
      for (i in 1:n_simulations) {
        # Calcul du hazard cumulé
        cumulative_hazard <- cumsum(mortality_matrix[i, ] * dt)
        
        # Probabilité de survie à chaque semaine
        survival_probs <- exp(-cumulative_hazard)
        
        # Contribution au payoff à chaque semaine
        # f(t) = μ(t) * S(t) (densité de décès au temps t)
        # Payoff = ∫ death_benefit * e^(-r*t) * f(t) dt
        
        for (t in 1:n_weeks) {
          death_density <- mortality_matrix[i, t] * survival_probs[t]
          discount_factor <- exp(-risk_free_rate * t * dt)
          payoffs[i] <- payoffs[i] + death_benefit * death_density * discount_factor * dt
        }
      }
    }
    
    time_pay <- as.numeric(difftime(Sys.time(), start_pay, units = "secs"))
    if (verbose) cat(sprintf("   ✓ Payoffs calculés en %.2fs\n", time_pay))
    
    total_time <- as.numeric(difftime(Sys.time(), total_start, units = "secs"))
    
    # Statistiques
    premium_value <- mean(payoffs)
    premium_se <- sd(payoffs) / sqrt(n_simulations)
    ci_lower <- premium_value - 1.96 * premium_se
    ci_upper <- premium_value + 1.96 * premium_se
    n_deaths <- sum(payoffs > 0)
    death_rate_empirical <- n_deaths / n_simulations
    
    if (verbose) {
      cat("\n📊 RÉSULTATS\n")
      cat(paste(rep("=", 60), collapse = ""), "\n")
      cat(sprintf("⏱️  Temps total: %.2f sec (%.0f sims/sec)\n", 
                  total_time, n_simulations / total_time))
      cat(sprintf("\n"))
      cat(sprintf("💰 Prime: %.2f € [IC 95%%: %.2f - %.2f]\n", 
                  premium_value, ci_lower, ci_upper))
      cat(sprintf("📈 Décès: %.2f%%\n", death_rate_empirical * 100))
      
      if (n_deaths > 0) {
        cat(sprintf("   Payoff moyen (si décès): %.2f €\n", mean(payoffs[payoffs > 0])))
      }
    }
    
    list(
      premium_value = premium_value,
      premium_se = premium_se,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      death_rate_empirical = death_rate_empirical,
      payoffs = payoffs,
      computation_time = total_time,
      breakdown_time = list(
        states = time_states,
        mu_sigma = time_mu,
        mortality = time_mort,
        payoffs = time_pay
      ),
      method = "monte_carlo_ultra_optimized",
      n_simulations = n_simulations,
      parameters = list(
        age_start = age_start,
        death_benefit = death_benefit,
        contract_duration = contract_duration,
        risk_free_rate = risk_free_rate,
        initial_state = initial_state
      )
    )
  }
  
  # === COMPARAISON ===
  
  compare_deterministic_vs_stochastic <- function(age_start,
                                                 death_benefit,
                                                 contract_duration,
                                                 initial_state = 1L,
                                                 n_simulations = 10000,
                                                 seed = NULL,
                                                 xi_H = NULL) {  # NOUVEAU
    
    cat("\n🔬 COMPARAISON: Déterministe vs Stochastique\n")
    cat(paste(rep("=", 60), collapse = ""), "\n\n")
    
    if (!is.null(seed)) set.seed(seed)
    
    # Déterministe
    cat("1️⃣  Déterministe...\n")
    start_det <- Sys.time()
    
    # Si xi_H est fourni, l'utiliser aussi pour le déterministe
    if (!is.null(xi_H)) {
      pricing_model$set_transition_adjustment(xi_H, homo = TRUE)
    }
    
    premium_det <- price_death_benefit_analytical(
      pricing_model = pricing_model,
      age_start = age_start,
      death_benefit = death_benefit,
      contract_duration = contract_duration,
      method = "continuous",
      initial_state = initial_state,
      verbose = FALSE
    )
    time_det <- as.numeric(difftime(Sys.time(), start_det, units = "secs"))
    cat(sprintf("   ✓ Prime: %.2f € (%.3fs)\n", premium_det$premium_value, time_det))
    
    # Stochastique
    cat("\n2️⃣  Stochastique...\n")
    result_mc <- price_monte_carlo(
      age_start = age_start,
      death_benefit = death_benefit,
      contract_duration = contract_duration,
      initial_state = initial_state,
      n_simulations = n_simulations,
      method = "discrete",  # Utiliser discrete
      verbose = FALSE,
      xi_H = xi_H  # Passer xi_H
    )
    cat(sprintf("   ✓ Prime: %.2f € [%.2f - %.2f] (%.2fs)\n", 
                result_mc$premium_value, result_mc$ci_lower, 
                result_mc$ci_upper, result_mc$computation_time))
    
    # Analyse
    cat("\n📊 ANALYSE\n")
    cat(paste(rep("-", 60), collapse = ""), "\n")
    
    diff_abs <- result_mc$premium_value - premium_det$premium_value
    diff_rel <- diff_abs / premium_det$premium_value * 100
    
    cat(sprintf("Différence: %.2f € (%.2f%%)\n", diff_abs, diff_rel))
    
    z_score <- diff_abs / result_mc$premium_se
    p_value <- 2 * (1 - pnorm(abs(z_score)))
    
    cat(sprintf("Significativité: z=%.2f, p=%.4f ", z_score, p_value))
    cat(ifelse(p_value < 0.05, "⚠️ Significatif\n", "✓ Non significatif\n"))
    
    list(
      deterministic = premium_det,
      stochastic = result_mc,
      comparison = list(
        premium_diff_abs = diff_abs,
        premium_diff_rel = diff_rel,
        z_score = z_score,
        p_value = p_value
      )
    )
  }
  
  # === VISUALISATION SIMPLIFIÉE ===
  
  plot_sample_trajectories <- function(age_start,
                                      contract_duration,
                                      initial_state = 1L,
                                      n_trajectories = 100,
                                      seed = NULL) {
    
    if (!is.null(seed)) set.seed(seed)
    
    cat(sprintf("\n📈 Simulation de %d trajectoires...\n", n_trajectories))
    
    result <- price_monte_carlo(
      age_start = age_start,
      death_benefit = 100000,
      contract_duration = contract_duration,
      initial_state = initial_state,
      n_simulations = n_trajectories,
      verbose = FALSE
    )
    
    cat("✓ Visualisation générée\n")
    
    invisible(result)
  }
  
  # Ajout des méthodes
  pricing_model$price_monte_carlo <- price_monte_carlo
  pricing_model$compare_methods <- compare_deterministic_vs_stochastic
  pricing_model$plot_trajectories <- plot_sample_trajectories
  
  cat("\n✅ Monte Carlo ultra-optimisé activé!\n")
  cat("   Performance cible: 2000-10000 sims/sec\n\n")
  
  return(pricing_model)
}

# =============================================================================
# UTILISATION
# =============================================================================

cat("\n📚 UTILISATION OPTIMALE\n")
cat("=======================\n\n")
cat("pricing_model <- add_monte_carlo_methods(pricing_model, resultats_extreme_2)\n")
cat("result <- pricing_model$price_monte_carlo(65, 100000, 40, n_simulations=10000)\n\n")
cat("Avec matrice personnalisée:\n")
cat("  xi_H <- matrix(c(0.95, 0.05, 0.05, 0.95), nrow=2)\n")
cat("  result <- pricing_model$price_monte_carlo(..., xi_H=xi_H)\n\n")
cat("Optimisations:\n")
cat("  • Transitions homogènes (matrice moyenne ou personnalisée)\n")
cat("  • Lookup table pour μ et σ\n")
cat("  • Vectorisation maximale\n\n")
}

# =============================================================================
# FONCTIONS POUR ANALYSE DE SENSIBILITÉ STOCHASTIQUE
# =============================================================================

# === FONCTION: Analyse pour 2 états ===

analyze_stochastic_2_states <- function(pricing_model,
                                       age_start = 65,
                                       death_benefit = 100000,
                                       contract_duration = 40,
                                       initial_state = 1,
                                       n_simulations = 10000,
                                       scenarios = list(
                                         "État 1 stable" = matrix(c(0.95, 0.05, 0.05, 0.95), 2, 2),
                                         "État 2 stable" = matrix(c(0.05, 0.95, 0.95, 0.05), 2, 2)
                                       ),
                                       seed = 123) {
  
  cat("\n📊 ANALYSE STOCHASTIQUE - 2 ÉTATS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  results_list <- list()
  
  for (scenario_name in names(scenarios)) {
    cat(sprintf("\n🎯 Scénario: %s\n", scenario_name))
    cat(paste(rep("-", 70), collapse = ""), "\n")
    
    xi_H <- scenarios[[scenario_name]]
    
    cat("Matrice de transition:\n")
    print(xi_H)
    cat("\n")
    
    # Tarification déterministe
    pricing_model$set_transition_adjustment(xi_H, homo = TRUE)
    det_result <- price_death_benefit_analytical(
      pricing_model = pricing_model,
      age_start = age_start,
      death_benefit = death_benefit,
      contract_duration = contract_duration,
      method = "continuous",
      initial_state = initial_state,
      verbose = FALSE
    )
    
    # Tarification stochastique
    set.seed(seed)
    mc_result <- pricing_model$price_monte_carlo(
      age_start = age_start,
      death_benefit = death_benefit,
      contract_duration = contract_duration,
      initial_state = initial_state,
      n_simulations = n_simulations,
      method = "discrete",
      xi_H = xi_H,
      verbose = FALSE
    )
    
    # Résumé
    cat(sprintf("Déterministe: %.2f €\n", det_result$premium_value))
    cat(sprintf("Stochastique: %.2f € [%.2f - %.2f]\n", 
                mc_result$premium_value, mc_result$ci_lower, mc_result$ci_upper))
    cat(sprintf("Différence:   %.2f € (%.2f%%)\n\n", 
                mc_result$premium_value - det_result$premium_value,
                (mc_result$premium_value - det_result$premium_value) / det_result$premium_value * 100))
    
    results_list[[scenario_name]] <- list(
      xi_H = xi_H,
      deterministic = det_result,
      stochastic = mc_result,
      difference = mc_result$premium_value - det_result$premium_value,
      difference_pct = (mc_result$premium_value - det_result$premium_value) / det_result$premium_value * 100
    )
  }
  
  # Graphique comparatif
  df_plot <- data.frame(
    scenario = rep(names(scenarios), each = 2),
    method = rep(c("Déterministe", "Stochastique"), length(scenarios)),
    premium = c(rbind(
      sapply(results_list, function(x) x$deterministic$premium_value),
      sapply(results_list, function(x) x$stochastic$premium_value)
    ))
  )
  
  p <- ggplot(df_plot, aes(x = scenario, y = premium, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Comparaison Déterministe vs Stochastique (2 états)",
         subtitle = sprintf("Âge=%d, Capital=%s, Durée=%d ans", 
                          age_start, format(death_benefit, big.mark=","), contract_duration),
         x = "Scénario",
         y = "Prime (€)",
         fill = "Méthode") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  
  invisible(results_list)
}

# === FONCTION: Analyse pour 3 états ===

analyze_stochastic_3_states <- function(pricing_model,
                                       age_start = 65,
                                       death_benefit = 100000,
                                       contract_duration = 40,
                                       initial_state = 2,
                                       n_simulations = 10000,
                                       scenarios = list(
                                         "État 1 stable" = matrix(c(0.95, 0.025, 0.025,
                                                                   0.05, 0.05, 0.90,
                                                                   0.05, 0.90, 0.05), 3, 3, byrow=TRUE),
                                         "État 2 stable" = matrix(c(0.05, 0.05, 0.90,
                                                                   0.025, 0.95, 0.025,
                                                                   0.90, 0.05, 0.05), 3, 3, byrow=TRUE),
                                         "État 3 stable" = matrix(c(0.05, 0.90, 0.05,
                                                                   0.90, 0.05, 0.05,
                                                                   0.025, 0.025, 0.95), 3, 3, byrow=TRUE)
                                       ),
                                       seed = 123) {
  
  cat("\n📊 ANALYSE STOCHASTIQUE - 3 ÉTATS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  results_list <- list()
  
  for (scenario_name in names(scenarios)) {
    cat(sprintf("\n🎯 Scénario: %s\n", scenario_name))
    cat(paste(rep("-", 70), collapse = ""), "\n")
    
    xi_H <- scenarios[[scenario_name]]
    
    cat("Matrice de transition:\n")
    print(round(xi_H, 3))
    cat("\n")
    
    # Tarification déterministe
    pricing_model$set_transition_adjustment(xi_H, homo = TRUE)
    det_result <- price_death_benefit_analytical(
      pricing_model = pricing_model,
      age_start = age_start,
      death_benefit = death_benefit,
      contract_duration = contract_duration,
      method = "continuous",
      initial_state = initial_state,
      verbose = FALSE
    )
    
    # Tarification stochastique
    set.seed(seed)
    mc_result <- pricing_model$price_monte_carlo(
      age_start = age_start,
      death_benefit = death_benefit,
      contract_duration = contract_duration,
      initial_state = initial_state,
      n_simulations = n_simulations,
      xi_H = xi_H,
      verbose = FALSE
    )
    
    # Résumé
    cat(sprintf("Déterministe: %.2f €\n", det_result$premium_value))
    cat(sprintf("Stochastique: %.2f € [%.2f - %.2f]\n", 
                mc_result$premium_value, mc_result$ci_lower, mc_result$ci_upper))
    cat(sprintf("Différence:   %.2f € (%.2f%%)\n\n", 
                mc_result$premium_value - det_result$premium_value,
                (mc_result$premium_value - det_result$premium_value) / det_result$premium_value * 100))
    
    results_list[[scenario_name]] <- list(
      xi_H = xi_H,
      deterministic = det_result,
      stochastic = mc_result,
      difference = mc_result$premium_value - det_result$premium_value,
      difference_pct = (mc_result$premium_value - det_result$premium_value) / det_result$premium_value * 100
    )
  }
  
  # Graphique comparatif
  df_plot <- data.frame(
    scenario = rep(names(scenarios), each = 2),
    method = rep(c("Déterministe", "Stochastique"), length(scenarios)),
    premium = c(rbind(
      sapply(results_list, function(x) x$deterministic$premium_value),
      sapply(results_list, function(x) x$stochastic$premium_value)
    ))
  )
  
  p <- ggplot(df_plot, aes(x = scenario, y = premium, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Comparaison Déterministe vs Stochastique (3 états)",
         subtitle = sprintf("Âge=%d, Capital=%s, Durée=%d ans", 
                          age_start, format(death_benefit, big.mark=","), contract_duration),
         x = "Scénario",
         y = "Prime (€)",
         fill = "Méthode") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  
  invisible(results_list)
}

# === FONCTION: Graphique des différences ===

plot_differences <- function(results_list, title = "Impact de la stochasticité") {
  
  df_diff <- data.frame(
    scenario = names(results_list),
    difference = sapply(results_list, function(x) x$difference),
    difference_pct = sapply(results_list, function(x) x$difference_pct)
  )
  
  p1 <- ggplot(df_diff, aes(x = scenario, y = difference)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste(title, "- Différence absolue"),
         x = "Scénario",
         y = "Différence (€)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(df_diff, aes(x = scenario, y = difference_pct)) +
    geom_bar(stat = "identity", fill = "orange") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste(title, "- Différence relative"),
         x = "Scénario",
         y = "Différence (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  gridExtra::grid.arrange(p1, p2, ncol = 2)
  
  invisible(df_diff)
}

# =============================================================================
# EXEMPLES D'UTILISATION
# =============================================================================

cat("\n📖 EXEMPLES - ANALYSE DE SENSIBILITÉ STOCHASTIQUE\n")
cat("==================================================\n\n")
cat("# Pour 2 états:\n")
cat("results_2 <- analyze_stochastic_2_states(pricing_model)\n")
cat("plot_differences(results_2, '2 états')\n\n")
cat("# Pour 3 états:\n")
cat("results_3 <- analyze_stochastic_3_states(pricing_model)\n")
cat("plot_differences(results_3, '3 états')\n\n")
cat("# Personnalisé:\n")
cat("xi_custom <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)\n")
cat("result <- pricing_model$price_monte_carlo(65, 100000, 40, xi_H=xi_custom)\n\n")

pricing_model <- add_monte_carlo_methods(pricing_model, resultats_2)

# Matrice personnalisée
xi_H <- matrix(c(0.05, 0.95, 0.05, 0.05), 2, 2, byrow = TRUE)

# Test
result <- pricing_model$price_monte_carlo(
  age_start = 65,
  death_benefit = 100000,
  contract_duration = 40,
  initial_state = 1,
  n_simulations = 1000,
  xi_H = xi_H,
  seed = 123
)

# La prime devrait maintenant être réaliste (quelques milliers d'euros, pas 100k!)
print(result$premium_value)




library(ggplot2)
library(dplyr)
library(scales)

# Créer les données pour chaque modèle
# On simule les distributions basées sur les quantiles fournis

# Fonction pour créer des données à partir des statistiques
create_boxplot_data <- function(min_val, q1, median, q3, max_val, n = 49) {
  # Génération approximative basée sur les quantiles
  data <- c(
    min_val,
    seq(q1, median, length.out = n/4),
    seq(median, q3, length.out = n/4),
    max_val
  )
  return(data)
}

# Données pour chaque modèle
univarie_data <- create_boxplot_data(0.8548, 0.8708, 0.8753, 0.8775, 0.8796)
bivarie2_data <- create_boxplot_data(0.8077, 0.8092, 0.8101, 0.8101, 0.8111)
bivarie3_data <- create_boxplot_data(0.8051, 0.8104, 0.8132, 0.8132, 0.8182)

# Créer le dataframe
df <- data.frame(
  Modele = factor(rep(c("Univarié\n2 états", "Bivarié\n2 états", "Bivarié\n3 états"), 
                      each = length(univarie_data)),
                  levels = c("Univarié\n2 états", "Bivarié\n2 états", "Bivarié\n3 états")),
  Taux = c(univarie_data, bivarie2_data, bivarie3_data)
)

# Statistiques pour les annotations
stats_df <- data.frame(
  Modele = factor(c("Univarié\n2 états", "Bivarié\n2 états", "Bivarié\n3 états"),
                  levels = c("Univarié\n2 états", "Bivarié\n2 états", "Bivarié\n3 états")),
  Mediane = c(0.8753, 0.8101, 0.8132),
  IQR = c(0.0041, 0.0009, 0.0028),
  Amplitude = c("2,90%", "0,42%", "1,63%"),
  x = c(1, 2, 3),
  y = c(0.805, 0.805, 0.805)
)

# Créer le graphique
p <- ggplot(df, aes(x = Modele, y = Taux, fill = Modele)) +
  
  # Boxplot
  geom_boxplot(
    width = 0.5,
    outlier.shape = 21,
    outlier.size = 3,
    outlier.fill = "red",
    outlier.color = "black",
    outlier.stroke = 1.5,
    color = "black",
    linewidth = 0.8,
    alpha = 0.7
  ) +
  
  # Couleurs personnalisées
  scale_fill_manual(values = c("#bbdefb", "#90caf9", "#64b5f6")) +
  
  # Échelle Y
  scale_y_continuous(
    breaks = seq(0.81, 0.88, by = 0.01),
    limits = c(0.80, 0.89),
    labels = number_format(accuracy = 0.01, decimal.mark = ",")
  ) +
  
  # Labels et titre
  labs(
    title = "Distribution des taux de prime (pour 1$ de prestation)",
    subtitle = "selon 49 scénarios de transition",
    y = "Taux de prime",
    x = NULL
  ) +
  
  # Annotations des statistiques
  geom_text(
    data = stats_df,
    aes(x = x, y = y, 
        label = paste0("Médiane: ", sprintf("%.4f", Mediane), "\n",
                      "IQR: ", sprintf("%.4f", IQR), "\n",
                      "Amplitude: ", Amplitude)),
    inherit.aes = FALSE,
    size = 3.5,
    vjust = 1,
    color = "black",
    fontface = "plain",
    lineheight = 0.9
  ) +
  
  # Annotations Min/Max
  annotate("text", x = 3.35, y = 0.88, label = "Max", 
           size = 3.5, color = "gray40", hjust = 0) +
  annotate("text", x = 3.35, y = 0.81, label = "Min", 
           size = 3.5, color = "gray40", hjust = 0) +
  
  # Thème
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, 
                             color = "#1a237e", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, 
                                 margin = margin(b = 15)),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
    axis.text.x = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.5),
    plot.margin = margin(20, 30, 20, 20)
  ) +
  
  # Légende personnalisée
  guides(fill = guide_legend(
    nrow = 1,
    label.position = "bottom",
    keywidth = 1.5,
    keyheight = 1
  ))

# Afficher le graphique
print(p)

# Sauvegarder en PDF haute qualité
ggsave(
  filename = "boxplot_stabilite_tarifaire.pdf",
  plot = p,
  width = 12,
  height = 10,
  units = "in",
  dpi = 300,
  device = cairo_pdf  # Pour meilleure qualité des polices
)

# Sauvegarder aussi en PNG haute résolution
ggsave(
  filename = "boxplot_stabilite_tarifaire.png",
  plot = p,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)

cat("✅ Graphiques générés avec succès!\n")
cat("   - boxplot_stabilite_tarifaire.pdf\n")
cat("   - boxplot_stabilite_tarifaire.png\n")