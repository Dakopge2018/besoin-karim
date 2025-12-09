# =============================================================================
# FICHIER DE CONFIGURATION GLOBALE (const.r)
# Ancre centralisée pour tous les paramètres et variables globales
# =============================================================================

# === LIBRAIRIES REQUISES ===
required_packages <- c(
  "readxl",
  "ggplot2",
  "scico",
  "hmmTMB",
  "gridExtra",
  "RColorBrewer",
  "grid",
  "dplyr",
  "pracma",
  "expm"
)

# Fonction d'initialisation
initialize_environment <- function() {
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

# === PARAMÈTRES DE DONNÉES ===
CONST$DATA <- list(
  path_to_data = "./data/main_database.xlsx",
  sheet_name = "database_4",
  seed = 123
)

# === PARAMÈTRES DE MODÈLE HMM ===
CONST$HMM <- list(
  # Nombre d'états à tester (2 et 3)
  n_states_to_fit = c(2, 3),
  
  # Paramètres de covariables
  degree_obs_pol = 1,          # Degré du polynôme observé
  degree_trans_pol = 1,         # Degré du polynôme de transition
  period = 52,                  # Période (semaines)
  
  # Paramètres d'ajustement
  maxit = 1000,                 # Nombre max d'itérations
  tol = 1e-6,                   # Tolérance de convergence
  n_ajustements = 150,          # Nombre de tentatives d'ajustement
  
  # Distributions observées
  dists = list(
    log_death_rate = "norm",
    temp_extreme = "norm"
  )
)

# === PARAMÈTRES DE TARIFICATION ===
CONST$PRICING <- list(
  # Taux sans risque
  risk_free_rate = 0.025,
  
  # Groupes d'âge (breaks et labels)
  age_breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, Inf),
  age_labels = c(
    "0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44",
    "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "+90"
  ),
  
  # Expositions moyennes par groupe d'âge (log-échelle)
  avg_log_exposure = c(
    "0-4" = 7.5, "5-9" = 7.6, "10-14" = 7.7, "15-19" = 7.8, "20-24" = 7.9,
    "25-29" = 8.0, "30-34" = 8.1, "35-39" = 8.2, "40-44" = 8.5,
    "45-49" = 8.79, "50-54" = 8.62, "55-59" = 8.48, "60-64" = 8.47,
    "65-69" = 8.31, "70-74" = 8.31, "75-79" = 8.07, "80-84" = 7.76,
    "85-89" = 7.20, "+90" = 6.35
  )
)

# === PARAMÈTRES DE MONTE CARLO ===
CONST$MONTE_CARLO <- list(
  n_simulations_default = 10000,
  n_steps_default = 52 * 40,  # 40 ans * 52 semaines
  seed_default = 123
)

# === CHEMINS DE SORTIE ===
CONST$PATHS <- list(
  # Répertoires de sortie pour graphiques et résultats
  output_images = "./results",
  output_results = "./results",
  output_models = "./models",
  output_pricing = "./results"
)

# === FORMULES DE MODÈLE ===
CONST$FORMULAS <- list(
  # Formule de transition
  transition = ~ cos_1 + sin_1,
  
  # Formules d'observation (log_death_rate et temp_extreme)
  observation = list(
    log_death_rate = list(
      mean = ~ sin_1 + cos_1 + trend + Age_factor,
      sd = ~ 1
    ),
    temp_extreme = list(
      mean = ~ sin_1 + cos_1 + trend,
      sd = ~ 1
    )
  )
)

# === OPTIONS DE VISUALISATION ===
CONST$VIZ <- list(
  # Palette de couleurs pour les états
  state_palette_name = "Set1",
  
  # Thème par défaut
  default_theme = "theme_bw",
  
  # Résolution des graphiques
  dpi_default = 300,
  width_pdf = 10,
  height_pdf = 6
)

# === CONSTANTES DE MODÈLE COMPARATIF ===
CONST$COMPARISON <- list(
  # Modèles à comparer
  models_to_fit = c("poisson_normal_2", "normal_normal_2", "normal_normal_3"),
  
  # Critères de sélection
  selection_criteria = c("AIC", "BIC", "LogLik"),
  
  # Paramètres d'ajustement des modèles
  n_init_strategies = 4,  # 4 stratégies pour paramètres initiaux
  
  # Groupes d'âge pour analyse détaillée
  age_groups_analysis = c("65-69", "+90")
)

# === PALETTES DE COULEURS ===
CONST$COLORS <- list(
  state_colors = RColorBrewer::brewer.pal(9, "Set1"),
  comparison_colors = RColorBrewer::brewer.pal(3, "Dark2"),
  model_colors = c(
    "poisson_normal_2" = "#1b9e77",
    "normal_normal_2" = "#d95f02",
    "normal_normal_3" = "#7570b3"
  )
)

# === FONCTION D'ACCÈS AUX CONSTANTES ===
# Usage: get_const("PRICING", "risk_free_rate")
get_const <- function(...) {
  args <- list(...)
  result <- CONST
  for (arg in args) {
    result <- result[[arg]]
    if (is.null(result)) {
      warning(sprintf("Constante non trouvée: %s", paste(args, collapse = " > ")))
      return(NULL)
    }
  }
  return(result)
}

# === FONCTION DE VALIDATION ===
validate_environment <- function() {
  cat("🔍 Validation de l'environnement...\n")
  
  checks <- list(
    "Chemin données valide" = file.exists(CONST$DATA$path_to_data),
    "Répertoires de sortie" = all(file.exists(c(
      CONST$PATHS$output_images,
      CONST$PATHS$output_results
    )))
  )
  
  for (check_name in names(checks)) {
    status <- ifelse(checks[[check_name]], "✓", "✗")
    cat(sprintf("%s %s\n", status, check_name))
  }
  
  return(all(unlist(checks)))
}

# === AFFICHAGE DES INFORMATIONS ===
print_constants <- function() {
  cat("\n📋 CONSTANTES GLOBALES CHARGÉES\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  
  cat("📊 HMM:\n")
  cat(sprintf("   • États à tester: %s\n", paste(CONST$HMM$n_states_to_fit, collapse = ", ")))
  cat(sprintf("   • Ajustements: %d tentatives\n", CONST$HMM$n_ajustements))
  
  cat("\n💰 TARIFICATION:\n")
  cat(sprintf("   • Taux sans risque: %.2f%%\n", CONST$PRICING$risk_free_rate * 100))
  cat(sprintf("   • Groupes d'âge: %d\n", length(CONST$PRICING$age_labels)))
  
  cat("\n🎲 MONTE CARLO:\n")
  cat(sprintf("   • Simulations par défaut: %s\n", format(CONST$MONTE_CARLO$n_simulations_default, big.mark = ",")))
  
  cat("\n📁 CHEMINS:\n")
  cat(sprintf("   • Images: %s\n", CONST$PATHS$output_images))
  cat(sprintf("   • Résultats: %s\n", CONST$PATHS$output_results))
}

# Initialisation automatique
CONST <- list()  # Conteneur global

# Remplissage du conteneur
CONST$DATA <- list(
  path_to_data = "C:\\Users\\samue\\OneDrive\\Documents\\cours\\Projet de mémoire\\git\\real_data_code\\data\\main_database.xlsx",
  sheet_name = "database_4",
  seed = 123
)

CONST$HMM <- list(
  n_states_to_fit = c(2, 3),
  degree_obs_pol = 1,
  degree_trans_pol = 1,
  period = 52,
  maxit = 1000,
  tol = 1e-6,
  n_ajustements = 150,
  dists = list(log_death_rate = "norm", temp_extreme = "norm")
)

CONST$PRICING <- list(
  risk_free_rate = 0.025,
  age_breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, Inf),
  age_labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44",
                 "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "+90"),
  avg_log_exposure = c("0-4" = 7.5, "5-9" = 7.6, "10-14" = 7.7, "15-19" = 7.8, "20-24" = 7.9,
                       "25-29" = 8.0, "30-34" = 8.1, "35-39" = 8.2, "40-44" = 8.5,
                       "45-49" = 8.79, "50-54" = 8.62, "55-59" = 8.48, "60-64" = 8.47,
                       "65-69" = 8.31, "70-74" = 8.31, "75-79" = 8.07, "80-84" = 7.76,
                       "85-89" = 7.20, "+90" = 6.35)
)

CONST$MONTE_CARLO <- list(
  n_simulations_default = 10000,
  n_steps_default = 52 * 40,
  seed_default = 123
)

CONST$PATHS <- list(
  output_images = "C:\\Users\\samue\\OneDrive\\Documents\\cours\\Projet de mémoire\\git\\real_data_code\\image\\nnnn\\hmm_tmb\\normal",
  output_results = "./results",
  output_models = "./models",
  output_pricing = "./pricing_results"
)

CONST$FORMULAS <- list(
  transition = ~ cos_1 + sin_1,
  observation = list(
    log_death_rate = list(mean = ~ sin_1 + cos_1 + trend + Age_factor, sd = ~ 1),
    temp_extreme = list(mean = ~ sin_1 + cos_1 + trend, sd = ~ 1)
  )
)

CONST$VIZ <- list(
  state_palette_name = "Set1",
  default_theme = "theme_bw",
  dpi_default = 300,
  width_pdf = 10,
  height_pdf = 6
)

CONST$COMPARISON <- list(
  models_to_fit = c("poisson_normal_2", "normal_normal_2", "normal_normal_3"),
  selection_criteria = c("AIC", "BIC", "LogLik"),
  age_groups_analysis = c("65-69", "+90")
)

CONST$COLORS <- list(
  state_colors = RColorBrewer::brewer.pal(9, "Set1"),
  comparison_colors = RColorBrewer::brewer.pal(3, "Dark2"),
  model_colors = c("poisson_normal_2" = "#1b9e77", "normal_normal_2" = "#d95f02", "normal_normal_3" = "#7570b3")
)

cat("✅ Constantes globales chargées avec succès!\n")
