# =============================================================================
# GLOBAL CONFIGURATION FILE (const.r)
# Centralized anchor for all parameters and global variables
# =============================================================================

# === REQUIRED LIBRARIES ===
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

# Initialization function
initialize_environment <- function() {
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

# === AUTOMATIC INITIALIZATION ===
CONST <- list()  # Global container

# === DATA PARAMETERS ===
CONST$DATA <- list(
  path_to_data = "./data/main_database.xlsx",
  sheet_name = "database_4",
  seed = 123,
  code_filter = "ITI43"  # Geographic code to analyze
)

# === HMM MODEL PARAMETERS ===
CONST$HMM <- list(
  # Number of states to test (2 and 3)
  n_states_to_fit = c(2, 3),
  
  # Covariate parameters
  degree_obs_pol = 1,          # Observation polynomial degree
  degree_trans_pol = 1,         # Transition polynomial degree
  period = 52,                  # Period (weeks)
  
  # Fitting parameters
  maxit = 1000,                 # Maximum iterations
  tol = 1e-6,                   # Convergence tolerance
  n_ajustements = 150,          # Number of fitting attempts
  
  # Observed distributions
  dists = list(
    log_death_rate = "norm",
    temp_extreme = "norm"
  )
)

# === PRICING PARAMETERS ===
CONST$PRICING <- list(
  # Risk-free rate
  risk_free_rate = 0.0225,
  
  # Age groups (breaks and labels)
  age_breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, Inf),
  age_labels = c(
    "0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44",
    "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "+90"
  ),
  
  # Average exposure by age group (log-scale)
  avg_log_exposure = c(
    "0-4" = 7.5, "5-9" = 7.6, "10-14" = 7.7, "15-19" = 7.8, "20-24" = 7.9,
    "25-29" = 8.0, "30-34" = 8.1, "35-39" = 8.2, "40-44" = 8.5,
    "45-49" = 8.79, "50-54" = 8.62, "55-59" = 8.48, "60-64" = 8.47,
    "65-69" = 8.31, "70-74" = 8.31, "75-79" = 8.07, "80-84" = 7.76,
    "85-89" = 7.20, "+90" = 6.35
  )
)

# === MONTE CARLO PARAMETERS ===
CONST$MONTE_CARLO <- list(
  n_simulations_default = 10000,
  n_steps_default = 52 * 40,  # 40 years * 52 weeks
  seed_default = 123
)

# === OUTPUT PATHS ===
CONST$PATHS <- list(
  # Output directories for plots and results
  output_images = "./results",
  output_results = "./results",
  output_models = "./models",
  output_pricing = "./results"
)

# === MODEL FORMULAS ===
CONST$FORMULAS <- list(
  # Transition formula
  transition = ~ cos_1 + sin_1,
  
  # Observation formulas (log_death_rate and temp_extreme)
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

# === VISUALIZATION OPTIONS ===
CONST$VIZ <- list(
  # Color palette for states
  state_palette_name = "Set1",
  
  # Default theme
  default_theme = "theme_bw",
  
  # Plot resolution
  dpi_default = 300,
  width_pdf = 10,
  height_pdf = 6
)

# === MODEL COMPARISON CONSTANTS ===
CONST$COMPARISON <- list(
  # Models to compare
  models_to_fit = c("poisson_normal_2", "normal_normal_2", "normal_normal_3"),
  
  # Selection criteria
  selection_criteria = c("AIC", "BIC", "LogLik"),
  
  # Model fitting parameters
  n_init_strategies = 4,  # 4 strategies for initial parameters
  
  # Age groups for detailed analysis
  age_groups_analysis = c("65-69", "+90")
)

# === INTERACTIVE PARAMETERS (defined at runtime) ===
CONST$INTERACTIVE <- list(
  model_type = NULL,           # "poisson_normal" or "normal_normal"
  temp_variable = NULL,         # "temp_norm" or "temp_extreme"
  plot_results = NULL,          # TRUE or FALSE
  n_simulations = 1000,
  n_ajustements = 100
)

# === COLOR PALETTES ===
CONST$COLORS <- list(
  state_colors = RColorBrewer::brewer.pal(9, "Set1"),
  comparison_colors = RColorBrewer::brewer.pal(3, "Dark2"),
  model_colors = c(
    "poisson_normal_2" = "#1b9e77",
    "normal_normal_2" = "#d95f02",
    "normal_normal_3" = "#7570b3"
  )
)

# === FUNCTION FOR ACCESSING CONSTANTS ===
# Usage: get_const("PRICING", "risk_free_rate")
get_const <- function(...) {
  args <- list(...)
  result <- CONST
  for (arg in args) {
    result <- result[[arg]]
    if (is.null(result)) {
      warning(sprintf("Constant not found: %s", paste(args, collapse = " > ")))
      return(NULL)
    }
  }
  return(result)
}

# === VALIDATION FUNCTION ===
validate_environment <- function() {
  cat("🔍 Validating environment...\n")
  
  checks <- list(
    "Data path valid" = file.exists(CONST$DATA$path_to_data),
    "Output directories" = all(file.exists(c(
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

# === DISPLAY OF INFORMATION ===
print_constants <- function() {
  cat("\n📋 GLOBAL CONSTANTS LOADED\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  
  cat("📊 HMM:\n")
  cat(sprintf("   • States to test: %s\n", paste(CONST$HMM$n_states_to_fit, collapse = ", ")))
  cat(sprintf("   • Fittings: %d attempts\n", CONST$HMM$n_ajustements))
  
  cat("\n💰 PRICING:\n")
  cat(sprintf("   • Risk-free rate: %.2f%%\n", CONST$PRICING$risk_free_rate * 100))
  cat(sprintf("   • Age groups: %d\n", length(CONST$PRICING$age_labels)))
  
  cat("\n🎲 MONTE CARLO:\n")
  cat(sprintf("   • Default simulations: %s\n", format(CONST$MONTE_CARLO$n_simulations_default, big.mark = ",")))
  
  cat("\n📁 PATHS:\n")
  cat(sprintf("   • Images: %s\n", CONST$PATHS$output_images))
  cat(sprintf("   • Results: %s\n", CONST$PATHS$output_results))
}

cat("✅ Global constants loaded successfully!\n")
