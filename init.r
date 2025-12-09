# ============================================================================
# Initialisation - Charger les modules essentiels
# ============================================================================

# Charger les constantes globales
source("./src/const.r")

# Charger les fonctions utilitaires
source("./src/utilities.r")

# Créer répertoires s'ils n'existent pas
dir.create("./models", showWarnings = FALSE)
dir.create("./results", showWarnings = FALSE)

cat("\n Environnement initialisé. Exécutez maintenant:\n")
cat("  source('./src/main_model_comparison.r')  # Ajuster modèles HMM\n")
cat("  source('./src/main_pricing.r')            # Calculer primes\n\n")
