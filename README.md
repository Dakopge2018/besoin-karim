# Modèle HMM d'Assurance Décès - Architecture Modulaire

## Structure

```
.
├── init.r                          # Initialisation: source(const.r, utilities.r)
├── data/                           # Données brutes (Excel)
├── src/
│   ├── const.r                     # Constantes globales - CONST$ 
│   ├── utilities.r                 # Fonctions réutilisables (7 sections)
│   ├── main_model_comparison.r    # Programme 1: ajuster & comparer 2 vs 3 états
│   ├── main_pricing.r              # Programme 2: calculer primes & sensibilité
│   └── README.md                   # Doc technique des modules
└── results/                        # Outputs (CSV, PDF, RDS)

```

## Démarrage rapide

```r
# 1. Initialiser
source("./init.r")

# 2. Ajuster modèles HMM (2 états vs 3 états)
source("./src/main_model_comparison.r")
# Génère: ./data/df_prepared.rds, ./models/hmm_model_*.rds

# 3. Calculer primes & sensibilité
source("./src/main_pricing.r")
# Génère: ./results/premiums_*.csv, ./results/sensitivity_*.pdf
```

## Modules principaux

### const.r (150 lignes)
Centralise tous les paramètres globaux dans `CONST$` :
- `CONST$DATA` : chemins, feuilles Excel, seed
- `CONST$HMM` : n_states, degrees, périodes, distribution
- `CONST$PRICING` : taux sans risque, groupes d'âge, expositions
- `CONST$PATHS` : répertoires de sortie
- `CONST$MONTE_CARLO`, `CONST$VIZ`, `CONST$COLORS`

**Modification facile** : changer `CONST$PRICING$risk_free_rate = 0.025` dans const.r et re-sourcer.

### utilities.r (400+ lignes)
Fonctions partagées entre les deux programmes, organisées en 7 sections :

1. **Data prep** : `load_and_prepare_data()`, `generate_trig_covariates()`, `prepare_hmm_dataframe()`
2. **Initialization** : `calculer_parametres_initiaux_diversifies_simple()` (4 stratégies)
3. **Visualization** : `sauvegarder_graphiques()`, `create_state_palette()`
4. **Model stats** : `calculate_aic()`, `calculate_bic()`, `extract_model_statistics()`
5. **Simulation** : `simulate_from_model()`
6. **Comparison** : `compare_models_summary()`
7. **Pricing** : `precompute_mortality_coefficients()`, `compute_mortality_intensity()`

### main_model_comparison.r (300+ lignes)

**Fonction auxiliaire** : `main_fit_hmm_model(data, nb_etats, n_simulations, plot_results, n_ajustements)`
- Boucle 150× avec paramètres diversifiés
- Sélectionne le meilleur modèle par AIC
- Retourne : modèle optimal, métriques, graphiques

**Main** : orchestration 5 étapes
1. Charger & préparer données → `./data/df_prepared.rds`
2. Ajuster modèles (2 et 3 états)
3. Comparer résultats (table AIC/BIC)
4. Générer visualisations convergence
5. Sauvegarder modèles → `./models/hmm_model_*.rds`

### main_pricing.r (350+ lignes)

**Fonction auxiliaire 1** : `main_init_pricing_model(hmm_model, n_states)`
- Interface retournée avec méthodes :
  - `compute_survival_probability(age, time_horizon, initial_state)`
  - `calculate_premium(age, benefit, duration)`
  - `set_xi_matrix(xi_H, homo)`
  - `get_info()`

**Fonction auxiliaire 2** : `main_sensitivity_analysis(pricing_model, age, benefit, duration, n_states)`
- Teste 100 combinaisons (10×10) de probabilités de transition
- Retourne matrice de primes

**Main** : orchestration 5 étapes
1. Charger modèles HMM
2. Initialiser moteurs tarifaires (2 et 3 états)
3. Calculer primes de base (4 scénarios)
4. Analyse sensibilité (100 combinaisons)
5. Générer heatmaps → `./results/sensitivity_*.pdf`

## Workflow typique

```r
# 1. Initialisation
source("./init.r")

# 2. Ajuster modèles (15-30 min)
source("./src/main_model_comparison.r")

# 3. Calculer primes avec sensibilité (5-10 min)
source("./src/main_pricing.r")

# 4. Modifier paramètres et réitérer
# Exemple: changer le taux sans risque
CONST$PRICING$risk_free_rate <- 0.03
source("./src/main_pricing.r")
```

## Modification des paramètres

Tous les paramètres sont dans `const.r`. Les plus courants :

```r
# Groupes d'âge
CONST$PRICING$age_breaks <- c(0, 20, 30, ..., 100)
CONST$PRICING$age_labels <- c("0-19", "20-29", ...)

# Taux sans risque
CONST$PRICING$risk_free_rate <- 0.025

# Robustesse HMM (nombre d'ajustements)
CONST$HMM$n_ajustements <- 150

# Nombre d'états à tester
CONST$HMM$n_states_to_fit <- c(2, 3)
```

Puis re-source le programme voulu.

## Outputs

**main_model_comparison.r** :
- `./data/df_prepared.rds` : données préparées
- `./models/hmm_model_2states.rds` : modèle 2 états
- `./models/hmm_model_3states.rds` : modèle 3 états
- Console : comparaison AIC/BIC

**main_pricing.r** :
- `./results/premiums_2states.csv` : primes base (2 états)
- `./results/premiums_3states.csv` : primes base (3 états)
- `./results/sensitivity_2states.csv` : matrice sensibilité (2 états)
- `./results/sensitivity_3states.csv` : matrice sensibilité (3 états)
- `./results/sensitivity_*.pdf` : heatmaps

## Documentation technique

Voir `src/README.md` pour les détails de chaque fonction.

## Dépannage

- **Erreur "const.r not found"** : exécuter `source("./init.r")` en premier
- **Modèles HMM non trouvés** : exécuter `source("./src/main_model_comparison.r")` d'abord
- **Primes vides** : vérifier que `./models/hmm_model_*.rds` existent
- **Graphiques non générés** : vérifier que `./results/` existe

---

**Créé** : Restructuration d'une monolith 2686 lignes en 5 modules clairs + paramètre centralisé.
