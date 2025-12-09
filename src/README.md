# STRUCTURE DU PROJET - TARIFICATION D'ASSURANCE VIE

## 📋 Vue d'ensemble

Ce projet modélise et tarife des contrats d'assurance vie décès basés sur les modèles de Markov cachés (HMM).

### Architecture du projet

```
besoin karim/
├── src/
│   ├── const.r                    # Variables globales et constantes
│   ├── utilities.r                # Fonctions réutilisables
│   ├── main_model_comparison.r    # Programme 1: Comparaison des modèles
│   └── main_pricing.r             # Programme 2: Tarification
│
├── data/                          # Données préparées (généré)
├── models/                        # Modèles HMM ajustés (généré)
├── results/                       # Résultats et graphiques (généré)
│
└── README.md                      # Ce fichier
```

---

## 🔧 FICHIERS DE CONFIGURATION

### 1. **const.r** - Constantes globales

**Rôle**: Centralise tous les paramètres constants du projet pour faciliter la maintenance.

**Contenu principal**:
- `CONST$DATA`: Chemins de données et seeds
- `CONST$HMM`: Paramètres d'ajustement HMM (états, itérations, tolérance, etc.)
- `CONST$PRICING`: Taux sans risque, groupes d'âge, expositions
- `CONST$MONTE_CARLO`: Paramètres de simulation Monte Carlo
- `CONST$PATHS`: Répertoires de sortie
- `CONST$FORMULAS`: Formules des modèles
- `CONST$VIZ`: Options de visualisation
- `CONST$COLORS`: Palettes de couleurs

**Utilisation**:
```r
source("./src/const.r")
taux_risque <- CONST$PRICING$risk_free_rate
states_a_tester <- CONST$HMM$n_states_to_fit
```

---

### 2. **utilities.r** - Fonctions réutilisables

**Rôle**: Bibliothèque de fonctions communes utilisées par les deux programmes principaux.

**Sections**:

#### Section 1: PRÉPARATION DES DONNÉES
- `load_and_prepare_data()`: Charge et nettoie les données brutes
- `generate_trig_covariates()`: Génère les covariables saisonnières (sin/cos)
- `prepare_hmm_dataframe()`: Prépare le dataframe pour HMM

#### Section 2: PARAMÈTRES INITIAUX
- `calculer_parametres_initiaux_diversifies_simple()`: Génère 4 stratégies de paramètres initiaux

#### Section 3: VISUALISATIONS
- `sauvegarder_graphiques()`: Exporte les graphiques en PDF/PNG
- `create_state_palette()`: Crée une palette adaptée aux états

#### Section 4: STATISTIQUES
- `calculate_aic()`, `calculate_bic()`: Critères d'information
- `extract_model_statistics()`: Extrait stats HMM

#### Section 5: SIMULATION
- `simulate_from_model()`: Simule depuis un HMM ajusté

#### Section 6: COMPARAISON
- `compare_models_summary()`: Résumé comparatif des modèles

#### Section 7: TARIFICATION
- `precompute_mortality_coefficients()`: Pré-calcule les coefficients
- `compute_mortality_intensity()`: Calcule l'intensité de mortalité

---

## 📊 PROGRAMMES PRINCIPAUX

### 3. **main_model_comparison.r** - Programme 1: Comparaison des modèles

**Objectif**: Estimer et comparer les modèles HMM avec 2 et 3 états.

**Flux d'exécution**:

```
MAIN()
├─ ÉTAPE 1: Chargement et préparation des données
│  ├─ load_and_prepare_data()
│  ├─ prepare_hmm_dataframe()
│  └─ Sauvegarder df dans ./data/df_prepared.rds
│
├─ ÉTAPE 2: Ajustement des modèles
│  ├─ main_fit_hmm_model() pour 2 états
│  │  ├─ Boucle: 150 tentatives d'ajustement (diversification)
│  │  ├─ Sélection du meilleur par AIC
│  │  └─ Génération des graphiques de convergence
│  │
│  └─ main_fit_hmm_model() pour 3 états
│     └─ Idem
│
├─ ÉTAPE 3: Comparaison des modèles
│  └─ compare_models_summary()
│     → Tableau AIC/BIC/LogLik
│
├─ ÉTAPE 4: Visualisations
│  └─ Graphiques de convergence et comparaison
│
└─ ÉTAPE 5: Sauvegarde des résultats
   ├─ ./models/hmm_model_2states.rds
   ├─ ./models/hmm_model_3states.rds
   ├─ ./results/model_comparison.csv
   └─ Graphiques en PDF
```

**Fonction auxiliaire**:
```r
main_fit_hmm_model(data, nb_etats, n_simulations, plot_results, n_ajustements)
├─ Ajuste un modèle HMM
├─ Sélectionne le meilleur
└─ Génère graphiques
```

**Sortie**:
- Modèles ajustés en `.rds`
- Tableau comparatif en `.csv`
- Graphiques de convergence en `.pdf`

**Utilisation**:
```bash
Rscript src/main_model_comparison.r
```

---

### 4. **main_pricing.r** - Programme 2: Tarification

**Objectif**: Tarifier les contrats d'assurance en fonction des probabilités de persistance.

**Flux d'exécution**:

```
MAIN()
├─ ÉTAPE 1: Charger les modèles HMM
│  └─ readRDS("./models/hmm_model_*.rds")
│
├─ ÉTAPE 2: Initialiser les moteurs de tarification
│  └─ main_init_pricing_model(hmm_model, n_states)
│     ├─ Pré-calcule les coefficients de mortalité
│     ├─ Initialise la matrice de transition par défaut
│     └─ Expose l'interface publique
│
├─ ÉTAPE 3: Calcul des primes (scénarios de base)
│  ├─ pricing_model$calculate_premium()
│  └─ Pour 4 scénarios d'âge et durée
│
├─ ÉTAPE 4: Analyse de sensibilité
│  └─ main_sensitivity_analysis()
│     ├─ Teste 100 combinaisons de probabilités
│     ├─ Crée une matrice n_p × n_p
│     └─ Exporte en CSV
│
├─ ÉTAPE 5: Visualisations
│  └─ Heatmap de sensibilité
│
└─ Sauvegarde: ./results/sensitivity_*.csv
```

**Fonction auxiliaire 1**:
```r
main_init_pricing_model(hmm_model, n_states)
├─ Initialise l'environnement privé
├─ Pré-calcule les coefficients
├─ Retourne l'interface publique avec:
│  ├─ compute_survival_probability(age, time_horizon, initial_state)
│  ├─ calculate_premium(age, benefit, duration, initial_state)
│  ├─ set_xi_matrix(xi_H, homo=FALSE)
│  └─ get_info()
```

**Fonction auxiliaire 2**:
```r
main_sensitivity_analysis(pricing_model, age, benefit, duration, n_states)
├─ Teste 10×10 = 100 combinaisons de probabilités
├─ Retourne une matrice de primes
└─ Exporte en CSV
```

**Interface de tarification**:
```r
# Charger après exécution
results <- source("./src/main_pricing.r")
pricing_model <- results$pricing_2states

# Utiliser
result <- pricing_model$calculate_premium(
  age_start = 65,
  death_benefit = 100000,
  contract_duration = 40,
  initial_state = 1,
  verbose = TRUE
)

# Modifier les transitions
xi_H <- matrix(c(0.95, 0.05, 0.05, 0.95), nrow=2)
pricing_model$set_xi_matrix(xi_H, homo=TRUE)
```

**Utilisation**:
```bash
Rscript src/main_pricing.r
```

---

## 🚀 FLUX DE TRAVAIL COMPLET

### Étape 1: Comparaison des modèles
```r
source("./src/main_model_comparison.r")
# Génère: ./models/hmm_model_*.rds
```

### Étape 2: Tarification
```r
source("./src/main_pricing.r")
# Génère: ./results/sensitivity_*.csv et graphiques
```

### Étape 3: Utiliser les modèles en R interactif
```r
source("./src/const.r")
source("./src/utilities.r")

# Charger les modèles
hmm_2 <- readRDS("./models/hmm_model_2states.rds")
df <- readRDS("./data/df_prepared.rds")

# Créer le moteur de tarification
source("./src/main_pricing.r")  # ou utiliser directly
```

---

## 📁 STRUCTURE DES FICHIERS GÉNÉRÉS

```
results/
├── model_comparison.csv          # Tableau comparatif
├── summary_2states.txt           # Résumé modèle 2 états
├── summary_3states.txt           # Résumé modèle 3 états
├── convergence_2states.pdf       # Graphique convergence
├── convergence_3states.pdf       # Graphique convergence
├── comparison_criteria.pdf       # Graphique AIC/BIC
├── sensitivity_2states.csv       # Matrice sensibilité 2 états
├── sensitivity_3states.csv       # Matrice sensibilité 3 états
└── heatmap_sensitivity_*.pdf     # Visualisation sensibilité

models/
├── hmm_model_2states.rds         # Modèle HMM 2 états
└── hmm_model_3states.rds         # Modèle HMM 3 états

data/
└── df_prepared.rds               # Dataframe préparé
```

---

## 🔑 CONCEPTS CLÉS

### 1. **Modèles HMM**
- **2 états**: État normal vs état perturbé
- **3 états**: Normal, Intermédiaire, Perturbé
- Distributions observées: Normal-Normal

### 2. **Tarification**
- **Méthode classique**: Discrète annuelle
- **Survie**: P(T > t) = exp(-∫μ(s)ds)
- **Prime**: E[D × v^T × 1_{T≤H}]

### 3. **Intensité de mortalité**
```
μ(t, age, state) = exp(β₀ + β₁ sin(t) + β₂ cos(t) + β₃ trend + β₄ age_factor)
```

### 4. **Probabilités de persistance**
- Matrice xi_H: Multiplicateur des transitions
- Transitions homogènes ou stochastiques

---

## 💡 CONSEILS D'UTILISATION

### Pour modifier les constantes:
```r
# Éditer const.r et changer CONST$PRICING$risk_free_rate = 0.03
```

### Pour ajouter une nouvelle fonction:
```r
# L'ajouter dans utilities.r dans la section appropriée
# et l'appeler dans main_model_comparison.r ou main_pricing.r
```

### Pour déboguer:
```r
# Ajouter verbose=TRUE dans les appels
# Activer le mode debug avec debugging_enabled <- TRUE
```

---

## 📈 PARAMÈTRES À AJUSTER

### Dans const.r:
- `CONST$HMM$n_ajustements`: Nombre de tentatives (défaut: 150)
- `CONST$HMM$maxit`: Itérations max (défaut: 1000)
- `CONST$PRICING$risk_free_rate`: Taux sans risque (défaut: 0.025)

### Dans les programs:
- `n_simulations`: Nombre de simulations Monte Carlo
- `age_ref`: Groupe d'âge de référence
- `p22_val`: Probabilités fixes (analyse sensibilité)

---

## 🐛 DÉPANNAGE

### Modèles ne convergent pas
→ Augmenter `CONST$HMM$n_ajustements`
→ Réduire `CONST$HMM$maxit`

### Primes incohérentes
→ Vérifier la matrice xi_H
→ Vérifier les risque-free-rate

### Fichiers manquants
→ Créer les répertoires: `./data`, `./models`, `./results`

---

**Version**: 1.0
**Dernière modification**: Décembre 2025
