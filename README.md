# ATP Competitive Balance

R code accompanying *An analysis of the evolution of competitive balance in elite professional male tennis over 2000–2025* (Guzmán Crespo-Vicente and Álvaro Briz-Redón).

The main analysis uses Bayesian logistic regression with time-varying effects represented by second-order random walks (RW2), fitted with **NIMBLE**.

## Data

The ATP match data are available on [Kaggle](https://www.kaggle.com/datasets/dissfya/atp-tennis-2000-2023daily-pull). Download `atp_tennis.csv` and place it in a local `Data/` folder; the dataset is not redistributed here.

**Data version:** The manuscript uses a 2000–2025 snapshot, retaining 40,988 matches between players ranked in the Top 100 and excluding carpet courts. The Kaggle dataset is updated periodically, so a newer download may not give exactly the same sample or numerical results.

## R scripts

- `R/01_fit_main_RW2.R`: fit the main RW2 model with four MCMC chains.
- `R/02_competitive_balance.R`: compute group-level and rank-specific competitive-balance summaries.
- `R/03_player_scores.R`: compute player-level descriptive scores and in-sample predictions.
- `R/04_heatmaps.R`: create effect-difference heatmaps with posterior-uncertainty markers.
- `R/05_exceedance.R`: create exceedance-probability heatmaps.

## NIMBLE model specifications

- `NIMBLE codes/basic_random_effects_model_revised_ncp.R`: **main model**; full regression with RW2 temporal effects, using a non-centred parameterisation.
- `NIMBLE codes/rankdiff_only_model.R`: **RankDiff-only benchmark**; logistic model using the ranking difference as its sole predictor.
- `NIMBLE codes/basic_static_full_model.R`: **static full model**; same fixed predictors as the main model, without temporal effects.
- `NIMBLE codes/basic_linear_time_model.R`: **linear-time model**; full fixed predictors with linear temporal variation instead of RW2 effects.
- `NIMBLE codes/model_player_rw2_ncp_v2.R`: **player-effects sensitivity model**; the main RW2 structure extended with player-specific random effects.

The prior-sensitivity analysis uses the same main RW2 model with `sigma_prior_sd = 0.1` rather than `1`. The Top 200 and COVID-period exclusions likewise use the main RW2 specification, but with different analysis samples. They therefore do not need separate NIMBLE code files.

## Use and scope

Run the R scripts from the **repository root**, for example:

```r
source("R/01_fit_main_RW2.R")
source("R/02_competitive_balance.R")
```

The main fitting script saves `Models/Basic_random_effects_model_revised_ncp.rda`; the remaining R scripts use this saved fit. They require packages including `nimble`, `dplyr`, `lubridate`, `ggplot2`, `tidyr`, and `scales`. Model fitting can be computationally demanding.

This repository contains R scripts for the main Bayesian analysis and the NIMBLE specifications of the main model, alternative models, and sensitivity analyses considered in the paper. Additional fitting scripts for the alternative models, scripts for supplementary sensitivity tables and figures, and saved MCMC samples are not included.

The optional RankDiff-only predictive comparison in `R/03_player_scores.R` runs only if a separately fitted `Models/RankDiff_only_model.rda` is available. The corresponding model **specification** is included above, but the fitted object is not distributed here.
