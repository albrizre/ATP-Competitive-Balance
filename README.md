# ATP Competitive Balance

R code accompanying *An analysis of the evolution of competitive balance in elite professional male tennis over 2000–2025* (Guzmán Crespo-Vicente and Álvaro Briz-Redón).

The main analysis uses a Bayesian logistic regression with time-varying effects represented by second-order random walks (RW2), fitted with **NIMBLE**.

## Data

The ATP match data are publicly available on Kaggle: **[ATP Tennis 2000–2026 Daily Update](https://www.kaggle.com/datasets/dissfya/atp-tennis-2000-2023daily-pull)**. Download `atp_tennis.csv` and place it in a local `Data/` folder. The dataset is not redistributed here.

**Data-version note:** The Kaggle dataset is updated periodically. The manuscript uses the 2000–2025 data snapshot and a subset of **40,988 matches** with both players ranked in the Top 100, excluding carpet courts. A newer Kaggle download may not reproduce exactly the same sample or results.

## Files

- `R/01_fit_main_RW2.R`: fit the main RW2 model (four MCMC chains).
- `R/02_competitive_balance.R`: model-based group and rank-specific summaries.
- `R/03_player_scores.R`: player-level descriptive scores and in-sample prediction summaries.
- `R/04_heatmaps.R`: effect-difference heatmaps with posterior uncertainty markers.
- `R/05_exceedance.R`: exceedance-probability heatmaps.
- `NIMBLE codes/basic_random_effects_model_revised_ncp.R`: main model specification.

## Use

Run scripts from the **repository root**. For example:

```r
source("R/01_fit_main_RW2.R")
source("R/02_competitive_balance.R")
```

The fitting script writes `Models/Basic_random_effects_model_revised_ncp.rda`. The remaining scripts use that saved fit and do not refit the model. The scripts use R packages such as `nimble`, `posterior`, `dplyr`, `lubridate`, `ggplot2`, `tidyr`, and `scales`; some calculations are computationally demanding.

The optional RankDiff-only benchmark in `R/03_player_scores.R` is run only if its separately fitted model file is available; that comparison model is not included in this minimal repository. Sensitivity analyses and large MCMC objects are also not distributed. The repository therefore provides the **main-model analysis code**, not a complete reproduction package for every comparison in the manuscript.
