library(dplyr)
library(nimble)
library(lubridate)
library(mgcv)
library(ggplot2)
library(parallel)


dir.create("Models", showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SETTINGS
# =============================================================================

# Main revised model: four chains, as requested by Reviewer 2.
n_chains <- 4
n_workers <- 2  # two chains run simultaneously
n_iter <- 15000
n_burnin <- 5000
thin_by <- 10

# Revised priors
beta_prior_sd <- 2.5
sigma_prior_sd <- 1

# Proper prior for the explicit linear component of each RW2.
# This is kept identical to the reparameterised player-effect sensitivity model.
time_slope_prior_sd <- 1

# Do not store pi[i] for every match during the four-chain fit.
monitor_pi <- FALSE

# =============================================================================
# DATA
# =============================================================================

atp_tennis <- read.csv(
  "Data/atp_tennis.csv",
  sep = ",",
  stringsAsFactors = FALSE
)

atp_tennis$Date <- as.Date(
  atp_tennis$Date,
  format = "%Y-%m-%d"
)

atp_tennis <- atp_tennis %>%
  filter(
    !is.na(Rank_1),
    !is.na(Rank_2),
    Rank_1 <= 100,
    Rank_2 <= 100,
    Rank_1 > 0,
    Rank_2 > 0
  ) %>%
  mutate(
    rank_best_raw = pmin(Rank_1, Rank_2),
    Y = if_else(
      Winner == if_else(Rank_1 < Rank_2, Player_1, Player_2),
      1,
      0
    ),
    rank_diff_raw = abs(Rank_1 - Rank_2)
  ) %>%
  filter(Surface != "Carpet") %>%
  mutate(
    Surface = droplevels(as.factor(Surface))
  )

atp_tennis$Best_of_5 <-
  ifelse(atp_tennis$Best.of == 5, 1, 0)

atp_tennis$Year <- year(atp_tennis$Date)

first_year <- min(atp_tennis$Year)
last_year <- max(atp_tennis$Year)

atp_tennis$Year_num <-
  atp_tennis$Year - first_year + 1

n_years <- length(
  unique(atp_tennis$Year_num)
)

# Years must form a complete consecutive sequence.
stopifnot(
  n_years == max(atp_tennis$Year_num)
)

atp_tennis$Clay <-
  as.numeric(atp_tennis$Surface == "Clay")

atp_tennis$Grass <-
  as.numeric(atp_tennis$Surface == "Grass")

atp_tennis$Surface_num <- 0
atp_tennis$Surface_num[atp_tennis$Surface == "Clay"] <- 1
atp_tennis$Surface_num[atp_tennis$Surface == "Grass"] <- 2
atp_tennis$Surface_num[atp_tennis$Surface == "Hard"] <- 3

# Same transformations as in the original revised model
atp_tennis$rank_best <-
  as.numeric(
    scale(
      log(atp_tennis$rank_best_raw + 1)
    )
  )

atp_tennis$rank_diff <-
  as.numeric(
    scale(
      log(atp_tennis$rank_diff_raw + 1)
    )
  )

cat("Matches:", nrow(atp_tennis), "\n")
cat(
  "Years:",
  first_year, "-", last_year,
  "(", n_years, "years )\n"
)

# =============================================================================
# WHITENED RW2 BASIS
# =============================================================================
#
# For T years:
#
#   Delta^2 f_t = f_t - 2 f_{t-1} + f_{t-2}.
#
# Let D2 be the second-difference matrix and Q = D2' D2.
# Q is the intrinsic RW2 precision matrix and has rank T - 2, with null
# space span{1, t}.
#
# On the positive-eigenvalue subspace, if Q = V Lambda V', define
#
#   B = V Lambda^{-1/2}.
#
# Then B' Q B = I. Hence sigma * B * z, z ~ N(0,I), retains the RW2
# quadratic penalty while using whitened non-centred coordinates.
# =============================================================================

T_year <- n_years

D2 <- matrix(
  0,
  nrow = T_year - 2,
  ncol = T_year
)

for(r in 1:(T_year - 2)) {
  D2[r, r] <- 1
  D2[r, r + 1] <- -2
  D2[r, r + 2] <- 1
}

Q_rw2 <- crossprod(D2)

eig <- eigen(
  Q_rw2,
  symmetric = TRUE
)

tol <- max(eig$values) * 1e-10

positive <- which(
  eig$values > tol
)

n_rw2 <- length(positive)

stopifnot(
  n_rw2 == T_year - 2
)

V_pos <- eig$vectors[
  , positive,
  drop = FALSE
]

lambda_pos <- eig$values[positive]

B_rw2 <- sweep(
  V_pos,
  2,
  sqrt(lambda_pos),
  FUN = "/"
)

# Explicit linear null-space component.
# Centering/scaling leaves it in the RW2 null space and makes the slope prior
# interpretable on a convenient time scale.
time_scaled <- as.numeric(
  scale(seq_len(T_year))
)

# -----------------------------------------------------------------------------
# Numerical checks
# -----------------------------------------------------------------------------

ones <- rep(1, T_year)

time_centered_raw <-
  seq_len(T_year) -
  mean(seq_len(T_year))

constant_orthogonality <-
  max(
    abs(
      crossprod(
        ones,
        B_rw2
      )
    )
  )

linear_orthogonality <-
  max(
    abs(
      crossprod(
        time_centered_raw,
        B_rw2
      )
    )
  )

penalty_check <-
  max(
    abs(
      crossprod(
        B_rw2,
        Q_rw2 %*% B_rw2
      ) -
      diag(n_rw2)
    )
  )

cat(
  "\nRW2 basis checks:\n",
  "rank(Q) =", n_rw2, "\n",
  "max |1'B| =", constant_orthogonality, "\n",
  "max |t'B| =", linear_orthogonality, "\n",
  "max |B'QB-I| =", penalty_check, "\n\n"
)

stopifnot(
  constant_orthogonality < 1e-7,
  linear_orthogonality < 1e-7,
  penalty_check < 1e-7
)

# =============================================================================
# NIMBLE INPUTS
# =============================================================================

data <- list(
  Y = atp_tennis$Y
)

constants <- list(
  rank_best = atp_tennis$rank_best,
  rank_diff = atp_tennis$rank_diff,

  Clay = atp_tennis$Clay,
  Grass = atp_tennis$Grass,
  Best_of_5 = atp_tennis$Best_of_5,

  # Kept for compatibility with downstream code if needed
  Surface_num = atp_tennis$Surface_num,

  n = nrow(atp_tennis),

  n_years = n_years,
  Year_num = atp_tennis$Year_num,

  n_rw2 = n_rw2,
  B_rw2 = B_rw2,
  time_scaled = time_scaled,

  beta_prior_sd = beta_prior_sd,
  sigma_prior_sd = sigma_prior_sd,
  time_slope_prior_sd = time_slope_prior_sd
)

# =============================================================================
# INITIAL VALUES
# =============================================================================

inits <- function() {

  list(
    beta0 = rnorm(1, 0, 0.1),
    beta_best = rnorm(1, 0, 0.1),
    beta_diff = rnorm(1, 0, 0.1),
    beta_inter = rnorm(1, 0, 0.1),

    beta_clay = rnorm(1, 0, 0.1),
    beta_grass = rnorm(1, 0, 0.1),
    beta_bestof = rnorm(1, 0, 0.1),

    beta_clay_best = rnorm(1, 0, 0.1),
    beta_grass_best = rnorm(1, 0, 0.1),
    beta_bestof_best = rnorm(1, 0, 0.1),

    beta_clay_diff = rnorm(1, 0, 0.1),
    beta_grass_diff = rnorm(1, 0, 0.1),
    beta_bestof_diff = rnorm(1, 0, 0.1),

    beta_clay_inter = rnorm(1, 0, 0.1),
    beta_grass_inter = rnorm(1, 0, 0.1),
    beta_bestof_inter = rnorm(1, 0, 0.1),

    # Explicit linear RW2 components
    slope_beta0 = rnorm(1, 0, 0.05),
    slope_beta_best = rnorm(1, 0, 0.05),
    slope_beta_diff = rnorm(1, 0, 0.05),
    slope_beta_inter = rnorm(1, 0, 0.05),

    # Whitened curvature coordinates
    z_beta0 = rnorm(n_rw2, 0, 0.10),
    z_beta_best = rnorm(n_rw2, 0, 0.10),
    z_beta_diff = rnorm(n_rw2, 0, 0.10),
    z_beta_inter = rnorm(n_rw2, 0, 0.10),

    # RW2 scale parameters
    sigma.beta0_year = runif(1, 0.02, 0.10),
    sigma.beta_year_best = runif(1, 0.01, 0.08),
    sigma.beta_year_diff = runif(1, 0.005, 0.05),
    sigma.beta_year_inter = runif(1, 0.005, 0.05)
  )
}

# =============================================================================
# MODEL
# =============================================================================

source(
  "NIMBLE codes/basic_random_effects_model_revised_ncp.R"
)

# Monitor the same substantive parameters as before, plus the explicit
# non-centred RW2 parameters needed for proper convergence diagnostics.
parameter_monitors <- c(
  "beta0",
  "beta_clay",
  "beta_grass",
  "beta_bestof",

  "beta_clay_best",
  "beta_grass_best",
  "beta_bestof_best",

  "beta_clay_diff",
  "beta_grass_diff",
  "beta_bestof_diff",

  "beta_clay_inter",
  "beta_grass_inter",
  "beta_bestof_inter",

  "beta_best",
  "beta_diff",
  "beta_inter",

  "slope_beta0",
  "slope_beta_best",
  "slope_beta_diff",
  "slope_beta_inter",

  "z_beta0",
  "z_beta_best",
  "z_beta_diff",
  "z_beta_inter",

  "sigma.beta0_year",
  "sigma.beta_year_best",
  "sigma.beta_year_diff",
  "sigma.beta_year_inter",

  # Deterministic trajectories, directly comparable with the previous model
  "beta0_year",
  "beta_year_best",
  "beta_year_diff",
  "beta_year_inter"
)

monitors <- if(monitor_pi) {
  c(
    parameter_monitors,
    "pi"
  )
} else {
  parameter_monitors
}


# =============================================================================
# PARALLEL MCMC
# =============================================================================
#
# Four chains are fitted in total, but only two R/NIMBLE processes run at the
# same time. Each worker fits two chains sequentially. This follows NIMBLE's
# recommendation that model construction/compilation occur inside the
# parallel worker, while avoiding four simultaneous heavy processes.
#
# Online WAIC is disabled here because the four chains are split across two
# independent workers. If needed, WAIC can be calculated post hoc from the
# combined posterior draws (see optional block below).
# =============================================================================

cat(
  "
Running main reparameterised RW2 model in parallel:
",
  "total chains =", n_chains, "
",
  "simultaneous workers =", n_workers, "
",
  "iterations per chain =", n_iter, "
",
  "burn-in per chain =", n_burnin, "
",
  "thin =", thin_by, "

"
)

stopifnot(n_chains == 4L, n_workers == 2L)

project_dir <- getwd()
model_file <- file.path(
  project_dir,
  "NIMBLE codes",
  "basic_random_effects_model_revised_ncp.R"
)
model_object_name <- "basic_random_effects_model_ncp"

if(!file.exists(model_file)) {
  stop("NIMBLE model file not found: ", model_file)
}

# Generate all initial values in the master R session before opening workers.
set.seed(20261008)
chain_inits <- lapply(seq_len(n_chains), function(i) inits())

# Reproducible, distinct seeds for the four chains.
chain_seeds <- 20261008 + seq_len(n_chains)

# Worker 1 runs chains 1-2; worker 2 runs chains 3-4.
chain_batches <- split(
  seq_len(n_chains),
  rep(seq_len(n_workers), each = n_chains / n_workers)
)

run_chain_batch <- function(chain_ids) {

  suppressPackageStartupMessages(library(nimble))

  setwd(project_dir)

  model_env <- new.env(parent = globalenv())
  sys.source(model_file, envir = model_env)
  model_code <- get(model_object_name, envir = model_env)

  nimbleMCMC(
    code = model_code,
    data = data,
    constants = constants,
    inits = chain_inits[chain_ids],
    monitors = monitors,
    niter = n_iter,
    nburnin = n_burnin,
    nchains = length(chain_ids),
    thin = thin_by,
    setSeed = chain_seeds[chain_ids],
    progressBar = FALSE,
    summary = TRUE,
    WAIC = FALSE
  )
}

start_time <- Sys.time()
print(start_time)

cl <- parallel::makeCluster(n_workers)

parallel::clusterExport(
  cl,
  varlist = c(
    "project_dir",
    "model_file",
    "model_object_name",
    "data",
    "constants",
    "chain_inits",
    "chain_seeds",
    "monitors",
    "n_iter",
    "n_burnin",
    "thin_by"
  ),
  envir = environment()
)

worker_results <- tryCatch(
  parallel::parLapply(
    cl,
    chain_batches,
    run_chain_batch
  ),
  finally = parallel::stopCluster(cl)
)

end_time <- Sys.time()
print(end_time)
print(end_time - start_time)

# Collect the four chain matrices in their original order.
samples_list <- unlist(
  lapply(
    worker_results,
    function(x) {
      s <- x$samples
      if(is.matrix(s)) list(s) else s
    }
  ),
  recursive = FALSE
)

stopifnot(length(samples_list) == n_chains)

samples_all <- do.call(
  rbind,
  samples_list
)

# Keep an object with the same component used by the downstream code.
# Worker-specific summaries are retained as well.
mcmc.output <- list(
  samples = samples_list,
  summary_by_worker = lapply(worker_results, function(x) x$summary)
)

save(
  mcmc.output,
  B_rw2,
  Q_rw2,
  time_scaled,
  file = "Models/Basic_random_effects_model_revised_ncp.rda"
)

# -----------------------------------------------------------------------------
# OPTIONAL: post-hoc WAIC
# -----------------------------------------------------------------------------
# Because the chains are split across two workers, online WAIC is not combined
# automatically. If WAIC is needed for model comparison, set this to TRUE.
# NIMBLE will then calculate conditional WAIC from the combined posterior draws.
calculate_waic_after <- FALSE

if(calculate_waic_after) {

  waic_model <- nimbleModel(
    basic_random_effects_model_ncp,
    data = data,
    constants = constants,
    inits = chain_inits[[1]],
    check = TRUE
  )

  WAIC_posthoc <- nimble::calculateWAIC(
    samples_all,
    waic_model
  )

  print(WAIC_posthoc)

  save(
    WAIC_posthoc,
    file = "Models/WAIC_revised_ncp.rda"
  )
}


# =============================================================================
# POSTERIOR SUMMARIES
# =============================================================================

beta0_year <- colMeans(
  samples_all[
    , grep(
      "^beta0_year\\[",
      colnames(samples_all)
    ),
    drop = FALSE
  ]
)

beta_year_best <- colMeans(
  samples_all[
    , grep(
      "^beta_year_best\\[",
      colnames(samples_all)
    ),
    drop = FALSE
  ]
)

beta_year_diff <- colMeans(
  samples_all[
    , grep(
      "^beta_year_diff\\[",
      colnames(samples_all)
    ),
    drop = FALSE
  ]
)

beta_year_inter <- colMeans(
  samples_all[
    , grep(
      "^beta_year_inter\\[",
      colnames(samples_all)
    ),
    drop = FALSE
  ]
)

par(mfrow = c(2, 2))

plot(
  first_year:last_year,
  beta0_year,
  type = "l",
  xlab = "Year",
  ylab = "",
  main = "beta0_year"
)

plot(
  first_year:last_year,
  beta_year_best,
  type = "l",
  xlab = "Year",
  ylab = "",
  main = "beta_year_best"
)

plot(
  first_year:last_year,
  beta_year_diff,
  type = "l",
  xlab = "Year",
  ylab = "",
  main = "beta_year_diff"
)

plot(
  first_year:last_year,
  beta_year_inter,
  type = "l",
  xlab = "Year",
  ylab = "",
  main = "beta_year_inter"
)

par(mfrow = c(1, 1))

bestof_parameters <- c(
  "beta_bestof",
  "beta_bestof_best",
  "beta_bestof_diff",
  "beta_bestof_inter"
)

print(
  summary(
    samples_all[
      , bestof_parameters,
      drop = FALSE
    ]
  )
)

# =============================================================================
# R-HAT, BULK-ESS AND TAIL-ESS
# =============================================================================

if(requireNamespace("posterior", quietly = TRUE)) {

  diagnostic_names <-
    colnames(samples_list[[1]])[
      !grepl(
        "^pi\\[",
        colnames(samples_list[[1]])
      )
    ]

  n_saved <- nrow(samples_list[[1]])
  n_chains_saved <- length(samples_list)
  n_variables <- length(diagnostic_names)

  draws_array <- array(
    NA_real_,
    dim = c(
      n_saved,
      n_chains_saved,
      n_variables
    ),
    dimnames = list(
      NULL,
      NULL,
      diagnostic_names
    )
  )

  for(ch in seq_len(n_chains_saved)) {
    draws_array[, ch, ] <-
      samples_list[[ch]][
        , diagnostic_names,
        drop = FALSE
      ]
  }

  draws <-
    posterior::as_draws_array(
      draws_array
    )

  diagnostics <-
    posterior::summarise_draws(
      draws,
      "rhat",
      "ess_bulk",
      "ess_tail"
    )

  diag_df <-
    as.data.frame(
      diagnostics
    )

  write.csv(
    diag_df,
    "Models/MCMC_diagnostics_revised_ncp.csv",
    row.names = FALSE
  )

  cat(
    "\nOVERALL DIAGNOSTICS\n"
  )

  cat(
    "Maximum R-hat:",
    max(
      diag_df$rhat,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Minimum bulk-ESS:",
    min(
      diag_df$ess_bulk,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Minimum tail-ESS:",
    min(
      diag_df$ess_tail,
      na.rm = TRUE
    ),
    "\n"
  )

  # Sampled structural parameters only:
  # exclude deterministic temporal trajectories.
  structural_diag <-
    diag_df[
      !grepl(
        "^beta0_year\\[|^beta_year_best\\[|^beta_year_diff\\[|^beta_year_inter\\[",
        diag_df$variable
      ),
      ,
      drop = FALSE
    ]

  cat(
    "\nSAMPLED STRUCTURAL PARAMETERS\n"
  )

  cat(
    "Maximum R-hat:",
    max(
      structural_diag$rhat,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Minimum bulk-ESS:",
    min(
      structural_diag$ess_bulk,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Minimum tail-ESS:",
    min(
      structural_diag$ess_tail,
      na.rm = TRUE
    ),
    "\n"
  )

  bad_structural <-
    structural_diag[
      !is.na(
        structural_diag$rhat
      ) &
      structural_diag$rhat > 1.01,
      ,
      drop = FALSE
    ]

  cat(
    "\nSampled structural parameters with R-hat > 1.01:\n"
  )

  if(nrow(bad_structural) > 0) {
    print(
      bad_structural[
        order(
          bad_structural$rhat,
          decreasing = TRUE
        ),
      ]
    )
  } else {
    cat("None.\n")
  }

  # Deterministic temporal trajectories
  temporal_diag <-
    diag_df[
      grepl(
        "^beta0_year\\[|^beta_year_best\\[|^beta_year_diff\\[|^beta_year_inter\\[",
        diag_df$variable
      ),
      ,
      drop = FALSE
    ]

  cat(
    "\nTEMPORAL TRAJECTORIES\n"
  )

  cat(
    "Maximum R-hat:",
    max(
      temporal_diag$rhat,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Minimum bulk-ESS:",
    min(
      temporal_diag$ess_bulk,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Minimum tail-ESS:",
    min(
      temporal_diag$ess_tail,
      na.rm = TRUE
    ),
    "\n"
  )

} else {

  message(
    "Install the 'posterior' package to compute R-hat, bulk-ESS and tail-ESS."
  )
}
