# =============================================================================
# SECTION 4.5: REVISED RW2 - PLAYER-LEVEL SCORES AND IN-SAMPLE PREDICTION
# =============================================================================
# Run from the Article/ project folder:
#   source("R/03_player_scores.R")
#
# Uses ONLY previously fitted models. Does NOT re-run MCMC.
#
# Required:
#   Data/atp_tennis.csv
#   Models/Basic_random_effects_model_revised_ncp.rda
#
# Optional for predictive benchmark (comparison model not included here):
#   Models/RankDiff_only_model.rda
#
# Outputs (never overwrite previous figures/models):
#   Models/Section_4_5_revised_RW2/
#      player_scores_all.csv
#      player_scores_eligible.csv
#      player_extremes.csv
#      player_extremes_wide.csv
#      Table_4_5_player_scores.tex
#      match_level_predictions.csv (optional)
#      predictive_performance.csv
#      confusion_matrices.csv
#      ROC_main_vs_RankDiff.png / .eps (only if benchmark is present)
#      data_and_checks.txt
#
# Methodological choices match the CURRENT manuscript (Section 4.5):
# * Top-100 against Top-100 matches, no Carpet; 2000-2025.
# * All 16 fixed coefficients and 4 temporal RW2 sequences are included.
# * Best_of_5 = 1 only when Best.of == 5 (NOT Best.of == 1).
# * Match probability = posterior mean E[invlogit(eta) | data].
#   IMPORTANT: it is NOT invlogit(E[eta | data]).
# * Favorite/underdog probabilities are oriented for EACH individual player.
# * Mean residual, surprise of wins, surprise of losses, net score.
# * Minimum 40 wins AND 40 losses as in the original player-score analysis.
# * The RankDiff-only benchmark, if present, is evaluated in sample.
# =============================================================================

# -------------------------- USER SETTINGS -----------------------------------
main_fit_path <- "Models/Basic_random_effects_model_revised_ncp.rda"
benchmark_fit_path <- "Models/RankDiff_only_model.rda"
data_path <- "Data/atp_tennis.csv"
out_dir <- "Models/Section_4_5_revised_RW2"

# 1,200 posterior draws (300 per chain when 4 chains are present) are usually
# sufficient for this initial descriptive table; set to Inf for ALL saved draws.
posterior_draws <- 1200L
seed <- 20261022L
block_size <- 600L
min_wins <- 40L
min_losses <- 40L
save_match_level_predictions <- TRUE

# If RStudio's working directory is not Article/, set it manually here:

if(!file.exists(main_fit_path) || !file.exists(data_path)) {
  stop("Model or data not found. Execute the script from the Article/ directory.\n",
       "Current directory: ", getwd(), "\n",
       "Expected: ", main_fit_path, " and ", data_path)
}
stopifnot(is.numeric(block_size), block_size >= 1,
          min_wins >= 1L, min_losses >= 1L)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

# ---------------------- DATA / REPRODUCE FILTER -----------------------------
raw <- read.csv(data_path, stringsAsFactors = FALSE)
required_fields <- c("Date", "Rank_1", "Rank_2", "Player_1", "Player_2",
                     "Winner", "Surface", "Best.of")
missing_fields <- setdiff(required_fields, names(raw))
if(length(missing_fields)) stop("Missing data columns: ", paste(missing_fields, collapse = ", "))

# Preserve EXACT filtering/order used to fit the principal model.
d <- raw[!is.na(raw$Rank_1) & !is.na(raw$Rank_2) &
           raw$Rank_1 > 0 & raw$Rank_2 > 0 &
           raw$Rank_1 <= 100 & raw$Rank_2 <= 100, , drop = FALSE]
d <- d[!is.na(d$Surface) & d$Surface != "Carpet", , drop = FALSE]
rownames(d) <- NULL

# Check the expected main-model sample: fail before predictions if it differs.
if(nrow(d) != 40988L) {
  stop("Sample has ", nrow(d), " matches, but main model used 40,988. ",
       "Check data version and filtering before calculating results.")
}

d$Date <- as.Date(d$Date, format = "%Y-%m-%d")
d$Year <- as.integer(format(d$Date, "%Y"))
if(anyNA(d$Date) || !identical(sort(unique(d$Year)), 2000:2025)) {
  stop("Unexpected or missing match dates / year range (expected 2000:2025).")
}
d$Year_num <- d$Year - 1999L

if(anyNA(d$Player_1) || anyNA(d$Player_2) || anyNA(d$Winner) ||
   any(!nzchar(trimws(d$Player_1))) || any(!nzchar(trimws(d$Player_2)))) {
  stop("Missing player names or outcomes in selected matches.")
}
# Like the original scripts, trim external spaces but do not rewrite identities.
d$Player_1 <- trimws(d$Player_1)
d$Player_2 <- trimws(d$Player_2)
d$Winner <- trimws(d$Winner)
if(any(!(d$Winner == d$Player_1 | d$Winner == d$Player_2))) {
  stop("At least one recorded winner does not match either player: check labels.")
}
if(any(d$Player_1 == d$Player_2)) stop("Same player identified on both sides of a match.")

# Same convention as original main-model likelihood: if ranks are equal, Player_2
# is classified as the 'favorite'. This is flagged for the manuscript audit.
d$favorite_is_P1 <- d$Rank_1 < d$Rank_2
d$Y <- as.integer(ifelse(d$favorite_is_P1,
                         d$Winner == d$Player_1,
                         d$Winner == d$Player_2))
n_tied_rank <- sum(d$Rank_1 == d$Rank_2)
if(n_tied_rank > 0L) {
  warning(n_tied_rank, " matches have tied rankings. ",
          "Following original model convention: Player_2 is treated as favorite.")
}

if(!all(d$Surface %in% c("Hard", "Clay", "Grass"))) {
  stop("Unexpected Surface value(s): ",
       paste(setdiff(unique(d$Surface), c("Hard", "Clay", "Grass")), collapse = ", "))
}
if(anyNA(d$Best.of)) stop("Missing Best.of values.")
d$Clay <- as.numeric(d$Surface == "Clay")
d$Grass <- as.numeric(d$Surface == "Grass")
d$Best_of_5 <- as.numeric(d$Best.of == 5)

# Same log transformation and SAME original Top-100 scaling, after filtering.
rb <- log(pmin(d$Rank_1, d$Rank_2) + 1)
rd <- log(abs(d$Rank_1 - d$Rank_2) + 1)
d$rank_best <- (rb - mean(rb)) / sd(rb)
d$rank_diff <- (rd - mean(rd)) / sd(rd)
stopifnot(all(is.finite(d$rank_best)), all(is.finite(d$rank_diff)))

cat("\nMain-model sample:", nrow(d), "matches; years 2000-2025.\n")
cat("Best-of-five matches:", sum(d$Best_of_5 == 1), "\n")
cat("Tied-ranking matches:", n_tied_rank, "\n")

# ---------------------- POSTERIOR LOADING / CHECKS ---------------------------
load_chains <- function(file) {
  e <- new.env(parent = baseenv())
  load(file, envir = e)
  if(!exists("mcmc.output", envir = e, inherits = FALSE)) {
    stop("Missing mcmc.output in ", file)
  }
  s <- e$mcmc.output$samples
  if(is.matrix(s) || inherits(s, "mcmc")) s <- list(s)
  if(!is.list(s) || !length(s)) stop("Unrecognized posterior chains in ", file)
  s <- lapply(s, as.matrix)
  nc <- colnames(s[[1]])
  if(is.null(nc) || anyDuplicated(nc)) stop("Missing/duplicate posterior parameter names.")
  for(j in seq_along(s)) {
    if(!identical(colnames(s[[j]]), nc)) stop("Chain parameter names differ in ", file)
    if(any(!is.finite(s[[j]]))) stop("Nonfinite posterior draws in chain ", j)
  }
  s
}

choose_draws <- function(chains, requested) {
  sizes <- vapply(chains, nrow, integer(1))
  n_total <- sum(sizes)
  if(is.infinite(requested) || requested >= n_total) return(do.call(rbind, chains))
  requested <- as.integer(requested)
  if(requested < length(chains)) stop("Request at least one draw per chain.")
  # Balanced allocation across chains, excluding accidental over-allocation.
  counts <- rep(requested %/% length(chains), length(chains))
  if(requested %% length(chains)) {
    counts[seq_len(requested %% length(chains))] <-
      counts[seq_len(requested %% length(chains))] + 1L
  }
  if(any(counts > sizes)) {
    stop("Requested too many draws from at least one chain; lower posterior_draws.")
  }
  do.call(rbind, Map(function(s, n) s[sample.int(nrow(s), n), , drop = FALSE],
                     chains, counts))
}

chains <- load_chains(main_fit_path)
if(length(chains) != 4L) {
  stop("Expected four main-model chains; found ", length(chains))
}
samples <- choose_draws(chains, posterior_draws)
cat("Selected posterior draws:", nrow(samples),
    "from", length(chains), "chains.\n")

fixed_names <- c(
  "beta0", "beta_best", "beta_diff", "beta_inter",
  "beta_bestof", "beta_clay", "beta_grass",
  "beta_clay_best", "beta_grass_best", "beta_bestof_best",
  "beta_clay_diff", "beta_grass_diff", "beta_bestof_diff",
  "beta_clay_inter", "beta_grass_inter", "beta_bestof_inter"
)

extract_year_matrix <- function(draws, prefix) {
  nm <- paste0(prefix, "[", seq_len(26L), "]")
  absent <- setdiff(nm, colnames(draws))
  if(length(absent)) stop("Missing time coefficients: ", paste(absent, collapse = ", "))
  draws[, nm, drop = FALSE]
}

missing_beta <- setdiff(fixed_names, colnames(samples))
if(length(missing_beta)) stop("Missing fixed betas: ", paste(missing_beta, collapse = ", "))
beta_draws <- samples[, fixed_names, drop = FALSE]
year_draws <- list(
  intercept = extract_year_matrix(samples, "beta0_year"),
  best = extract_year_matrix(samples, "beta_year_best"),
  diff = extract_year_matrix(samples, "beta_year_diff"),
  inter = extract_year_matrix(samples, "beta_year_inter")
)
rm(samples, chains)
gc(verbose = FALSE)

# -------------------- COMPLETE LINEAR PREDICTOR ------------------------------
# These are the 16 covariate terms appearing in
#   NIMBLE codes/basic_random_effects_model_revised_ncp.R
# Column order MUST equal 'fixed_names'.
X <- cbind(
  beta0 = 1,
  beta_best = d$rank_best,
  beta_diff = d$rank_diff,
  beta_inter = d$rank_best * d$rank_diff,
  beta_bestof = d$Best_of_5,
  beta_clay = d$Clay,
  beta_grass = d$Grass,
  beta_clay_best = d$Clay * d$rank_best,
  beta_grass_best = d$Grass * d$rank_best,
  beta_bestof_best = d$Best_of_5 * d$rank_best,
  beta_clay_diff = d$Clay * d$rank_diff,
  beta_grass_diff = d$Grass * d$rank_diff,
  beta_bestof_diff = d$Best_of_5 * d$rank_diff,
  beta_clay_inter = d$Clay * d$rank_best * d$rank_diff,
  beta_grass_inter = d$Grass * d$rank_best * d$rank_diff,
  beta_bestof_inter = d$Best_of_5 * d$rank_best * d$rank_diff
)
stopifnot(identical(colnames(X), fixed_names),
          ncol(X) == 16L, nrow(X) == nrow(d))

# Matrix block: posterior draws x match indices. Peak memory is controlled
# by block_size instead of constructing all draws x 40,988 matches at once.
posterior_mean_predictions <- function(beta, tcoefs, design, years, size = 600L) {
  n <- nrow(design)
  p <- numeric(n)
  n_blocks <- ceiling(n / size)
  for(k in seq_len(n_blocks)) {
    lo <- (k - 1L) * size + 1L
    hi <- min(k * size, n)
    ids <- lo:hi
    yy <- years[ids]
    rb <- design[ids, "beta_best"]
    rd <- design[ids, "beta_diff"]
    rbrd <- rb * rd
    eta <- beta %*% t(design[ids, , drop = FALSE])
    eta <- eta + tcoefs$intercept[, yy, drop = FALSE]
    eta <- eta + sweep(tcoefs$best[, yy, drop = FALSE], 2, rb, `*`)
    eta <- eta + sweep(tcoefs$diff[, yy, drop = FALSE], 2, rd, `*`)
    eta <- eta + sweep(tcoefs$inter[, yy, drop = FALSE], 2, rbrd, `*`)
    p[ids] <- colMeans(plogis(eta))
    if(k %% 10L == 0L || k == n_blocks) {
      cat(sprintf("  RW2 predictions: block %d/%d (%d matches)\n", k, n_blocks, hi))
    }
  }
  p
}

cat("\nComputing RW2 match probabilities using the COMPLETE model...\n")
p_favorite <- posterior_mean_predictions(beta_draws, year_draws, X,
                                          d$Year_num, block_size)
stopifnot(length(p_favorite) == nrow(d),
          all(is.finite(p_favorite)), all(p_favorite > 0), all(p_favorite < 1))

# ------------------ IN-SAMPLE PREDICTIVE PERFORMANCE ------------------------
# All results here are IN-SAMPLE; do not present as held-out model validation.
auc_from_ranks <- function(y, p) {
  if(length(unique(y)) != 2L) stop("AUC requires both classes.")
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  (sum(rank(p, ties.method = "average")[y == 1L]) - n1 * (n1 + 1) / 2) /
    (n1 * n0)
}

get_performance <- function(y, p, label) {
  pred <- as.integer(p >= 0.5)
  tp <- sum(pred == 1L & y == 1L)
  fp <- sum(pred == 1L & y == 0L)
  tn <- sum(pred == 0L & y == 0L)
  fn <- sum(pred == 0L & y == 1L)
  precision <- if(tp + fp > 0) tp / (tp + fp) else NA_real_
  recall <- if(tp + fn > 0) tp / (tp + fn) else NA_real_
  f1 <- if(is.finite(precision) && is.finite(recall) &&
           precision + recall > 0) 2 * precision * recall / (precision + recall)
        else NA_real_
  list(
    metrics = data.frame(Model = label, N = length(y),
                         Accuracy = mean(pred == y),
                         AUC = auc_from_ranks(y, p),
                         Brier = mean((y - p)^2),
                         Precision = precision, Recall = recall, F1 = f1),
    confusion = data.frame(Model = label, TP = tp, FP = fp, TN = tn, FN = fn)
  )
}

performance <- list(get_performance(d$Y, p_favorite, "Main RW2"))
all_predictions <- list(`Main RW2` = p_favorite)

if(file.exists(benchmark_fit_path)) {
  cat("\nComputing the RankDiff-only in-sample benchmark...\n")
  # This benchmark uses the SAME Top-100 sample and its standardised log rank gap.
  benchmark_chains <- load_chains(benchmark_fit_path)
  benchmark_draws <- choose_draws(benchmark_chains, posterior_draws)
  req <- c("beta0", "beta_diff")
  if(!all(req %in% colnames(benchmark_draws))) {
    stop("The RankDiff-only fit does not contain beta0 and beta_diff.")
  }
  bb <- benchmark_draws[, req, drop = FALSE]
  p_bench <- numeric(nrow(d))
  nb <- ceiling(nrow(d) / block_size)
  for(k in seq_len(nb)) {
    ids <- ((k - 1L) * block_size + 1L):min(k * block_size, nrow(d))
    eta <- outer(bb[, "beta0"], rep(1, length(ids))) +
      tcrossprod(bb[, "beta_diff", drop = FALSE], d$rank_diff[ids])
    p_bench[ids] <- colMeans(plogis(eta))
  }
  performance[[2L]] <- get_performance(d$Y, p_bench, "RankDiff-only")
  all_predictions[["RankDiff-only"]] <- p_bench
  rm(benchmark_chains, benchmark_draws, bb)
} else {
  message("Optional benchmark not found: ", benchmark_fit_path,
          ". Player scores and main RW2 performance will still be exported.")
}

perf_df <- do.call(rbind, lapply(performance, `[[`, "metrics"))
conf_df <- do.call(rbind, lapply(performance, `[[`, "confusion"))
write.csv(perf_df, file.path(out_dir, "predictive_performance.csv"), row.names = FALSE)
write.csv(conf_df, file.path(out_dir, "confusion_matrices.csv"), row.names = FALSE)
cat("\nIN-SAMPLE PERFORMANCE (favorite win = positive class):\n")
print(perf_df, row.names = FALSE, digits = 5)

# ROC plotting (only if both models available). AUC is computed without
# dependencies; curves use a direct threshold-based construction.
roc_curve <- function(y, p) {
  o <- order(p, decreasing = TRUE)
  pred <- p[o]
  yy <- y[o]
  last_in_tie <- !duplicated(pred, fromLast = TRUE)
  data.frame(
    FPR = c(0, cumsum(yy == 0L)[last_in_tie] / sum(yy == 0L)),
    TPR = c(0, cumsum(yy == 1L)[last_in_tie] / sum(yy == 1L))
  )
}

if(length(all_predictions) == 2L) {
  curve_main <- roc_curve(d$Y, p_favorite)
  curve_base <- roc_curve(d$Y, all_predictions[["RankDiff-only"]])
  draw_roc <- function() {
    plot(curve_main$FPR, curve_main$TPR, type = "l", xlim = c(0, 1), ylim = c(0, 1),
         col = "#2C7FB8", lwd = 2.5,
         xlab = "False positive rate", ylab = "True positive rate",
         main = "In-sample ROC: revised RW2 vs RankDiff-only")
    lines(curve_base, col = "#238B45", lwd = 2, lty = 2)
    abline(0, 1, col = "grey60", lty = 3)
    legend("bottomright", bty = "n", lwd = c(2.5, 2), lty = c(1, 2),
           col = c("#2C7FB8", "#238B45"),
           legend = sprintf("%s (AUC %.3f)", perf_df$Model, perf_df$AUC))
  }
  png(file.path(out_dir, "ROC_main_vs_RankDiff.png"), width = 1500, height = 1200,
      res = 180)
  draw_roc()
  dev.off()
  postscript(file.path(out_dir, "ROC_main_vs_RankDiff.eps"),
             width = 7, height = 5.5, horizontal = FALSE, paper = "special")
  draw_roc()
  dev.off()
}

# -------------------- FAVORITE / UNDERDOG ORIENTATION ------------------------
# Player_1's probability of winning flips when Player_1 is underdog.
p_P1 <- ifelse(d$favorite_is_P1, p_favorite, 1 - p_favorite)
p_P2 <- 1 - p_P1
win_P1 <- as.integer(d$Winner == d$Player_1)
win_P2 <- as.integer(d$Winner == d$Player_2)
stopifnot(all(win_P1 + win_P2 == 1L),
          max(abs(p_P1 + p_P2 - 1)) < 1e-12,
          all(d$Y == ifelse(d$favorite_is_P1, win_P1, win_P2)))

if(save_match_level_predictions) {
  match_out <- data.frame(
    Match_id = seq_len(nrow(d)), Year = d$Year, Date = d$Date,
    Player_1 = d$Player_1, Player_2 = d$Player_2,
    Rank_1 = d$Rank_1, Rank_2 = d$Rank_2,
    Surface = d$Surface, Best.of = d$Best.of, Winner = d$Winner,
    P_favorite = p_favorite,
    P_player_1 = p_P1, P_player_2 = p_P2,
    Y_favorite = d$Y, stringsAsFactors = FALSE
  )
  write.csv(match_out, file.path(out_dir, "match_level_predictions.csv"), row.names = FALSE)
}

# ---------------------- FOUR PLAYER-LEVEL SCORES -----------------------------
# NOTE: All individual scores are averages of posterior mean match probabilities.
# They are DESCRIPTIVE summaries conditional on the already-observed match set,
# not independent estimates of ability, and not out-of-sample evaluation.

pm <- data.frame(
  Player = c(d$Player_1, d$Player_2),
  Won = c(win_P1, win_P2),
  P_win = c(p_P1, p_P2),
  stringsAsFactors = FALSE
)

ix_by_player <- split(seq_len(nrow(pm)), pm$Player)
score_one_player <- function(nm, ids) {
  a <- pm[ids, , drop = FALSE]
  w <- a$Won == 1L
  l <- !w
  nW <- sum(w)
  nL <- sum(l)
  SW <- if(nW) mean(1 - a$P_win[w]) else NA_real_
  SL <- if(nL) mean(a$P_win[l]) else NA_real_
  data.frame(
    Player = nm,
    Total_Matches = nrow(a),
    n_wins = nW,
    n_losses = nL,
    Mean_Residual = mean(a$Won - a$P_win),
    Avg_Surprise_Win = SW,
    Avg_Surprise_Loss = SL,
    Net_Score = if(nW && nL) SW - SL else NA_real_,
    Mean_Rank_Observed = NA_real_, # populated separately below
    stringsAsFactors = FALSE
  )
}

all_scores <- do.call(rbind, Map(score_one_player, names(ix_by_player), ix_by_player))
rownames(all_scores) <- NULL
# Mean observed ranking is a descriptive aid, NOT a model covariate summary.
mean_player_rank <- tapply(c(d$Rank_1, d$Rank_2),
                           c(d$Player_1, d$Player_2), mean)
all_scores$Mean_Rank_Observed <- as.numeric(mean_player_rank[all_scores$Player])
all_scores <- all_scores[order(all_scores$Player), , drop = FALSE]
write.csv(all_scores, file.path(out_dir, "player_scores_all.csv"), row.names = FALSE)

eligible <- all_scores[all_scores$n_wins >= min_wins &
                         all_scores$n_losses >= min_losses, , drop = FALSE]
if(nrow(eligible) < 20L) {
  stop("Only ", nrow(eligible), " players meet the original >=40 wins, >=40 losses rule. ",
       "Cannot form Top/Bottom 10 lists safely.")
}
write.csv(eligible, file.path(out_dir, "player_scores_eligible.csv"), row.names = FALSE)

score_columns <- c("Mean_Residual", "Avg_Surprise_Win",
                   "Avg_Surprise_Loss", "Net_Score")
score_titles <- c("Mean residual score", "Average surprise score",
                  "Average loss score", "Net score")

ordered_extremes <- function(data, metric) {
  # Alphabetical tie-break ensures reproducible order.
  ix_top <- order(-data[[metric]], data$Player)[seq_len(10L)]
  ix_bot <- order( data[[metric]], data$Player)[seq_len(10L)]
  data.frame(
    Position = c(paste0("Top ", 1:10), paste0("Bottom ", 1:10)),
    Player = c(data$Player[ix_top], data$Player[ix_bot]),
    Metric = metric,
    Value = c(data[[metric]][ix_top], data[[metric]][ix_bot]),
    Total_Matches = c(data$Total_Matches[ix_top], data$Total_Matches[ix_bot]),
    n_wins = c(data$n_wins[ix_top], data$n_wins[ix_bot]),
    n_losses = c(data$n_losses[ix_top], data$n_losses[ix_bot]),
    Mean_Rank_Observed = c(data$Mean_Rank_Observed[ix_top],
                           data$Mean_Rank_Observed[ix_bot]),
    stringsAsFactors = FALSE
  )
}

extremes_long <- do.call(rbind, lapply(score_columns, function(m) {
  ordered_extremes(eligible, m)
}))
rownames(extremes_long) <- NULL
write.csv(extremes_long, file.path(out_dir, "player_extremes.csv"), row.names = FALSE)

wide <- data.frame(Position = c(paste0("Top ", 1:10),
                                paste0("Bottom ", 1:10)), stringsAsFactors = FALSE)
format_score <- function(nm, x, metric) {
  decimals <- if(metric == "Mean_Residual") 4L else 3L
  paste0(nm, " (", sprintf(paste0("%.", decimals, "f"), x), ")")
}
for(k in seq_along(score_columns)) {
  sub <- extremes_long[extremes_long$Metric == score_columns[k], ]
  stopifnot(identical(sub$Position, wide$Position))
  wide[[score_titles[k]]] <- format_score(sub$Player, sub$Value, score_columns[k])
}
write.csv(wide, file.path(out_dir, "player_extremes_wide.csv"), row.names = FALSE)

# Explicitly preserve the original 20-row x 4-score layout, with blue data.
tex_escape <- function(x) {
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x <- gsub("#", "\\\\#", x, fixed = TRUE)
  x
}
blue_cell <- function(x) paste0("\\textcolor{blue}{", tex_escape(x), "}")
tex_rows <- vapply(seq_len(nrow(wide)), function(i) {
  paste0(wide$Position[i], " & ",
         paste(vapply(wide[i, -1, drop = TRUE], blue_cell, character(1)),
               collapse = " & "),
         " \\\\")
}, character(1))

tex_lines <- c(
  "% Automatically generated from the REVISED four-chain RW2 model.",
  "% Eligible players: at least 40 wins AND 40 losses in the Top-100 sample.",
  "% Values use posterior mean match probabilities (all model interactions).",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "Position & \\textcolor{blue}{Mean residual score} & Average surprise score & Average loss score & Net score \\\\",
  "\\midrule",
  tex_rows[1:10],
  "\\midrule",
  tex_rows[11:20],
  "\\bottomrule",
  "\\end{tabular}}",
  "\\caption{Top and bottom players for the mean residual score, average surprise score, average loss score, and net score. Players with at least 40 wins and 40 losses in the sample are included.}",
  "\\label{tab:wide_player_metrics}",
  "\\end{table}"
)
writeLines(tex_lines, con = file.path(out_dir, "Table_4_5_player_scores.tex"),
           useBytes = TRUE)

# ----------------------- DIAGNOSTIC REPORT ----------------------------------
# Identity: for each match, both player-oriented residuals add to zero.
residual_balance <- sum(win_P1 - p_P1) + sum(win_P2 - p_P2)
if(abs(residual_balance) > 1e-8) stop("Player residual symmetry check failed.")
report <- c(
  paste0("Script date: ", Sys.Date()),
  paste0("Model: ", main_fit_path),
  paste0("Matched data: ", data_path),
  paste0("Matches: ", nrow(d)),
  paste0("Posterior draws used: ", nrow(beta_draws)),
  paste0("Years: ", min(d$Year), "-", max(d$Year)),
  paste0("Best-of-five matches: ", sum(d$Best_of_5)),
  paste0("Tied rankings (Player_2 favorite by original convention): ", n_tied_rank),
  paste0("Players in selected matches: ", nrow(all_scores)),
  paste0("Players with >=40 wins AND >=40 losses: ", nrow(eligible)),
  paste0("Overall observed favorite win rate: ", sprintf("%.6f", mean(d$Y))),
  paste0("Overall fitted favorite win probability: ", sprintf("%.6f", mean(p_favorite))),
  paste0("Sum of both-sided player residuals (should be zero): ",
         format(residual_balance, digits = 7)),
  "All predictive performance numbers are IN-SAMPLE.",
  "The scores are descriptive; they do not identify intrinsic player ability."
)
writeLines(report, file.path(out_dir, "data_and_checks.txt"))

cat("\nCOMPLETED SECTION 4.5.\n")
cat("Players eligible for the main table:", nrow(eligible), "\n")
cat("Leading/bottom players for the four metrics:\n")
print(wide, row.names = FALSE, right = FALSE)
cat("\nOutputs saved in: ", normalizePath(out_dir, winslash = "/"), "\n")
cat("Main files: player_scores_eligible.csv; player_extremes_wide.csv; ",
    "Table_4_5_player_scores.tex; predictive_performance.csv\n", sep = "")
