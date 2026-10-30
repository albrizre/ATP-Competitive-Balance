# =============================================================================
# REGENERATE FIGURES 2, 3, 6 (+ clay heatmap) WITH POSTERIOR-UNCERTAINTY DOTS
# =============================================================================
# Run from the Article/ folder.
# Required:
#   Data/atp_tennis.csv
#   Models/Basic_random_effects_model_revised_ncp.rda
# Outputs (same names used in main.tex):
#   Figures/Heatmap_2_Format.eps / .pdf
#   Figures/Heatmap_1_Surface.eps / .pdf
#   Figures/Heatmap_3_Era.eps / .pdf
# Additional:
#   Figures/Heatmap_4_Clay.eps / .pdf
#   Models/Heatmap_uncertainty/heatmap_difference_summaries.csv
#   Models/Heatmap_uncertainty/heatmap_difference_draws.rds
#
# Each heatmap still displays the posterior mean difference, but now a black dot
# is added to cells whose 95% posterior credible interval excludes zero.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(tidyr)
})

model_file <- "Models/Basic_random_effects_model_revised_ncp.rda"
data_file  <- "Data/atp_tennis.csv"
figure_dir <- "Figures"
out_dir    <- "Models/Heatmap_uncertainty"

# Balanced posterior subsample to keep the script fast and reproducible.
n_draws <- 1200L
seed <- 20261027L

for (path in c(figure_dir, out_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(model_file) || !file.exists(data_file)) {
  stop("Run from Article/. Missing: ",
       paste(c(model_file, data_file)[!file.exists(c(model_file, data_file))], collapse = ", "))
}

cat("Loading and filtering ATP data...\n")
raw <- read.csv(data_file, stringsAsFactors = FALSE)
raw$Date <- as.Date(raw$Date, format = "%Y-%m-%d")

d <- raw %>%
  filter(!is.na(Rank_1), !is.na(Rank_2),
         Rank_1 > 0, Rank_2 > 0,
         Rank_1 <= 100, Rank_2 <= 100,
         !is.na(Surface), Surface != "Carpet") %>%
  mutate(
    rank_best_orig = pmin(Rank_1, Rank_2),
    rank_diff_orig = abs(Rank_1 - Rank_2),
    Year = as.integer(format(Date, "%Y"))
  )

stopifnot(nrow(d) == 40988L)
base_year <- min(d$Year)
years <- sort(unique(d$Year))

mean_log_best <- mean(log(d$rank_best_orig + 1))
sd_log_best   <- sd(log(d$rank_best_orig + 1))
mean_log_diff <- mean(log(d$rank_diff_orig + 1))
sd_log_diff   <- sd(log(d$rank_diff_orig + 1))

cat("Loading main RW2 model...\n")
env <- new.env(parent = baseenv())
load(model_file, envir = env)
if (!exists("mcmc.output", envir = env, inherits = FALSE)) {
  stop("The file does not contain 'mcmc.output'.")
}
chains <- env$mcmc.output$samples
if (is.matrix(chains)) chains <- list(chains)
if (!is.list(chains) || length(chains) != 4L) {
  stop("Expected four chains in mcmc.output$samples.")
}
chains <- lapply(chains, as.matrix)
if (!all(vapply(chains, function(x) identical(colnames(x), colnames(chains[[1]])), logical(1)))) {
  stop("Posterior parameter names differ between chains.")
}

set.seed(seed)
n_each <- as.integer(n_draws / length(chains))
stopifnot(n_each * length(chains) == n_draws)
post <- do.call(rbind, lapply(chains, function(ch) {
  if (nrow(ch) < n_each) stop("Chain has fewer samples than requested.")
  ch[sort(sample.int(nrow(ch), n_each)), , drop = FALSE]
}))

fixed_names <- c(
  "beta0", "beta_best", "beta_diff", "beta_inter",
  "beta_clay", "beta_grass", "beta_bestof",
  "beta_clay_best", "beta_grass_best", "beta_bestof_best",
  "beta_clay_diff", "beta_grass_diff", "beta_bestof_diff",
  "beta_clay_inter", "beta_grass_inter", "beta_bestof_inter"
)
temporal_names <- unlist(lapply(
  c("beta0_year", "beta_year_best", "beta_year_diff", "beta_year_inter"),
  function(nm) paste0(nm, "[", seq_along(years), "]")
))
needed <- c(fixed_names, temporal_names)
missing_cols <- setdiff(needed, colnames(post))
if (length(missing_cols)) stop("Missing posterior columns: ", paste(missing_cols, collapse = ", "))

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
make_grid <- function(rank_best_seq = seq(5, 50, by = 5),
                      rank_diff_seq = seq(5, 25, by = 5),
                      year = 2025L,
                      surface = "Hard",
                      bo5 = 0L) {
  expand.grid(
    rank_best_orig = rank_best_seq,
    rank_diff_orig = rank_diff_seq,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      rb_std = (log(rank_best_orig + 1) - mean_log_best) / sd_log_best,
      rd_std = (log(rank_diff_orig + 1) - mean_log_diff) / sd_log_diff,
      is_clay = as.integer(surface == "Clay"),
      is_grass = as.integer(surface == "Grass"),
      Best_of_5 = as.integer(bo5),
      Year = as.integer(year),
      Year_num = match(Year, years)
    )
}

eta_from_grid <- function(grid) {
  rb  <- grid$rb_std
  rd  <- grid$rd_std
  # Every grid in this script represents ONE calendar year. Using the entire
  # 50-element Year_num vector to index posterior columns would create
  # matrices rather than one posterior vector per temporal coefficient.
  yr_values <- unique(as.integer(grid$Year_num))
  if (length(yr_values) != 1L || is.na(yr_values)) {
    stop("eta_from_grid(): every grid must contain exactly one valid year.")
  }
  yr <- yr_values[[1L]]
  icl <- grid$is_clay
  igr <- grid$is_grass
  bo5 <- grid$Best_of_5
  inter <- rb * rd

  beta0      <- post[, "beta0"]      + post[, paste0("beta0_year[", yr, "]")]
  beta_best  <- post[, "beta_best"]  + post[, paste0("beta_year_best[", yr, "]")]
  beta_diff  <- post[, "beta_diff"]  + post[, paste0("beta_year_diff[", yr, "]")]
  beta_inter <- post[, "beta_inter"] + post[, paste0("beta_year_inter[", yr, "]")]

  eta <- outer(rep(1, nrow(grid)), beta0) +
    outer(rb, beta_best) +
    outer(rd, beta_diff) +
    outer(inter, beta_inter) +
    outer(icl, post[, "beta_clay"]) +
    outer(igr, post[, "beta_grass"]) +
    outer(bo5, post[, "beta_bestof"]) +
    outer(icl * rb, post[, "beta_clay_best"]) +
    outer(igr * rb, post[, "beta_grass_best"]) +
    outer(bo5 * rb, post[, "beta_bestof_best"]) +
    outer(icl * rd, post[, "beta_clay_diff"]) +
    outer(igr * rd, post[, "beta_grass_diff"]) +
    outer(bo5 * rd, post[, "beta_bestof_diff"]) +
    outer(icl * inter, post[, "beta_clay_inter"]) +
    outer(igr * inter, post[, "beta_grass_inter"]) +
    outer(bo5 * inter, post[, "beta_bestof_inter"])

  # Defensive shape check to catch posterior-indexing mistakes immediately.
  if (!is.matrix(eta) ||
      !identical(dim(eta), c(nrow(grid), nrow(post)))) {
    stop("Unexpected eta dimensions: ", paste(dim(eta), collapse = " x "),
         "; expected ", nrow(grid), " x ", nrow(post))
  }
  plogis(eta)
}

summarise_delta <- function(delta_draws, grid, label) {
  out <- grid[, c("rank_best_orig", "rank_diff_orig")]
  out$Figure <- label
  out$Delta <- rowMeans(delta_draws)
  qs <- t(apply(delta_draws, 1, quantile, probs = c(0.025, 0.975)))
  out$Lower <- qs[, 1]
  out$Upper <- qs[, 2]
  out$Dot <- (out$Lower > 0) | (out$Upper < 0)
  out
}

plot_heatmap <- function(df, title, filename_base, width = 10, height = 7) {
  # A symmetric colour scale ensures that zero is exactly white; not simply
  # the numerical midpoint between the observed minimum and maximum.
  limit <- max(abs(df$Delta), na.rm = TRUE)
  if (!is.finite(limit) || limit == 0) limit <- 1e-8

  p <- ggplot(df, aes(x = rank_diff_orig, y = rank_best_orig, fill = Delta)) +
    geom_tile(color = NA) +
    geom_point(data = df %>% filter(Dot),
               aes(x = rank_diff_orig, y = rank_best_orig),
               inherit.aes = FALSE,
               shape = 16, size = 3.2, colour = "black") +
    scale_y_reverse(breaks = seq(5, 50, by = 5), expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(5, 25, by = 5), expand = c(0, 0)) +
    scale_fill_gradient2(
      low = "#8c510a", mid = "#FFFFFF", high = "#01665e", midpoint = 0,
      limits = c(-limit, limit),
      labels = percent_format(accuracy = 1), name = "Difference"
    ) +
    labs(
      title = title,
      x = "Ranking difference", y = "Favorite's rank"
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 20),
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

  pdf_file <- file.path(figure_dir, paste0(filename_base, ".pdf"))
  eps_file <- file.path(figure_dir, paste0(filename_base, ".eps"))
  ggsave(pdf_file, p, device = cairo_pdf, width = width, height = height, units = "in")
  postscript(eps_file, family = "Helvetica", width = width, height = height)
  print(p)
  dev.off()
  cat("Saved:", pdf_file, "and", eps_file, "\n")
}

# -----------------------------------------------------------------------------
# Surface: Grass vs Hard
# -----------------------------------------------------------------------------
cat("Computing surface contrast: Grass vs Hard...\n")
grid_hard_ref  <- make_grid(year = 2025L, surface = "Hard",  bo5 = 0L)
grid_grass_ref <- make_grid(year = 2025L, surface = "Grass", bo5 = 0L)
grid_clay_ref  <- make_grid(year = 2025L, surface = "Clay",  bo5 = 0L)

p_hard  <- eta_from_grid(grid_hard_ref)
p_grass <- eta_from_grid(grid_grass_ref)
p_clay  <- eta_from_grid(grid_clay_ref)

delta_surface <- p_grass - p_hard
delta_clay    <- p_clay  - p_hard
sum_surface   <- summarise_delta(delta_surface, grid_hard_ref, "Grass vs Hard")
sum_clay      <- summarise_delta(delta_clay,    grid_hard_ref, "Clay vs Hard")

plot_heatmap(sum_surface, "Surface effect: Grass vs. Hard", "Heatmap_1_Surface")
plot_heatmap(sum_clay,    "Surface effect: Clay vs. Hard",  "Heatmap_4_Clay")

# -----------------------------------------------------------------------------
# Format: Best-of-5 vs Best-of-3
# -----------------------------------------------------------------------------
cat("Computing format contrast: Bo5 vs Bo3...\n")
grid_bo3 <- make_grid(year = 2025L, surface = "Hard", bo5 = 0L)
grid_bo5 <- make_grid(year = 2025L, surface = "Hard", bo5 = 1L)

p_bo3 <- eta_from_grid(grid_bo3)
p_bo5 <- eta_from_grid(grid_bo5)

delta_format <- p_bo5 - p_bo3
sum_format   <- summarise_delta(delta_format, grid_bo3, "Bo5 vs Bo3")

plot_heatmap(sum_format, "Match format effect: 5 sets vs. 3 sets", "Heatmap_2_Format")

# -----------------------------------------------------------------------------
# Era: 2015 vs 2000
# -----------------------------------------------------------------------------
cat("Computing era contrast: 2015 vs 2000...\n")
grid_2000 <- make_grid(year = 2000L, surface = "Hard", bo5 = 0L)
grid_2015 <- make_grid(year = 2015L, surface = "Hard", bo5 = 0L)

p_2000 <- eta_from_grid(grid_2000)
p_2015 <- eta_from_grid(grid_2015)

delta_era <- p_2015 - p_2000
sum_era   <- summarise_delta(delta_era, grid_2000, "2015 vs 2000")

plot_heatmap(sum_era, "Era effect: 2015 vs. 2000", "Heatmap_3_Era")

# Save numerical summaries too.
all_summaries <- bind_rows(sum_surface, sum_clay, sum_format, sum_era)
write.csv(all_summaries,
          file.path(out_dir, "heatmap_difference_summaries.csv"),
          row.names = FALSE)
saveRDS(list(
  surface = delta_surface,
  clay = delta_clay,
  format = delta_format,
  era = delta_era
), file.path(out_dir, "heatmap_difference_draws.rds"))

cat("\nDone. Black dots mark cells whose 95% credible interval excludes zero.\n")
cat("Summary CSV:", file.path(out_dir, "heatmap_difference_summaries.csv"), "\n")
