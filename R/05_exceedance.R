# =============================================================================
# FIGURE 9 + SUPPLEMENTARY EXCEEDANCE FIGURES -- MAIN RW2 (NO REFITTING)
# =============================================================================
# Run this script from the Article/ directory, which contains Data/ and Models/.
# Required: Data/atp_tennis.csv
#           Models/Basic_random_effects_model_revised_ncp.rda
# Produces:
#   Figures/Metric_Exceedance_Prob_2000_2015.pdf        [MAIN Figure 9]
#   Figures/Metric_Exceedance_Prob_2000_2015.eps        [optional EPS]
#   Supplementary_Overleaf/Figure_S5_Exceedance_065_full.pdf
#   Supplementary_Overleaf/Figure_S6_Exceedance_thresholds.pdf
#   Models/Exceedance_Figure9/exceedance_posterior_summary.csv
# No NIMBLE fitting. Every posterior probability is obtained from main-RW2
# posterior samples under hard court and best-of-three reference conditions.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

model_file <- "Models/Basic_random_effects_model_revised_ncp.rda"
data_file <- "Data/atp_tennis.csv"
main_figure_dir <- "Figures"
supplement_figure_dir <- "Supplementary_Overleaf"
results_dir <- "Models/Exceedance_Figure9"

# 4 chains x 300 draws; balanced subsampling with fixed seed.
n_draws <- 1200L
seed <- 20261026L
ranking_grid <- seq(5, 50, by = 5)
main_years <- c(2000L, 2015L)
full_years <- c(2000L, 2008L, 2016L, 2024L)
thresholds <- c(0.60, 0.65, 0.70)

if (!file.exists(data_file) || !file.exists(model_file)) {
  stop("Work from your Article/ folder. Missing: ",
       paste(c(data_file, model_file)[!file.exists(c(data_file, model_file))],
             collapse = ", "))
}
for (path in c(main_figure_dir, supplement_figure_dir, results_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Use precisely the sample and transformations of the original main fit.
raw <- read.csv(data_file, stringsAsFactors = FALSE)
raw$Date <- as.Date(raw$Date, format = "%Y-%m-%d")
keep <- !is.na(raw$Rank_1) & !is.na(raw$Rank_2) &
  raw$Rank_1 > 0 & raw$Rank_2 > 0 &
  raw$Rank_1 <= 100 & raw$Rank_2 <= 100 &
  !is.na(raw$Surface) & raw$Surface != "Carpet"
d <- raw[keep, , drop = FALSE]
stopifnot(nrow(d) == 40988L, !anyNA(d$Date))
years <- sort(unique(as.integer(format(d$Date, "%Y"))))
stopifnot(identical(years, 2000:2025))

best_ref <- log(pmin(d$Rank_1, d$Rank_2) + 1)
diff_ref <- log(abs(d$Rank_1 - d$Rank_2) + 1)
best_center <- mean(best_ref); best_sd <- sd(best_ref)
diff_center <- mean(diff_ref); diff_sd <- sd(diff_ref)

# Read four chains from a separate environment to avoid overwriting objects.
env <- new.env(parent = baseenv())
loaded_names <- load(model_file, envir = env)
if (!"mcmc.output" %in% loaded_names) stop("Model lacks 'mcmc.output'")
chains <- env$mcmc.output$samples
if (is.matrix(chains)) chains <- list(chains)
if (!is.list(chains) || length(chains) != 4L) {
  stop("Expected four MCMC chains, found: ", length(chains))
}
chains <- lapply(chains, as.matrix)
if (!all(vapply(chains, function(x) identical(colnames(x), colnames(chains[[1]])), logical(1)))) {
  stop("Posterior parameter names differ between chains")
}
set.seed(seed)
n_each <- as.integer(n_draws / length(chains))
stopifnot(n_each * length(chains) == n_draws)
post <- do.call(rbind, lapply(chains, function(ch) {
  if (nrow(ch) < n_each) stop("Chain has fewer samples than required")
  ch[sort(sample.int(nrow(ch), n_each)), , drop = FALSE]
}))

necessary <- c("beta0", "beta_best", "beta_diff", "beta_inter",
  unlist(lapply(c("beta0_year", "beta_year_best", "beta_year_diff", "beta_year_inter"),
                function(nm) paste0(nm, "[", seq_along(years), "]"))))
if (!all(necessary %in% colnames(post))) {
  stop("Missing posterior quantities: ",
       paste(setdiff(necessary, colnames(post)), collapse = ", "))
}

# Same 45 possible pairs among the ten grid positions 5,10,...,50 as Fig 9.
grid <- expand.grid(Fav_Rank = ranking_grid, Und_Rank = ranking_grid)
grid <- grid[grid$Fav_Rank < grid$Und_Rank, , drop = FALSE]
stopifnot(nrow(grid) == 45L)
rb <- (log(grid$Fav_Rank + 1) - best_center) / best_sd
rd <- (log(grid$Und_Rank - grid$Fav_Rank + 1) - diff_center) / diff_sd
int <- rb * rd

# Main and supplementary scenarios are hard-court / best-of-three. All 12
# surface- and format-related columns vanish in the *reference* linear predictor.
all_years <- sort(unique(c(main_years, full_years)))
results <- vector("list", length(all_years) * length(thresholds))
idx <- 0L
for (year in all_years) {
  tt <- match(year, years)
  get_beta <- function(base, temporal) {
    post[, base] + post[, paste0(temporal, "[", tt, "]")]
  }
  eta <- outer(rep(1, nrow(grid)), get_beta("beta0", "beta0_year")) +
    outer(rb, get_beta("beta_best", "beta_year_best")) +
    outer(rd, get_beta("beta_diff", "beta_year_diff")) +
    outer(int, get_beta("beta_inter", "beta_year_inter"))
  for (threshold in thresholds) {
    idx <- idx + 1L
    results[[idx]] <- data.frame(
      grid, Year = year, Threshold = threshold,
      exceedance_prob = rowMeans(eta > qlogis(threshold))
    )
  }
}
result <- do.call(rbind, results)
write.csv(result, file.path(results_dir, "exceedance_posterior_summary.csv"),
          row.names = FALSE)

# Fixed 0--1 sequential colour scale across ALL panels and thresholds.
prob_colors <- c("#F4F9FA", "#D3EAE7", "#92CAC5", "#4E9DBB", "#244A78")
base_plot <- function(x) {
  ggplot(x, aes(x = Und_Rank, y = Fav_Rank, fill = exceedance_prob)) +
    geom_tile(color = "white", linewidth = 0.35) +
    scale_fill_gradientn(colours = prob_colors, limits = c(0, 1),
                         oob = scales::squish,
                         breaks = c(0, .25, .5, .75, 1),
                         labels = scales::percent_format(accuracy = 1),
                         name = "Exceedance\nprobability") +
    scale_x_continuous(breaks = ranking_grid, expand = c(0, 0)) +
    scale_y_reverse(breaks = ranking_grid, expand = c(0, 0)) +
    coord_fixed() +
    labs(x = "Lower-ranked player's rank", y = "Favorite's rank") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"),
          legend.position = "right", plot.margin = margin(5, 8, 5, 5))
}

# MAIN Figure 9: two large panels, no numbers in cells.
d_main <- subset(result, Year %in% main_years & abs(Threshold - .65) < 1e-12)
p_main <- base_plot(d_main) + facet_wrap(~ Year, nrow = 1) +
  ggtitle("Exceedance probabilities: 2000 and 2015")
main_pdf <- file.path(main_figure_dir, "Metric_Exceedance_Prob_2000_2015.pdf")
ggsave(main_pdf, p_main, device = "pdf", width = 11, height = 5.2,
       units = "in", useDingbats = FALSE)

# Optional EPS for projects that still use EPS graphics.
main_eps <- file.path(main_figure_dir, "Metric_Exceedance_Prob_2000_2015.eps")
tryCatch({
  postscript(main_eps, width = 11, height = 5.2, paper = "special",
             horizontal = FALSE, onefile = FALSE, family = "Helvetica")
  print(p_main)
  dev.off()
}, error = function(e) {
  message("Could not save optional EPS (PDF was saved): ", conditionMessage(e))
  if (dev.cur() != 1L) dev.off()
})

# SUPPLEMENT S5: the original four-year display, but no overlapping numbers.
d_full <- subset(result, Year %in% full_years & abs(Threshold - .65) < 1e-12)
p_full <- base_plot(d_full) + facet_wrap(~ Year, nrow = 2, ncol = 2) +
  ggtitle("Exceedance probabilities across four selected years")
ggsave(file.path(supplement_figure_dir, "Figure_S5_Exceedance_065_full.pdf"),
       p_full, device = "pdf", width = 9, height = 8.5,
       units = "in", useDingbats = FALSE)

# SUPPLEMENT S6: threshold check, holding years fixed at 2000 and 2015.
d_threshold <- subset(result, Year %in% main_years &
                                  (abs(Threshold - .60) < 1e-12 |
                                   abs(Threshold - .70) < 1e-12))
d_threshold$Threshold_label <- factor(d_threshold$Threshold,
                                     levels = c(.60, .70),
                                     labels = c("Threshold 0.60", "Threshold 0.70"))
p_threshold <- base_plot(d_threshold) +
  facet_grid(Threshold_label ~ Year) +
  ggtitle("Sensitivity to exceedance thresholds")
ggsave(file.path(supplement_figure_dir, "Figure_S6_Exceedance_thresholds.pdf"),
       p_threshold, device = "pdf", width = 9, height = 8.5,
       units = "in", useDingbats = FALSE)

cat("\nDONE. No new model was fitted.\n")
cat("Main figure:", main_pdf, "\n")
cat("Supplement figures: Figure_S5 and Figure_S6 in",
    supplement_figure_dir, "\n")
cat("All cell probabilities:",
    file.path(results_dir, "exceedance_posterior_summary.csv"), "\n")
