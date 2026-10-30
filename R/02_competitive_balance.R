# =============================================================================
# Competitive-balance summaries from the main RW2 model (no new MCMC fitting)
# =============================================================================
# Run from the repository root, after R/01_fit_main_RW2.R.
# Input : Data/atp_tennis.csv; Models/Basic_random_effects_model_revised_ncp.rda
# Output: Figures/Metric_1_Predictability_Index.eps,
#         Figures/Metric_2_Rank_k_Expected.eps,
#         Results/competitive_balance_main_RW2.csv
# All hypothetical ranking pairs are equally weighted. Match conditions are
# fixed at hard court and best-of-three, as in the manuscript.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

data_file <- "Data/atp_tennis.csv"
model_file <- "Models/Basic_random_effects_model_revised_ncp.rda"
if (!file.exists(data_file) || !file.exists(model_file)) {
  stop("Run this script from the repository root after fitting the main RW2 model.")
}
dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("Results", showWarnings = FALSE, recursive = TRUE)

raw <- read.csv(data_file, stringsAsFactors = FALSE)
d <- raw[!is.na(raw$Rank_1) & !is.na(raw$Rank_2) &
           raw$Rank_1 > 0 & raw$Rank_2 > 0 &
           raw$Rank_1 <= 100 & raw$Rank_2 <= 100 &
           !is.na(raw$Surface) & raw$Surface != "Carpet", , drop = FALSE]
if (nrow(d) != 40988L) stop("Expected 40,988 matches; check the Kaggle dataset version.")
d$Date <- as.Date(d$Date, format = "%Y-%m-%d")
years <- sort(unique(as.integer(format(d$Date, "%Y"))))
stopifnot(identical(years, 2000:2025))

best_center <- mean(log(pmin(d$Rank_1, d$Rank_2) + 1))
best_scale <- sd(log(pmin(d$Rank_1, d$Rank_2) + 1))
diff_center <- mean(log(abs(d$Rank_1 - d$Rank_2) + 1))
diff_scale <- sd(log(abs(d$Rank_1 - d$Rank_2) + 1))

model_env <- new.env(parent = baseenv())
load(model_file, envir = model_env)
chains <- model_env$mcmc.output$samples
if (is.matrix(chains)) chains <- list(chains)
if (!is.list(chains) || length(chains) != 4L) stop("Expected four MCMC chains")
chains <- lapply(chains, as.matrix)

# Draws selected equally across four chains, reproducibly.
set.seed(20261013L)
n_each <- 250L
post <- do.call(rbind, lapply(chains, function(x) {
  if (nrow(x) < n_each) stop("Insufficient draws in a chain")
  x[sort(sample.int(nrow(x), n_each)), , drop = FALSE]
}))

pairs <- expand.grid(r1 = 1:100, r2 = 1:100)
pairs <- pairs[pairs$r1 < pairs$r2, , drop = FALSE]
stopifnot(nrow(pairs) == 4950L)
rb <- (log(pairs$r1 + 1) - best_center) / best_scale
rd <- (log(pairs$r2 - pairs$r1 + 1) - diff_center) / diff_scale

regions <- list(
  "Top 100 pairs" = seq_len(nrow(pairs)),
  "Top 10 pairs" = which(pairs$r1 <= 10 & pairs$r2 <= 10),
  "Rank 1" = which(pairs$r1 == 1),
  "Rank 10" = which(pairs$r1 == 10),
  "Rank 20" = which(pairs$r1 == 20),
  "Rank 50" = which(pairs$r1 == 50)
)
stopifnot(length(regions[["Top 10 pairs"]]) == 45L)

summaries <- vector("list", length(years) * length(regions))
i <- 0L
for (year in years) {
  tt <- match(year, years)
  b <- function(fixed, temporal) {
    post[, fixed] + post[, paste0(temporal, "[", tt, "]")]
  }
  eta <- outer(rep(1, nrow(pairs)), b("beta0", "beta0_year")) +
    outer(rb, b("beta_best", "beta_year_best")) +
    outer(rd, b("beta_diff", "beta_year_diff")) +
    outer(rb * rd, b("beta_inter", "beta_year_inter"))
  probability <- plogis(eta)
  for (name in names(regions)) {
    i <- i + 1L
    draw_mean <- colMeans(probability[regions[[name]], , drop = FALSE])
    summaries[[i]] <- data.frame(
      Year = year,
      Group = name,
      Mean = mean(draw_mean),
      Lower = unname(quantile(draw_mean, 0.025)),
      Upper = unname(quantile(draw_mean, 0.975))
    )
  }
  message("Processed year: ", year)
}
results <- do.call(rbind, summaries)
write.csv(results, "Results/competitive_balance_main_RW2.csv", row.names = FALSE)

plot_groups <- function(groups, title, filename, width = 9, height = 6) {
  x <- results[results$Group %in% groups, , drop = FALSE]
  x$Group <- factor(x$Group, levels = groups)
  p <- ggplot(x, aes(x = Year, y = Mean, colour = Group, fill = Group)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.18,
                colour = NA, show.legend = FALSE) +
    geom_line(linewidth = 0.9) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = "Year", y = "Favorite's winning probability", title = title) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  # EPS is used by the manuscript project.
  postscript(filename, family = "Helvetica", width = width, height = height)
  print(p)
  dev.off()
}

plot_groups(
  c("Top 100 pairs", "Top 10 pairs"),
  "Average winning probability for selected groups",
  "Figures/Metric_1_Predictability_Index.eps"
)
plot_groups(
  c("Rank 1", "Rank 10", "Rank 20", "Rank 50"),
  "Rank-k winning probability for selected ranks",
  "Figures/Metric_2_Rank_k_Expected.eps"
)
message("Saved figures and Results/competitive_balance_main_RW2.csv")
