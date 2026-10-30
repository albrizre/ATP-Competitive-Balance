basic_linear_time_model <- nimbleCode({

  for(i in 1:n) {

    Y[i] ~ dbern(pi[i])

    logit(pi[i]) <-
      beta0 +
      slope_beta0 * time_scaled[Year_num[i]] +

      beta_bestof * Best_of_5[i] +

      (beta_best +
         slope_beta_best * time_scaled[Year_num[i]] +
         beta_clay_best * Clay[i] +
         beta_grass_best * Grass[i] +
         beta_bestof_best * Best_of_5[i]) *
         rank_best[i] +

      (beta_diff +
         slope_beta_diff * time_scaled[Year_num[i]] +
         beta_clay_diff * Clay[i] +
         beta_grass_diff * Grass[i] +
         beta_bestof_diff * Best_of_5[i]) *
         rank_diff[i] +

      (beta_inter +
         slope_beta_inter * time_scaled[Year_num[i]] +
         beta_clay_inter * Clay[i] +
         beta_grass_inter * Grass[i] +
         beta_bestof_inter * Best_of_5[i]) *
         rank_best[i] * rank_diff[i] +

      beta_clay * Clay[i] +
      beta_grass * Grass[i]
  }

  beta0 ~ dnorm(0, sd = beta_prior_sd)
  beta_best ~ dnorm(0, sd = beta_prior_sd)
  beta_diff ~ dnorm(0, sd = beta_prior_sd)
  beta_inter ~ dnorm(0, sd = beta_prior_sd)

  beta_clay ~ dnorm(0, sd = beta_prior_sd)
  beta_grass ~ dnorm(0, sd = beta_prior_sd)
  beta_bestof ~ dnorm(0, sd = beta_prior_sd)

  beta_clay_best ~ dnorm(0, sd = beta_prior_sd)
  beta_grass_best ~ dnorm(0, sd = beta_prior_sd)
  beta_bestof_best ~ dnorm(0, sd = beta_prior_sd)

  beta_clay_diff ~ dnorm(0, sd = beta_prior_sd)
  beta_grass_diff ~ dnorm(0, sd = beta_prior_sd)
  beta_bestof_diff ~ dnorm(0, sd = beta_prior_sd)

  beta_clay_inter ~ dnorm(0, sd = beta_prior_sd)
  beta_grass_inter ~ dnorm(0, sd = beta_prior_sd)
  beta_bestof_inter ~ dnorm(0, sd = beta_prior_sd)

  slope_beta0 ~ dnorm(0, sd = time_slope_prior_sd)
  slope_beta_best ~ dnorm(0, sd = time_slope_prior_sd)
  slope_beta_diff ~ dnorm(0, sd = time_slope_prior_sd)
  slope_beta_inter ~ dnorm(0, sd = time_slope_prior_sd)

  for(t in 1:n_years) {
    beta0_year[t] <- slope_beta0 * time_scaled[t]
    beta_year_best[t] <- slope_beta_best * time_scaled[t]
    beta_year_diff[t] <- slope_beta_diff * time_scaled[t]
    beta_year_inter[t] <- slope_beta_inter * time_scaled[t]
  }
})
