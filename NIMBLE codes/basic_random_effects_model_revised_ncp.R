basic_random_effects_model_ncp <- nimbleCode({

  # ---------------------------------------------------------------------------
  # Likelihood
  # ---------------------------------------------------------------------------
  for(i in 1:n) {
    Y[i] ~ dbern(pi[i])

    logit(pi[i]) <-
      beta0 + beta0_year[Year_num[i]] +

      # Main effect for match format
      beta_bestof * Best_of_5[i] +

      # Best Rank
      (beta_best +
         beta_year_best[Year_num[i]] +
         beta_clay_best * Clay[i] +
         beta_grass_best * Grass[i] +
         beta_bestof_best * Best_of_5[i]) * rank_best[i] +

      # Rank Difference
      (beta_diff +
         beta_year_diff[Year_num[i]] +
         beta_clay_diff * Clay[i] +
         beta_grass_diff * Grass[i] +
         beta_bestof_diff * Best_of_5[i]) * rank_diff[i] +

      # Best Rank x Rank Difference
      (beta_inter +
         beta_year_inter[Year_num[i]] +
         beta_clay_inter * Clay[i] +
         beta_grass_inter * Grass[i] +
         beta_bestof_inter * Best_of_5[i]) *
         rank_best[i] * rank_diff[i] +

      # Main effects for surface
      beta_clay * Clay[i] +
      beta_grass * Grass[i]
  }

  # ---------------------------------------------------------------------------
  # Static regression coefficients
  # ---------------------------------------------------------------------------
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

  # ---------------------------------------------------------------------------
  # Non-centred / whitened RW2 representation
  #
  # Q = D2' D2 has rank n_years - 2.  B_rw2 is built in R so that
  # B_rw2' Q B_rw2 = I and its columns are orthogonal to the constant and
  # linear null-space directions.
  #
  # Each time-varying coefficient is represented as
  #
  #   beta_k(t) = slope_k * time_scaled(t)
  #               + sigma_k * B_rw2(t,.) z_k,
  #
  # with z_k ~ N(0, I).
  # ---------------------------------------------------------------------------

  # Explicit linear component of each RW2 null space
  slope_beta0 ~ dnorm(0, sd = time_slope_prior_sd)
  slope_beta_best ~ dnorm(0, sd = time_slope_prior_sd)
  slope_beta_diff ~ dnorm(0, sd = time_slope_prior_sd)
  slope_beta_inter ~ dnorm(0, sd = time_slope_prior_sd)

  # Whitened curvature coefficients
  for(k in 1:n_rw2) {
    z_beta0[k] ~ dnorm(0, 1)
    z_beta_best[k] ~ dnorm(0, 1)
    z_beta_diff[k] ~ dnorm(0, 1)
    z_beta_inter[k] ~ dnorm(0, 1)
  }

  # RW2 scale parameters
  sigma.beta0_year ~ T(dnorm(0, sd = sigma_prior_sd), 0, Inf)
  sigma.beta_year_best ~ T(dnorm(0, sd = sigma_prior_sd), 0, Inf)
  sigma.beta_year_diff ~ T(dnorm(0, sd = sigma_prior_sd), 0, Inf)
  sigma.beta_year_inter ~ T(dnorm(0, sd = sigma_prior_sd), 0, Inf)

  # Deterministic temporal trajectories
  for(t in 1:n_years) {

    rw_beta0[t] <-
      inprod(B_rw2[t, 1:n_rw2], z_beta0[1:n_rw2])

    rw_beta_best[t] <-
      inprod(B_rw2[t, 1:n_rw2], z_beta_best[1:n_rw2])

    rw_beta_diff[t] <-
      inprod(B_rw2[t, 1:n_rw2], z_beta_diff[1:n_rw2])

    rw_beta_inter[t] <-
      inprod(B_rw2[t, 1:n_rw2], z_beta_inter[1:n_rw2])

    beta0_year[t] <-
      slope_beta0 * time_scaled[t] +
      sigma.beta0_year * rw_beta0[t]

    beta_year_best[t] <-
      slope_beta_best * time_scaled[t] +
      sigma.beta_year_best * rw_beta_best[t]

    beta_year_diff[t] <-
      slope_beta_diff * time_scaled[t] +
      sigma.beta_year_diff * rw_beta_diff[t]

    beta_year_inter[t] <-
      slope_beta_inter * time_scaled[t] +
      sigma.beta_year_inter * rw_beta_inter[t]
  }
})
