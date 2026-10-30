basic_static_full_model <- nimbleCode({

  for(i in 1:n) {

    Y[i] ~ dbern(pi[i])

    logit(pi[i]) <-
      beta0 +

      beta_bestof * Best_of_5[i] +

      (beta_best +
         beta_clay_best * Clay[i] +
         beta_grass_best * Grass[i] +
         beta_bestof_best * Best_of_5[i]) *
         rank_best[i] +

      (beta_diff +
         beta_clay_diff * Clay[i] +
         beta_grass_diff * Grass[i] +
         beta_bestof_diff * Best_of_5[i]) *
         rank_diff[i] +

      (beta_inter +
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
})
