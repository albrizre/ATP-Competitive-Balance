rankdiff_only_model <- nimbleCode({

  for(i in 1:n) {
    Y[i] ~ dbern(pi[i])

    logit(pi[i]) <-
      beta0 +
      beta_diff * rank_diff[i]
  }

  beta0 ~ dnorm(0, sd = beta_prior_sd)
  beta_diff ~ dnorm(0, sd = beta_prior_sd)
})
