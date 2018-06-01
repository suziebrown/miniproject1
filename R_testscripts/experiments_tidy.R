## EXPERIMENTS (new tidier source!)

rm(list = ls())

## Required packages ---------------------

library(foreach)
library(doParallel)
library(future)
library(doRNG)


## Function dependencies -----------------

standard_SMC <- function(n_particles, observations, emission_density, transition_sam, initial_sam, ...){
  n_obs <- length(observations)
  ancestors <- matrix(NA, nrow = n_obs - 1, ncol = n_particles)
  positions <- matrix(NA, nrow = n_obs, ncol = n_particles)
  weights <- matrix(NA, nrow = n_obs, ncol = n_particles)

  ## Initial states
  positions[1, ] <- initial_sam(n_particles, ...)
  weights_tmp <- emission_density(observations[1], positions[1, ], ...)
  weights[1, ] <- weights_tmp / sum(weights_tmp)

  for (t in 2:(n_obs)){
    ## Resample
    ancestors[t-1, ] <- sample(1:n_particles, n_particles, prob=weights_tmp, replace=TRUE)
    positions_tmp <- positions[t-1, ancestors[t-1, ]]

    ## Propagate
    positions[t, ] <- transition_sam(positions_tmp, ...)

    ## Calculate importance weights
    weights_tmp <- emission_density(observations[t], positions[t, ], ...)
    weights[t, ] <- weights_tmp / sum(weights_tmp)
  }

  list(positions = positions, weights = weights, ancestors = ancestors)
}

standard_SMC_anc <- function(n_particles, observations, emission_density, transition_sam, initial_sam, ...){
  n_obs <- length(observations)
  ancestors <- matrix(NA, nrow = n_obs - 1, ncol = n_particles)
  positions <- numeric(n_particles)

  ## Initial states
  positions <- initial_sam(n_particles, ...)
  weights_tmp <- emission_density(observations[1], positions, ...)

  for (t in 2:(n_obs)){
    ## Resample
    ancestors[t-1, ] <- sample(1:n_particles, n_particles, prob=weights_tmp, replace=TRUE)
    positions_tmp <- positions[ancestors[t-1, ]]

    ## Propagate
    positions <- transition_sam(positions_tmp, ...)

    ## Calculate importance weights
    weights_tmp <- emission_density(observations[t], positions, ...)
  }

  ancestors
}

conditional_SMC_anc <- function(n_particles, observations, cond_trajectory, emission_density, transition_sam, initial_sam, ...){
  n_obs <- length(observations)
  ancestors <- matrix(NA, nrow = n_obs - 1, ncol = n_particles)
  positions <- numeric(n_particles)

  ## Initial states
  positions[1] <- cond_trajectory[1]
  positions[2:n_particles] <- initial_sam(n_particles - 1, ...)
  weights_tmp <- emission_density(observations[1], positions, ...)

  for (t in 2:(n_obs)){
    ## Resample
    ancestors[t-1, ] <- sample(1:n_particles, n_particles, prob=weights_tmp, replace=TRUE)
    ancestors[t-1, 1] <- 1
    positions_tmp <- positions[ancestors[t-1, ]]

    ## Propagate
    positions[2:n_particles] <- transition_sam(positions_tmp[2:n_particles], ...)
    positions[1] <- cond_trajectory[t]

    ## Calculate importance weights
    weights_tmp <- emission_density(observations[t], positions, ...)
  }

  ancestors
}

ancestrySize_treeht <- function(history, sampl=NULL, maxgen=NULL){
  N <- attr(history, 'N')
  N.gen <- (attr(history, 'Ngen'))
  min.gen <- max(1, N.gen-maxgen)
  if (is.null(sampl)){
    sampl <- 1:N
    n <- N
  }
  else{
    n <- length(sampl)
  }

  Size <- rep(NA, N.gen)
  Size[N.gen] <- n
  ancestry <- sampl
  MRCA <- NA

  for (t in ((N.gen-1):min.gen)){
    ancestry <- history[t,unique(ancestry)] ## find ancestors of current sample
    Size[t] <- length(unique(ancestry))
    if (is.na(MRCA)){ ## once MRCA is reached, can exit loop & set all previous generations to 1
      if (Size[t]==1){
        MRCA <- t
        break
      }
    }
  }
  if (!is.na(MRCA)){ ## set size to 1 for all generations above MRCA
    Size[min.gen:(MRCA-1)] <- rep(1, MRCA-1)
  }

  N.gen-MRCA
}

traceAncestry <- function(history, sampl, maxgen=NULL){
  n <- length(sampl)
  N <- attr(history, 'N')
  N.gen <- (attr(history, 'Ngen'))
  min.gen <- max(1, N.gen-maxgen)

  ancestry <- matrix(NA, nrow=N.gen, ncol=n)
  ancestry[N.gen,] <- sort(sampl)
  MRCA <- NA

  for (t in ((N.gen-1):min.gen)){ ## try eliminating the for loops
    for (i in 1:n){
      ancestry[t,i] <- history[t,ancestry[t+1,i]]
    }
    if (is.na(MRCA)){
      if (length(unique(ancestry[t,]))==1){
        MRCA <- t ## could leave NAs above once MRCA is found?
      }
    }
  }
  class(ancestry) <- 'samplegenealogy'
  ancestry@samplesize <- n
  ancestry@sample <- sampl
  ancestry@MRCA <- as.integer(MRCA)
  ancestry@history <- unclass(history)
  ancestry
}

ou_emission_density <- function(observation, position, delta, sigma) {
  (2 * pi) ^ (-0.5) * sigma ^ (-1) * exp(-(observation - position) ^ 2 / (2 * sigma ^ 2))
}

ou_transition_sam <- function(pos_old, delta, sigma) {
  rnorm(length(pos_old), (1 - delta) * pos_old, delta ^ 0.5)
}

ou_initial_sam <- function(n_particles, delta, sigma) {
  rnorm(n_particles)
}

ou_rts_smooth <- function(observations, delta, sigma, prior_mean = 0, prior_var = 1) {
  n_obs <- length(observations) - 1

  ## initialise variables
  kal_mean <- numeric(n_obs + 1) ## states x(t|t) after time & measurement updates
  kal_var <- numeric(n_obs + 1)
  kal_mean_tmp <- numeric(n_obs + 1) ## 'tmp' = states x(t|t-1) after time update and before measurement update
  kal_var_tmp <- numeric(n_obs + 1)

  ## Kalman filter recursion
  kal_mean_tmp[1] <- prior_mean
  kal_var_tmp[1] <- prior_var
  for (t in 1:(n_obs)) {
    ## measurement update:
    kal_mean[t] <- kal_mean_tmp[t] + kal_var_tmp[t] * (observations[t] - kal_mean_tmp[t]) / (kal_var_tmp[t] + sigma ^ 2)
    kal_var[t] <- kal_var_tmp[t] - kal_var_tmp[t] ^ 2 / (kal_var_tmp[t] + sigma ^ 2)
    ## time update:
    kal_mean_tmp[t + 1] <- (1 - delta) * kal_mean[t]
    kal_var_tmp[t + 1] <- delta + kal_var[t] * (1 - delta) ^ 2
  }
  kal_mean[n_obs + 1] <- kal_mean_tmp[n_obs + 1] + kal_var_tmp[n_obs + 1] * (observations[n_obs + 1] - kal_mean_tmp[n_obs + 1]) / (kal_var_tmp[n_obs + 1] + sigma ^ 2)
  kal_var[n_obs + 1] <- kal_var_tmp[n_obs + 1] - kal_var_tmp[n_obs + 1] ^ 2 / (kal_var_tmp[n_obs + 1] + sigma ^ 2)

  ## initialise variables
  rts_mean <- numeric(n_obs + 1)
  rts_var <- numeric(n_obs + 1)
  pbar <- numeric(n_obs + 1)
  g <- numeric(n_obs)

  ## final state
  rts_mean[n_obs + 1] <- kal_mean[n_obs + 1]
  rts_var[n_obs + 1] <- kal_var[n_obs + 1]

  ## RTS smoother recursion
  for (t in n_obs:1) {
    pbar[t + 1] <- ((1 - delta) ^ 2) * kal_var[t] + delta
    g[t] <- kal_var[t] * (1 - delta) / pbar[t + 1]
    rts_mean[t] <- kal_mean[t] + g[t] * (rts_mean[t + 1] - (1 - delta) * kal_mean[t])
    rts_var[t] <- kal_var[t] + (g[t] ^ 2) * (rts_var[t + 1] - pbar[t + 1])
  }

  list(mean = rts_mean, variance = rts_var)
}

setClass("genealogy",representation(history="matrix", N="integer", Ngen="integer", model="character"))
setClass("samplegenealogy", representation(ancestry="matrix", samplesize="integer", sample="integer", MRCA="integer"),contains="genealogy")


## Initialise variables ------------------

## logarithmic series:
# log_n_particles_vals <- 4:12
# n_particles_vals <- 2 ^ (log_n_particles_vals)
# n_leaves <- n_particles_vals[1]

## linear series:
n_leaves <- 2
n_particles_vals <- seq(from = 256, to = 4096, by = 256)

n_reps <- 1000


## Generate (now read in) synthetic data ---------------

delta <- 0.1
sigma <- 0.1
n_obs <- 6 * n_particles_vals[length(n_particles_vals)]

ou_data <- read.csv("oudata.csv", header = FALSE)[[1]]

## Run simulations -----------------------

no_cores <- future::availableCores()
registerDoParallel(makeCluster(no_cores, type='FORK', outfile="debug_file.txt"))

for (n_sd_away in 0:3) {
  ## Generate immortal trajectory (using the Kalman filter solution)
  n_sd_away <- 3 ## conditioned trajectory should be how many SDs away from the MAP smoothed trajectory?
  rts <- ou_rts_smooth(ou_data, delta, sigma)
  rts_mean <- rts$mean
  rts_sd <- (rts$variance) ^ (0.5)
  imtl_states <- rts_mean + n_sd_away * rts_sd

  treeht_iters <- function(data, j, n_particles, imtl_states){
    anc <- conditional_SMC_anc(n_particles, data, imtl_states, ou_emission_density, ou_transition_sam, ou_initial_sam, sigma = sigma, delta = delta)
    class(anc) <- 'genealogy'
    anc@N <- as.integer(n_particles)
    anc@Ngen <- as.integer(n_obs - 1)
    ancestrySize_treeht(anc, sampl=sample(1:n_particles, n_leaves, replace=F))
  }

  for (i in 1:length(n_particles_vals)) {
    tree_height <- foreach(j = 1:n_reps, .combine = c)  %dorng% treeht_iters(ou_data, j, n_particles_vals[i], imtl_states)
    write.table(t(tree_height), file="treeht_out_n2.csv", sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
  }
}

stopImplicitCluster()


# ### Standard SMC ======================
# ## Run simulations -----------------------
#

# no_cores <- future::availableCores()
# registerDoParallel(makeCluster(no_cores, type='FORK', outfile="debug_file.txt"))
#
# treeht_iters <- function(data, j, n_particles){
#   anc <- standard_SMC_anc(n_particles, data, ou_emission_density, ou_transition_sam, ou_initial_sam, sigma = sigma, delta = delta)
#
#   class(anc) <- 'genealogy'
#   anc@N <- as.integer(n_particles)
#   anc@Ngen <- as.integer(n_obs - 1)
#
#   ancestrySize_treeht(anc, sampl=sample(1:n_particles, n_leaves, replace=F))
# }
#
# for (i in 1:length(n_particles_vals)) {
#   tree_height <- foreach(j = 1:n_reps, .combine = c)  %dorng% treeht_iters(ou_data, j, n_particles_vals[i])
#   write.table(t(tree_height), file="treeht_out_std.csv", sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
# }
#
#
# stopImplicitCluster()

