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

ou_sim <- function(n_obs, delta, sigma, return_states = FALSE) {
  states <- numeric(n_obs+1)
  obs <- numeric(n_obs+1)

  states[1] <- rnorm(1)
  obs[1] <- rnorm(1, states[1], sigma)
  for (i in 2:(n_obs+1)) {
    states[i] <- rnorm(1, (1 - delta) * states[i-1], delta ^ 0.5)
    obs[i] <- rnorm(1, states[i], sigma)
  }

  if (return_states) {
    return(list(states = states, observations = obs))
  }else{
    return(obs)
  }
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

setClass("genealogy",representation(history="matrix", N="integer", Ngen="integer", model="character"))
setClass("samplegenealogy", representation(ancestry="matrix", samplesize="integer", sample="integer", MRCA="integer"),contains="genealogy")


## Initialise variables ------------------

## logarithmic series:
# log_n_particles_vals <- 4:12
# n_particles_vals <- 2 ^ (log_n_particles_vals)
# n_leaves <- n_particles_vals[1]

## linear series:
n_leaves <- 16
n_particles_vals <- seq(from = 256, to = 4096, by = 256)

n_reps <- 100

## Generate synthetic data ---------------

delta <- 0.1
sigma <- 0.1
n_obs <- 2 * n_particles_vals[length(n_particles_vals)]

ou_data <- ou_sim(n_obs, delta, sigma)

### Conditional SMC ======================
## Generate immortal trajectory ----------

std_n_particles <- 100

std_smc <- standard_SMC(std_n_particles, ou_data, ou_emission_density, ou_transition_sam, ou_initial_sam, sigma = sigma, delta = delta)
imtl_index <- sample(1:std_n_particles, 1, prob = std_smc$weights[nrow(std_smc$weights), ])

std_anc <- std_smc$ancestors
class(std_anc) <- 'genealogy'
std_anc@N <- as.integer(std_n_particles)
std_anc@Ngen <- as.integer(n_obs - 1)

imtl_anc <- traceAncestry(std_anc, sampl = imtl_index)

imtl_states <- numeric(n_obs)
for (t in 1:(n_obs - 1)) {
  imtl_states[t] <- std_smc$positions[t, imtl_anc[t]]
}
imtl_states[n_obs] <- std_smc$positions[n_obs, imtl_index]

## Run simulations -----------------------

treeht_iters <- function(data, j, n_particles, imtl_states){
  anc <- conditional_SMC_anc(n_particles, data, imtl_states, ou_emission_density, ou_transition_sam, ou_initial_sam, sigma = sigma, delta = delta)

  class(anc) <- 'genealogy'
  anc@N <- as.integer(n_particles)
  anc@Ngen <- as.integer(n_obs - 1)

  ancestrySize_treeht(anc, sampl=sample(1:n_particles, n_leaves, replace=F))
}

no_cores <- future::availableCores()
registerDoParallel(makeCluster(no_cores, type='FORK', outfile="debug_file.txt"))

for (i in 1:length(n_particles_vals)) {
  tree_height <- foreach(j = 1:n_reps, .combine = c)  %dorng% treeht_iters(ou_data, j, n_particles_vals[i], imtl_states)
  write.table(t(tree_height), file="treeht_out.csv", sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
}

stopImplicitCluster()

# ## Process results -----------------------
#
# tree_height <- read.csv("treeht_out.csv", header=FALSE)
#
# mean_tree_ht <- apply(tree_height, 1, mean)
# var_tree_ht <- apply(tree_height, 1, var)
# se_tree_ht <- (var_tree_ht / n_reps) ^ 0.5
#
# plot(n_particles_vals, mean_tree_ht / n_particles_vals, ylim = c(0.25, 0.45), type = 'b', pch = 16, col = 2, lwd = 2, xlab = "number of particles", ylab = "mean tree height / no. particles", main = paste("Tree height profile: conditional SMC, no. leaves=", n_leaves))
# lines(n_particles_vals, (mean_tree_ht - se_tree_ht) / n_particles_vals, col = 2, lty = 2)
# lines(n_particles_vals, (mean_tree_ht + se_tree_ht) / n_particles_vals, col = 2, lty = 2)
# legend("topright", c("mean", "+/- 1 std err"), lty = c(1, 2), lwd = c(2,1), col = 2, pch = c(16, NA))
#

# ### Standard SMC ======================
# ## Run simulations -----------------------
#
# treeht_iters <- function(data, j, n_particles, imtl_states){
#   anc <- standard_SMC_anc(n_particles, data, ou_emission_density, ou_transition_sam, ou_initial_sam, sigma = sigma, delta = delta)
#
#   class(anc) <- 'genealogy'
#   anc@N <- as.integer(n_particles)
#   anc@Ngen <- as.integer(n_obs - 1)
#
#   ancestrySize_treeht(anc, sampl=sample(1:n_particles, n_leaves, replace=F))
# }
#
# no_cores <- future::availableCores()
# registerDoParallel(makeCluster(no_cores, type='FORK', outfile="debug_file.txt"))
#
# for (i in 1:length(n_particles_vals)) {
#   tree_height <- foreach(j = 1:n_reps, .combine = c)  %dorng% treeht_iters(ou_data, j, n_particles_vals[i], imtl_states)
#   write.table(t(tree_height), file="treeht_out.csv", sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
# }
#
# stopImplicitCluster()
#
## Process results -----------------------

# tree_height <- read.csv("treeht_out.csv", header=FALSE)
#
# mean_tree_ht <- apply(tree_height, 1, mean)
# var_tree_ht <- apply(tree_height, 1, var)
# se_tree_ht <- (var_tree_ht / n_reps) ^ 0.5
#
# plot(n_particles_vals, mean_tree_ht / n_particles_vals, ylim=c(0.1, 0.3), type = 'b', pch = 16, col = 2, lwd = 2, xlab = "number of particles", ylab = "mean tree height / no. particles", main = paste("Tree height profile: standard SMC, no. leaves=", n_leaves))
# lines(n_particles_vals, (mean_tree_ht - se_tree_ht) / n_particles_vals, col = 2, lty = 2)
# lines(n_particles_vals, (mean_tree_ht + se_tree_ht) / n_particles_vals, col = 2, lty = 2)
# legend("topright", c("mean", "+/- 1 std err"), lty = c(1, 2), lwd = c(2,1), col = 2, pch = c(16, NA))
