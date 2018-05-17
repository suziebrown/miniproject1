#' Standard SMC algorithm
#'
#' Standard SMC implementation with multinomial resampling for a generic HMM
#'
#' @param n_particles number of particles
#' @param observations a vector of observations on which to perform inference
#' @param emission_density a function of the form f(y,x) giving the emission density of y|x
#' @param transition_sam a function of the form f(x) sampling a new position, i.e. propagating the particle
#' @param initial_sam a function of the form f() sampling the initial positions
#' @param ... arguments to be passed to emission_density, transition_sam and initial_sam
#' (e.g. for Ornstein-Uhlenbeck process, provide values for sigma and delta)
#'
#' @return a list containing matrices of positions and weights
#' (rows represent generations/time steps, columns represent individuals/particles)
#' and an ancestry matrix (where element (i,j) is the parent in generation i of individual j in generation i+1)
#'
#' @author Suzie Brown
#'
#' @export standard_SMC
#'

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


#' Conditional SMC algorithm
#'
#' Conditional SMC implementation with multinomial resampling for a generic (exchangeable) HMM
#'
#' @param n_particles number of particles
#' @param observations a vector of observations on which to perform inference
#' @param cond_trajectory a vector of states representing the trajectory on which to condition
#' @param emission_density a function of the form f(y,x) giving the emission density of y|x
#' @param transition_sam a function of the form f(x) sampling a new position, i.e. propagating the particle
#' @param initial_sam a function of the form f() sampling the initial positions
#' @param ... arguments to be passed to emission_density, transition_sam and initial_sam
#' (e.g. for Ornstein-Uhlenbeck process, provide values for sigma and delta)
#'
#' @return a list containing matrices of positions and weights
#' (rows represent generations/time steps, columns represent individuals/particles)
#' and an ancestry matrix (where element (i,j) is the parent in generation i of individual j in generation i+1)
#'
#' @author Suzie Brown
#'
#' @export conditional_SMC
#'

conditional_SMC <- function(n_particles, observations, cond_trajectory, emission_density, transition_sam, initial_sam, ...){
  n_obs <- length(observations)
  ancestors <- matrix(NA, nrow = n_obs - 1, ncol = n_particles)
  positions <- matrix(NA, nrow = n_obs, ncol = n_particles)
  weights <- matrix(NA, nrow = n_obs, ncol = n_particles)

  ## Initial states
  positions[1, 1] <- cond_trajectory[1]
  positions[1, 2:n_particles] <- initial_sam(n_particles - 1, ...)
  weights_tmp <- emission_density(observations[1], positions[1, ], ...)
  weights[1, ] <- weights_tmp / sum(weights_tmp)

  for (t in 2:(n_obs)){
    ## Resample
    ancestors[t-1, ] <- sample(1:n_particles, n_particles, prob=weights_tmp, replace=TRUE)
    ancestors[t-1, 1] <- 1
    positions_tmp <- positions[t-1, ancestors[t-1, ]]

    ## Propagate
    positions[t, 2:n_particles] <- transition_sam(positions_tmp[2:n_particles], ...)
    positions[t, 1] <- cond_trajectory[t]

    ## Calculate importance weights
    weights_tmp <- emission_density(observations[t], positions[t, ], ...)
    weights[t, ] <- weights_tmp / sum(weights_tmp)
  }

  list(positions = positions, weights = weights, ancestors = ancestors)
}
