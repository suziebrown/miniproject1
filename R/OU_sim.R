#' Ornstein-Uhlenbeck simulation
#'
#' Simulate synthetic data from the Ornstien-Uhlenbeck model
#'
#' @param n_obs the number of observations (excluding initial state) to produce, i.e. the number of time steps
#' @param delta drift parameter of the Ornstein-Uhlenbeck process
#' @param sigma noise parameter of the Ornstein-Uhlenbeck process
#' @param return_states return the hidden states as well as the noisy observations?
#'
#' @return a vector of (n_obs + 1) observations; or if return_states = TRUE, a list containing a vector of hidden
#' states and a corresponding vector of observations
#'
#' @author Suzie Brown
#'
#' @export ou_sim
#'

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
