#' RTS smoother for OU model
#'
#' Calculate the Rauch-Tung-Striebel smoothing distribution for observations from an Ornstein-Uhlenbeck model
#'
#' @param observations vector of observations from the OU model
#' @param delta drift parameter of the Ornstein-Uhlenbeck process
#' @param sigma noise parameter of the Ornstein-Uhlenbeck process
#'
#' @return list containing the smoothed means and variances
#'
#' @author Suzie Brown
#'
#' @export ou_rts_smooth
#'

ou_rts_smooth <- function(observations, delta, sigma) {
  n_obs <- length(observations) - 1

  ## run Kalman filter
  kalman <- KalmanSmooth(observations, mod=list(T = 1 - delta, Z = 1, h = sigma ^ 2, V = delta, a = 0, P = 1, Pn = 1))
  kalman_mean <- kalman$smooth
  kalman_var <- kalman$var[,,1]

  ## initialise variables
  rts_mean <- numeric(n_obs + 1)
  rts_var <- numeric(n_obs + 1)
  pbar <- numeric(n_obs + 1)
  g <- numeric(n_obs)

  ## final state
  rts_mean[n_obs + 1] <- kalman_mean[n_obs + 1]
  rts_var[n_obs + 1] <- kalman_var[n_obs + 1]

  ## RTS smoother recursion
  for (t in n_obs:1) {
    pbar[t + 1] <- ((1 - delta) ^ 2) * kalman_var[t] + delta
    g[t] <- kalman_var[t] * (1 - delta) / pbar[t + 1]
    rts_mean[t] <- kalman_mean[t] + g[t] * (rts_mean[t + 1] - (1 - delta) * kalman_mean[t])
    rts_var[t] <- kalman_var[t] + (g[t] ^ 2) * (rts_var[t + 1] - pbar[t + 1])
  }

  list(mean = rts_mean, variance = rts_var)
}
