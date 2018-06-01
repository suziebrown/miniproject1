#' RTS smoother for OU model
#'
#' Calculate the Rauch-Tung-Striebel smoothing distribution for observations from an Ornstein-Uhlenbeck model
#'
#' @param observations vector of observations from the OU model
#' @param delta drift parameter of the Ornstein-Uhlenbeck process
#' @param sigma noise parameter of the Ornstein-Uhlenbeck process
#' @param prior_mean prior mean
#' @param prior_var prior variance
#'
#' @return list containing the smoothed means and variances
#'
#' @author Suzie Brown
#'
#' @export ou_rts_smooth
#'

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
