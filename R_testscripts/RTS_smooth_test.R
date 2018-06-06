## Smoothing solution for Ornstein-Uhlenbeck model
## This is the test script: the approved version is OU_RTS_smooth.R in R folder.

## Test values ---------------------------

delta <- 0.1
sigma <- 0.1
n_obs <- 20
ou_data <- ou_sim(n_obs, delta, sigma)

## compute Kalman filter solution --------

kalman <- KalmanSmooth(ou_data, mod=list(T = 1 - delta, Z = 1, h = sigma ^ 2, V = delta, a = 0, P = 1, Pn = 1))
kalman_mean <- kalman$smooth
kalman_var <- kalman$var[,,1]


## Rauch-Tung-Striebel smoother ----------
## https://users.aalto.fi/~ssarkka/course_k2011/pdf/handout7.pdf (p.18)

rts_mean <- numeric(n_obs + 1)
rts_var <- numeric(n_obs + 1)
pbar <- numeric(n_obs + 1)
g <- numeric(n_obs)

rts_mean[n_obs + 1] <- kalman_mean[n_obs + 1]
rts_var[n_obs + 1] <- kalman_var[n_obs + 1]

for (t in n_obs:1) {
  pbar[t + 1] <- ((1 - delta) ^ 2) * kalman_var[t] + delta
  g[t] <- kalman_var[t] * (1 - delta) / pbar[t + 1]
  rts_mean[t] <- kalman_mean[t] + g[t] * (rts_mean[t + 1] - (1 - delta) * kalman_mean[t])
  rts_var[t] <- kalman_var[t] + (g[t] ^ 2) * (rts_var[t + 1] - pbar[t + 1])
}


## Test plot -----------------------------

plot(ou_data, type='l', col=2, lwd=2)
points(rts_mean, pch='+')
lines(rts_mean + rts_var^(0.5), lty=2)
lines(rts_mean - rts_var^(0.5), lty=2)
