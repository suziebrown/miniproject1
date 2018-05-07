Delta <- 0.1
sigma <- 0.1
T <- 1000

## discretised Ornstein-Uhlenbeck process

sim.data <- function(T, Delta, sigma){
  x <- numeric(T+1)
  y <- numeric(T+1)
  x[1] <- rnorm(1)
  y[1] <- rnorm(1,x[1], sigma^2)
  for (t in 2:(T+1)){
    x[t] <- (1-Delta)*x[t-1] + Delta^(-0.5)*rnorm(1)
    y[t] <- rnorm(1,x[t],sigma^2)
  }
  list(x=x,y=y)
}

p <- function(x,xprime){
  (2*pi()*Delta)^(-0.5) * exp(-(xprime - (1-Delta)*x)^2/(2*Delta))
}

phi <- function(xprime, yt){
  (2*pi()*Delta)^(-0.5) * exp(-(yt-xprime)^2/(2*sigma^2))
}


data <- sim.data(T, Delta, sigma)
datax <- data$x
datay <- data$y
plot(datax, type='l', col=2)
lines(datay, col=1)

