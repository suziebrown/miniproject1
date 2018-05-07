Delta <- 0.1
sigma <- 0.1
tmax <- 15

## discretised Ornstein-Uhlenbeck process

sim.data <- function(tmax, Delta, sigma){
  x <- numeric(tmax+1)
  y <- numeric(tmax+1)
  x[1] <- rnorm(1)
  y[1] <- rnorm(1,x[1], sigma)
  for (t in 2:(tmax+1)){
    x[t] <- (1-Delta)*x[t-1] + Delta^(-0.5)*rnorm(1)
    y[t] <- rnorm(1,x[t],sigma)
  }
  list(x=x,y=y)
}

K <- function(x){ ## transition function of HMM
  (1-Delta)*x + Delta^(0.5)*rnorm(1)
}

p <- function(x,xprime){
  (2*pi*Delta)^(-0.5) * exp(-(xprime - (1-Delta)*x)^2/(2*Delta))
}
q <- function(x,xprime){ ## proposal distribution for sampling p
  p(x,xprime)
}

phi <- function(x, yt){
  (2*pi*Delta)^(-0.5) * exp(-(yt-x)^2/(2*sigma^2))
}
g <- function(x,yt){
  phi(x, yt)
}

data <- sim.data(tmax, Delta, sigma)
x <- data$x
y <- data$y
plot(x, type='l', col=2)
lines(y, col=1)

stdSMC <- function(N, n.step=length(y), y){
  A <- matrix(NA, nrow=n.step, ncol=N)
  x <- rnorm(N) ## initial state
  w <- g(x, y[1]) ## calculate weights
  w <- w/sum(w) ## normalise weights
  for (t in 1:n.step){
    A[t,] <- sample(1:N, N, prob=w, replace=T) ## choose the parent of each offspring
    x <- x[a] ## resample the particles
    x <- K(x) ## propagate particles
    w <- g(x,y[t+1]) ## caluculate weights
    w <- w/sum(w) ## normalise weights
  }
  list(x=x,w=w,A=A)
}

foo <- stdSMC(N,y)

