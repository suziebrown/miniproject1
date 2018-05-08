Delta <- 0.1
sigma <- 0.1
tmax <- 10
N <- 200

## discretised Ornstein-Uhlenbeck process

sim.data <- function(tmax, Delta, sigma){
  x <- numeric(tmax+1)
  y <- numeric(tmax+1)
  x[1] <- rnorm(1)
  y[1] <- rnorm(1,x[1], sigma)
  for (t in 2:(tmax+1)){
    x[t] <- rnorm(1, (1-Delta)*x[t-1], Delta^(0.5))
    y[t] <- rnorm(1,x[t],sigma)
  }
  list(x=x,y=y)
}

# K <- function(x){ ## transition function of HMM
#   (1-Delta)*x + Delta^(0.5)*rnorm(length(x))
# }
#
# g <- function(x,xprime){ ## transition density x[t+1] | x[t]
#   (2*pi*Delta)^(-0.5) * exp(-(xprime - (1-Delta)*x)^2/(2*Delta))
# }

p <- function(y,x){ ## emission density y[t] | x[t]
  (2*pi)^(-0.5)*sigma^(-1) * exp(-(y-x)^2/(2*sigma^2))
}

data <- sim.data(tmax, Delta, sigma)
x <- data$x
y <- data$y
plot(x, type='l', col=2)
lines(y, col=1)

stdSMC <- function(N, y, n.step=length(y)-1){
  A <- matrix(NA, nrow=n.step, ncol=N)
  X <- matrix(NA, nrow=n.step+1, ncol=N)
  W <- matrix(NA, nrow=n.step+1, ncol=N)
  X[1,] <- rnorm(N) ## initial state
  w <- p(y[1], X[1,]) ## calculate weights
  W[1,] <- w/sum(w) ## normalise weights
  for (t in 2:(n.step+1)){
    A[t-1,] <- sample(1:N, N, prob=w, replace=T) ## choose the parent of each offspring
    x <- X[t-1,A[t-1,]] ## resample the particles
    X[t,] <- rnorm(N, (1-Delta)*X[t-1,], Delta^(0.5)) ## propagate particles
    w <- p(y[t], X[t,]) ## caluculate weights
    W[t,] <- w/sum(w) ## normalise weights
  }
  list(positions=X,weights=W,ancestry=A)
}

foo <- stdSMC(N,y)
meanx <- rowSums(foo$positions)/N
varx <- apply(foo$positions,1,var)/(N^(0.5))
points(meanx, pch=16, col=2)
lines(meanx+varx^(0.5), col=2, lty=2)
lines(meanx-varx^(0.5), col=2, lty=2)

# A <- foo$ancestry
# class(A) <- 'genealogy'
# A@N <-as.integer(N)
# A@Ngen <- as.integer(tmax)
# plot(A, highlight.sample = T)
