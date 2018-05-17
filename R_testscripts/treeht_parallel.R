library(foreach)
library(doParallel)
library(future)
library(doRNG)

## function dependencies --------------

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

condSMC <- function(N, y, cond_states, n_step=length(y)-1){
  A <- matrix(NA, nrow=n_step, ncol=N)
  X <- matrix(NA, nrow=n_step+1, ncol=N)
  W <- matrix(NA, nrow=n_step+1, ncol=N)

  X[1,2:N] <- rnorm(N-1) ## initial state of free particles
  X[1,1] <- cond_states[1] ## known initial state of immortal particle

  w <- p(y[1], X[1,]) ## calculate weights
  W[1,] <- w/sum(w) ## normalise weights

  for (t in 2:(n_step+1)){
    A[t-1,] <- sample(1:N, N, prob=w, replace=T) ## choose the parent of each offspring
    A[t-1,1] <- 1 ## insist on continuing the immortal line
    x <- X[t-1,A[t-1,]] ## resample the particles

    X[t,2:N] <- rnorm(N-1, (1-Delta)*X[t-1,2:N], Delta^(0.5)) ## propagate free particles
    X[t,1] <- cond_states[t] ## known next state of immortal particle

    w <- p(y[t], X[t,]) ## caluculate weights
    W[t,] <- w/sum(w) ## normalise weights
  }
  list(positions=X,weights=W,ancestry=A)
}


p <- function(y,x){ ## emission density y[t] | x[t]
  (2*pi)^(-0.5)*sigma^(-1) * exp(-(y-x)^2/(2*sigma^2))
}

ancestrySize <- function(history, sampl=NULL, maxgen=NULL){
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

  list(familySize=Size, treeHeight=N.gen-MRCA)
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

##  variable initialisation -----------------
setClass("genealogy",representation(history="matrix", N="integer", Ngen="integer", model="character"))
setClass("samplegenealogy", representation(ancestry="matrix", samplesize="integer", sample="integer", MRCA="integer"),contains="genealogy")

Delta <- 0.1
sigma <- 0.1
max_sam_size <- 13
min_sam_size <- 4
n_sam_size <- max_sam_size-min_sam_size+1
Nvals <- 2^(min_sam_size:max_sam_size)
tmax <- max(Nvals) ## number of observations
n.reps <- 100 ## number of repetitions each
n <- 2^min_sam_size
Nmax <- max(Nvals)

data <- sim.data(tmax, Delta, sigma)
x <- data$x
y <- data$y

treeHeight <- matrix(NA, nrow=n_sam_size, ncol=n.reps)


## parallel execution (standard SMC) --------------------

# SMC_treeht_reps <- function(y,i,j, Nvals){
#   samSMC <- stdSMC(Nvals[i],y)
#
#   anc <- samSMC$ancestry
#   class(anc) <- 'genealogy'
#   anc@N <-as.integer(Nvals[i])
#   anc@Ngen <- as.integer(tmax)
#
#   ancestrySize(anc, sampl=sample(1:Nvals[i],n,replace=F))$treeHeight
# }
#
#
# no_cores <- future::availableCores() # -1  # if using desktop
# registerDoParallel(makeCluster(no_cores, type='FORK', outfile="debug_file.txt"))
#
# for (i in 1:n_sam_size){
#   treeHeight_i <- foreach(j=1:n.reps, .combine = c)  %dorng% SMC_treeht_reps(y,i,j,Nvals)
#   write.table(t(treeHeight_i), file="treeht_out.csv", sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
# }
#
# stopImplicitCluster()


## parallel execution (conditional SMC) --------------------

SMC_treeht_reps <- function(y, i, j, Nvals, cond_states){
  samSMC <- condSMC(Nvals[i], y, cond_states)

  anc <- samSMC$ancestry
  class(anc) <- 'genealogy'
  anc@N <- as.integer(Nvals[i])
  anc@Ngen <- as.integer(tmax)

  ancestrySize(anc, sampl=sample(1:Nvals[i], n, replace=F))$treeHeight
}

## find states to condition on
cond_run <- stdSMC(Nmax,y)
cond_ind_sam <- sample(1:Nmax, 1, prob=cond_run$weights[nrow(cond_run$weights),])
cond_anc <- cond_run$ancestry
class(cond_anc) <- 'genealogy'
cond_anc@N <- as.integer(Nmax)
cond_anc@Ngen <- as.integer(tmax)
cond_anc_sam <- traceAncestry(cond_anc, sampl=cond_ind_sam)
cond_states <- numeric(tmax+1)
for (t in 1:tmax){
  cond_states[t] <- cond_run$positions[t,cond_anc_sam[t]]
}
cond_states[tmax+1] <- cond_run$positions[tmax+1,cond_ind_sam]


no_cores <- future::availableCores() # -1  # if using desktop
registerDoParallel(makeCluster(no_cores, type='FORK', outfile="debug_file.txt"))

for (i in 1:n_sam_size){
  treeHeight_i <- foreach(j=1:n.reps, .combine = c)  %dorng% SMC_treeht_reps(y,i,j,Nvals,cond_states)
  write.table(t(treeHeight_i), file="treeht_out.csv", sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
}

stopImplicitCluster()


## process results -----------------------

# treeHeight <- read.csv("treeht_out.csv", header=FALSE)
#
# meanTreeHt <- apply(treeHeight,1,mean)
# varTreeHt <- apply(treeHeight,1,var)
# seTreeHt <- (varTreeHt/n.reps)^0.5
#
# plot(2^(min_sam_size:max_sam_size), meanTreeHt/Nvals, type='b', pch=16, col=2, xlab='N', ylab='mean tree height / N', main=paste('Tree height profile: conditional SMC, n=',n))
# lines(2^(min_sam_size:max_sam_size), (meanTreeHt-seTreeHt)/Nvals, col=2, lty=2)
# lines(2^(min_sam_size:max_sam_size), (meanTreeHt+seTreeHt)/Nvals, col=2, lty=2)
# legend("topright", c("mean","+/- 1 std error"), lty=c(1,2), col=2, pch=c(16,NA))
