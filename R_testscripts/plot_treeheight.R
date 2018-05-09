## make plots of tree height, following KJJS p.29
## standard SMC with exchangeable multinomial resampling
## data from Ornstein-Uhlenbeck process


## simulate observed data ----------------------

Delta <- 0.1
sigma <- 0.1
n_sam_size <- 12
N <- 2^n_sam_size
tmax <- 7*N
n.reps <- 1000

data <- sim.data(tmax, Delta, sigma)
x <- data$x
y <- data$y


## run particle filter and store tree height -----------

treeHeight <- matrix(NA, nrow=n_sam_size, ncol=n.reps)

for (j in 1:n.reps){
  samSMC <- condSMC(N,y) ## change here to switch between standard and conditional SMC

  anc <- samSMC$ancestry
  class(anc) <- 'genealogy'
  anc@N <-as.integer(N)
  anc@Ngen <- as.integer(tmax)

  for (i in 1:n_sam_size){
    treeHeight[i,j] <- ancestrySize(anc, sampl=sample(1:N,2^i,replace=F))$treeHeight
  }
}


## make plot -----------------------

meanTreeHt <- apply(treeHeight,1,mean) ## na.rm=T => mean and var are actually HIGHER than the calculated ones
varTreeHt <- apply(treeHeight,1,var)
seTreeHt <- (varTreeHt/n.reps)^0.5

plot(2^(1:n_sam_size), meanTreeHt/N, type='b', pch=16, col=2, xlab='n', ylab='mean tree height / N', main=paste('Tree height profile: conditional SMC, N=',N))
lines(2^(1:n_sam_size), (meanTreeHt-seTreeHt)/N, col=2, lty=2)
lines(2^(1:n_sam_size), (meanTreeHt+seTreeHt)/N, col=2, lty=2)


## investigate treeHeight -----------

treeHeight_temp <- treeHeight

summary(t(treeHeight))
par(mfrow=c(4,3))
for (i in 1:10){
  hist(treeHeight[i,], xlab=c(0,max(treeHeight)), main=paste('Histogram of treeHeight[',i,',]'))
}
par(mfrow=c(1,1))
boxplot(t(treeHeight),main='Distribution of tree height: condSMC', ylab='tree height', xlab='log2(n)')
points(2:10, rep(5200,9), pch='+', col=2)

numNA <- apply(treeHeight,1,function(x) sum(is.na(x)))








