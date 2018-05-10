## make plots of tree height: with n fixed, N varying
## standard SMC with exchangeable multinomial resampling
## data from Ornstein-Uhlenbeck process


## simulate observed data ----------------------

Delta <- 0.1
sigma <- 0.1
max_sam_size <- 12
min_sam_size <- 4
n_sam_size <- max_sam_size-min_sam_size+1
Nvals <- 2^(min_sam_size:max_sam_size)
tmax <- 7*max(Nvals)
n.reps <- 100
n <- 2^min_sam_size

data <- sim.data(tmax, Delta, sigma)
x <- data$x
y <- data$y


## run particle filter and store tree height -----------

treeHeight <- matrix(NA, nrow=n_sam_size, ncol=n.reps)

for (i in 1:n_sam_size){
  for (j in 1:n.reps){

    samSMC <- condSMC(Nvals[i],y) ## change here to switch between standard and conditional SMC

    anc <- samSMC$ancestry
    class(anc) <- 'genealogy'
    anc@N <-as.integer(Nvals[i])
    anc@Ngen <- as.integer(tmax)

    treeHeight[i,j] <- ancestrySize(anc, sampl=sample(1:Nvals[i],n,replace=F))$treeHeight
  }
}


## make plot -----------------------

meanTreeHt <- apply(treeHeight,1,mean)
varTreeHt <- apply(treeHeight,1,var)
seTreeHt <- (varTreeHt/n.reps)^0.5

plot(2^(1:n_sam_size), meanTreeHt/Nvals, type='b', pch=16, col=2,ylim=c(0,2), xlab='N', ylab='mean tree height / N', main=paste('Tree height profile: conditional SMC, n=',n))
lines(2^(1:n_sam_size), (meanTreeHt-seTreeHt)/Nvals, col=2, lty=2)
lines(2^(1:n_sam_size), (meanTreeHt+seTreeHt)/Nvals, col=2, lty=2)
legend("topright", c("mean","+/- 1 std error"), lty=c(1,2), col=2, pch=c(16,NA))

