## Absorbing Markov chain analysis for Moran model coalescence

N <- 20 ## absorbing state is c=N-1

## for WF model:
Q <- matrix(0, nrow=N-2, ncol=N-2)
for (i in 1:(N-2)){
  for (j in i:(N-2)){
    Q[i,j] <- Stirling2(N-i,N-j)*factorial(N)/(factorial(j)*N^(N-i))
  }
}

## for Moran model:
# Q <- matrix(0, nrow=N-2, ncol=N-2) ## build transition matrix over transient states (0,...,N-2)
# for (i in 1:(N-3)){ ## c=i-1
#   Q[i,i] <- 1-(N-i+1)*(N-i)/N^2
#   Q[i,i+1] <- (N-i+1)*(N-i)/N^2
# }
# Q[N-2,N-2] <- 1-(N-i+1)*(N-i)/N^2

I <- diag(N-2)
R <- solve(I-Q) ## R is the matrix called N on wikipedia

tt <- rowSums(R) ## called t on wikipedia
Et <- tt[1] ## expected number of steps before absorption (full coalescence)
Vt <- (2*R-I)%*%tt - tt^2
Vart <- Vt[1] ## variance of number of steps before absorption (full coalescence)


#~~~ PLOTS

minN <- 4
maxN <- 500
Et <- numeric(maxN-minN+1)
Vart <- numeric(maxN-minN+1)

for (N in minN:maxN){
  ## MORAN:
  # Q <- matrix(0, nrow=N-2, ncol=N-2) ## build transition matrix over transient states (0,...,N-2)
  # for (i in 1:(N-3)){ ## c=i-1
  #   Q[i,i] <- 1-(N-i+1)*(N-i)/N^2
  #   Q[i,i+1] <- (N-i+1)*(N-i)/N^2
  # }
  # Q[N-2,N-2] <- 1-(N-i+1)*(N-i)/N^2
  ##~~

  ## WRIGHT-FISHER:
  Q <- matrix(0, nrow=N-2, ncol=N-2)
  for (i in 1:(N-2)){
    for (j in i:(N-2)){
      Q[i,j] <- Stirling2(N-i,N-j)*factorial(N)/(factorial(j)*N^(N-i))
    }
  }
  ##~~

  I <- diag(N-2)
  R <- solve(I-Q) ## R is the matrix called N on wikipedia

  tt <- rowSums(R) ## called t on wikipedia
  Et[N-minN+1] <- tt[1] ## expected number of steps before absorption (full coalescence)
  Vt <- (2*R-I)%*%tt - tt^2
  Vart[N-minN+1] <- Vt[1] ## variance of number of steps before absorption (full coalescence)
}

## WF plots:
plot(minN:maxN, Et, xlab='population size', ylab='expected number of generations to full coalescence', main='Wright-Fisher model coalescence')
plot(Et/(minN:maxN), type='l', lwd=3,col='purple', xlab='population size N', ylab='E(tau)*1/N', main='Wright-Fisher model coalescence')

## moran plots:
plot(minN:maxN, Et/(minN:maxN), xlab='population size', ylab='expected no. rescaled generations to full coalescence', main='N-rescaled Moran model coalescence')
plot(2*Et/(minN:maxN)^2, type='l', lwd=3,col='purple', xlab='population size N', ylab='E(tau)*2/N^2', main='Moran model coalescence')

## plot transition matirx (originally for debugging):
N <- 70
P <- matrix(0, nrow=N-1, ncol=N-1)
for (i in 1:(N-2)){
  for (j in i:(N-1)){
    P[i,j] <- Stirling2(N-i,N-j)*factorial(N)/(factorial(j)*N^(N-i))
  }
}
P[N-1,N-1] <- 1
image(P, ylim=c(1,0), col=suzie.colors(1000, rev=F), xaxt='n', yaxt='n', main=paste('transition matrix elements P_ij with N=',N, sep=''))

plot(diag(P), type='l')
