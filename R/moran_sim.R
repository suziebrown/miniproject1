#' Moran Sample
#'
#' Generate a realisation of the Moran population model
#'
#' @param N constant population size
#' @param N.gen number of generations (single reproductions) to simulate
#' @param rescale return only every Nth generation? (for comparability with W-F model)
#'
#' @details The attribute 'history' contains an (N.gen-1)xN matrix , where the (i,j)th entry specifies the parent of individual j of generation i+1.
#' For comparability to Wright-Fisher model, multiply N.gen by N: one generation in W-F model comprises N reproduction events.
#' Note that when rescale=TRUE, the returned value of N.gen is the scaled number of generations, and is not equal to the inputted N.gen.
#'
#' @return an object of class 'genealogy'
#'
#' @author Suzie Brown
#'
#' @export moranSim
#'

moranSim <- function(N, N.gen, rescale=FALSE){
  history <- matrix(rep(1:N,N.gen-1), ncol=N, nrow=N.gen-1, byrow=T)
  for (t in 1:(N.gen-1)){ ## can probably eliminate this loop
    inds <- sample(1:N,2,replace=T) ## inds1=birth individual; inds2=death individual
    if (inds[1] != inds[2]){
      history[t,inds[2]] <- inds[1]
      history[t,] <- sort(history[t,])
    }
  }

  if (rescale){
    if(N.gen<2*N){stop("must have N.gen>=2*N for N-rescaling")} ## check arguments are valid
    N.regen <- floor((N.gen-1)/N)+1 ## how many rescaled generations will there be
    rehistory <- matrix(NA, ncol=N, nrow=N.regen-1)
    ancestry <- 1:N

    for (T in (N.regen):2){
      for (t in 1:N){
        ancestry <- history[N.gen+N*(T-N.regen)-t, ancestry]
      }
      rehistory[T-1,] <- ancestry
    }
    history <- rehistory
  }

  class(history) <- 'genealogy'
  history@N <- as.integer(N)
  if(rescale){
    history@model <- 'N-rescaled Moran'
    history@Ngen <- as.integer(N.regen)
  }
  else{
    history@model <- 'Moran'
    history@Ngen <- as.integer(N.gen)
  }
  history
}
