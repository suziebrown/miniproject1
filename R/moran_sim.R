#' Moran Sample
#'
#' Generate a realisation of the Moran population model
#'
#' @param N constant population size
#' @param N.gen number of generations (single reproductions) to simulate
#'
#' @details The attribute 'history' contains an (N.gen-1)xN matrix , where the (i,j)th entry specifies the parent of individual j of generation i+1. For comparability to Wright-Fisher model, multiply N.gen by N: one generation in W-F model comprises N reproduction events.
#'
#' @return an object of class 'genealogy'
#'
#' @author Suzie Brown
#'
#' @export moranSim
#'

moranSim <- function(N, N.gen){
  history <- matrix(rep(1:N,N.gen-1), ncol=N, nrow=N.gen-1, byrow=T)
  for (t in 1:(N.gen-1)){ ## can probably eliminate this loop
    inds <- sample(1:N,2,replace=T) ## inds1=birth individual; inds2=death individual
    if (inds[1] != inds[2]){
      history[t,inds[2]] <- inds[1]
      history[t,] <- sort(history[t,])
    }
  }
  class(history) <- 'genealogy'
  history@model <- 'Moran'
  history@N <- N
  history@Ngen <- N.gen
  history
}
