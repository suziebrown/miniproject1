#' Wright-Fisher Sample
#'
#' Generate a realisation of the Wright-Fisher population model
#'
#' @param N constant population size
#' @param N.gen number of generations to simulate
#' @param prob function providing probabilities for (non-uniform) multinomial sampling of parents
#'
#' @return an object of class 'genealogy'
#'
#' @details The attribute 'history' contains an (N.gen-1)xN matrix , where the (i,j)th entry specifies the parent of individual j of generation i+1.
#' The function prob should take a vector argument specifying a set of individuals; and should return a vector of (unnormalised) probability weights.
#'
#' @author Suzie Brown
#'
#' @export WFsim
#'

WFsim <- function(N, N.gen, prob=NULL){
  history <- matrix(NA, ncol=N, nrow=N.gen-1)
  if (is.null(prob)){
    history <- t(replicate(N.gen-1,sort(sample(1:N, N, replace=T))))
  }
  else{
    history <- t(replicate(N.gen-1,sort(sample(1:N, N, replace=T, prob=prob(1:N)))))
  }
  class(history) <- 'genealogy'
  if (is.null(prob)){
    history@model <- "Wright-Fisher"
  }
  else{
    history@model <- "Generalised Wright-Fisher"
  }
  history@N <- as.integer(N)
  history@Ngen <- as.integer(N.gen)
  history
}
