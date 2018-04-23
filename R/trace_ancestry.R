#' Trace sample ancestry
#'
#' Trace the ancestry of individuals in a sample and find their most recent common ancestor
#'
#' @param history an ocject of class genealogy, as returned from WFsim or moranSim
#' @param sampl a vector of integers indicating which individuals from the most recent generation are to be included in the sample
#' @param maxgen the maximum number of generations back in time to investigate
#'
#' @return an object of class samplegenealogy
#'
#' @author Suzie Brown
#'
#' @export traceAncestry
#'

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

