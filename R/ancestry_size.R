#' Sample ancestry size
#'
#' Find the number of ancestors of a given sample in each generation
#'
#' @param history an ocject of class genealogy, as returned from WFsim or moranSim
#' @param sampl a vector of integers indicating which individuals from the most recent generation are to be included in the sample
#' @param maxgen the maximum number of generations back in time to investigate
#'
#' @return a vector indicating the number of ancestors in each generation 1,...,N.gen
#'
#' @author Suzie Brown
#'
#' @export ancestrySize
#'

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
