#' Wright-Fisher Sample
#'
#' Generate a realisation of the Wright-Fisher population model
#'
#' @param N constant population size
#' @param N.gen number of generations to simulate
#'
#' @return an (N.gen-1)xN matrix , where the (i,j)th entry specifies the parent of individual j of generation i+1
#'
#' @author Suzie Brown
#'
#' @export WFsim
#'

WFsim <- function(N, N.gen){
  history <- matrix(NA, ncol=N, nrow=N.gen-1)
  for (t in 1:(N.gen-1)){ ## can probably eliminate this loop e.g. using replicate?
    history[t,] <- sort(sample(1:N, N, replace=T))
  }
  history
}

#' Wright-Fisher Plot
#'
#' Plot a realisation of the Wright-Fisher population model
#'
#' @param history a matrix as returned from WFsim, where the (i,j)th entry specifies the parent of individual j of generation i+1
#' @param col.line the colour used to plot parent-offspring lines
#' @param shade.dead use a different colour to shade individuals with no offspring (dead lineages)?
#' @param col.alive if shade.dead=TRUE, the colour to shade individuals with offspring
#' @param col.dead if shade.dead=TRUE, the colour to shade individuals without offspring
#' @param highlight.sample highlight the ancestry of a sample from the most recent generation?
#' @param sampl if highlight.sample=TRUE, which individuals from the most recent generation to include in the sample
#' @param col.highlight if highlight.sample=TRUE, the colour used to plot the highlighted ancestry
#'
#' @author Suzie Brown
#'
#' @export WFplot
#'

WFplot <- function(history, col.line='darkgrey', shade.dead=!highlight.sample, col.alive='black', col.dead='grey', highlight.sample=FALSE, sampl=1:N, col.highlight='purple'){
  N <- ncol(history)
  N.gen <- nrow(history)+1

  x <- rep(1:N, N.gen)
  y <- rep(1:N.gen, each=N)

  plot(NA, xlim=c(1,N), ylim=c(1,N.gen), bty='n',xaxt='n', yaxt='n', xlab='<-- individuals -->', ylab='<-- generations', main='Wright-Fisher model')

  for (t in 2:N.gen){
   for (i in 1:N){
     lines(c(i,history[t-1,i]),c(N.gen+1-t,N.gen+1-(t-1)), lwd=4, col=col.line)
   }
  }

  ## find which individuals have no offspring and shade them appropriately:
  if (shade.dead){
    dead <- matrix(NA, ncol=N, nrow=N.gen)
    deadcol <- matrix(NA, ncol=N, nrow=N.gen)
    for (t in 1:(N.gen-1)){
      for (i in 1:N){
        dead[t,i] <- !(i %in% history[t,])
      }
      ## convert T/F to colour for plotting
      deadcol[t,dead[t,]] <- col.dead
      deadcol[t,!dead[t,]] <- col.alive
    }
    dead[N.gen,] <- rep(0,N)
    deadcol[N.gen,]<- rep(col.alive, N)

    deadcol <- t(apply(deadcol,2,rev)) ## reverse generations for plotting

    points(x,y, pch=21, col=1, bg=deadcol, cex=4)
  }
  ## or colour all individuals as if alive:
  else {
    points(x,y, pch=21, col=1, bg=col.alive, cex=4)
  }

  ## highlight the ancestry of a sample:
  if (highlight.sample){
    ii <- 1:N
    hl.history <- matrix(NA, ncol=N, nrow=N.gen-1) ## matrix of parent-child links
    hl.nodes <- matrix(NA, ncol=N, nrow=N.gen) ## matrix of T/F: individual active in sample ancestry?

    hl.nodes[N.gen,] <- 1:N %in% sampl
    points(ii[hl.nodes[N.gen,]],rep(1, sum(hl.nodes[N.gen,])), pch=21, col=1, bg=col.highlight, cex=4)

    for (t in (N.gen-1):1){
      for (i in ii[hl.nodes[t+1,]]){
        hl.history[t,i] <- history[t,i]

        lines(c(i, hl.history[t,i]), c(N.gen+1-(t+1),N.gen+1-t), col=col.highlight, lwd=4)
      }

      hl.nodes[t,] <- ii %in% hl.history[t,]
      points(ii[hl.nodes[t+1,]], rep(N.gen+1-(t+1), sum(hl.nodes[t+1,])), pch=21, col=1, bg=col.highlight, cex=4)
    }
    points(ii[hl.nodes[1,]],rep(N.gen, sum(hl.nodes[1,])), pch=21, col=1, bg=col.highlight, cex=4)
  }
}

