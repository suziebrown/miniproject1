Stirling2 <- function(n,m)
{
  ## Purpose:  Stirling Numbers of the 2-nd kind
  ## 		S^{(m)}_n = number of ways of partitioning a set of
  ##                      $n$ elements into $m$ non-empty subsets
  ## Author: Martin Maechler, Date:  May 28 1992, 23:42
  ## ----------------------------------------------------------------
  ## Abramowitz/Stegun: 24,1,4 (p. 824-5 ; Table 24.4, p.835)
  ## Closed Form : p.824 "C."
  ## ----------------------------------------------------------------

  if (0 > m || m > n) stop("'m' must be in 0..n !")
  k <- 0:m
  sig <- rep(c(1,-1)*(-1)^m, length= m+1)# 1 for m=0; -1 1 (m=1)
  ## The following gives rounding errors for (25,5) :
  ## r <- sum( sig * k^n /(gamma(k+1)*gamma(m+1-k)) )
  ga <- gamma(k+1)
  round(sum( sig * k^n /(ga * rev(ga))))
}
