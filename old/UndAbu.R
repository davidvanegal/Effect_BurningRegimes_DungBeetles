#' UndAbu(x) is a function of estimating undetected species relative abundance.
#' @param x a vector of species abundance frequency
#' @return a numerical vector
UndAbu <- function(x){
  x <- unlist(x)
  n <- sum(x)  
  f1 <- sum(x==1)
  f2 <- sum(x==2)
  f3 <- sum(x==3)
  f4 <- max(sum(x == 4), 1)
  f0.hat <- ceiling(ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2))  # Chao1 estimator
  if(f2==0){
    f1 <- max(f1 - 1, 0)
    f2 <- 1
  }
  A1 <- f1 / n * ((n-1)*f1 / ((n-1)*f1 + 2*max(f2,1)))
  A2 <- f2 / choose(n, 2) * ((n-2)*f2 / ((n-2)*f2 + 3*max(f3,1)))^2
  R <- A1^2/A2
  j <- 1:f0.hat
  f.solve <- function(x){ 
    out <- sum(x^j)^2 / sum((x^j)^2) - R
    abs(out)
  }
  b <-  tryCatch(optimize(f.solve, lower=(R-1)/(R+1), upper=1, tol=1e-5)$min, error = function(e) {(R-1)/(R+1)})
  a <- A1 / sum(b^j)
  p <- a * b^j
  if(f0.hat == 1) p <- A1
  p
}
