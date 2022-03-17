#==========================================================================#
# Fit the AH model
#==========================================================================#

AH <- function(yobs,delta,X){
  N <- length(yobs)
  p <- ncol(X)
  y.sort <- sort( yobs )
  y.rank <- rank( yobs, ties.method = 'min')
  X.bar.sort <- array(0, dim=c(N, p))
  for( j in 1:N){
    Y <- (yobs >= y.sort[j])
    X.bar.sort[j,] <- apply( X * Y, 2, sum ) / sum(Y)
  }
  y.sort.diff <- diff( c(0,y.sort) )
  # calculate A, B and d
  A0 <- B0 <- array(0, dim=c(p, p))
  d0 <- rep(0,p)
  for( i in 1:N ){
    Ki <- y.rank[i]
    Xi.aug <- matrix(rep(X[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Xi.aug - X.bar.sort[1:Ki,]
    di <- y.sort.diff[1:Ki]
    I2i <- t(Ri) %*% di %*% X[i,]
    I1i <- ( X[i,] - X.bar.sort[y.rank[i],] ) * delta[i]
    A0 <- A0 + I2i
    d0 <- d0 + I1i
    B0 <- B0 + delta[i]*( Ri[Ki,]%*%t(Ri[Ki,]) )
  }
  A <- A0/N; B <- B0/N; d <- d0/N
  # calculate the estimate of beta and SE
  Est    <- solve(A,d)
  Sigma  <- solve(A) %*% B %*% solve(A) # asymptotic var-cov matrix
  SE     <- sqrt( diag(Sigma)/N )
  zvalue <- Est/SE
  pvalue <- 2*(1-pnorm(abs(zvalue)))
  coef <- data.frame(Est=Est, SE=SE, zvalue=zvalue, pvalue=pvalue,
                     row.names=colnames(X))
  # output
  out <- list(coef = coef,
              var_cov = Sigma/N
              )
  return(out)
}
