# Lin and Ying (1994)'s estimating equation: in individual level details
AH_EE <- function(yobs, delta, z, theta){
  N <- length(yobs)
  p <- length(theta)
  y.sort <- sort( yobs )
  y.rank <- rank( yobs, ties.method = 'min')
  Z.bar.sort <- array(0, dim=c(N, p))
  for( j in 1:N){
    Y <- (yobs >= y.sort[j])
    Z.bar.sort[j,] <- apply( z * Y, 2, sum ) / sum(Y)
  }
  y.sort.diff <- c(y.sort[1], diff( y.sort ) )
  # calculate phi
  phi.all <- array(0, dim=c(N, p))
  for( i in 1:N ){
    Ki <- y.rank[i]
    Zi.aug <- matrix(rep(z[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Zi.aug - Z.bar.sort[1:Ki,]
    di <- y.sort.diff[1:Ki]
    I2i <- t(Ri) %*% di %*% z[i,]
    I1i <- ( z[i,] - Z.bar.sort[y.rank[i],] ) * delta[i]
    phi.all[i,] <- I1i - I2i %*% theta
  }
  return(phi.all)
}


