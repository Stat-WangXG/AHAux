AH_Aux_GMM <- function(yobs,delta,Z,beta.star,
                            maxit, rounding.digits = 6){
  # Basic elements
  N <- length(yobs)
  p <- ncol(Z)
  N0 <- sum(delta)
  W<- Z[,1:length(beta.star)]
  # Get the initial value for theta
  ahfit <- AH(yobs, delta, Z)
  theta.LY <- ahfit$coef[,1]
  theta0 <- theta.LY

  # return A and d matrix of Additive hazard model
  y.sort <- sort( yobs )
  y.rank <- rank( yobs, ties.method = 'min')
  Z.bar.sort <- array(0, dim=c(N, p))
  for( j in 1:N){
    Y <- (yobs >= y.sort[j])
    Z.bar.sort[j,] <- apply( Z * Y, 2, sum ) / sum(Y)
  }
  y.sort.diff <- diff( c(0,y.sort) )
  # calculate A and d
  A0 <- array(0, dim=c(p, p))
  d0 <- rep(0,p)
  for( i in 1:N ){ # i
    Ki <- y.rank[i]
    Zi.aug <- matrix(rep(Z[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Zi.aug - Z.bar.sort[1:Ki,]
    di <- y.sort.diff[1:Ki]
    I2i <- t(Ri) %*% di %*% Z[i,]
    I1i <- ( Z[i,] - Z.bar.sort[y.rank[i],] ) * delta[i]
    A0 <- A0 + I2i
    d0 <- d0 + I1i
  }
  A0 <- A0/N; d0 <- d0/N
  Psin <- NULL
  Psin <- cbind(Psin, AH_EE(yobs,delta,W,beta.star))
  Psii <- apply(Psin,2,mean)
  d <- c(d0,Psii)
  A <- rbind(A0,matrix(0,nrow=length(Psii),ncol=p))

  ### circle to get optimal W ###
  bet.all <- theta0
  tol <- 1e-6
  bet.old <- theta0
  numit <- 0
  repeat{
    # specify the new W (The Asymptotic Covariance Matrix)
    Sigma.hat <- gmm_weight(yobs,delta,Z,bet.old,beta.star)$sigma
    #### Solve using explict formula
    bet <- as.vector( solve(t(A)%*%solve(Sigma.hat)%*%A) %*% t(A) %*% solve(Sigma.hat) %*% d )

    ### whether to stop
    if(sqrt(mean((bet.old-bet)^2)) > tol & numit < maxit){
      bet.old <- bet
      bet.all <- rbind(bet.all,bet)
      numit <- numit + 1
    }else{
      bet.all <- rbind(bet.all,bet)
      break
    }
  }
  # SE of GMM without iteration
  weight1 <- gmm_weight(yobs,delta,Z,bet.all[2,],beta.star)
  gamma1 <- weight1$gamma
  b1 <- weight1$b
  omiga1 <- weight1$omiga
  v1 <- solve(A0)%*%(gamma1 - b1%*%solve(omiga1)%*%t(b1))%*%solve(A0)
  s1 <- sqrt(diag(v1)/N)
  # SE of GMM with iteration
  weight2 <- gmm_weight(yobs,delta,Z,bet.all[nrow(bet.all),],beta.star)
  gamma2 <- weight2$gamma
  b2 <- weight2$b
  omiga2 <- weight2$omiga
  v2 <- solve(A0)%*%(gamma2 - b2%*%solve(omiga2)%*%t(b2))%*%solve(A0)
  s2 <- sqrt(diag(v2)/N)

  list(coef_GMM = round(bet.all[2,],rounding.digits),
       coef_IGMM = round(bet.all[nrow(bet.all),],rounding.digits),
       SE_GMM = round(s1,rounding.digits),
       SE_IGMM = round(s2,rounding.digits)
       )
}

