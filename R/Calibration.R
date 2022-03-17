Calibration <- function(yobs,delta,Z,beta.star){
  LY <- AH(yobs,delta,Z)
  thetaEst.LY <- LY$coef[,1]
  GMM <- AH_Aux_GMM(yobs,delta,Z,beta.star,100, rounding.digits = 6)
  thetaEst.IGMM <- GMM$coef_IGMM
  sigma1 <- LY$var_cov
  reduce.ah <- AH(yobs, delta, Z[,1:length(beta.star)])
  betaEst.LY <- reduce.ah$coef[,1]
  IESD <- (betaEst.LY[1:length(beta.star)] - beta.star)%*%t(betaEst.LY[1:length(beta.star)] - beta.star)
  theta.delta1 <- AH_Aux_GMM(yobs,delta,Z,beta.star+c(10^(-6),0),maxit = 100, rounding.digits = 6)
  theta.delta2 <- AH_Aux_GMM(yobs,delta,Z,beta.star+c(0,10^(-6)),maxit = 100, rounding.digits = 6)
  Delta <- (cbind(theta.delta1$coef_IGMM,theta.delta2$coef_IGMM) - thetaEst.IGMM)/10^(-6)
  sigma2 <- Delta %*% IESD %*% t(Delta)
  thetaEst.SIGMM <- t(sigma2%*%solve(sigma1+sigma2)%*%thetaEst.LY
                      + sigma1%*%solve(sigma1+sigma2)%*%thetaEst.IGMM)
  return(thetaEst.SIGMM)
}
