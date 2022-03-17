#==========================================================================#
# over-identified test (Sargan-Hansen-type-J-test)
#==========================================================================#
J_test <- function(yobs,delta,Z,theta0,beta.star){
  #theta0 : IGMM estimator
  #beta.star : auxiliary information
  #Z : covariance of full model
  N <- length(yobs)
  W <- Z[,1:length(beta.star)]
  Sigma.hat <- gmm_weight(yobs,delta,Z,theta0,beta.star)$sigma
  U <- c(apply(AH_EE(yobs, delta, Z, theta0),2,mean),
         apply(AH_EE(yobs, delta, W, beta.star),2,mean))
  J <- N * t(U) %*% solve(Sigma.hat) %*% U #J test statistic
  p <- 1 - pchisq(J,2) #pvalue
  out <- list(J_statistic = J, pvalue = p)
  return(out)
}
