#Estimation of the weight matrix Sigma
gmm_weight <- function(yobs,delta,Z,theta.hat,beta.star){
  N <- length(yobs)
  W <- Z[,1:length(beta.star)]
  phi.hat <- AH_EE(yobs,delta,Z,theta.hat)
  phi.bar <- apply(phi.hat,2,mean)
  psi.hat <- AH_EE(yobs,delta,W,beta.star)
  psi.bar <- apply(psi.hat,2,mean)
  #Estimation of gamma
  l <- array(0, dim=c(ncol(Z),ncol(Z)))
  for(i in 1:N){
    k<-(phi.hat[i,]-phi.bar)%*%t(phi.hat[i,]-phi.bar)
    l<-l+k
  }
  gamma <- 1/N*l
  #Estimation of omiga
  l <- array(0, dim=c(ncol(W),2))
  for(i in 1:N){
    k<-(psi.hat[i,]-psi.bar)%*%t(psi.hat[i,]-psi.bar)
    l<-l+k
  }
  omiga <- 1/N*l
  #Estimation of B
  l <- array(0, dim=c(ncol(Z),2))
  for(i in 1:N){
    k<-(phi.hat[i,]-phi.bar)%*%t(psi.hat[i,]-psi.bar)
    l<-l+k
  }
  b <- 1/N*l

  #Estimation of the weight matrix Sigma
  sigma <- rbind(cbind(gamma,b),cbind(t(b),omiga))

  list(gamma = gamma, omiga = omiga, b = b, sigma = sigma)
}
