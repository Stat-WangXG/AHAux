\name{J_test}
\alias{J_test}
\title{Sargan-Hansen-type-J-test
}
\description{The Sargan-Hansen-type-J-test is to test the over-identified restrictions. And the null hypothesis H0 is that the auxiliary information is consistent with the individual-level data.
}
\usage{
J_test(yobs,delta,Z,theta0,beta.star)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yobs}{observed failure time.
  }
  \item{delta}{censoring indicator (right-censored).
  }
  \item{Z}{covariates of full model (should be a matrix and have colnames).
  }
  \item{theta0}{GMM estimator.
  }
  \item{beta.star}{auxiliary information that is the population_level covariate effect information coming from external sources.
  }
}
\value{
  \item{J_statistic}{the value of J statistic.
  }
  \item{p_value}{p-value for the test.
  }
}

\examples{
data(simdata)

# get observed survival time and censoring indicator
yobs  <- simdata[,1]
delta <- simdata[,2]
Z <- as.matrix( simdata[,3:6] )

# auxiliary information
beta.star <- c(0.47,0.76)

#get the GMM estimator of theta
theta.hat <- AH_Aux_GMM(yobs,delta,Z,beta.star,maxit=100, rounding.digits = 6)$coef_IGMM

#J-test
J_test(yobs,delta,Z,theta.hat,beta.star)
}
