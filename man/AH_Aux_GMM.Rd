\name{AH_Aux_GMM}
\alias{AH_Aux_GMM}
\title{Fit additive hazard model with auxiliary covariate effect information  using GMM
}
\description{This is a function used to fit additive hazard model with auxiliary covariate effect information based on GMM.
}
\usage{
AH_Aux_GMM(yobs,delta,Z,beta.star,maxit, rounding.digits = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yobs}{observed failure time.
  }
  \item{delta}{censoring indicator (right-censored).
  }
  \item{Z}{covariates (should be a matrix and have colnames).
  }
  \item{beta.star}{the population_level covariate effect information coming from external sources.
  }
  \item{maxit}{number of iterations.
  }
  \item{rounding.digits}{the number of decimal points that you wish to have for the number.
  }
}

\value{
  \item{coef_GMM }{GMM estimates without iteration}
  \item{coef_IGMM }{GMM estimates with iteration}
  \item{SE_GMM }{asymptotic standard error for GMM estimates without iteration}
  \item{SE_IGMM }{asymptotic standard error for GMM estimates with iteration}
}

\examples{
data(simdata)

# get observed survival time and censoring indicator
yobs  <- simdata[,1]
delta <- simdata[,2]
Z <- as.matrix( simdata[,3:6] )

# auxiliary information
beta.star <- c(0.47,0.76)

#fit the modle
AH_Aux_GMM(yobs,delta,Z,beta.star,maxit=100, rounding.digits = 6)
}

