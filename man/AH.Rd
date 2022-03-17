\name{AH}
\alias{AH}
\title{Fit Additive Hazard Model
}
\description{This is a function used to fit the additive hazard model based on the estimation procedure proposed by Lin and Ying (1994).
}
\usage{
AH(yobs,delta,X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yobs}{observed failure time.
  }
  \item{delta}{censoring indicator (right-censored).
  }
  \item{X}{covariates (should be a matrix and have colnames).
  }
}
\value{
  \item{coef}{for estimated coefficients
  }
  \item{var_cov}{asymptotic var-cov matrix coefficient estimates
  }
}
\references{
Lin, D. Y. and Ying, Z. (1994). Semiparametric analysis of the additive risk model. Biometrika, 81(1):61–71.
}
\examples{
data(simdata)

# get observed survival time and censoring indicator
yobs  <- simdata[,1]
delta <- simdata[,2]
Z <- as.matrix( simdata[,3:6] )

#fit the modle
AH(yobs,delta,Z)
}
