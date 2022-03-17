\name{Calibration}
\alias{Calibration}
\title{Calibration method
}
\description{When the auxiliary survival information is not consistent with the individual-level data, we can use the calibration method to get new estimates, thereby alleviating the influence of the heterogeneity between internal and external data on the estimation. It is a convex combination of LY estimates and GMM estimates.
}
\usage{
Calibration(yobs,delta,Z,beta.star)
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
}
\value{Returns the parameter estimates for the additive hazard model based on the calibration method.
}

\examples{
data(simdata)

# get observed survival time and censoring indicator
yobs  <- simdata[,1]
delta <- simdata[,2]
Z <- as.matrix( simdata[,3:6] )

# auxiliary information
beta.star <- c(0.47,0.76)

#Calibration method
Calibration(yobs,delta,Z,beta.star)
}
