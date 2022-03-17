\name{gmm_weight}
\alias{gmm_weight}
\title{Estimation of Weight Matrix
}
\description{This is a function used to calculate estimation of weight matrix in objective function of generalized moments method.
}
\usage{
gmm_weight(yobs,delta,Z,theta.hat,beta.star)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yobs}{observed failure time.
  }
  \item{delta}{censoring indicator (right-censored).
  }
  \item{Z}{covariates (should be a matrix and have colnames).
  }
  \item{theta.hat}{estimation of covariate coefficients based on additive hazard model (ncol(Z) = length(theta.hat)).
  }
  \item{beta.star}{the population_level covariate effect information coming from external sources (length(beta.star) <= length(theta.hat)).
  }
}
\value{
  \item{gamma}{estimation of the variance of the estimating equation at the individual level based on additive hazard model.
  }
  \item{omiga}{estimation of the variance of the estimating equation at the individual level based on additive hazard model with auxiliary information.
  }
  \item{b}{estimation of the covariance of two estimating equations based on GMM framework }
  \item{sigma}{estimation of weight matrix in objective function of generalized moments method}
}


