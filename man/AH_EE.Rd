\name{AH_EE}
\alias{AH_EE}
\title{Estimating equation of Lin and Ying(1994)
}
\description{This is a function used to calculate the details of estimating equation (Lin and Ying, 1994) based on the additive hazard model at the individual level
}
\usage{
AH_EE(yobs, delta, z, theta)
}
\arguments{
  \item{yobs}{observed failure time.
  }
  \item{delta}{censoring indicator (right-censored).
  }
  \item{z}{covariates (should be a matrix and have colnames).
  }
  \item{theta}{coefficients of covariates based on additive hazard model.
  }
}

\references{
Lin, D. Y. and Ying, Z. (1994). Semiparametric analysis of the additive risk model. Biometrika, 81(1):61–71.
}
