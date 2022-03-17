\name{simdata}
\alias{simdata}
\docType{data}
\title{a simulated dataset from additive risk model}

\description{Simdata is a simulated dataset from additive hazard model with four covariates. The censoring rate is approximately 25\%.
}
\usage{data(simdata)}
\format{The variables represented in the data set are as follows:
  \describe{
  \item{time}{survival time (month)(time to death or end of study).}
  \item{delta}{death indicator: 1 - Dead, 0 - Alive.}
  \item{X1}{covariate X1 comes from U(0,1) distribution.}
  \item{X2}{covariate X2 comes from b(1,0.5) distribution.}
  \item{X3}{(covariate X3, covariate X4) come from a mean zero bivariate normal distribution with var(X3) = 1,var(X4) = 2 and a correlation coefficient of 0.5 between X3 and X4.}
  \item{X4}{see description of X3.}
  }
}
\keyword{datasets}
