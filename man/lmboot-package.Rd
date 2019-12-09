\name{lmboot-package}
\alias{lmboot-package}
\alias{lmboot}
\docType{package}

\title{
\packageTitle{lmboot}
}

\description{
\packageDescription{lmboot}
}

\details{
\packageDESCRIPTION{lmboot}

\packageIndices{lmboot}

This package is useful to users who wish to perform bootstrap in linear models.  The package contains functions to create
the sampling distributions for linear model parameters using either efficient or robust bootstrap methods.  

As classified by
Liu and Singh (1992), efficient bootstrap types include the residual bootstrap (\code{residual.boot()}).  These types of
bootstrap are useful when it is not reasonable to assume that errors come from a normal distribution, but you may make other
classical assumptions:  errors are independent, have mean 0, and have constant variance.

Robust bootstrap types include the paired bootstrap (\code{paired.boot}), wild bootstrap (\code{wild.boot}), and the jackknife (\code{jackknife}).
These types of bootstrap are useful when it is not reasonable to assumet that errors have constant variance, but you may make other
classical assumptions:  errors are independent and have mean 0.

The package also contains a function for Bayesian bootstrap (\code{bayesian.boot} and a function to perform bootstrap in the
ANOVA hypothesis test (\code{ANOVA.boot}).  The ANOVA bootstrap function has options to use the wild or residual bootstrap techniques
and has been tested to work in 2-way ANOVA.  Its functionality allows K-way ANOVA, however those capabilities have not been fully tested.

Currently, the user must manipulate the output of the function to conduct hypothesis tests and create 
  confidence intervals for the predictor coefficients.  More convenient/streamlined output is expected
  in future package versions.
}


\author{
\packageAuthor{lmboot}

Maintainer: \packageMaintainer{lmboot}
}


\references{
Efron, B.  (1979). "Bootstrap methods: Another look at the jackknife."  \emph{Annals of Statistics}.  Vol. 7, pp.1-26.

Liu, R. Y. and Singh, K. (1992). "Efficiency and Robustness in Resampling."  \emph{Annals of Statistics}. Vol. 20, No. 1, pp.370-384.

Rubin, D. B.  (1981). "The Bayesian Bootstrap."  \emph{Annals of Statistics}.  Vol. 9, No. 1,  pp.130-134.

Wu, C.F.J.  (1986). "Jackknife, Bootstrap, and Other Resampling Methods in Regression Analysis."  \emph{Annals of Statistics}.  Vol. 14, No. 4, pp.1261 - 1295.

}


\keyword{ package }

\examples{
Seed <- 14
set.seed(Seed)
y <- rnorm(20) #randomly generated response
x <- rnorm(20) #randomly generated predictor

ResidObj <- residual.boot(y~x, B=100, seed=Seed) #perform the residual bootstrap
WildObj <- wild.boot(y~x, B=100, seed=Seed) #perform the wild bootstrap

#residual bootstrap 95% CI for slope parameter (percentile method)
quantile(ResidObj$bootEstParam[,2], probs=c(.025, .975))

#bootstrap 95% CI for slope parameter (percentile method)
quantile(WildObj$bootEstParam[,2], probs=c(.025, .975))
}
