\name{bayesian.boot}
\alias{bayesian.boot}



\title{
Bayesian Bootstrap in Linear Models
}



\description{
This function performs the bayesian bootstrap in linear models as described by Rubin (1981) <doi:10.1214/aos/1176345338>.  
Linear models incorporating categorical and/or quantitative predictor variables with a quantitative response are allowed.  
The function output creates the bootstrap sampling distribution for each coefficient.  
Estimation is performed via least squares.
}



\usage{
bayesian.boot(formula, B = 1000, seed = NULL, data = NULL)
}



\arguments{
  \item{formula}{
  input a linear model formula of the form \code{response}~\code{predictors} as you would in the \code{lm()} 
  function.  All variables must contain non-missing entries.
  }
  \item{B}{
  number of bootstrap samples.  This should be a large, positive integer value.
  }
  \item{seed}{
  optionally, set a value for the seed for the bootstrap sample generation.  The default \code{NULL} will
  pick a random value for the seed.
  }
  \item{data}{
  optionally, input the name of the dataset where variables appearing in the model are stored.
  }
}



\details{
  Currently, the user must manipulate the output of the function to conduct hypothesis tests and create 
  confidence intervals for the predictor coefficients.  More convenient/streamlined output is expected
  in future package versions.
}



\value{
  \item{bootEstParam}{matrix containing the bootstrap parameter estimates.  Each column corresponds to a 
                      coefficient.  There are \code{B} rows, each corresponding to a bootstrap sample.}
  \item{origEstParam}{vector containing the least squares parameter estimates.  These are the same as
                      estimates obtained from \code{lm}.}
  \item{seed}{numerical value set for the seed.  This is associated with the set of bootstrap parameter
              estimates and helps the process to be reproducible.}
}



\references{
Rubin, D. B.  (1981). "The Bayesian Bootstrap."  \emph{Annals of Statistics}.  Vol. 9, No. 1,  pp.130-134.
}


\author{
  Megan Heyman, heyman@rose-hulman.edu
}




\examples{
Seed <- 14
set.seed(Seed)
y <- rnorm(20) #randomly generated response
x <- rnorm(20) #randomly generated predictor
BayesObj <- bayesian.boot(y~x, B=100, seed=Seed) #perform the Bayesian bootstrap

#plot the sampling distribution of the slope coefficient
hist(BayesObj$bootEstParam[,2], main="Bayesian Bootstrap Sampling Distn.",
     xlab="Slope Estimate") 

#bootstrap 95% CI for slope parameter (percentile method)
quantile(BayesObj$bootEstParam[,2], probs=c(.025, .975))

}

\keyword{ models regression }
\keyword{ nonparametric }
\keyword{ distribution }
