\name{jackknife}
\alias{jackknife}


\title{
Delete-1 Jackknife in Linear Models
}


\description{
This function performs the delete-1 jackknife in linear models as described by Quenouille (1956) <doi:10.2307/2332914>.  
Linear models incorporating categorical and/or quantitative predictor variables with a quantitative response are allowed.  
The function output creates the jackknife sampling distribution for each coefficient.  
Estimation is performed via least squares.
}


\usage{
jackknife(formula, data = NULL)
}


\arguments{
  \item{formula}{
  input a linear model formula of the form \code{response}~\code{predictors} as you would in the \code{lm()} 
  function.  All variables must contain non-missing entries.
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
  \item{bootEstParam}{matrix containing the jackknife parameter estimates.  Each column corresponds to a 
                      coefficient.  There are \code{n-1} rows, each corresponding to a jackknife sample.}
  \item{origEstParam}{vector containing the least squares parameter estimates.  These are the same as
                      estimates obtained from \code{lm}.}
}


\references{
Quenouille, M. (1956). "Notes on bias in estimation."  \emph{Biometrika}.  Vol. 61,  pp.1-15
}



\author{
  Megan Heyman, heyman@rose-hulman.edu
}




\examples{
Seed <- 14
set.seed(Seed)
y <- rnorm(20) #randomly generated response
x <- rnorm(20) #randomly generated predictor
JackObj <- jackknife(y~x) #perform the jackknife

#plot the sampling distribution of the slope coefficient
hist(JackObj$bootEstParam[,2], main="Jackknife Sampling Distn.",
     xlab="Slope Estimate") 

#jackknife 95% CI for slope parameter (percentile method)
quantile(JackObj$bootEstParam[,2], probs=c(.025, .975))
}



\keyword{ models regression }
\keyword{ nonparametric }
