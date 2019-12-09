\name{wild.boot}
\alias{wild.boot}

\title{Wild Bootstrap in Linear Models
}



\description{
This function performs the wild/external bootstrap in linear models as described by Wu (1986) <doi:10.1214/aos/1176350142>.  
Linear models incorporating categorical and/or quantitative predictor variables with a quantitative response are allowed.  
The function output creates the bootstrap sampling distribution for each coefficient.  
Estimation is performed via least squares.
}



\usage{
wild.boot(formula, B = 1000, data = NULL, seed = NULL, bootDistn = "normal")
}



\arguments{
  \item{formula}{
  input a linear model formula of the form \code{response}~\code{predictors} as you would in the \code{lm()} 
  function.  All variables must contain non-missing entries.
  }
  \item{B}{
  number of bootstrap samples.  This should be a large, positive integer value.
  }
  \item{data}{
  optionally, input the name of the dataset where variables appearing in the model are stored.
  }
  \item{seed}{
  optionally, set a value for the seed for the bootstrap sample generation.  The default \code{NULL} will
  pick a random value for the seed.
  }
  \item{bootDistn}{
  distribution used to create the wild bootstrap weights for the residuals.  Allowed distributions include
  \code{"normal"}, \code{"uniform"}, \code{"exponential"}, \code{"laplace"}, \code{"lognormal"}, 
  \code{"gumbel"}, \code{"t5"}, \code{"t8"}, and \code{"t14"}.  The numbers after the t-distributions
  indicate the degrees of freedom.  Any selected distribution creates weights with mean 0 and variance 1 from
  the named distribution.
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
  \item{bootDistn}{type of distribution used to generate the wild bootstrap weights for the residuals}
}



\references{
Wu, C.F.J.  (1986). "Jackknife, Bootstrap, and Other Resampling Methods in Regression Analysis."  \emph{Annals of Statistics}.  Vol. 14, No. 4, pp.1261 - 1295.
}



\author{
  Megan Heyman, heyman@rose-hulman.edu
}




\examples{
Seed <- 14
set.seed(Seed)
y <- rnorm(20) #randomly generated response
x <- rnorm(20) #randomly generated predictor
WildObj <- wild.boot(y~x, B=100, seed=Seed) #perform the wild bootstrap

#plot the sampling distribution of the slope coefficient
hist(WildObj$bootEstParam[,2], main="Wild Bootstrap Sampling Distn.",
     xlab="Slope Estimate") 

#bootstrap 95% CI for slope parameter (percentile method)
quantile(WildObj$bootEstParam[,2], probs=c(.025, .975))
}


\keyword{ distribution }
\keyword{ models regression }
\keyword{ nonparametric }
