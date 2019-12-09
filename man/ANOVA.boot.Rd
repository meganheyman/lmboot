\name{ANOVA.boot}
\alias{ANOVA.boot}


\title{
Residual and wild bootstrap in 1-way and 2-way ANOVA
}


\description{
This function performs the residual bootstrap as described by Efron (1979) and wild bootstrap as described by Wu (1986)
for ANOVA hypothesis testing.  Linear models 
incorporating categorical and/or quantitative predictor variables with a quantitative response are allowed.  
The function output creates the bootstrap null distribution for each term to be tested.  
Estimation is performed via least squares and only Type I sum of squares are calculated.
}


\usage{
ANOVA.boot(formula, B = 1000, type = "residual", wild.dist = "normal", 
            seed = NULL, data = NULL, keep.boot.resp = FALSE)
}


\arguments{
   \item{formula}{
  input a linear model formula of the form \code{response}~\code{predictors} as you would in the \code{lm()} 
  function.  All variables must contain non-missing entries.
  }

  \item{B}{
  number of bootstrap samples.  This should be a large, positive integer value.
  }

  \item{type}{
  type of bootstrap to perform.  Select either "residual" for residual bootstrap or "wild" for wild bootstrap.
  }

  \item{wild.dist}{
  distribution used to create the wild bootstrap weights for the residuals.  Allowed distributions include
  \code{"normal"}, \code{"uniform"}, \code{"exponential"}, \code{"laplace"}, \code{"lognormal"}, 
  \code{"gumbel"}, \code{"t5"}, \code{"t8"}, and \code{"t14"}.  The numbers after the t-distributions
  indicate the degrees of freedom.  Any selected distribution creates weights with mean 0 and variance 1 from
  the named distribution.
  }

  \item{seed}{
  optionally, set a value for the seed for the bootstrap sample generation.  The default \code{NULL} will
  pick a random value for the seed.
  }

  \item{data}{
  optionally, input the name of the dataset where variables appearing in the model are stored.
  }

  \item{keep.boot.resp}{
  a boolean indicating whether the list of returns includes raw bootstrap responses.  Setting this to TRUE
  may not be possible for larger datasets or too many bootstrap samples due to memory usage.
  }
}



\details{
  Currently, the user must manipulate the output of the function manually to view the bootstrap ANOVA table components
  and visualize the null distribution.  More convenient/streamlined output is expected
  in future package versions.

  Thanks to Bochuan Lyu who helped to coding to this function.
}


\value{
  \item{terms}{names of the terms/rows of the ANOVA table.  These correspond to each predictor variable input to the formula.}

  \item{df}{degrees of freedom associated with each term/row in the ANOVA table.  These correspond to the number of categories in
          each predictor variable (or are 1 for quantitative predictors)}

  \item{origFStats}{original F-statistic value.  Same value as obtained by \code{aov()} using type I sum of squares.}

  \item{origSSE}{original sum of squares, error.  Same value as obtained by \code{aov()} using type I sum of squares.}

  \item{origSSTr}{original sum of squares, treatment.  Vector containing the sum of squares for each term in the ANOVA model.
                These are the same values as obtained by \code{aov()} using type I sum of squares.}

  \item{bootFStats}{matrix containing the bootstrap F statistics.  Each column corresponds to a term in the ANOVA table.  There
                  are \code{B} rows.}

  \item{bootSSE}{matrix containing the bootstrap sum of squares, error.  Each column corresponds to a term in the ANOVA table.  There
               are \code{B} rows.  These are calculated using type I sum of squares.}

  \item{bootSSTr}{matrix containing the bootstrap sum of squares, treatment.  Each column corresponds to a term in the ANOVA table.  There
               are \code{B} rows.  These are calculated using type I sum of squares.}

  \item{`p-values`}{vector containing the bootstrap p-values for each predictor term in the ANOVA model.  These are calculated by
                  counting the number of bootstrap test statistics which are greater than the original observed test statistic and
                  dividing by \code{B} }
}


\references{
  Efron, B.  (1979). "Bootstrap methods: Another look at the jackknife."  \emph{Annals of Statistics}.  Vol. 7, pp.1-26.

  Wu, C.F.J.  (1986). "Jackknife, Bootstrap, and Other Resampling Methods in Regression Analysis."  \emph{Annals of Statistics}.  Vol. 14, No. 4, pp.1261 - 1295.
}


\author{
  Megan Heyman, heyman@rose-hulman.edu
}



\seealso{
  \code{\link{wild.boot}}, \code{\link{residual.boot}}
}


\examples{

data(mtcars)         #load an example dataset
myANOVA2 <- ANOVA.boot(mpg~as.factor(cyl)*as.factor(am), data=mtcars)
myANOVA2$`p-values`  #bootstrap p-values for 2-way interactions model

myANOVA1 <- ANOVA.boot(mpg~as.factor(cyl), data=mtcars)
myANOVA1$`p-values` #bootstrap p-values for 1-way model

myANOVA2a <- ANOVA.boot(mpg~as.factor(cyl)+as.factor(am), data=mtcars)
myANOVA2a$`p-values` #bootstrap p-values for 1-way additive model

}


\keyword{ models regression }
\keyword{ nonparametric }
