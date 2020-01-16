\name{inns}
\alias{inns}
\title{
Measure of dependence based on multivariate spatial ranks of nearest neighbors.
}
\description{
It measures the dependence among several random vectors of arbitrary dimensions.
}
\usage{
inns(Xlist = NULL , Dlist = NULL )
}
\arguments{
  \item{Xlist}{
A list matrices where each matrix corresponds to each random vector. Each row of a particular matrix is an observation on the corresponding random vector.
}
  \item{Dlist}{
A list of "dist" objects where each distance object corresponds to each random vector. It is required if \code{Xlist} is \code{NULL}.
}
}
\details{
See reference.
}
\value{
A numeric vector containing sum and max values of the spatial statistics respectively.
}
\references{
Roy et al. (2019)
}
\author{
Angshuman Roy
}
\examples{
set.seed(0)
x <- rnorm(10)
y <- rnorm(10)
z <- matrix(x+y+rnorm(10,sd=0.1),ncol=1)
xlist <- list(x,y,z)
inns(xlist)

##
set.seed(0)
x <- rnorm(10)
y <- rnorm(10)
z <- matrix(x+y+rnorm(10,sd=0.1),ncol=1)
xlist <- list(x,y,z)
dlist <- lapply(xlist,dist)
inns(Dlist=dlist)
}
\keyword{nonparametric}
