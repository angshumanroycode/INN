\name{inns.test}
\alias{inns.test}
\title{
Test of independence based on spatial ranks of nearest neighbors.
}
\description{
Two tests of independence among several random vectors of arbitrary dimension. The null hypothesis is that the random vectors are jointly independent. These two tests are based on spatial ranks of nearest neighbors.
}
\usage{
inns.test(Xlist = NULL , Dlist = NULL , B = 100, alpha = 0.05)
}
\arguments{
  \item{Xlist}{
A list matrices where each matrix corresponds to each random vector. Each row of a particular matrix is an observation on the corresponding random vector.
}
  \item{Dlist}{
A list of "dist" objects where each distance object corresponds to each random vector. It is required if \code{Xlist} is \code{NULL}.
}
  \item{B}{
An integer value specifying the number of Monte-Carlo iterations made in the permutation test.
}
  \item{alpha}{
A numeric value between 0 and 1 specifying the level of the test.
}
}
\details{
See reference.
}
\value{
A list of the following components:
  \item{TSsum.stat}{
The value of the test statistic TSsum.
}
  \item{TSsum.max}{
The value of the test statistic TSmax.
}
  \item{TSsum.cutoff}{
The critical value for the test statistic TSsum.
}
  \item{TSmax.cutoff}{
The critical value for the test statistic TSmax.
}
  \item{TSsum.pvalue}{
P-value of the hypothesis test based on TSsum.
}
  \item{TSmax.pvalue}{
P-value of the hypothesis test based on TSmax.
}
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
set.seed(0)
inns.test(xlist,B=100,alpha=0.05)

#
set.seed(0)
x <- rnorm(10)
y <- rnorm(10)
z <- matrix(x+y+rnorm(10,sd=0.1),ncol=1)
xlist <- list(x,y,z)
dlist <- lapply(xlist,dist)
set.seed(0)
inns.test(Dlist=dlist,B=100,alpha=0.05)
}
\keyword{nonparametric}
\keyword{htest}
