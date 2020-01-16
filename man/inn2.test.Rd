\name{inn2.test}
\alias{inn2.test}
\title{
Test of independence based on pair-wise ranks of nearest neighbors.
}
\description{
Two tests of independence among several random vectors of arbitrary dimension. The null hypothesis is that the random vectors are jointly independent. These two tests are pair-wise ranks of nearest neighbors.
}
\usage{
inn2.test(Xlist = NULL , Dlist = NULL , B = 100, alpha = 0.05)
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
  \item{T2sum.stat}{
The value of the test statistic T2sum.
}
  \item{T2sum.max}{
The value of the test statistic T2max.
}
  \item{T2sum.cutoff}{
The critical value for the test statistic T2sum.
}
  \item{T2max.cutoff}{
The critical value for the test statistic T2max.
}
  \item{T2sum.pvalue}{
P-value of the hypothesis test based on T2sum.
}
  \item{T2max.pvalue}{
P-value of the hypothesis test based on T2max.
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
inn2.test(xlist,B=100,alpha=0.05)

#
set.seed(0)
x <- rnorm(10)
y <- rnorm(10)
z <- matrix(x+y+rnorm(10,sd=0.1),ncol=1)
xlist <- list(x,y,z)
dlist <- lapply(xlist,dist)
set.seed(0)
inn2.test(Dlist=dlist,B=100,alpha=0.05)
}
\keyword{nonparametric}
\keyword{htest}
