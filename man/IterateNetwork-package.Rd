\name{iterateNetwork-package}
\alias{iterateNetwork-package}
\alias{iterateNetwork}
\docType{package}
\title{
Simulation and iteration over network metrics
}
\description{
iterateNetwork takes an igraph or network object and calculates network metrics over subsamples of the network.
}
\details{
\tabular{ll}{
Package: \tab iterateNetwork\cr
Type: \tab Package\cr
Version: \tab 0.1.8\cr
Date: \tab 2015-10-30\cr
License: \tab GPL (>= 2)\cr
}
iterateNetwork takes an igraph or network object and calculates network metrics over subsamples of the network. Simulations can be run over random samples, but also degree, betweenness, and closeness distribution. Iterative sampling can also be defined by node and edge attributes or by node intearaction.
}
\author{
Marco T. Bastos
Maintainer: Marco T. Bastos <marco@toledobastos.com>
}
\keyword{ package }
\examples{
iterateNetwork(barabasi.game(50))
}
