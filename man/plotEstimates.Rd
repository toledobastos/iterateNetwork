\name{plotEstimates}
\alias{plotEstimates}
\title{
Plot list of network attributes generate by iterateComplete
}
\description{
Description: The function plotEstimates plots the estiamtes in a list of network attributes generate by iterateComplete.
}
\usage{
plotEstimates(list.complete, return.estimates, attribute, col, plot.type)
}
\arguments{
  \item{list.complete}{
List of network estimates.
}
  \item{return.estimates}{
Selection of estiamtes to be plotted. Defaults to ALL.
}
  \item{attribute}{
Attribute used to generate the estimates.
}
  \item{col}{
Colors to be asigned to each of the network attributes.
}
  \item{plot.type}{
Type of plot. Accepts p for points or l for lines.
}
}
\author{
Marco T. Bastos
}
\examples{
# Not run:
plotEstimates(list.complete, return.estimates = c(1:5, 9:15, 17), attribute = "org_type", plot.type = "p")
}
\keyword{ iteration }
\keyword{ network }
\keyword{ sample }
\keyword{ plot }
