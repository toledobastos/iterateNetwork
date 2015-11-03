\name{iterateNetwork}
\alias{iterateNetwork}
\title{
Simulation and iteration over network metrics
}
\description{
Description: The function iterateNetwork takes an igraph or network object and calculates network metrics over subsamples of the network. The argument net.samples defaults to a numeric vector of length 100 giving the fractions of the network from 1% to 100%, but it can be set to any other fraction higher or lower than 1:100. The argument iteration.type defines the type of iteration and defaults to random, but it also accepts degree, betweenness, closeness, and attribute. If iteration.type='attribute' you have to define which vertex attribute to run the iteration with attribute=vertex.attribute.name. The argument net.iterate defaults to 10 and specifies the number of iterations to calculate. Lastly, the argument plot.estimators defaults to TRUE to include a plot of the estimates generated during the iterations.
}
\usage{
net.object <- barabasi.game(50)
iterateNetwork(net.object, net.iterate = 10, plot.estimators = TRUE)
}
\arguments{
  \item{net.object}{
The network object. Accepts both igraph and network objects
}
  \item{net.samples}{
Defaults to a numeric vector of length 100 giving the fractions of the network from 1% to 100%. It can be set to any other fraction higher or lower than 1:100.
}
  \item{net.iterate}{
Defaults to 10 and specifies the number of iterations to calculate.
}
  \item{iteration.type}{
Defaults to random iteration. It also accepts degree, betweenness, closeness, and attribute. If iteration.type='attribute' you have to define which vertex attribute to run the iteration with attribute=vertex.attribute.name.
}
  \item{attribute}{
Required for running iterations over vertex attribute. Defaults to NULL.
}
  \item{plot.estimators}{
Defaults to TRUE to include a plot of the estimates generated during the iterations.
}
  \item{return.estimates}{
Defines which network metrics to iterate thru. It defaults to a selection of 15 metrics but can be set to ALL to return all metrics available (see below). You can manually select or remove any of the following metrics: 
# 01. nodes
# 02. edges
# 03. centralization
# 04. diameter
# 05. eigenvector
# 06. permutation
# 07. transitivity
# 08. articulations
# 09. cluster
# 10. path.length
# 11. degree
# 12. closeness
# 13. page.rank
# 14. betweenness
# 15. density
# 16. largest.component
# 17. small.world
}
}
\author{
Marco T. Batos
}

\examples{
# Not run:
net.object <- barabasi.game(50) # generate sample network
iterateNetwork(net.object)
}
\keyword{ iteration }
\keyword{ network }
\keyword{ sample }
