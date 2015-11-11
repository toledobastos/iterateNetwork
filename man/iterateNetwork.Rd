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
  \item{removal}{
Defaults to node and defines whether the iteration should be sampled over nodes, edges, or node.interaction. Removal by edge cannot be performed over node metrics (i.e. degree, betweenness, or closeness), but node.interaction accepts node attributes and will generate edge attributes based on node interactions. 
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
  \item{stepwise.removal}{
Defines the sample of size of targeted nodes to be removed each iteration. Required for running iterations over vertex attribute. Defaults to auto which calculates the ideal fraction of stepwise.removal and net.iterate based on the maximum number of nodes across subgroups of an attribute.
}
  \item{plot.estimators}{
Defaults to TRUE to include a plot of the estimates generated during iterations.
}
}
  \item{plot.type}{
Defines the plot type used in R base graphics. Defaults to 'p' for points, but lines ('l') or lines and points ('o') might be more suitable for targeted node removal.
}
  \item{return.estimates}{
Defines which network metrics to iterate over. It defaults to a selection of 15 metrics but can be set to ALL to return all metrics available (see below). You can manually select or remove any of the following metrics: 
# 01. nodes
# 02. edges
# 03. degree
# 04. eigenvector
# 05. local.clustering
# 06. centralization
# 07. diameter
# 08. permutation
# 09. transitivity
# 10. articulations
# 11. cluster
# 12. path.length
# 13. closeness
# 14. page.rank
# 15. betweenness
# 16. density
# 17. largest.component
# 18. small.world
by setting return.estimates=c(1:10,12:17), with numbers indicating the above variables.
}
}
\author{
Marco T. Batos
}

\examples{
# Not run:
net.object <- barabasi.game(50) # generates sample network
iterateNetwork(net.object) # iterate over network with default arguments
}
\keyword{ iteration }
\keyword{ network }
\keyword{ sample }
