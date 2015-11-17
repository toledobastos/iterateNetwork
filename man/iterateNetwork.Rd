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
  \item{plot.type}{
Defines the plot type used in R base graphics. Defaults to 'p' for points, but lines ('l') or lines and points ('o') might be more suitable for targeted node removal.
}
  \item{return.estimates}{
Defines which network metrics to iterate over. It defaults to a selection of 15 metrics but can be set to ALL to return all metrics available (see below). You can manually select or remove any of the following metrics: 
# 01. nodes
# 02. edges
# 03. degree
# 04. eigenvector
# 05. centralization
# 06. diameter
# 07. permutation
# 08. transitivity
# 09. articulations
# 10. average local clustering coefficient
# 11. number of clusters
# 12. average cluster size
# 13. median cluster size
# 14. average path length
# 15. closeness
# 16. page rank
# 17. betweenness
# 18. density
# 19. largest component
# 20. largest component as a fraction of total nodes in network
# 21. size of second largest component
# 22. size of third largest component
# 23. triangles
# 24. small world quotient
# 25. small world quotient w/o isolates
# 26. small world quotient w/o isolates or isolated clusters
# 27. isolates
# 28. isolates as a fraction of total nodes in network
by setting return.estimates=c(1:10,12:17), with numbers indicating the above variables.
}
}
\author{
Marco T. Bastos
}
\examples{
# Not run:
net.object <- barabasi.game(50) # generates sample network
iterateNetwork(net.object) # iterate over network with default arguments
}
\keyword{ iteration }
\keyword{ network }
\keyword{ sample }
