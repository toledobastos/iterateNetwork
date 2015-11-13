\name{iterateComplete}
\alias{iterateComplete}
\title{
Complete iteration with node attributes, random removal, and degree removal.
}
\description{
Description: The function iterateComplete takes an igraph or network object and proceeds to complete node removal based on node attributes, random removal, and degree removal.
}
\usage{
iterateComplete(net.object, attribute="color", net.iterate = "min", plot.estimators = TRUE)
}
\arguments{
  \item{net.object}{
The network object. Accepts both igraph and network objects.
}
  \item{attribute}{
Required for running iterations over vertex/node attribute.
}
  \item{net.iterate}{
Defaults to low and defines the number of interations. Also accepts mid for half the potential iterations and max for all possible iterations and permutations. 
}
  \item{collapse.minor}{
Defaults to 1 and defines the minimum number of nodes per subgroup in the selected attribute. Subgroups with fewer nodes than defined by collapse.minor will be grouped together with the label ETAL. 
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
Marco T. Bastos
}

\examples{
# Not run:
net.object <- barabasi.game(50)
V(net.object)$color <- ifelse(V(net.object) %% 2==1, "darkblue", "darkred") # create random node attribute
iterateComplete(net.object, attribute="color", net.iterate = "min", plot.estimators = TRUE)
}
\keyword{ iteration }
\keyword{ network }
\keyword{ sample }
