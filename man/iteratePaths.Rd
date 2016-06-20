\name{iteratePaths}
\alias{iteratePaths}
\title{
Calculates the shortest path between two actor types
}
\description{
Description: The function iteratePaths computes the shortest path between two groups of actors by iteratively removing the other groups.
}
\usage{
iteratePaths(net.object, attribute="orgsector", from.node.group="EDU", to.node.group="IGOV", stepwise.removal=1)
}

\arguments{
  \item{net.object}{
The network object. Accepts both igraph and network objects
}
  \item{iteration.type}{
Defaults to "shortest.path"
}
  \item{from.node.group}{
Starting group from which to calculate the shortest path  
}
  \item{to.node.group}{
Ending group from which to calculate the shortest path  
}
  \item{plot.type}{
Select variable to plot. Accepts shortest.path or harmonic.path 
}
}
\author{
Marco T. Bastos
}
\examples{
# Not run:
net.object <- watts.strogatz.game(1, 300, 4, 0.01)
V(net.object)$group <- sample(rep(LETTERS[c(1,4,6,8,4,20)],1000),100)
iteratePaths(net.object, attribute="group", from.node.group="A", to.node.group="T", stepwise.removal=1)
}
\keyword{ iteration }
\keyword{ network }
\keyword{ path length }

