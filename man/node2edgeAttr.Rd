\name{node2edgeAttr}
\alias{node2edgeAttr}
\title{
Assigns edge attributes to network based on node attributes
}
\description{
Description: The function node2edgeAttr assigns edge attributes to a network object based on the node attribute specified in the argument.
}
\usage{
node2edgeAttr(net.object, vertex.attribute)
}
\arguments{
  \item{net.object}{
The network object. Accepts both igraph and network objects
}
  \item{net.samples}{
Node or vertex.attribute to be used to create the edge attribute.
}
  \item{directed}{
Defaults to FALSE. Indicates whether the network is directed or not. If undirected the function aggregates edges nodeA->nodeB and nodeB->nodeA as the same .
}
}
\author{
Marco T. Bastos
}
\examples{
# Not run:
net.object <- barabasi.game(50)
V(net.object)$color <- ifelse(V(net.object) %% 2==1, "darkblue", "darkred") # create random node attribute
E(net.object)$color <- node2edgeAttr(net.object, vertex.attribute = "color") # create edge attribute based on node attribute
}
\keyword{ iteration }
\keyword{ network }
\keyword{ sample }
\keyword{ nodes }
\keyword{ edges }
