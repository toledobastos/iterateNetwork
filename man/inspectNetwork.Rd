\name{inspectNetwork}
\alias{inspectNetwork}
\title{
Returns sample network for further inspection
}
\description{
Description: The function inspectNetwork returns a sample igraph network for further inspection.
}
\usage{
inspectNetwork(net.object, iteration.type = "random", sample.size = .5)
}
\arguments{
  \item{net.object}{
The network object. Accepts both igraph and network objects
}
  \item{iteration.type}{
Type of network iteration to be inspected. It accepts degree, betweenness, closeness, and random (default).
}
  \item{sample.size}{
Percentage of the network to be returned.
}
}
\author{
Marco T. Bastos
}
\examples{
# Not run:
net.object <- barabasi.game(50)
inspectNetwork(net.object, iteration.type = "random", sample.size = .8)
}
\keyword{ iteration }
\keyword{ network }
\keyword{ sample }
\keyword{ nodes }
