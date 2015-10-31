\name{iterateNetwork}
\alias{iterateNetwork}
\title{
Sample and Iteration for Network Metrics
}
\description{
Description: The function iterateNetwork takes an igraph or network object and calculates network metrics over subsamples of the network. The argument directed.net defaults to FALSE and specifies if the network is directed. The argument net.samples defaults to a numeric vector of length 100 giving the fractions of the network from 1% to 100%, but it can be set to any other fraction higher or lower than 1:100. The argument net.iterate defaults to 10 and specifies the number of iterations to calculate. Lastly, the argument plot.estimators defaults to TRUE to include a plot of the estimates generated during the iterations.
}
\usage{
net.object <- barabasi.game(50)
iterateNetwork(net.object, directed.net = FALSE, net.samples = rev(seq(0.01:1, by = 0.01)), net.iterate = 10, plot.estimators = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net.object}{
The network object. Accepts both igraph and network objects
}
  \item{directed.net}{
Defaults to FALSE and specifies if the network is directed.
}
  \item{net.samples}{
Defaults to a numeric vector of length 100 giving the fractions of the network from 1% to 100%. It can be set to any other fraction higher or lower than 1:100.
}
  \item{net.iterate}{
Defaults to 10 and specifies the number of iterations to calculate.
}
  \item{plot.estimators}{
Defaults to TRUE to include a plot of the estimates generated during the iterations.
}
}
\author{
Marco T. Batos
}

\examples{
net.object <- barabasi.game(50) # generate sample network
iterateNetwork(net.object)
}
\keyword{ iteration }
\keyword{ network }
\keyword{ sample }
