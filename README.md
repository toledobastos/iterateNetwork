# IterateNetworkMetrics

The function iterateNetwork takes an igraph or network object and calculates network metrics over subsamples of the network. The argument directed.net defaults to FALSE and specifies if the network is directed. The argument net.samples defaults to a numeric vector of length 100 giving the fractions of the network from 1% to 100%, but it can be set to any other fraction higher or lower than 1:100. The argument net.iterate defaults to 10 and specifies the number of iterations to calculate. Lastly, the argument plot.estimators defaults to TRUE to include a plot of the estimates generated during the iterations.

