iterateNetwork <- function(net.object, 
                           directed.net = FALSE, 
                           net.samples = rev(seq(0.01:1,by=0.01)),
                           net.iterate = 10,
                           plot.estimators=TRUE) {
    
    # load dependencies
    require(network)
    require(igraph)
    require(sna)
    
    # generate network & igraph objects
    if(class(net.object)=="igraph") { corenet <- as.network(as.matrix(get.adjacency(net.object)), directed = directed.net) }
    if(class(net.object)=="network") { corenet <- net.object }
    if(directed.net==FALSE & class(net.object)=="network") { corenet.g <- graph.adjacency(as.sociomatrix(net.object), mode="undirected") }
    if(directed.net==TRUE & class(net.object)=="network") { corenet.g <- graph.adjacency(as.sociomatrix(net.object), mode="directed") }
    
    # set seed 
    seed <- gsub("-","",as.character(Sys.Date()))
    set.seed(as.numeric(seed))
    print(paste0("Setting seed to ",seed))
    
    # prepare for loop
    estimates.df <- data.frame()
    net.size <- length(V(corenet.g))
    nodes.num.list <- list()
    edges.num.list <- list()
    centralization.list <- list()
    diameter.list <- list()
    eigenvector.list <- list()
    permutation.list <- list()
    transitivity.list <- list()
    articulations.list <- list()
    clusters.list <- list()
    avr.pathlength.list <- list()
    avr.degree.list <- list()
    avr.closeness.list <- list()
    page.rank.list <- list()
    betweenness.list <- list()
    density.list <- list()
    largest.component.list <- list()
    
    for(u in 1:length(net.samples)) {
        # set graph sample size
        graph.size <- round(net.size*net.samples[u], digits = 0)
        # reset estimates
        nodes.num.vec <- as.numeric()
        edges.num.vec <- as.numeric()
        centralization.vec <- as.numeric()
        diameter.vec <- as.numeric()
        eigenvector.vec <- as.numeric()
        permutation.vec <- as.numeric()
        transitivity.vec <- as.numeric()
        articulations.vec <- as.numeric()
        clusters.vec <- as.numeric()
        avr.pathlength.vec <- as.numeric()
        avr.degree.vec <- as.numeric()
        avr.closeness.vec <- as.numeric()
        page.rank.vec <- as.numeric()
        betweenness.vec <- as.numeric()
        density.vec <- as.numeric()
        largest.component.vec <- as.numeric()
        # start iteration
        for(j in 1:net.iterate) {
            corenet.gx <- induced.subgraph(corenet.g, which(V(corenet.g)$name %in% V(corenet.g)$name[sample(1:vcount(corenet.g), graph.size)]))
            nodes.num.vec <- c(nodes.num.vec,vcount(corenet.gx))
            edges.num.vec <- c(edges.num.vec,ecount(corenet.gx))
            centralization.vec <- c(centralization.vec,igraph::centralization.degree(corenet.gx)$centralization)
            diameter.vec <- c(diameter.vec,igraph::diameter(corenet.gx))
            eigenvector.vec <- c(eigenvector.vec,igraph::evcent(corenet.gx)$value)
            permutation.vec <- c(permutation.vec,igraph::canonical.permutation(corenet.gx)$info$nof_nodes)
            transitivity.vec <- c(transitivity.vec,igraph::transitivity(corenet.gx))
            articulations.vec <- c(articulations.vec,length(igraph::articulation.points(corenet.gx)))
            clusters.vec <- c(clusters.vec,igraph::clusters(corenet.gx)$csize[1])
            avr.pathlength.vec <- c(avr.pathlength.vec,igraph::average.path.length(corenet.gx))
            avr.degree.vec <- c(avr.degree.vec,mean(igraph::degree(corenet.gx)))
            avr.closeness.vec <- c(avr.closeness.vec,mean(igraph::closeness(corenet.gx)))
            page.rank.vec <- c(page.rank.vec,mean(igraph::page.rank(corenet.g)$vector))
            betweenness.vec <- c(betweenness.vec,igraph::centralization.betweenness(corenet.gx)$centralization)
            density.vec <- c(density.vec,igraph::graph.density(corenet.gx))
            largest.component.vec <- c(largest.component.vec,sum(sna::component.largest(as.network(as.matrix(get.adjacency(corenet.gx)), directed = directed.net), connected=c("strong"))))
        }
        nodes.num.list[[u]] <- as.list(nodes.num.vec)
        edges.num.list[[u]] <- as.list(edges.num.vec)
        centralization.list[[u]] <- as.list(centralization.vec)
        diameter.list[[u]] <- as.list(diameter.vec)
        eigenvector.list[[u]] <- as.list(eigenvector.vec)
        permutation.list[[u]] <- as.list(permutation.vec)
        transitivity.list[[u]] <- as.list(transitivity.vec)
        articulations.list[[u]] <- as.list(articulations.vec)
        clusters.list[[u]] <- as.list(clusters.vec)
        avr.pathlength.list[[u]] <- as.list(avr.pathlength.vec)
        avr.degree.list[[u]] <- as.list(avr.degree.vec)
        avr.closeness.list[[u]] <- as.list(avr.closeness.vec)
        page.rank.list[[u]] <- as.list(page.rank.vec)
        betweenness.list[[u]] <- as.list(betweenness.vec)
        density.list[[u]] <- as.list(density.vec)
        largest.component.list[[u]] <- as.list(largest.component.vec)
        # clear sample network
        rm(corenet.gx)
        # print process
        print(paste0("Interation ",u," of ", length(net.samples)," complete."))
    }

# agregate results
estimates.df <- data.frame(sample = rep(net.samples, each = net.iterate),
                           nodes=unlist(nodes.num.list),
                           edges=unlist(edges.num.list),
                           centralization=unlist(centralization.list),
                           diameter=unlist(diameter.list),
                           eigenvector=unlist(eigenvector.list),
                           permutation=unlist(permutation.list),
                           transitivity=unlist(transitivity.list),
                           articulations=unlist(articulations.list),
                           cluster=unlist(clusters.list),
                           path.length=unlist(avr.pathlength.list),
                           degree=unlist(avr.degree.list),
                           closeness=unlist(avr.closeness.list),
                           page.rank=unlist(page.rank.list),
                           betweenness=unlist(betweenness.list),
                           density=unlist(density.list),
                           largest.component=unlist(largest.component.list))

if(plot.estimators==TRUE) {
    # plot observed estimators
    colorsmetric <- rainbow(ncol(estimates.df))
    png(paste0("network_estimates_",net.iterate,"_iterations.png"), type='cairo', width=20,height=12, units='in', res=200)
    par(mfrow=c(4,4))
    for(i in 2:ncol(estimates.df)) {
        plot(as.numeric(estimates.df[,i]), xlab="", ylab="", col=colorsmetric[i], cex=.5, xaxt="n",
             main=paste(colnames(estimates.df)[i]), type="l",lwd=3,cex.lab=1.6, cex.axis=1.6, cex.main=2.5, cex.sub=2)
        axis(1, at=1:length(estimates.df$sample), labels=paste0((rev(seq(from=1, to=length(estimates.df$sample), by=1))/net.iterate),"%"))
    }
    dev.off()
    }
return(estimates.df)
}
