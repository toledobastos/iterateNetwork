iterateNetwork <- function(net.object,
                           net.samples = rev(seq(0.01:1,by=0.01)),
                           net.iterate = 10,
                           iteration.type="random",
                           attribute=NULL,
                           stepwise.removal=3,
                           return.estimates="selected",
                           plot.estimators=TRUE) {
    
    # load dependencies
    require(intergraph)
    require(network)
    require(igraph)
    require(sna)

    # check for request error
    if(iteration.type=="attribute" && is.null(attribute)) { stop(print(paste0("iteration.type by attribute requires specifying vertex attribute.")))}
    
    # check vector names in network object
    if(class(net.object)=="igraph" && is.null(V(net.object)$name)) { V(net.object)$name <- 1:vcount(net.object)}
    if(class(net.object)=="network") { if(any(is.na(network::get.vertex.attribute(net.object, "name")))) { network::set.vertex.attribute(net.object, "name", value = 1:length(network::get.vertex.attribute(net.object, "name"))) } }

    # generate network & igraph objects
    if(class(net.object)=="igraph") { corenet <- intergraph::asNetwork(net.object) }
    if(class(net.object)=="igraph") { corenet.g <- net.object }
    if(class(net.object)=="network") { corenet.g <- intergraph::asIgraph(net.object) }
    if(class(net.object)=="network") { corenet <- net.object }
    
    # check vector names for igraph object again
    if(is.null(V(corenet.g)$name)) { V(corenet.g)$name <- V(corenet.g)$vertex.names }
    
    # set seed 
    seed <- gsub("-","",as.character(Sys.Date()))
    set.seed(as.numeric(seed))
    print(paste0("Setting seed to ",seed))
    
    # prepare for loop
    estimates.df <- data.frame()
    net.size <- vcount(corenet.g)
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
    small.world.list <- list()
    
    # index node for attribute iteration
    if(iteration.type=="attribute") {
        attribute.index <- data.frame(nodes=igraph::V(corenet.g)$name, attribute=igraph::get.vertex.attribute(corenet.g, attribute, index=V(corenet.g)))
        attribute.index$attribute <- as.factor(as.character(attribute.index$attribute))
        attribute.index <- attribute.index[order(attribute.index$attribute),]
        net.samples.list <- list()
        attribute.unique <- unique(as.character(attribute.index$attribute))
        for(x in 1:length(unique(attribute.index$attribute))) {
            net.samples.list[[x]] <- attribute.index$nodes[as.character(attribute.index$attribute)==attribute.unique[x] ] }
        net.samples <- as.character(sort(unique(attribute.index$attribute)))
        module.sizes <- unlist(lapply(net.samples.list, length))
        min.stepwise <- min(module.sizes)-1
        print(paste0("Maximum stepwise removal for this network attribute is ", min.stepwise))
        if(stepwise.removal>min.stepwise) { 
            stop(print(paste0("Maximum stepwise removal for this network attribute is ", min.stepwise,"! Set stepwise.removal accordingly.")))}
    }
    
    # start network slicing
    for(u in 1:length(net.samples)) {
        # set graph sample size
        if(iteration.type=="attribute") { graph.size <- vcount(corenet.g)-stepwise.removal } 
        if(iteration.type!="attribute") { graph.size <- round(net.size*net.samples[u]+.5, digits = 0) }
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
        small.world.vec <- as.numeric()
        # start iteration
        for(j in 1:net.iterate) {
            if(iteration.type=="random") { 
                cat("\r","Starting random iteration",j,"of",net.iterate)
                corenet.gx <- induced.subgraph(corenet.g, which(V(corenet.g)$name %in% V(corenet.g)$name[sample(1:vcount(corenet.g), graph.size)])) }
            if(iteration.type=="degree") { 
                cat("\r","Starting degree iteration",j,"of",net.iterate)
                nodes.select <- names(sort(degree(corenet.g), decreasing=T)[(vcount(corenet.g)-graph.size):vcount(corenet.g)])
                corenet.gx <- induced.subgraph(corenet.g, which(V(corenet.g)$name %in% nodes.select)) }
            if(iteration.type=="betweenness") { 
                cat("\r","Starting betweenness iteration",j,"of",net.iterate)
                nodes.select <- names(sort(betweenness(corenet.g), decreasing=T)[(vcount(corenet.g)-graph.size):vcount(corenet.g)])
                corenet.gx <- induced.subgraph(corenet.g, which(V(corenet.g)$name %in% nodes.select)) }
            if(iteration.type=="closeness") { 
                cat("\r","Starting closeness iteration",j,"of",net.iterate)
                nodes.select <- names(sort(closeness(corenet.g), decreasing=T)[(vcount(corenet.g)-graph.size):vcount(corenet.g)])
                corenet.gx <- induced.subgraph(corenet.g, which(V(corenet.g)$name %in% nodes.select)) }
            if(iteration.type=="attribute") { 
                cat("\r","Starting attribute iteration",j,"of",net.iterate) 
                nodes.deselect <- sample(net.samples.list[[u]], stepwise.removal)
                nodes.select <- V(corenet.g)$name[!V(corenet.g)$name %in% nodes.deselect]
                corenet.gx <- induced.subgraph(corenet.g, which(V(corenet.g)$name %in% nodes.select)) }
            # collect metrics per iteration
            nodes.num.vec <- c(nodes.num.vec,igraph::vcount(corenet.gx))
            edges.num.vec <- c(edges.num.vec,igraph::ecount(corenet.gx))
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
            largest.component.vec <- c(largest.component.vec,sum(sna::component.largest(as.network(as.matrix(igraph::get.adjacency(corenet.gx)), directed = igraph::is.directed(corenet.gx)), connected=c("strong"))))
            small.world.vec <- c(small.world.vec, (mean(igraph::transitivity(corenet.gx, type=c("localundirected"), isolates=c("zero")))/igraph::graph.density(corenet.gx)/igraph::vcount(corenet.gx))/(igraph::average.path.length(corenet.gx)/log(igraph::vcount(corenet.gx))/log(igraph::graph.density(corenet.gx)*(igraph::vcount(corenet.gx)-1))))
        }
        # aggregate estimates        
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
        small.world.list[[u]] <- as.list(small.world.vec)
        # clear sample network
        rm(corenet.gx)
        # print process
        cat("\n")
        print(paste0("Completed iteration ",u," of ", length(net.samples),".")) 
    }
    estimates.df <- data.frame(nodes=unlist(nodes.num.list),
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
                               largest.component=unlist(largest.component.list),
                               small.world=unlist(small.world.list))

    # select output
    if(return.estimates!="ALL") {
        if(return.estimates=="selected") { estimates.df <- estimates.df[,c(1:6,8:12,14:17)] } else {
            estimates.df <- estimates.df[,c(return.estimates)] } }
    estimates.total <- ncol(estimates.df)
    
    # add sample marker
    estimates.df <- cbind(data.frame(sample=rep(net.samples, each = net.iterate)), estimates.df)
    
    # plot data
    if(plot.estimators==TRUE) {
        # calculate divisor
        divisors <- function(x) { y <- seq_len(x); y[ x%%y == 0 ] }
        # plot observed estimators
        if(net.iterate<50) { lwd.by.iteration <- 2}
        if(net.iterate>50 && net.iterate<100) { lwd.by.iteration <- 1}
        if(net.iterate>100 && net.iterate<500) { lwd.by.iteration <- 0.5}
        if(net.iterate>500) { lwd.by.iteration <- 0.15}
        colorsmetric <- rainbow(estimates.total+1)
        # set plot window
        if(estimates.total<6) { plot.panels <- c(estimates.total,1) }
        if(estimates.total==6) { plot.panels <- c(3,2) }
        if(estimates.total>6 && estimates.total<9) { plot.panels <- c(4,4) }
        if(estimates.total==9) { plot.panels <- c(3,3) }
        if(estimates.total==10) { plot.panels <- c(5,2) }
        if(estimates.total>10 && estimates.total<17) { plot.panels <- c(4,4) }
        if(estimates.total==15) { plot.panels <- c(3,5) }
        if(estimates.total>16 && estimates.total<20) { plot.panels <- c(4,5) }
        if(estimates.total>20) { plot.panels <- c(5,round(median(divisors(estimates.total)))) }        
        png(paste0("network_estimates_",net.iterate,"_iterations_type_l.png"), type='cairo', width=20,height=12, units='in', res=200)
        par(mfrow=plot.panels)
        if(iteration.type!="attribute") {
            labels.plot1 <- 1:length(estimates.df$sample)
            labels.plot2 <- paste0(rev(seq(from=100/length(estimates.df$sample), to=100, by=(100/length(estimates.df$sample)))),"%") }
        if(iteration.type=="attribute") { 
            labels.plot1 <- 1:length(estimates.df$sample)
            labels.plot2 <- paste0(estimates.df$sample) }
        for(i in 2:ncol(estimates.df)) {
            plot(as.numeric(estimates.df[,i]), xlab="", ylab="", col=colorsmetric[i], cex=.5, xaxt="n",
                 main=paste(colnames(estimates.df)[i]), type="p", lwd=lwd.by.iteration,cex.lab=1.6, cex.axis=1.6, cex.main=2.5, cex.sub=2)
<<<<<<< HEAD
            lines(as.numeric(estimates.df[,i]), col="black", lwd = lwd.by.iteration)
=======
            lines(as.numeric(estimates.df[,i]), col=black, lwd = lwd.by.iteration)
>>>>>>> fc250b33c2df95ab75c1fc1458e0412597a2174c
            axis(1, at=labels.plot1, labels=labels.plot2)
        }
        dev.off()
    }
    return(estimates.df)
}
