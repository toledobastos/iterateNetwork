iterateComplete <- function(net.object,
                            attribute = NULL,
                            net.iterate = "low",
                            collapse.minor = 1,
                            return.estimates = "ALL",
                            plot.estimators = TRUE,
                            plot.type = "p") {
    
    # load dependencies
    require(intergraph)
    require(network)
    require(igraph)
    require(sna)

    # check node & edge names in network object
    if(class(net.object)=="igraph" && is.null(V(net.object)$name)) { V(net.object)$name <- 1:vcount(net.object) }
    if(class(net.object)=="igraph" && is.null(E(net.object)$name)) { E(net.object)$name <- 1:ecount(net.object) }
    if(class(net.object)=="network") { if(any(is.na(network::get.vertex.attribute(net.object, "name")))) { network::set.vertex.attribute(net.object, "name", value = 1:length(network::get.vertex.attribute(net.object, "name"))) } }
    if(class(net.object)=="network") { if(is.null(network::get.edge.attribute(net.object, "name"))) { network::set.edge.attribute(net.object, attrname = "name", value = 1:network::network.edgecount(net.object)) } }

    # generate network & igraph objects
    if(class(net.object)=="igraph") { corenet <- intergraph::asNetwork(net.object) }
    if(class(net.object)=="igraph") { corenet.g <- net.object }
    if(class(net.object)=="network") { corenet.g <- intergraph::asIgraph(net.object) }
    if(class(net.object)=="network") { corenet <- net.object }
    
    # check again node & edge names for igraph and network objects
    if(is.null(igraph::V(corenet.g)$name)) { igraph::V(corenet.g)$name <- igraph::V(corenet.g)$vertex.names }
    if(is.null(igraph::E(corenet.g)$name)) { igraph::E(corenet.g)$name <- 1:igraph::ecount(net.object) }
    if(any(is.na(network::get.vertex.attribute(corenet, "name")))) { network::set.vertex.attribute(corenet, "name", value = 1:length(network::get.vertex.attribute(corenet, "name"))) }
    if(is.null(network::get.edge.attribute(corenet, "name"))) { network::set.edge.attribute(corenet, attrname = "name", value = 1:network::network.edgecount(corenet)) }
    
    # check node & edge attribute
    if(is.null(igraph::get.vertex.attribute(corenet.g, attribute))) { stop(print(paste0("Node attribute ",attribute," not found."))) }
    
    # set seed 
    seed <- gsub("-","",as.character(Sys.Date()))
    set.seed(as.numeric(seed))
    print(paste0("Setting seed to ",seed))
    
    # load small-word quotient calculation
    small.wordness <- function(net.object) {
        n <- igraph::vcount(net.object)
        p <- igraph::graph.density(net.object)
        k <- p*(n-1)
        cc.rand <- k/n
        apl.rand <- log(n)/log(k)
        cc <- igraph::transitivity(net.object, type=c("localundirected"), isolates=c("zero"))
        cc.obs <- mean(cc)
        apl.obs <- igraph::average.path.length(net.object, directed=F, unconnected=T)
        Q <- (cc.obs/cc.rand)/(apl.obs/apl.rand) 
        return(Q) }
    
    # load divisors function
    divisors <- function(x) { y <- seq_len(x); y[ x%%y == 0 ] }
                
    # define network sample by nodes
    net.size <- vcount(corenet.g)
    
    # load loop reset function
    loop.vec.reset <- function() {
        sample.num.vec <<- as.numeric()
        nodes.num.vec <<- as.numeric()
        edges.num.vec <<- as.numeric()
        centralization.vec <<- as.numeric()
        diameter.vec <<- as.numeric()
        eigenvector.vec <<- as.numeric()
        permutation.vec <<- as.numeric()
        transitivity.vec <<- as.numeric()
        local.clustering.vec <<- as.numeric()
        articulations.vec <<- as.numeric()
        clusters.vec <<- as.numeric()
        avr.pathlength.vec <<- as.numeric()
        avr.degree.vec <<- as.numeric()
        avr.closeness.vec <<- as.numeric()
        page.rank.vec <<- as.numeric()
        betweenness.vec <<- as.numeric()
        density.vec <<- as.numeric()
        largest.component.vec <<- as.numeric()
        small.world.vec <<- as.numeric()
    }
    loop.list.reset <- function() {        
        sample.num.list <<- list()
        nodes.num.list <<- list()
        edges.num.list <<- list()
        centralization.list <<- list()
        diameter.list <<- list()
        eigenvector.list <<- list()
        permutation.list <<- list()
        transitivity.list <<- list()
        articulations.list <<- list()
        clusters.list <<- list()
        avr.pathlength.list <<- list()
        avr.degree.list <<- list()
        avr.closeness.list <<- list()
        page.rank.list <<- list()
        betweenness.list <<- list()
        density.list <<- list()
        largest.component.list <<- list()
        small.world.list <<- list()
        local.clustering.list <<- list() 
    }
    
    # index network for attribute iteration
    attribute.index <- data.frame(nodes=igraph::V(corenet.g)$name, attribute=igraph::get.vertex.attribute(corenet.g, attribute, index=igraph::V(corenet.g)))
    attribute.index$attribute <- as.factor(as.character(attribute.index$attribute))
    attribute.index <- attribute.index[order(attribute.index$attribute),]
    # check and fix singletons
    if(collapse.minor>1) {
        attribute.index.table <- as.data.frame(table(as.character(attribute.index$attribute)))
        if(any(attribute.index.table$Freq<collapse.minor)) {
            attribute.index <- merge(attribute.index, attribute.index.table, by.x="attribute", by.y="Var1", all.x=T)
            attribute.index$attribute <- as.character(attribute.index$attribute)
            attribute.index$attribute[attribute.index$Freq<collapse.minor] <- "ETAL"
            attribute.index$Freq <- NULL
        }}
    
    # calculate max interations
    net.samples.list <- list()
    attribute.unique <- unique(as.character(attribute.index$attribute))
    for(x in 1:length(unique(attribute.index$attribute))) {
        net.samples.list[[x]] <- attribute.index$nodes[as.character(attribute.index$attribute)==attribute.unique[x] ] }
    module.sizes <- unlist(lapply(net.samples.list, length))
    min.stepwise <- min(module.sizes)
    print(paste0("Max node removal for ",attribute, " is ", min(module.sizes), ". Possible stepwise removal: ",paste(divisors(min.stepwise), collapse = ", ")))        
    if(net.iterate=="low") { net.iterate <- min.stepwise }
    if(net.iterate=="med") { net.iterate <- min.stepwise*(min.stepwise/2) }
    if(net.iterate=="max") { net.iterate <- min.stepwise*min.stepwise }
    stepwise.removal <- 1
    
    # check for small subgroups
    if(min.stepwise/net.size<0.1) { warning(print(paste0("The attribute ",attribute," includes subgroups with few nodes preventing more comprehensive simulations. Increase collapse.minor to combine them and allow for larger iterations"))) }
    
    # calculate net.sample
    net.size.0 <- vcount(corenet.g)
    net.removed <- net.size.0-min.stepwise
    net.int2 <- net.removed/net.size.0
    net.samples <- rev(seq(net.int2:1, by=net.int2/min.stepwise))

    # prepare for loop
    selected.variables <- c("random", "degree", attribute.unique)
    list.complete <- vector("list", length(selected.variables))
    names(list.complete) <- selected.variables
    estimates.df <- data.frame()
    
    # select output
    if(!is.character(return.estimates)) { estimates.df <- estimates.df[,c(return.estimates)] }
    if(is.character(return.estimates)) {
        if(return.estimates=="selected") { estimates.df <- estimates.df[,c(1:6,8:12,14:17)] }
        if(return.estimates=="ALL") { estimates.df <- estimates.df } }
    
    ###
    ### step 1: iterate random 
    ###
    
    # generate output list
    loop.list.reset()
    
    # start network slicing
        for(u in 1:length(net.samples)) {
            
            # set graph sample size
            if(net.samples[u]==1) { graph.size <- net.size } 
            else { graph.size <- round(net.size*net.samples[u]+.5, digits = 0) }
            
            # reset estimates
            loop.vec.reset()
            
            # random iteration
            for(j in 1:net.iterate) {
                cat("\r","Starting random iteration",j,"of",net.iterate)
                corenet.gx <- igraph::induced.subgraph(corenet.g, which(igraph::V(corenet.g)$name %in% igraph::V(corenet.g)$name[sample(1:igraph::vcount(corenet.g), graph.size)]))
                # collect metrics per iteration
                sample.num.vec <- c(sample.num.vec, net.samples[u])
                nodes.num.vec <- c(nodes.num.vec,igraph::vcount(corenet.gx))
                edges.num.vec <- c(edges.num.vec,igraph::ecount(corenet.gx))
                centralization.vec <- c(centralization.vec,igraph::centralization.degree(corenet.gx)$centralization)
                diameter.vec <- c(diameter.vec,igraph::diameter(corenet.gx))
                eigenvector.vec <- c(eigenvector.vec,igraph::evcent(corenet.gx)$value)
                permutation.vec <- c(permutation.vec,igraph::canonical.permutation(corenet.gx)$info$nof_nodes)
                transitivity.vec <- c(transitivity.vec,igraph::transitivity(corenet.gx,type=c("globalundirected"),isolates=c("zero")))
                local.clustering.vec <- c(local.clustering.vec,mean(igraph::transitivity(corenet.gx,type=c("localundirected"),isolates=c("zero"))))
                articulations.vec <- c(articulations.vec,length(igraph::articulation.points(corenet.gx)))
                clusters.vec <- c(clusters.vec,igraph::no.clusters(corenet.gx))
                avr.pathlength.vec <- c(avr.pathlength.vec,igraph::average.path.length(corenet.gx))
                avr.degree.vec <- c(avr.degree.vec,mean(igraph::degree(corenet.gx)))
                avr.closeness.vec <- c(avr.closeness.vec,mean(igraph::closeness(corenet.gx)))
                page.rank.vec <- c(page.rank.vec,mean(igraph::page.rank(corenet.g)$vector))
                betweenness.vec <- c(betweenness.vec,igraph::centralization.betweenness(corenet.gx)$centralization)
                density.vec <- c(density.vec,igraph::graph.density(corenet.gx))
                small.world.vec <- c(small.world.vec, small.wordness(corenet.gx))
                largest.component.vec <- c(largest.component.vec,base::max(igraph::clusters(corenet.gx)$csize))
            }
            # aggregate estimates
            sample.num.list[[u]] <- as.list(sample.num.vec)
            nodes.num.list[[u]] <- as.list(nodes.num.vec)
            edges.num.list[[u]] <- as.list(edges.num.vec)
            centralization.list[[u]] <- as.list(centralization.vec)
            diameter.list[[u]] <- as.list(diameter.vec)
            eigenvector.list[[u]] <- as.list(eigenvector.vec)
            permutation.list[[u]] <- as.list(permutation.vec)
            transitivity.list[[u]] <- as.list(transitivity.vec)
            local.clustering.list[[u]] <- as.list(local.clustering.vec)
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
    # parse results
    estimates.df <- data.frame(sample=unlist(sample.num.list),nodes=unlist(nodes.num.list),
                               edges=unlist(edges.num.list),degree=unlist(avr.degree.list),
                               eigenvector=unlist(eigenvector.list),local.clustering=unlist(local.clustering.list),
                               centralization=unlist(centralization.list),diameter=unlist(diameter.list),permutation=unlist(permutation.list),
                               transitivity=unlist(transitivity.list),articulations=unlist(articulations.list),
                               cluster=unlist(clusters.list),path.length=unlist(avr.pathlength.list),
                               closeness=unlist(avr.closeness.list),page.rank=unlist(page.rank.list),
                               betweenness=unlist(betweenness.list),density=unlist(density.list),
                               largest.component=unlist(largest.component.list),small.world=unlist(small.world.list))
    list.complete[[1]] <- estimates.df
    
    ###
    ### step 2: iterate degree 
    ###
    
    # generate output list
    loop.list.reset()
    
    # start network slicing
    for(u in 1:length(net.samples)) {
        
        # set graph sample size
        if(net.samples[u]==1) { graph.size <- net.size } 
        else { graph.size <- round(net.size*net.samples[u]+.5, digits = 0) }
        
        # reset estimates
        loop.vec.reset()
        
        # degree iteration
        for(j in 1:net.iterate) {
            cat("\r","Starting degree iteration",j,"of",net.iterate)
            nodes.select <- names(sort(igraph::degree(corenet.g), decreasing=T)[(igraph::vcount(corenet.g)-graph.size):igraph::vcount(corenet.g)])
            corenet.gx <- igraph::induced.subgraph(corenet.g, which(V(corenet.g)$name %in% nodes.select))
            # collect metrics per iteration
            sample.num.vec <- c(sample.num.vec, net.samples[u])
            nodes.num.vec <- c(nodes.num.vec,igraph::vcount(corenet.gx))
            edges.num.vec <- c(edges.num.vec,igraph::ecount(corenet.gx))
            centralization.vec <- c(centralization.vec,igraph::centralization.degree(corenet.gx)$centralization)
            diameter.vec <- c(diameter.vec,igraph::diameter(corenet.gx))
            eigenvector.vec <- c(eigenvector.vec,igraph::evcent(corenet.gx)$value)
            permutation.vec <- c(permutation.vec,igraph::canonical.permutation(corenet.gx)$info$nof_nodes)
            transitivity.vec <- c(transitivity.vec,igraph::transitivity(corenet.gx,type=c("globalundirected"),isolates=c("zero")))
            local.clustering.vec <- c(local.clustering.vec,mean(igraph::transitivity(corenet.gx,type=c("localundirected"),isolates=c("zero"))))
            articulations.vec <- c(articulations.vec,length(igraph::articulation.points(corenet.gx)))
            clusters.vec <- c(clusters.vec,igraph::no.clusters(corenet.gx))
            avr.pathlength.vec <- c(avr.pathlength.vec,igraph::average.path.length(corenet.gx))
            avr.degree.vec <- c(avr.degree.vec,mean(igraph::degree(corenet.gx)))
            avr.closeness.vec <- c(avr.closeness.vec,mean(igraph::closeness(corenet.gx)))
            page.rank.vec <- c(page.rank.vec,mean(igraph::page.rank(corenet.g)$vector))
            betweenness.vec <- c(betweenness.vec,igraph::centralization.betweenness(corenet.gx)$centralization)
            density.vec <- c(density.vec,igraph::graph.density(corenet.gx))
            small.world.vec <- c(small.world.vec, small.wordness(corenet.gx))
            largest.component.vec <- c(largest.component.vec,base::max(igraph::clusters(corenet.gx)$csize))
        }
        # aggregate estimates
        sample.num.list[[u]] <- as.list(sample.num.vec)
        nodes.num.list[[u]] <- as.list(nodes.num.vec)
        edges.num.list[[u]] <- as.list(edges.num.vec)
        centralization.list[[u]] <- as.list(centralization.vec)
        diameter.list[[u]] <- as.list(diameter.vec)
        eigenvector.list[[u]] <- as.list(eigenvector.vec)
        permutation.list[[u]] <- as.list(permutation.vec)
        transitivity.list[[u]] <- as.list(transitivity.vec)
        local.clustering.list[[u]] <- as.list(local.clustering.vec)
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
    # parse results
    estimates.df <- data.frame(sample=unlist(sample.num.list),nodes=unlist(nodes.num.list),
                               edges=unlist(edges.num.list),degree=unlist(avr.degree.list),
                               eigenvector=unlist(eigenvector.list),local.clustering=unlist(local.clustering.list),
                               centralization=unlist(centralization.list),diameter=unlist(diameter.list),permutation=unlist(permutation.list),
                               transitivity=unlist(transitivity.list),articulations=unlist(articulations.list),
                               cluster=unlist(clusters.list),path.length=unlist(avr.pathlength.list),
                               closeness=unlist(avr.closeness.list),page.rank=unlist(page.rank.list),
                               betweenness=unlist(betweenness.list),density=unlist(density.list),
                               largest.component=unlist(largest.component.list),small.world=unlist(small.world.list))
    list.complete[[2]] <- estimates.df
    
    ###
    ### step 3: iterate attribute 
    ###
    
    for(a in 3:length(list.complete)) {
        # print process
        cat("\n")
        print(paste0("Targetting nodes of ",tolower(attribute)," type ",a,"."))
        
        # generate output list
        loop.list.reset()
        
        # start network slicing
        for(u in 1:length(net.samples)) {
            
            # set graph sample size
            if(net.samples[u]==1) { graph.size <- net.size } 
            else { graph.size <- round(net.size*net.samples[u]+.5, digits = 0) }
            
            # reset estimates
            loop.vec.reset()
            
            # degree iteration
            for(j in 1:net.iterate) {
                cat("\r","Iterative removal of targeted nodes",j,"of",net.iterate)
                nodes.deselect <- sample(net.samples.list[[a-2]], length(net.samples.list[[a-2]])*net.samples[u])
                nodes.deselect <- net.samples.list[[a-2]][!net.samples.list[[a-2]] %in% nodes.deselect]
                nodes.select <- V(corenet.g)$name[!V(corenet.g)$name %in% nodes.deselect]
                corenet.gx <- igraph::induced.subgraph(corenet.g, which(V(corenet.g)$name %in% nodes.select))
                # collect metrics per iteration
                sample.num.vec <- c(sample.num.vec, net.samples[u])
                nodes.num.vec <- c(nodes.num.vec,igraph::vcount(corenet.gx))
                edges.num.vec <- c(edges.num.vec,igraph::ecount(corenet.gx))
                centralization.vec <- c(centralization.vec,igraph::centralization.degree(corenet.gx)$centralization)
                diameter.vec <- c(diameter.vec,igraph::diameter(corenet.gx))
                eigenvector.vec <- c(eigenvector.vec,igraph::evcent(corenet.gx)$value)
                permutation.vec <- c(permutation.vec,igraph::canonical.permutation(corenet.gx)$info$nof_nodes)
                transitivity.vec <- c(transitivity.vec,igraph::transitivity(corenet.gx,type=c("globalundirected"),isolates=c("zero")))
                local.clustering.vec <- c(local.clustering.vec,mean(igraph::transitivity(corenet.gx,type=c("localundirected"),isolates=c("zero"))))
                articulations.vec <- c(articulations.vec,length(igraph::articulation.points(corenet.gx)))
                clusters.vec <- c(clusters.vec,igraph::no.clusters(corenet.gx))
                avr.pathlength.vec <- c(avr.pathlength.vec,igraph::average.path.length(corenet.gx))
                avr.degree.vec <- c(avr.degree.vec,mean(igraph::degree(corenet.gx)))
                avr.closeness.vec <- c(avr.closeness.vec,mean(igraph::closeness(corenet.gx)))
                page.rank.vec <- c(page.rank.vec,mean(igraph::page.rank(corenet.g)$vector))
                betweenness.vec <- c(betweenness.vec,igraph::centralization.betweenness(corenet.gx)$centralization)
                density.vec <- c(density.vec,igraph::graph.density(corenet.gx))
                small.world.vec <- c(small.world.vec, small.wordness(corenet.gx))
                largest.component.vec <- c(largest.component.vec,base::max(igraph::clusters(corenet.gx)$csize))
            }
            # aggregate estimates
            sample.num.list[[u]] <- as.list(sample.num.vec)
            nodes.num.list[[u]] <- as.list(nodes.num.vec)
            edges.num.list[[u]] <- as.list(edges.num.vec)
            centralization.list[[u]] <- as.list(centralization.vec)
            diameter.list[[u]] <- as.list(diameter.vec)
            eigenvector.list[[u]] <- as.list(eigenvector.vec)
            permutation.list[[u]] <- as.list(permutation.vec)
            transitivity.list[[u]] <- as.list(transitivity.vec)
            local.clustering.list[[u]] <- as.list(local.clustering.vec)
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
        # parse results
        estimates.df <- data.frame(sample=unlist(sample.num.list),nodes=unlist(nodes.num.list),
                                   edges=unlist(edges.num.list),degree=unlist(avr.degree.list),
                                   eigenvector=unlist(eigenvector.list),local.clustering=unlist(local.clustering.list),
                                   centralization=unlist(centralization.list),diameter=unlist(diameter.list),permutation=unlist(permutation.list),
                                   transitivity=unlist(transitivity.list),articulations=unlist(articulations.list),
                                   cluster=unlist(clusters.list),path.length=unlist(avr.pathlength.list),
                                   closeness=unlist(avr.closeness.list),page.rank=unlist(page.rank.list),
                                   betweenness=unlist(betweenness.list),density=unlist(density.list),
                                   largest.component=unlist(largest.component.list),small.world=unlist(small.world.list))
        
        list.complete[[a]] <- estimates.df
    }
  
    # prepare for plotting
    estimates.total <- ncol(list.complete[[1]])-1
    estimates.count <- nrow(list.complete[[1]])
    
    # plot data
    if(plot.estimators==TRUE) {
        # plot observed estimators
        if(estimates.count<=500) { lwd.by.iteration <- 1}
        if(estimates.count>500 && estimates.count<=1000) { lwd.by.iteration <- .7}
        if(estimates.count>1000 && estimates.count<=2000) { lwd.by.iteration <- .5}
        if(estimates.count>2000 && estimates.count<=4000) { lwd.by.iteration <- .3}
        if(estimates.count>4000) { lwd.by.iteration <- 0.1}
        colorsmetric <- rainbow(length(list.complete))
        # set plot window
        if(estimates.total<6) { plot.panels <- c(estimates.total,1) }
        if(estimates.total==6) { plot.panels <- c(3,2) }
        if(estimates.total>6 && estimates.total<9) { plot.panels <- c(4,4) }
        if(estimates.total==9) { plot.panels <- c(3,3) }
        if(estimates.total==10) { plot.panels <- c(5,2) }
        if(estimates.total>10 && estimates.total<17) { plot.panels <- c(4,4) }
        if(estimates.total==15) { plot.panels <- c(3,5) }
        if(estimates.total>16 && estimates.total<18) { plot.panels <- c(4,5) }
        if(estimates.total==18) { plot.panels <- c(3,6) }
        if(estimates.total>18) { plot.panels <- c(5,5) }
        png(paste0("network_estimates_complete_by_",tolower(attribute),"_with_",net.iterate,"_iterations_over_",length(net.samples),"_projections.png"), type='cairo', width=plot.panels[2]*4,height=plot.panels[1]*4, units='in', res=200)
        par(mfrow=plot.panels, oma = c(4, 1, 1, 1))
        labels.plot1 <- 1:length(list.complete[[1]]$sample)
        labels.plot2 <- paste0(round(list.complete[[1]]$sample, 2)*100,"%")
        for(i in 2:ncol(list.complete[[1]])) { 
            plot(as.numeric(list.complete[[1]][,i]), xlab="", ylab="", col=colorsmetric[u], cex=0.5, xaxt="n", main=paste(colnames(list.complete[[1]])[i]), type=plot.type, lwd=lwd.by.iteration,cex.lab=1.6, cex.axis=1.6, cex.main=2.5, cex.sub=2)
            for(u in 1:length(list.complete)) {
                if(plot.type=="l") { lines(as.numeric(list.complete[[u]][,i]), xlab="", ylab="", col=colorsmetric[u], cex=0.5, xaxt="n", lwd=lwd.by.iteration,cex.lab=1.6, cex.axis=1.6, cex.main=2.5, cex.sub=2, type="l") }
                if(plot.type=="p") { points(as.numeric(list.complete[[u]][,i]), xlab="", ylab="", col=colorsmetric[u], cex=0.5, xaxt="n", lwd=lwd.by.iteration,cex.lab=1.6, cex.axis=1.6, cex.main=2.5, cex.sub=2, type="p") }
                
            }
            axis(1, at=labels.plot1, labels=labels.plot2)
        }
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
        legend("bottom", names(list.complete), xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty="n", pch=19, col=colorsmetric, cex=3)
        dev.off() 
    }    
    print(paste0("Iteration completed. Plot saved at ",getwd(),"/network_estimates_complete_by_",tolower(attribute),"_with_",net.iterate,"_iterations_over_",length(net.samples),"_projections.png"))
    return(list.complete)
}

