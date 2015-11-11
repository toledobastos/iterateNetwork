iterateNetwork <- function(net.object,
                           iteration.type ="random",
                           net.samples = rev(seq(0.01:1,by=0.01)),
                           removal = "node",
                           net.iterate = 10,
                           attribute=NULL,
                           stepwise.removal = "auto",
                           return.estimates = "ALL",
                           plot.estimators = TRUE,
                           plot.type = "p") {
    
    # load dependencies
    require(intergraph)
    require(network)
    require(igraph)
    require(sna)

    # check for request error
    if(iteration.type=="attribute" && is.null(attribute)) { stop(print(paste0("iteration.type by attribute requires specifying vertex attribute.")))}
    if(iteration.type=="node.interaction" && is.null(attribute)) { stop(print(paste0("iteration.type by attribute requires specifying node.interaction attribute.")))}
    if(iteration.type=="edge" && is.null(attribute)) { stop(print(paste0("iteration.type by edge requires specifying edge attribute.")))}

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
    
    # transform node to edge attribute for node.interaction
    if(removal=="node.interaction") {
        df1 <- as.data.frame(igraph::get.edgelist(corenet.g))
        df1$V1 <- as.character(df1$V1)
        df1$V2 <- as.character(df1$V2)
        user.from <- merge(df1, data.frame(V1=as.character(igraph::V(corenet.g)$name), attr=igraph::get.vertex.attribute(corenet.g, attribute)), by="V1", all.x=T)[3]$attr
        user.to <- merge(df1, data.frame(V2=as.character(igraph::V(corenet.g)$name), attr=igraph::get.vertex.attribute(corenet.g, attribute)), by="V2", all.x=T)[3]$attr
        edge.attr <- paste(user.from, user.to, sep="+")
        if(!igraph::is.directed(corenet.g)) {
            edge.attr.df <- as.data.frame(igraph::get.edgelist(igraph::graph.data.frame(data.frame(from=user.from, to=user.to), directed=F)))
            edge.attr <- paste(edge.attr.df$V1, edge.attr.df$V2, sep="+") }
        corenet.g <- igraph::set.edge.attribute(graph=corenet.g, name=attribute, value=edge.attr)
        removal="edge"
    }
    
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
    divisors <- function(x) { 
        y <- seq_len(x); y[ x%%y == 0 ] }
    
    # define network sample by nodes or edges
    if(removal=="node") { net.size <- vcount(corenet.g) }
    if(removal=="edge") { net.size <- ecount(corenet.g) }
    
    # prepare for loop
    estimates.df <- data.frame()
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
    local.clustering.list <- list()
    
    # index network for attribute iteration
    if(iteration.type=="attribute") {
        if(removal=="node") { 
        attribute.index <- data.frame(nodes=igraph::V(corenet.g)$name, attribute=igraph::get.vertex.attribute(corenet.g, attribute, index=igraph::V(corenet.g)))
        attribute.index$attribute <- as.factor(as.character(attribute.index$attribute))
        attribute.index <- attribute.index[order(attribute.index$attribute),]
        # check for singletons and fix if necessary
        attribute.index.table <- as.data.frame(table(as.character(attribute.index$attribute)))
        if(any(attribute.index.table$Freq<4)) {
            attribute.index <- merge(attribute.index, attribute.index.table, by.x="attribute", by.y="Var1", all.x=T)
            attribute.index$attribute <- as.character(attribute.index$attribute)
            attribute.index$attribute[attribute.index$Fre<round(median(attribute.index.table$Freq))] <- "ETAL"
            attribute.index$Freq <- NULL
        }
        net.samples.list <- list()
        attribute.unique <- unique(as.character(attribute.index$attribute))
        for(x in 1:length(unique(attribute.index$attribute))) {
            net.samples.list[[x]] <- attribute.index$nodes[as.character(attribute.index$attribute)==attribute.unique[x] ] }
        module.sizes <- unlist(lapply(net.samples.list, length))
        min.stepwise <- min(module.sizes)
        print(paste0("Max node removal for ",attribute, " is ", min(module.sizes), ". Possible stepwise removal: ",paste(divisors(min.stepwise), collapse = ", ")))        
        if(stepwise.removal=="auto") { stepwise.removal <- net.iterate <- round(sqrt(min.stepwise)-.5) }
        else { net.iterate <- round(min(module.sizes)/stepwise.removal) }
        net.samples <- attribute.unique
        if(stepwise.removal>min.stepwise | as.integer(min(module.sizes)/stepwise.removal*stepwise.removal)>min(module.sizes)) { 
            stop(print(paste0("Maximum stepwise node removal for this network attribute is ", round(min(module.sizes)/2), " per iteration. Decrease stepwise.removal to match this threshold.")))}
    } 
    if(removal=="edge") { 
        attribute.index <- data.frame(edges=igraph::E(corenet.g)$name, attribute=igraph::get.edge.attribute(corenet.g, attribute, index=igraph::E(corenet.g)))
        attribute.index$attribute <- as.factor(as.character(attribute.index$attribute))
        attribute.index <- attribute.index[order(attribute.index$attribute),]
        # check for singletons and fix if necessary
        attribute.index.table <- as.data.frame(table(as.character(attribute.index$attribute)))
        if(any(attribute.index.table$Freq<4)) {
            attribute.index <- merge(attribute.index, attribute.index.table, by.x="attribute", by.y="Var1", all.x=T)
            attribute.index$attribute <- as.character(attribute.index$attribute)
            attribute.index$attribute[attribute.index$Fre<round(median(attribute.index.table$Freq))] <- "ETAL"
            attribute.index$Freq <- NULL
        }
        net.samples.list <- list()
        attribute.unique <- unique(as.character(attribute.index$attribute))
        for(x in 1:length(unique(attribute.index$attribute))) {
            net.samples.list[[x]] <- attribute.index$edges[as.character(attribute.index$attribute)==attribute.unique[x] ] }
        module.sizes <- unlist(lapply(net.samples.list, length))
        min.stepwise <- min(module.sizes)
        print(paste0("Max edge removal for ",attribute, " is ", min(module.sizes), ". Possible stepwise removal: ",paste(divisors(min.stepwise), collapse = ", ")))        
        if(stepwise.removal=="auto") { stepwise.removal <- net.iterate <- round(sqrt(min.stepwise)-.5) }
        else { net.iterate <- round(min(module.sizes)/stepwise.removal) }
        net.samples <- attribute.unique
        if(stepwise.removal>min.stepwise | as.integer(min(module.sizes)/stepwise.removal*stepwise.removal)>min(module.sizes)) { 
            stop(print(paste0("Maximum stepwise edge removal for this network attribute is ", round(min(module.sizes)/2), " per iteration. Decrease stepwise.removal to match this threshold.")))}
    }
    }
        
    # start network slicing
    for(u in 1:length(net.samples)) {
        
        # set graph sample size
        if(iteration.type!="attribute") { 
            if(removal=="node") {
                if(net.samples[u]==1) { graph.size <- net.size } 
                else { graph.size <- round(net.size*net.samples[u]+.5, digits = 0) } } 
            if(removal=="edge") {
                if(net.samples[u]==1) { graph.size <- net.size } 
                else { graph.size <- round(net.size*net.samples[u]+.5, digits = 0) } } 
        }
        
        # reset estimates
        nodes.num.vec <- as.numeric()
        edges.num.vec <- as.numeric()
        centralization.vec <- as.numeric()
        diameter.vec <- as.numeric()
        eigenvector.vec <- as.numeric()
        permutation.vec <- as.numeric()
        transitivity.vec <- as.numeric()
        local.clustering.vec <- as.numeric()
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
                if(removal=="node") {
                    corenet.gx <- igraph::induced.subgraph(corenet.g, which(igraph::V(corenet.g)$name %in% igraph::V(corenet.g)$name[sample(1:igraph::vcount(corenet.g), graph.size)])) }
                if(removal=="edge") {
                    corenet.gx <- igraph::subgraph.edges(corenet.g, eids=which(igraph::E(corenet.g)$name %in% igraph::E(corenet.g)$name[sample(1:igraph::ecount(corenet.g), graph.size)]), delete.vertices=TRUE) }
            }
            if(iteration.type=="degree") {
                if(removal=="node") {
                    cat("\r","Starting degree iteration",j,"of",net.iterate)
                    nodes.select <- names(sort(igraph::degree(corenet.g), decreasing=T)[(igraph::vcount(corenet.g)-graph.size):igraph::vcount(corenet.g)])
                    corenet.gx <- igraph::induced.subgraph(corenet.g, which(V(corenet.g)$name %in% nodes.select)) }
                if(removal=="edge") { stop(print(paste0("iteration.type ",iteration.type, " is not a valid attribute for edge iterations."))) }
            }
            if(iteration.type=="betweenness") {
                if(removal=="node") {
                    cat("\r","Starting betweenness iteration",j,"of",net.iterate)
                    nodes.select <- names(sort(igraph::betweenness(corenet.g), decreasing=T)[(igraph::vcount(corenet.g)-graph.size):igraph::vcount(corenet.g)])
                    corenet.gx <- igraph::induced.subgraph(corenet.g, which(V(corenet.g)$name %in% nodes.select)) }
                if(removal=="edge") { stop(print(paste0("iteration.type ",iteration.type, " is not a valid attribute for edge iterations."))) }
            }
            if(iteration.type=="closeness") { 
                if(removal=="node") {
                    cat("\r","Starting closeness iteration",j,"of",net.iterate)
                    nodes.select <- names(sort(igraph::closeness(corenet.g), decreasing=T)[(igraph::vcount(corenet.g)-graph.size):igraph::vcount(corenet.g)])
                    corenet.gx <- igraph::induced.subgraph(corenet.g, which(V(corenet.g)$name %in% nodes.select)) }
                if(removal=="edge") { stop(print(paste0("iteration.type ",iteration.type, " is not a valid attribute for edge iterations."))) }
            }
            if(iteration.type=="attribute") {
                if(removal=="node") {
                    cat("\r","Iterative removal of targeted nodes",j,"of",net.iterate)
                    nodes.deselect <- sample(net.samples.list[[u]], stepwise.removal*j)
                    nodes.select <- V(corenet.g)$name[!V(corenet.g)$name %in% nodes.deselect]
                    corenet.gx <- igraph::induced.subgraph(corenet.g, which(V(corenet.g)$name %in% nodes.select)) }
                if(removal=="edge") { 
                    cat("\r","Iterative removal of targeted edges",j,"of",net.iterate)
                    edges.deselect <- sample(net.samples.list[[u]], stepwise.removal*j)
                    edges.select <- E(corenet.g)$name[!E(corenet.g)$name %in% edges.deselect]
                    corenet.gx <- igraph::subgraph.edges(corenet.g, eids=which(igraph::E(corenet.g)$name %in% edges.select), delete.vertices=TRUE) }
            }
            
            # collect metrics per iteration
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
            #         largest.component.vec <- c(largest.component.vec,sum(sna::component.largest(as.network(as.matrix(igraph::get.adjacency(corenet.gx)), directed = igraph::is.directed(corenet.gx)), connected=c("strong"))))
            largest.component.vec <- c(largest.component.vec,base::max(igraph::clusters(corenet.gx)$csize))
        }
        # aggregate estimates        
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
    estimates.df <- data.frame(nodes=unlist(nodes.num.list),
                               edges=unlist(edges.num.list),
                               degree=unlist(avr.degree.list),
                               eigenvector=unlist(eigenvector.list),
                               local.clustering=unlist(local.clustering.list),
                               centralization=unlist(centralization.list),
                               diameter=unlist(diameter.list),
                               permutation=unlist(permutation.list),
                               transitivity=unlist(transitivity.list),
                               articulations=unlist(articulations.list),
                               cluster=unlist(clusters.list),
                               path.length=unlist(avr.pathlength.list),
                               closeness=unlist(avr.closeness.list),
                               page.rank=unlist(page.rank.list),
                               betweenness=unlist(betweenness.list),
                               density=unlist(density.list),
                               largest.component=unlist(largest.component.list),
                               small.world=unlist(small.world.list))

    # select output
    if(!is.character(return.estimates)) { estimates.df <- estimates.df[,c(return.estimates)] }
    if(is.character(return.estimates)) {
        if(return.estimates=="selected") { estimates.df <- estimates.df[,c(1:6,8:12,14:17)] }
        if(return.estimates=="ALL") { estimates.df <- estimates.df } }
    estimates.total <- ncol(estimates.df)
    
    # add identifier for each network projection
    if(iteration.type=="attribute") { identifier <- rep(attribute.unique, each = net.iterate) }
    if(iteration.type!="attribute") { identifier <- rep(net.samples, each = net.iterate) }
    estimates.df <- cbind(data.frame(sample=identifier), estimates.df)
    
    # plot data
    if(plot.estimators==TRUE) {
        # plot observed estimators
        if(nrow(estimates.df)<=500) { lwd.by.iteration <- 4}
        if(nrow(estimates.df)>500 && nrow(estimates.df)<=1000) { lwd.by.iteration <- 3}
        if(nrow(estimates.df)>1000 && nrow(estimates.df)<=2000) { lwd.by.iteration <- 2}
        if(nrow(estimates.df)>2000 && nrow(estimates.df)<=4000) { lwd.by.iteration <- 1}
        if(nrow(estimates.df)>4000) { lwd.by.iteration <- 0.3}
        colorsmetric <- rainbow(estimates.total+1)
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
        las.plot <- 0
        png(paste0("network_estimates_by_",removal,"_",net.iterate,"_iterations_over_",length(net.samples),"_projections_",iteration.type,"_",tolower(attribute),".png"), type='cairo', width=plot.panels[2]*4,height=plot.panels[1]*4, units='in', res=200)
        par(mfrow=plot.panels)
        if(iteration.type!="attribute") { 
            labels.plot1 <- 1:length(estimates.df$sample)
            labels.plot2 <- labels.plot2 <- paste0(estimates.df$sample*100,"%") }
        if(iteration.type=="attribute") { 
            labels.plot1 <- 1:length(estimates.df$sample)
            labels.plot2 <- paste0(estimates.df$sample)
            if(max(nchar(labels.plot2))>5) {
                labels.plot2[which(duplicated(estimates.df$sample))] <- ""
                las.plot <- 2 } 
        }
        for(i in 2:ncol(estimates.df)) { 
            plot(as.numeric(estimates.df[,i]), xlab="", ylab="", col=colorsmetric[i], cex=0.5, xaxt="n",
                 main=paste(colnames(estimates.df)[i]), type=plot.type, lwd=lwd.by.iteration,cex.lab=1.6, cex.axis=1.6, cex.main=2.5, cex.sub=2)
            axis(1, at=labels.plot1, labels=labels.plot2, las=las.plot)
        }
        dev.off()
    }
    print(paste0("Iteration completed. Plot saved at ",getwd(),"/network_estimates_by_",removal,"_",net.iterate,"_iterations_over_",length(net.samples),"_projections_",iteration.type,"_",tolower(attribute),".png"))
    return(estimates.df)
}
