# # Not run:
# # generate random data
# net.object <- watts.strogatz.game(1, 300, 4, 0.01)
# V(net.object)$group <- sample(rep(LETTERS[c(1,4,6,8,4,20)],1000),100)
# table(V(net.object)$group)

<<<<<<< HEAD
=======
temp <- iteratePaths(net.object, attribute="orgsector", from.node.group="EDU", to.node.group="IGOV", stepwise.removal=1)

>>>>>>> 020257b15807a4091e441498baf822258b018820
iteratePaths <- function(net.object,
                         attribute=NULL,
                         from.node.group=NULL,
                         to.node.group=NULL,
                         net.iterate = 10,
                         collapse.minor = 1,
                         iteration.type ="shortest.path",
                         net.samples = rev(seq(0.01:1,by=0.01)),
                         removal = "node",
                         stepwise.removal = "auto",
                         plot.estimators = TRUE,
                         plot.type = "p") {
    
    # check from & to
    if(is.null(from.node.group) | is.null(to.node.group)) { stop(print("iteratePaths require start and end node group")) }
    
    # load dependencies
    require(intergraph)
    require(network)
    require(lattice)
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
    
    # check node & edge attribute
    if(removal=="node" && is.null(igraph::get.vertex.attribute(corenet.g, attribute))) { stop(print(paste0("Node attribute ",attribute," not found."))) }
    if(removal=="node.interaction" && is.null(igraph::get.vertex.attribute(corenet.g, attribute))) { stop(print(paste0("Node attribute ",attribute," not found."))) }
    if(removal=="edge" && is.null(igraph::get.edge.attribute(corenet.g, attribute))) { stop(print(paste0("Edge attribute ",attribute," not found."))) } 
    
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
    
    # load divisors function
    divisors <- function(x) { 
        y <- seq_len(x); y[ x%%y == 0 ] }
    
    # define network sample by nodes or edges
    if(removal=="node") { net.size <- vcount(corenet.g) }
    if(removal=="edge") { net.size <- ecount(corenet.g) }
    
    # prepare for loop
    estimates.df <- data.frame()
    shortest.paths.list <- list()
    nodes.num.list <- list()
    edges.num.list <- list()
    
    # index network for attribute iteration
    if(iteration.type=="shortest.path") {
        
        # set attribute group
        V(corenet.g)$attr.to.iterate <- V(corenet.g)$name <- igraph::get.vertex.attribute(corenet.g, attribute)
        V(corenet.g)$index.num <- 1:vcount(corenet.g)
        
        if(removal=="node") { 
            attribute.index <- data.frame(nodes=igraph::V(corenet.g)$name, attribute=igraph::get.vertex.attribute(corenet.g, attribute, index=igraph::V(corenet.g)))
            attribute.index$attribute <- as.factor(as.character(attribute.index$attribute))
            attribute.index <- attribute.index[order(attribute.index$attribute),]
            # check and fix singletons
            if(collapse.minor>1) {
                attribute.index.table <- as.data.frame(table(as.character(attribute.index$attribute)))
                if(any(attribute.index.table$Freq<collapse.minor)) {
                    attribute.index <- merge(attribute.index, attribute.index.table, by.x="attribute", by.y="Var1", all.x=T)
                    attribute.index$attribute <- as.character(attribute.index$attribute)
                    attribute.index$attribute[attribute.index$Freq<collapse.minor & attribute.index$attribute!=from.node.group & attribute.index$attribute!=to.node.group] <- "ETAL"
                    attribute.index$Freq <- NULL
                }}
            net.samples.list <- list()
            attribute.unique <- unique(as.character(attribute.index$attribute))
            attribute.unique <- attribute.unique[!attribute.unique==from.node.group]
            attribute.unique <- attribute.unique[!attribute.unique==to.node.group]
            for(x in 1:length(unique(attribute.unique))) {
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
            # check and fix singletons
            if(collapse.minor>1) {
                attribute.index.table <- as.data.frame(table(as.character(attribute.index$attribute)))
                if(any(attribute.index.table$Freq<collapse.minor)) {
                    attribute.index <- merge(attribute.index, attribute.index.table, by.x="attribute", by.y="Var1", all.x=T)
                    attribute.index$attribute <- as.character(attribute.index$attribute)
                    attribute.index$attribute[attribute.index$Freq<collapse.minor & attribute.index$attribute!=from.node.group & attribute.index$attribute!=to.node.group] <- "ETAL"
                    attribute.index$Freq <- NULL
                }}
            net.samples.list <- list()
            attribute.unique <- unique(as.character(attribute.index$attribute))
            attribute.unique <- attribute.unique[!attribute.unique==from.node.group]
            attribute.unique <- attribute.unique[!attribute.unique==to.node.group]
            
            for(x in 1:length(unique(attribute.unique))) {
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
        
        # reset estimates
        shortest.paths.vec <- as.numeric()
        nodes.num.vec <- as.numeric()
        edges.num.vec <- as.numeric()

        # start iteration
        for(j in 1:net.iterate) {
            if(iteration.type=="shortest.path") {
                if(removal=="node") {
                    cat("\r","Iterative calculation",j,"of",net.iterate,"of shortest path between node types",from.node.group, "and", to.node.group, "by targetting the remainder node types")
                    nodes.deselect <- sample(which(V(corenet.g)$name==net.samples[u]), stepwise.removal*j)
                    corenet.gx <- igraph::induced.subgraph(corenet.g, which(!V(corenet.g)$index.num %in% nodes.deselect)) }
                if(removal=="edge") { 
                    cat("\r","Iterative removal of targeted edges",j,"of",net.iterate)
                    edges.deselect <- sample(net.samples.list[[u]], stepwise.removal*j)
                    edges.select <- E(corenet.g)$name[!E(corenet.g)$name %in% edges.deselect]
                    corenet.gx <- igraph::subgraph.edges(corenet.g, eids=which(igraph::E(corenet.g)$name %in% edges.select), delete.vertices=TRUE) }
            
            # collect metrics per iteration
            nodes.num.vec <- c(nodes.num.vec, vcount(corenet.gx))
            edges.num.vec <- c(edges.num.vec, ecount(corenet.gx))
            shortest.paths.vec <- c(shortest.paths.vec, mean(unlist(lapply(igraph::get.shortest.paths(corenet.gx, 
                                    from=V(corenet.gx)$name[V(corenet.gx)$attr.to.iterate==from.node.group], 
                                    to=V(corenet.gx)$name[V(corenet.gx)$attr.to.iterate==to.node.group], 
                                    mode = "all")$vpath, length))))
            }
        
        # aggregate estimates        
        shortest.paths.list[[u]] <- as.list(shortest.paths.vec)
        nodes.num.list[[u]] <- as.list(nodes.num.vec)
        edges.num.list[[u]] <- as.list(edges.num.vec)
        
        }
        
        # clear sample network
        rm(corenet.gx)
        # print process
        cat("\n")
        print(paste0("Completed iteration ",u," of ", length(net.samples),".")) 
    }

    # create response data frame
    if(iteration.type=="shortest.path") { identifier <- rep(attribute.unique, each = net.iterate) }
    if(iteration.type!="shortest.path") { identifier <- rep(net.samples, each = net.iterate) }
    estimates.df <- data.frame(removal=rep(1:net.iterate, length(net.samples)), 
                               percent=round(rep(1:net.iterate, length(net.samples))/vcount(corenet.g), digits = 2), 
                               group=identifier, nodes=unlist(nodes.num.vec), edges=unlist(edges.num.vec), 
                               shortest.path=unlist(shortest.paths.list))
    
    # plot data
    if(plot.estimators==TRUE) {
        lattice::xyplot(shortest.path~removal|group, data=estimates.df) }
    
    print(paste0("Iteration completed"))
    return(estimates.df)
}
