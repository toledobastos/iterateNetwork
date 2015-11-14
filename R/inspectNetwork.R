inspectNetwork <- function(net.object,
                           iteration.type="random",
                           sample.size=.8) {
    
    # load dependencies
    require(intergraph)
    require(network)
    require(igraph)
    require(sna)
    
    # convert to igraph if object is network
    if(class(net.object)=="igraph") { corenet.g <- net.object }
    if(class(net.object)=="network") { corenet.g <- intergraph::asIgraph(net.object) }
    
    # check node names in network object
    if(is.null(V(corenet.g)$name)) { V(corenet.g)$name <- as.character(1:vcount(corenet.g)) }
    if(!is.character(V(corenet.g)$name)) { V(corenet.g)$name <- as.character(V(corenet.g)$name) }
    
    # define network size
    graph.size <- round((vcount(corenet.g)*sample.size)+.5, digits = 0)
    
    # sample network for inspection
    if(iteration.type=="degree"){
        nodes.select <- names(sort(igraph::degree(corenet.g), decreasing=T)[(igraph::vcount(corenet.g)-graph.size):igraph::vcount(corenet.g)])
        corenet.gx <- igraph::induced.subgraph(corenet.g, which(V(corenet.g)$name %in% nodes.select))
        print(paste0("Network sample includes ",vcount(corenet.gx)," nodes and ",ecount(corenet.gx)," edges.")) 
    }
    if(iteration.type=="random") { 
        corenet.gx <- igraph::induced.subgraph(corenet.g, which(igraph::V(corenet.g)$name %in% igraph::V(corenet.g)$name[sample(1:igraph::vcount(corenet.g), graph.size)]))
        print(paste0("Network sample includes ",vcount(corenet.gx)," nodes and ",ecount(corenet.gx)," edges.")) 
    }    
    if(iteration.type=="betweenness") {
        nodes.select <- names(sort(igraph::betweenness(corenet.g), decreasing=T)[(igraph::vcount(corenet.g)-graph.size):igraph::vcount(corenet.g)])
        corenet.gx <- igraph::induced.subgraph(corenet.g, which(V(corenet.g)$name %in% nodes.select))
        print(paste0("Network sample includes ",vcount(corenet.gx)," nodes and ",ecount(corenet.gx)," edges."))
    }
    if(iteration.type=="closeness") {
        nodes.select <- names(sort(igraph::closeness(corenet.g), decreasing=T)[(igraph::vcount(corenet.g)-graph.size):igraph::vcount(corenet.g)])
        corenet.gx <- igraph::induced.subgraph(corenet.g, which(V(corenet.g)$name %in% nodes.select))
        print(paste0("Network sample includes ",vcount(corenet.gx)," nodes and ",ecount(corenet.gx)," edges."))
    }
    return(corenet.gx)
}
