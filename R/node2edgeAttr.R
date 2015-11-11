node2edgeAttr <- function(net.object, vertex.attribute, directed=FALSE) {
    require(intergraph)
    require(igraph)
    require(network)
    if(class(net.object)=="network") { net.object <- intergraph::asIgraph(net.object) }
    df1 <- as.data.frame(igraph::get.edgelist(net.object))
    df1$V1 <- as.character(df1$V1)
    df1$V2 <- as.character(df1$V2)
    user.from <- merge(df1, data.frame(V1=as.character(V(net.object)), attr=igraph::get.vertex.attribute(net.object, vertex.attribute)), by="V1", all.x=T)[3]$attr
    user.to <- merge(df1, data.frame(V2=as.character(V(net.object)), attr=igraph::get.vertex.attribute(net.object, vertex.attribute)), by="V2", all.x=T)[3]$attr
    edge.attr <- paste(user.from, user.to, sep="+")
    if(directed==FALSE) {
        edge.attr.df <- as.data.frame(igraph::get.edgelist(igraph::graph.data.frame(data.frame(from=user.from, to=user.to), directed=F)))
        edge.attr <- paste(edge.attr.df$V1, edge.attr.df$V2, sep="+")       
    }
    return(edge.attr) }
