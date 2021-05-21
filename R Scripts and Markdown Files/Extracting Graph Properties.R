library(igraph)

#Graph Basics
g<- graph.formula(1-2, 1-3, 2-3, 2-4, 3-5, 4-5, 4-6, 4-7, 5-6, 6-7)
plot(g)
V(g)
get.adjacency(g)

## Subsetting Graphs
h<-induced.subgraph(g, 1:5)
plot(h)

##Multi-graphs and weighting
mg <- g + edge(2,3)
plot(mg)
is.simple(mg)
str(mg)
E(mg)$weight<-1
wg2<-simplify(mg)
plot(wg2)
degree(wg2)

##special types of graphs
g.full <- graph.full(20)
plot(g.full)
g.tree <- graph.tree(7, children=2, mode="undirected")
g.tree2 <- graph.tree(7, children=2, mode="in")
plot(g.tree)
plot(g.tree2)