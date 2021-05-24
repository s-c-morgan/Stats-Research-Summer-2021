library(igraph)
library(dplyr)
library(stringr)
library(readr)


authors <- read_table2("ca-sandi_auths.csv")
authlist<- unique(c(authors$node1, authors$node2))
g.authors<- graph.data.frame(authors[1:2], vertices = authlist, directed = FALSE)
E(g.authors)$weight<-authors$weight
plot(g.authors)
A<-get.adjacency(g.authors)

degree(g.authors)
weighted_degree<-E(g.authors)$weight*A







