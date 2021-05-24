library(igraph)
library(dplyr)
library(stringr)
library(readr)


authors <- read_table2("ca-sandi_auths.csv")
auth_temp<-authors$node1
authors$node1<-authors$node2
authors$node2<-auth_temp
authlist<- unique(c(authors$node1, authors$node2))



g.authors<- graph.data.frame(authors[1:2], vertices = authlist, directed = FALSE)
E(g.authors)$weight<-authors$weight
plot(g.authors)
A<-get.adjacency(g.authors)
A<-as.matrix(A)
D<-diag(strength(g.authors))
U<-D-A







