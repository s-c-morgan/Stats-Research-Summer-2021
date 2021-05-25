library(igraph)
library(dplyr)
library(stringr)
library(readr)
library(stats)


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

W<-A
for (i in 1:nrow(authors)){
  W[as.character(authors$node1[i]), as.character(authors$node2[i])]<-authors$weight[i]
  W[as.character(authors$node2[i]), as.character(authors$node1[i])]<-authors$weight[i]
}

L_unnorm<-D-W

plot(1:10,rev(evL$values)[1:10], log="y")
evL <- eigen(L_unnorm, symmetric=TRUE)

k   <- 2
Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]

km <- kmeans(Z, centers=k, nstart=5)
V(g.authors)$cluster<-km$cluster
V(g.authors)[cluster == 1]$color<-"red"
V(g.authors)[cluster == 2]$color<-"dodgerblue"
plot(g.authors)

k   <- 3
Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]

km <- kmeans(Z, centers=k, nstart=5)
V(g.authors)$cluster<-km$cluster
V(g.authors)[cluster == 1]$color<-"red"
V(g.authors)[cluster == 2]$color<-"dodgerblue"
V(g.authors)[cluster == 3]$color<-"yellow"
plot(g.authors)


D_neghalf<-1/sqrt(D)
D_neghalf[D_neghalf==Inf]<-0
L_normsym<-D_neghalf %*% L_unnorm %*% D_neghalf

evL_sym <- eigen(L_normsym, symmetric=TRUE)

k   <- 2
Z   <- evL_sym$vectors[,(ncol(evL_sym$vectors)-k+1):ncol(evL_sym$vectors)]
km <- kmeans(Z, centers=k, nstart=5)
V(g.authors)$cluster<-km$cluster
V(g.authors)[cluster == 1]$color<-"red"
V(g.authors)[cluster == 2]$color<-"dodgerblue"
plot(g.authors)

k   <- 3
Z   <- evL_sym$vectors[,(ncol(evL_sym$vectors)-k+1):ncol(evL_sym$vectors)]

km <- kmeans(Z, centers=k, nstart=5)
V(g.authors)$cluster<-km$cluster
V(g.authors)[cluster == 1]$color<-"red"
V(g.authors)[cluster == 2]$color<-"dodgerblue"
V(g.authors)[cluster == 3]$color<-"yellow"
setseed(123)
plot(g.authors)









