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
D_un<-diag(degree(g.authors))

W<-A
for (i in 1:nrow(authors)){
  W[as.character(authors$node1[i]), as.character(authors$node2[i])]<-authors$weight[i]
  W[as.character(authors$node2[i]), as.character(authors$node1[i])]<-authors$weight[i]
}

L_unnorm_w<-D-W
L_unnorm_u<-D_un-A

##Weighted unnormalized spectral clustering
evL_w <- eigen(L_unnorm_w, symmetric=TRUE)
plot(1:10,rev(evL_w$values)[1:10], log="y")

k   <- 2
Z   <- evL_w$vectors[,(ncol(evL_w$vectors)-k+1):ncol(evL_w$vectors)]

set.seed(123)
km <- kmeans(Z, centers=k, nstart=5)
unw2<-km$cluster
V(g.authors)$cluster<-km$cluster
V(g.authors)[cluster == 1]$color<-"red"
V(g.authors)[cluster == 2]$color<-"dodgerblue"
set.seed(123)
plot(g.authors)

k   <- 3
Z   <- evL_w$vectors[,(ncol(evL_w$vectors)-k+1):ncol(evL_w$vectors)]

set.seed(123)
km <- kmeans(Z, centers=k, nstart=5)
unw3<-km$cluster
V(g.authors)$cluster<-km$cluster
V(g.authors)[cluster == 1]$color<-"red"
V(g.authors)[cluster == 2]$color<-"dodgerblue"
V(g.authors)[cluster == 3]$color<-"yellow"
set.seed(123)
plot(g.authors)
##Unweighted unnormalized spectral clustering
evL_u <- eigen(L_unnorm_u, symmetric=TRUE)
plot(1:10,rev(evL_u$values)[1:10], log="y")

k   <- 2
Z   <- evL_u$vectors[,(ncol(evL_u$vectors)-k+1):ncol(evL_u$vectors)]

set.seed(123)
km <- kmeans(Z, centers=k, nstart=5)
unu2<-km$cluster
V(g.authors)$cluster<-km$cluster
V(g.authors)[cluster == 1]$color<-"red"
V(g.authors)[cluster == 2]$color<-"dodgerblue"
set.seed(123)
plot(g.authors)

k   <- 3
Z   <- evL_u$vectors[,(ncol(evL_u$vectors)-k+1):ncol(evL_u$vectors)]

set.seed(123)
km <- kmeans(Z, centers=k, nstart=5)
unu3<-km$cluster
V(g.authors)$cluster<-km$cluster
V(g.authors)[cluster == 1]$color<-"red"
V(g.authors)[cluster == 2]$color<-"dodgerblue"
V(g.authors)[cluster == 3]$color<-"yellow"
set.seed(123)
plot(g.authors)


D_neghalf<-1/sqrt(D)
D_neghalf[D_neghalf==Inf]<-0
L_sym_w<-D_neghalf %*% L_unnorm_w %*% D_neghalf

D_un_neghalf<-1/sqrt(D_un)
D_un_neghalf[D_un_neghalf==Inf]<-0
L_sym_u<-D_un_neghalf %*% L_unnorm_u %*% D_un_neghalf

##Weighted normalized symmetric spectral clustering
evL_sym_w <- eigen(L_sym_w, symmetric=TRUE)

k   <- 2
Z   <- evL_sym_w$vectors[,(ncol(evL_sym_w$vectors)-k+1):ncol(evL_sym_w$vectors)]

set.seed(123)
km <- kmeans(Z, centers=k, nstart=5)
V(g.authors)$cluster<-km$cluster
nsymw2<-km$cluster
V(g.authors)[cluster == 1]$color<-"red"
V(g.authors)[cluster == 2]$color<-"dodgerblue"
set.seed(123)
plot(g.authors)

k   <- 3
Z   <- evL_sym_w$vectors[,(ncol(evL_sym_w$vectors)-k+1):ncol(evL_sym_w$vectors)]

set.seed(123)
km <- kmeans(Z, centers=k, nstart=5)
V(g.authors)$cluster<-km$cluster
nsymw3<-km$cluster
V(g.authors)[cluster == 1]$color<-"red"
V(g.authors)[cluster == 2]$color<-"dodgerblue"
V(g.authors)[cluster == 3]$color<-"yellow"
set.seed(123)
plot(g.authors)

##Unweighted normalized symmetric spectral clustering

evL_sym_u <- eigen(L_sym_u, symmetric=TRUE)

k   <- 2
Z   <- evL_sym_u$vectors[,(ncol(evL_sym_u$vectors)-k+1):ncol(evL_sym_u$vectors)]

set.seed(123)
km <- kmeans(Z, centers=k, nstart=5)
V(g.authors)$cluster<-km$cluster
nsymu2<-km$cluster
V(g.authors)[cluster == 1]$color<-"red"
V(g.authors)[cluster == 2]$color<-"dodgerblue"
set.seed(123)
plot(g.authors)

k   <- 3
Z   <- evL_sym_u$vectors[,(ncol(evL_sym_u$vectors)-k+1):ncol(evL_sym_u$vectors)]

set.seed(123)
km <- kmeans(Z, centers=k, nstart=5)
V(g.authors)$cluster<-km$cluster
nsymu3<-km$cluster
V(g.authors)[cluster == 1]$color<-"red"
V(g.authors)[cluster == 2]$color<-"dodgerblue"
V(g.authors)[cluster == 3]$color<-"yellow"
set.seed(123)
plot(g.authors)


D_negone_w<-1/D
D_negone_w[D_negone_w==Inf]<-0
L_rw_w<-D_negone_w %*% L_unnorm_w

D_negone_u<-1/D_un
D_negone_u[D_negone_u==Inf]<-0
L_rw_w<-D_negone_u %*% L_unnorm_u

evL_rw_w <- eigen(L_rw_w, symmetric=FALSE)

k   <- 2
Z   <- evL_rw_w$vectors[,(ncol(evL_rw_w$vectors)-k+1):ncol(evL_rw_w$vectors)]

set.seed(123)
km <- kmeans(Z, centers=k, nstart=5)
V(g.authors)$cluster<-km$cluster
V(g.authors)[cluster == 1]$color<-"red"
V(g.authors)[cluster == 2]$color<-"dodgerblue"
set.seed(123)
plot(g.authors)

k   <- 3
Z   <- evL_rw_w$vectors[,(ncol(evL_rw_w$vectors)-k+1):ncol(evL_rw_w$vectors)]

set.seed(123)
km <- kmeans(Z, centers=k, nstart=5)
V(g.authors)$cluster<-km$cluster

V(g.authors)[cluster == 1]$color<-"red"
V(g.authors)[cluster == 2]$color<-"dodgerblue"
V(g.authors)[cluster == 3]$color<-"yellow"
set.seed(123)
plot(g.authors)







