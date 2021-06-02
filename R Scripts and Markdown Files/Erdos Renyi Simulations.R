library(igraph)
library(dplyr)
library(stringr)
library(readr)
library(stats)

er.sim<-erdos.renyi.game(15,0.2, directed = TRUE)
plot(er.sim)

g.full1<-graph.full(50)
g.full2<-graph.full(30)
g.full3<-graph.full(20)
g.full<-g.full1+g.full2+g.full3


plot(g.full)

edgelist<-as.data.frame(get.edgelist(g.full))
pairlist<-as.data.frame(t(combn(vcount(g.full),2)))

for (i in 1:nrow(pairlist)){
  seed<-runif(1)
  if (seed<=0.05){
    if (tail(duplicated(rbind(edgelist,c(pairlist[i,][,1],pairlist[i,][,2]))),1) 
        | tail(duplicated(rbind(edgelist,c(pairlist[i,][,2],pairlist[i,][,1]))),1) == TRUE){  
      g.full<-g.full-edges(c(pairlist[i,][,1],pairlist[i,][,2]))
    } else {
      g.full<-g.full+edges(c(pairlist[i,][,1],pairlist[i,][,2]))
    }
  }
}


erdosrenyi_perturb<-function(graph, p, directed){
  edgelist<-as.data.frame(get.edgelist(graph))
  pairlist<-as.data.frame(t(combn(vcount(graph),2)))
  seed_list<<-c(0)
  for (i in 1:nrow(pairlist)){
    seed<-runif(1)
    seed_list<<-append(seed_list,seed)
    if (seed<p){
      if (tail(duplicated(rbind(edgelist,c(pairlist[i,][,1],pairlist[i,][,2]))),1)
          | tail(duplicated(rbind(edgelist,c(pairlist[i,][,2],pairlist[i,][,1]))),1) == TRUE){  
        graph<-graph-edges(c(pairlist[i,][,1],pairlist[i,][,2]))
      } else {
        graph<-graph+edges(c(pairlist[i,][,1],pairlist[i,][,2]))
      }
    }
  }
  return(graph)
}