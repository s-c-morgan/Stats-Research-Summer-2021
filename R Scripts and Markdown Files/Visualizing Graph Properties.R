library(igraph)
library(dplyr)
library(stringr)

rh<-read.csv("reddithyperlinks.tsv", sep = "\t")

g.rh<- graph.data.frame(rh)
is.simple(g.rh)
E(g.rh)$weight<-1
g.rh.weighted<-simplify(g.rh)

degree(g.rh)

hist(graph.strength(g.rh.weighted), col="pink", xlab="Vertex Strength", ylab="Frequency", xlim = c(0,200), breaks = seq(0,10000,by=10))
hist(degree(g.rh.weighted), col="blue", xlab="Vertex Degree", ylab="Frequency")

d.rh <- degree(g.rh.weighted)
dd.rh <- degree.distribution(g.rh.weighted)
d <- 1:max(d.rh)-1
ind <- (dd.rh != 0)
plot(d[ind], dd.rh[ind], log="xy", col="red", xlab=c("Log-Degree"), ylab=c("Log-Intensity"), main="Log-Log Degree Distribution")

head(sort(graph.strength(g.rh.weighted), decreasing=TRUE),25)
