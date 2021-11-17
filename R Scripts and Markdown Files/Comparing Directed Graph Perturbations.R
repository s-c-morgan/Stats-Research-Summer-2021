library(igraph)
library(dplyr)
library(stringr)
library(readr)
library(stats)
library(ggplot2)
library(rARPACK)
library(aricode)
# Functions
# Generate block model with N nodes, assigned to J blocks randomly, with connection probabilities within and between blocks given by B
genSBM <- function(N, J, B, groupProbs=NULL){
  gs <- table(sample(1:J, N, replace=T, prob=groupProbs))
  trueClust = rep(1:J, gs)
  
  A <- do.call(rbind, lapply(1:J, function(i){
    do.call(cbind,
            lapply(1:J, function(j){
              matrix(rbinom(gs[i]*gs[j], 1, B[i,j]), nrow = gs[i])
            }))
  }))
  
  diag(A) <- 0
  
  g <- graph_from_adjacency_matrix(A, mode="directed")
  isConnected <- components(g, mode = "strong")$no == 1
  
  return(list(
    A=A,
    trueClust=trueClust,
    isConnected=isConnected
  ))
}
#Reed's directed clustering algorithm
#get directed Laplacian
getL <- function(A){
  n <- nrow(A)
  
  d <- apply(A,1,sum) #weigthed out-degrees/strength in largest component
  P <- diag(1/d) %*% A #transition matrix
  
  #get stationary distrubution
  #note this is unique (have irreducibility, assume aperiodic)
  p <- Re(eigs(t(P), 1, "LM")$vectors[,1]) #unnormalized
  p.s <- p/sum(p) #stationary distrubution
  
  #create Laplacian and compute eigenvectors
  L <- diag(n) - (diag(sqrt(p.s)) %*% P %*% diag(1/sqrt(p.s)) + 
                    diag(1/sqrt(p.s)) %*% t(P) %*% diag(sqrt(p.s)))/2
  
  return(list(L=L, p.s=p.s))
}

#calculate laplacian from adjacency matrix and return eigenvectors and eigenvalues
getEigenDir <- function(A){
  return(eigen(getL(A)$L))
}

#return eigenvalues and eigenvector ratios
#entrywise ratios of each eigenvector divided by the first eigenvector
getEigsRatiosAndValues <- function(A){
  n <- nrow(A)
  L_p.s <- getL(A)
  eigenDir <- eigen(L_p.s$L)
  E.score <- (diag(1/sqrt(L_p.s$p.s)) %*% eigenDir$vectors)[,-n]
  eigs <- (E.score - matrix(1/n,nrow=n,ncol=n) %*% E.score) %*% diag(1/apply(E.score,2,sd))
  return(list(
    eigs=eigs,
    values=eigenDir$values
  ))
}

#return results of kmeans clustering on set of eigenvectors
getClust <- function(eigens, k){
  if (is.null(dim(eigens))){
    n <- length(eigens$values)
    km <- kmeans(eigens$vectors[,(n-(k-1)):(n)], centers=k, nstart=10)
  }else{
    n <- ncol(eigens)
    km <- kmeans(eigens[,(n-(k-2)):(n)], centers=k, nstart=10)
  }
  membership <- km$cluster
  return(membership)
}

#Simulation
##Create Ground Truth SBM
set.seed(123)
crossp<-rep(.005,16)
B_crossp<-matrix(crossp, nrow = 4)
diag(B_crossp)<- 0.95
A_sbm<-genSBM(60,4,B_crossp)
g_truth<-graph_from_adjacency_matrix(A_sbm$A, mode = "directed")
V(g_truth)$cluster<-A_sbm$trueClust
V(g_truth)[cluster == 1]$color<-"firebrick2"
V(g_truth)[cluster == 2]$color<-"dodgerblue"
V(g_truth)[cluster == 3]$color<-"gold1"
V(g_truth)[cluster == 4]$color<-"forestgreen"
coords <- layout_in_circle(g_truth)
plot(g_truth, layout = coords)
##Create Initial Clustering 
eigens<-getEigsRatiosAndValues(A_sbm$A)
g_init<-g_truth
V(g_init)$cluster<-getClust(eigens$eigs, 4)
V(g_init)[cluster == 1]$color<-"firebrick2"
V(g_init)[cluster == 2]$color<-"dodgerblue"
V(g_init)[cluster == 3]$color<-"gold1"
V(g_init)[cluster == 4]$color<-"forestgreen"
coords <- layout_in_circle(g_init)
plot(g_init, layout = coords)
##Light Perturbation
### Degree-Preserving Rewiring
g_depresv_lt<-g_init %>% rewire(keeping_degseq(niter = vcount(g_init)*5))
eigens<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_depresv_lt)))
V(g_depresv_lt)$cluster<-getClust(eigens$eigs, 4)
V(g_depresv_lt)[cluster == 1]$color<-"firebrick2"
V(g_depresv_lt)[cluster == 2]$color<-"dodgerblue"
V(g_depresv_lt)[cluster == 3]$color<-"gold1"
V(g_depresv_lt)[cluster == 4]$color<-"forestgreen"
coords <- layout_in_circle(g_depresv_lt)
plot(g_depresv_lt, layout = coords)
### Edge-Preserving Rewiring
g_edgepresv_lt<-g_init %>% rewire(each_edge(p = .1))
eigens<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_edgepresv_lt)))
V(g_edgepresv_lt)$cluster<-getClust(eigens$eigs, 4)
V(g_edgepresv_lt)[cluster == 1]$color<-"firebrick2"
V(g_edgepresv_lt)[cluster == 2]$color<-"dodgerblue"
V(g_edgepresv_lt)[cluster == 3]$color<-"gold1"
V(g_edgepresv_lt)[cluster == 4]$color<-"forestgreen"
coords <- layout_in_circle(g_edgepresv_lt)
plot(g_edgepresv_lt, layout = coords)
### Random Rewiring
g_Y_lt<-erdos.renyi.game(60, .1, directed = TRUE)
Y_lt<-as.matrix(get.adjacency(g_Y_lt))
random_lt<-abs(as.matrix(get.adjacency(g_init))-Y_lt)
g_random_lt<-graph_from_adjacency_matrix(random_lt, mode = "directed")
eigens<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_random_lt)))
V(g_random_lt)$cluster<-getClust(eigens$eigs, 4)
V(g_random_lt)[cluster == 1]$color<-"firebrick2"
V(g_random_lt)[cluster == 2]$color<-"dodgerblue"
V(g_random_lt)[cluster == 3]$color<-"gold1"
V(g_random_lt)[cluster == 4]$color<-"forestgreen"
coords <- layout_in_circle(g_random_lt)
plot(g_random_lt, layout = coords)
##Medium Pertubation
### Degree-Preserving Rewiring
g_depresv_med<-g_init %>% rewire(keeping_degseq(niter = vcount(g_init)*10))
eigens<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_depresv_med)))
V(g_depresv_med)$cluster<-getClust(eigens$eigs, 4)
V(g_depresv_med)[cluster == 1]$color<-"firebrick2"
V(g_depresv_med)[cluster == 2]$color<-"dodgerblue"
V(g_depresv_med)[cluster == 3]$color<-"gold1"
V(g_depresv_med)[cluster == 4]$color<-"forestgreen"
coords <- layout_in_circle(g_depresv_med)
plot(g_depresv_med, layout = coords)
### Edge-Preserving Rewiring
g_edgeprev_med<-g_init %>% rewire(each_edge(p = .25))
eigens<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_edgeprev_med)))
V(g_edgeprev_med)$cluster<-getClust(eigens$eigs, 4)
V(g_edgeprev_med)[cluster == 1]$color<-"firebrick2"
V(g_edgeprev_med)[cluster == 2]$color<-"dodgerblue"
V(g_edgeprev_med)[cluster == 3]$color<-"gold1"
V(g_edgeprev_med)[cluster == 4]$color<-"forestgreen"
coords <- layout_in_circle(g_edgeprev_med)
plot(g_edgeprev_med, layout = coords)
### Random Rewiring
g_Y_med<-erdos.renyi.game(60, .25, directed = TRUE)
Y_med<-as.matrix(get.adjacency(g_Y_med))
random_med<-abs(as.matrix(get.adjacency(g_init))-Y_med)
g_random_med<-graph_from_adjacency_matrix(random_med, mode = "directed")
eigens<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_random_med)))
V(g_random_med)$cluster<-getClust(eigens$eigs, 4)
V(g_random_med)[cluster == 1]$color<-"firebrick2"
V(g_random_med)[cluster == 2]$color<-"dodgerblue"
V(g_random_med)[cluster == 3]$color<-"gold1"
V(g_random_med)[cluster == 4]$color<-"forestgreen"
coords <- layout_in_circle(g_random_med)
plot(g_random_med, layout = coords)
##Heavy Pertubation
### Degree-Preserving Rewiring
g_depresv_hev<-g_init %>% rewire(keeping_degseq(niter = vcount(g_init)*15))
eigens<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_depresv_hev)))
V(g_depresv_hev)$cluster<-getClust(eigens$eigs, 4)
V(g_depresv_hev)[cluster == 1]$color<-"firebrick2"
V(g_depresv_hev)[cluster == 2]$color<-"dodgerblue"
V(g_depresv_hev)[cluster == 3]$color<-"gold1"
V(g_depresv_hev)[cluster == 4]$color<-"forestgreen"
coords <- layout_in_circle(g_depresv_hev)
plot(g_depresv_hev, layout = coords)
### Edge-Preserving Rewiring
g_edgeprev_hev<-g_init %>% rewire(each_edge(p = .35))
eigens<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_edgeprev_hev)))
V(g_edgeprev_hev)$cluster<-getClust(eigens$eigs, 4)
V(g_edgeprev_hev)[cluster == 1]$color<-"firebrick2"
V(g_edgeprev_hev)[cluster == 2]$color<-"dodgerblue"
V(g_edgeprev_hev)[cluster == 3]$color<-"gold1"
V(g_edgeprev_hev)[cluster == 4]$color<-"forestgreen"
coords <- layout_in_circle(g_edgeprev_hev)
plot(g_edgeprev_hev, layout = coords)
### Random Rewiring
g_Y_hev<-erdos.renyi.game(60, .35, directed = TRUE)
Y_hev<-as.matrix(get.adjacency(g_Y_hev))
random_hev<-abs(as.matrix(get.adjacency(g_init))-Y_hev)
g_random_hev<-graph_from_adjacency_matrix(random_hev, mode = "directed")
eigens<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_random_hev)))
V(g_random_hev)$cluster<-getClust(eigens$eigs, 4)
V(g_random_hev)[cluster == 1]$color<-"firebrick2"
V(g_random_hev)[cluster == 2]$color<-"dodgerblue"
V(g_random_hev)[cluster == 3]$color<-"gold1"
V(g_random_hev)[cluster == 4]$color<-"forestgreen"
coords <- layout_in_circle(g_random_hev)
plot(g_random_hev, layout = coords)

calc_NMI_after_perturbation<-function(graph, clusters, mode, prob, nrewire){
  if (mode == "degprev"){
    g_degprev <- graph %>% rewire(keeping_degseq(niter = nrewire))
    eigens<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_degprev)))
    V(g_degprev)$cluster<-getClust(eigens$eigs, clusters)
    return(NMI(V(graph)$cluster,V(g_degprev)$cluster))
  }
  if (mode == "outdegprev"){
    g_outdegprev <- graph %>% rewire(each_edge(p = prob))
    eigens<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_outdegprev)))
    V(g_outdegprev)$cluster<-getClust(eigens$eigs, clusters)
    return(NMI(V(graph)$cluster,V(g_outdegprev)$cluster))
  }
  if (mode == "randomer"){
    g_Y<-erdos.renyi.game(vcount(graph), prob, directed = TRUE)
    Y<-as.matrix(get.adjacency(g_Y))
    random<-abs(as.matrix(get.adjacency(graph))-Y)
    g_randomer<-graph_from_adjacency_matrix(random, mode = "directed")
    eigens<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_randomer)))
    V(g_randomer)$cluster<-getClust(eigens$eigs, clusters)
    return(NMI(V(graph)$cluster,V(g_randomer)$cluster))
  }
  if (mode == "balanceder"){
    g_balanceder <- balanced_erdos_renyi_perturbation(graph, prob)
    eigens<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_balanceder)))
    V(g_balanceder)$cluster<-getClust(eigens$eigs, clusters)
    return(NMI(V(graph)$cluster,V(g_balanceder)$cluster))
  }
}



f1<-function(x) calc_NMI_after_perturbation(g_truth,4,mode = "randomer",prob = x,1)
f2<-function(x) calc_NMI_after_perturbation(g_truth,4,mode = "outdegprev",prob = x,1)
f3<-function(x) calc_NMI_after_perturbation(g_truth,4,"degprev",0.1,x)
f4<-function(x) calc_NMI_after_perturbation(g_truth,4,"balanceder",x,1)

tibble(p = seq(0,1,0.01), NMI = sapply(seq(0,1,0.01),f1))%>%
  ggplot(mapping = aes(x=p, y=NMI))+
  geom_line()+
  xlim(0,1)

tibble(p = seq(0,1,0.01), NMI = sapply(seq(0,1,0.01),f2))%>%
  ggplot(mapping = aes(x=p, y=NMI))+
  geom_line()+
  xlim(0,1)

tibble(p = seq(1,100,10), NMI = sapply(seq(1,100,1),f3))%>%
  ggplot(mapping = aes(x=p, y=NMI))+
  geom_line()

f1_sim_results<-as.data.frame(sapply(seq(0,1,0.01), function(x) sapply(rep(x,50),f1)))
f2_sim_results<-as.data.frame(sapply(seq(0,1,0.01), function(x) sapply(rep(x,50),f2)))
f4_sim_results<-as.data.frame(sapply(seq(0,1,0.01), function(x) sapply(rep(x,10),f4)))


tibble(p = seq(0,1,0.01), median = apply(f1_sim_results,2,median), q95 = apply(f1_sim_results, 2, function(x) quantile(x,0.95)),
       q5 = apply(f1_sim_results, 2, function(x) quantile(x,0.05))) %>%
  ggplot() + 
  geom_line(mapping = aes(x=p, y=median, color = "Median across 50 simulations")) +
  geom_line(mapping = aes(x=p, y=q95, color = "95th Percentile across 50 simulations")) +
  geom_line(mapping = aes(x=p, y=q5, color = "5th Percentile across 50 simulations")) +
  labs(title = "Effect of Erdos-Renyi Perturbations on Normalized Mutual Information Scores",
       y = "Normalized Mutual Information",
       x = "Probability")+
  theme(legend.position = "bottom")

tibble(p = seq(0,1,0.01), median = apply(f2_sim_results,2,median), q95 = apply(f2_sim_results, 2, function(x) quantile(x,0.95)),
       q5 = apply(f2_sim_results, 2, function(x) quantile(x,0.05))) %>%
  ggplot() + 
  geom_line(mapping = aes(x=p, y=median, color = "Median across 50 simulations")) +
  geom_line(mapping = aes(x=p, y=q95, color = "95th Percentile across 50 simulations")) +
  geom_line(mapping = aes(x=p, y=q5, color = "5th Percentile across 50 simulations")) +
  labs(title = "Effect of Out-Degree-Preserving Perturbations on Normalized Mutual Information Scores",
       y = "Normalized Mutual Information",
       x = "Probability")+
  theme(legend.position = "bottom")

tibble(p = seq(0,1,0.01), median = apply(f4_sim_results,2,median), q95 = apply(f4_sim_results, 2, function(x) quantile(x,0.95)),
       q5 = apply(f4_sim_results, 2, function(x) quantile(x,0.05))) %>%
  ggplot() + 
  geom_line(mapping = aes(x=p, y=median, color = "Median across 50 simulations")) +
  geom_line(mapping = aes(x=p, y=q95, color = "95th Percentile across 50 simulations")) +
  geom_line(mapping = aes(x=p, y=q5, color = "5th Percentile across 50 simulations")) +
  labs(title = "Effect of Edge-Preserving Perturbations on Normalized Mutual Information Scores",
       y = "Normalized Mutual Information",
       x = "Probability")+
  theme(legend.position = "bottom")
  



##How are eigenvalues perturbed by the three types of perturbation techniques
set.seed(123)
crossp<-rep(.01,9)
B_crossp<-matrix(crossp, nrow = 3)
diag(B_crossp)<- 0.90
A_sbm<-genSBM(100,3,B_crossp)
g_truth<-graph_from_adjacency_matrix(A_sbm$A, mode = "directed")
V(g_truth)$cluster<-A_sbm$trueClust
V(g_truth)[cluster == 1]$color<-"firebrick2"
V(g_truth)[cluster == 2]$color<-"dodgerblue"
V(g_truth)[cluster == 3]$color<-"gold1"
layout <- layout.fruchterman.reingold(g_truth)
plot(g_truth, layout=layout)

eigens1<-getEigsRatiosAndValues(A_sbm$A)
g_init<-g_truth

g_Y_low<-erdos.renyi.game(100, .25, directed = TRUE)
Y_low<-as.matrix(get.adjacency(g_Y_low))
random_low<-abs(as.matrix(get.adjacency(g_init))-Y_low)
g_random_low<-graph_from_adjacency_matrix(random_low, mode = "directed")
reed_eigens_random_low<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_random_low)))

g_Y_mid<-erdos.renyi.game(100, .35, directed = TRUE)
Y_mid<-as.matrix(get.adjacency(g_Y_mid))
random_mid<-abs(as.matrix(get.adjacency(g_init))-Y_mid)
g_random_mid<-graph_from_adjacency_matrix(random_mid, mode = "directed")
reed_eigens_random_mid<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_random_mid)))

g_Y_hev<-erdos.renyi.game(100, .45, directed = TRUE)
Y_hev<-as.matrix(get.adjacency(g_Y_hev))
random_hev<-abs(as.matrix(get.adjacency(g_init))-Y_hev)
g_random_hev<-graph_from_adjacency_matrix(random_hev, mode = "directed")
reed_eigens_random_hev<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_random_hev)))

g_depresv_lt<-g_init %>% rewire(keeping_degseq(niter = vcount(g_init)*5))
reed_eigens_dp_low<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_depresv_lt)))

g_depresv_med<-g_init %>% rewire(keeping_degseq(niter = vcount(g_init)*10))
reed_eigens_dp_med<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_depresv_med)))

g_depresv_high<-g_init %>% rewire(keeping_degseq(niter = 2000))
reed_eigens_dp_high<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_depresv_high)))

g_edgeprev_low_low<-g_init %>% rewire(each_edge(p = .05))
reed_eigens_ep_low_low<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_edgeprev_low_low)))

g_edgeprev_low<-g_init %>% rewire(each_edge(p = .25))
reed_eigens_ep_low<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_edgeprev_low)))

g_edgeprev_med<-g_init %>% rewire(each_edge(p = .35))
reed_eigens_ep_med<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_edgeprev_med)))

g_edgeprev_hev<-g_init %>% rewire(each_edge(p = .45))
reed_eigens_ep_high<-getEigsRatiosAndValues(as.matrix(get.adjacency(g_edgeprev_hev)))

ggplot()+
  geom_line(mapping = aes(x=1:10, y=rev(reed_eigens_random_hev$values)[1:10], color = "p =.45"))+
  geom_line(mapping = aes(x=1:10, y=rev(reed_eigens_random_mid$values)[1:10], color = "p = .35"))+
  geom_line(mapping = aes(x=1:10, y=rev(reed_eigens_random_low$values)[1:10], color = "p = .25"))+
  geom_line(mapping = aes(x=1:10, y=rev(eigens1$values)[1:10], color = "No Perturbations"))+
  labs(title = "Effect of Erdos-Renyi Random Perturbations on Eigenvalues of Reed's Laplacian",
       y = "Eigenvalue magnitude",
       x = "nth smallest Eigenvalue")+
  theme(legend.position = "bottom")


ggplot()+
  geom_line(mapping = aes(x=1:10, y=rev(reed_eigens_dp_high$values)[1:10], color = "2000 Perturbations"))+
  geom_line(mapping = aes(x=1:10, y=rev(reed_eigens_dp_med$values)[1:10], color = "1000 Perturbations"))+
  geom_line(mapping = aes(x=1:10, y=rev(reed_eigens_dp_low$values)[1:10], color = "500 Perturbation"))+
  geom_line(mapping = aes(x=1:10, y=rev(eigens1$values)[1:10], color = "No Perturbations"))+
  labs(title = "Effect of Degree-Preserving Perturbations on Eigenvalues of Reed's Laplacian",
       y = "Eigenvalue magnitude",
       x = "nth smallest Eigenvalue")+
  theme(legend.position = "bottom")

ggplot()+
  geom_line(mapping = aes(x=1:10, y=rev(reed_eigens_ep_high$values)[1:10], color = "p = .45"))+
  geom_line(mapping = aes(x=1:10, y=rev(reed_eigens_ep_med$values)[1:10], color = "p = .35"))+
  geom_line(mapping = aes(x=1:10, y=rev(reed_eigens_ep_low$values)[1:10], color = "p = .25"))+
  geom_line(mapping = aes(x=1:10, y=rev(reed_eigens_ep_low_low$values)[1:10], color = "p = .05"))+
  geom_line(mapping = aes(x=1:10, y=rev(eigens1$values)[1:10], color = "No Perturbations"))+
  labs(title = "Effect of Edge-Preserving Perturbations on Eigenvalues of Reed's Laplacian",
       y = "Eigenvalue magnitude",
       x = "nth smallest Eigenvalue")+
  theme(legend.position = "bottom")


balanced_erdos_renyi_perturbation<-function(initial_graph, p){
  ogm<-as.vector(get.adjacency(initial_graph))
  onepos<-which(ogm==1)
  zeropos<-which(ogm==0)
  zeropos<-ifelse(zeropos %in% seq(1,vcount(initial_graph)^2,vcount(initial_graph)+1),NA,zeropos)
  one_flipper<-sample(onepos, p*ecount(initial_graph))
  zero_flipper<-sample(zeropos[!is.na(zeropos)], p*ecount(initial_graph))
  for (n in one_flipper){ogm[n]<-0}
  for (m in zero_flipper){ogm[m]<-1}
  updated_ogm<-matrix(ogm,nrow=vcount(initial_graph))
  g_perturb<-graph_from_adjacency_matrix(updated_ogm)
  g_perturb
}

f4_sims<-c()
for (i in 1:100){
  f4_sims<-append(f4_sims,f4(.6))
}