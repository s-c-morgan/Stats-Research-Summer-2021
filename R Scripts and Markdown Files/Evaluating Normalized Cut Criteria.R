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

#evaluate normalized cut criteria
#or return block matrix g associated with membership and adjacency
evalObj <- function(membership, A=NULL, L=NULL, p.s=NULL, return_g=FALSE){
  groups <- sort(unique(membership))
  
  memMat <- do.call(cbind,
                    lapply(groups, function(group){
                      ifelse(membership==group, 1,0)
                    })
  )
  
  if (is.null(L)){
    L_p.s <- getL(A)
    p.s <- L_p.s$p.s
    L <- L_p.s$L
  }
  
  alphas <- t(p.s) %*% memMat
  
  #h <- memMat %*% diag(as.numeric(1/sqrt(alphas)))
  #g <- sapply(1:length(alphas), function(i) memMat[,i]*sqrt(p.s)/sqrt(alphas[i]))
  
  h <- sapply(1:length(alphas), function(i) memMat[,i]/sqrt(alphas[i]))
  g <- diag(sqrt(p.s)) %*% h
  
  if (return_g){
    return(g)
  }else{
    return(sum(diag(t(g) %*% L %*% g)))
  }
}

##Perform Balanced Erdos-Reyni Perturbation
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

##Quickly return the normalized cut criterion after a given perturbation
evalObj_after_perturbation<-function(graph, true_clusters, mode, prob, nrewire){
  if (mode == "degprev"){
    g_degprev <- graph %>% rewire(keeping_degseq(niter = nrewire))
    A_degprev<-as.matrix(get.adjacency(g_degprev))
    return(evalObj(true_clusters,A_degprev))
  }
  if (mode == "edgeprev"){
    g_edgeprev <- graph %>% rewire(each_edge(p = prob))
    A_edgeprev<-as.matrix(get.adjacency(g_edgeprev))
    return(evalObj(true_clusters,A_edgeprev))
  }
  if (mode == "randomer"){
    g_Y<-erdos.renyi.game(vcount(graph), prob, directed = TRUE)
    Y<-as.matrix(get.adjacency(g_Y))
    A_random<-abs(as.matrix(get.adjacency(graph))-Y)
    return(evalObj(true_clusters,A_random))
  }
  if (mode == "balanceder"){
    g_balanceder <- balanced_erdos_renyi_perturbation(graph, prob)
    A_balanceder<-as.matrix(get.adjacency(g_balanceder))
    return(evalObj(true_clusters,A_balanceder))
  }
}
#Simulation
##Create Ground Truth SBM
set.seed(123)
crossp<-rep(0,16)
B_crossp<-matrix(crossp, nrow = 4)
diag(B_crossp)<- .6
A_sbm<-genSBM(100,2,B_crossp)
g_truth<-graph_from_adjacency_matrix(A_sbm$A, mode = "directed")
evalObj(A_sbm$trueClust,A_sbm$A)


evalObj_after_perturbation(g_truth,A_sbm$trueClust,"randomer",.7)

f1<-function(x) evalObj_after_perturbation(g_truth,A_sbm$trueClust,mode = "randomer",prob = x,1)
f2<-function(x) evalObj_after_perturbation(g_truth,A_sbm$trueClust,mode = "edgeprev",prob = x,1)
f3<-function(x) evalObj_after_perturbation(g_truth,A_sbm$trueClust,"degprev",0.1,x)
f4<-function(x) evalObj_after_perturbation(g_truth,A_sbm$trueClust,"balanceder",x,1)

f1_sim_results<-as.data.frame(sapply(seq(0,1,0.01), function(x) sapply(rep(x,50),f1)))
f2_sim_results<-as.data.frame(sapply(seq(0,1,0.01), function(x) sapply(rep(x,50),f2)))
f4_sim_results<-as.data.frame(sapply(seq(0,1,0.01), function(x) sapply(rep(x,10),f4)))

tibble(p = seq(0,1,0.01), f1median = apply(f1_sim_results,2,median), f1q95 = apply(f1_sim_results, 2, function(x) quantile(x,0.95)),
       f1q5 = apply(f1_sim_results, 2, function(x) quantile(x,0.05)), f2median = apply(f2_sim_results,2,median), f2q95 = apply(f2_sim_results, 2, function(x) quantile(x,0.95)),
       f2q5 = apply(f2_sim_results, 2, function(x) quantile(x,0.05)),  f4median = apply(f4_sim_results,2,median), f4q95 = apply(f4_sim_results, 2, function(x) quantile(x,0.95)),
       f4q5 = apply(f4_sim_results, 2, function(x) quantile(x,0.05)))%>%
  ggplot() + 
  geom_ribbon(mapping = aes(x=p, ymin=f1q5, ymax=f1q95, fill = "Random ER Perturbations")) +
  geom_line(mapping = aes(x=p, y=f1median))+
  geom_line(mapping = aes(x=p, y=f2median))+
  geom_line(mapping = aes(x=p, y=f4median))+
  geom_ribbon(mapping = aes(x=p, ymin=f2q5, ymax=f2q95, fill = "Out Degree Preserving Perturbations")) +
  geom_ribbon(mapping = aes(x=p, ymin=f4q5, ymax=f4q95, fill = "Balanced ER Perturbations")) +
  labs(title = "Effect of Erdos-Renyi Perturbations on  Objective Cut Criterion",
       y = "Cut Criterion",
       x = "Probability")+
  ylim(0,4)+
  theme(legend.position = "bottom")

tibble(p = seq(0,1,0.01), median = apply(f2_sim_results,2,median), q95 = apply(f2_sim_results, 2, function(x) quantile(x,0.95)),
       q5 = apply(f2_sim_results, 2, function(x) quantile(x,0.05))) %>%
  ggplot() + 
  geom_line(mapping = aes(x=p, y=median, color = "Median across 50 simulations")) +
  geom_line(mapping = aes(x=p, y=q95, color = "95th Percentile across 50 simulations")) +
  geom_line(mapping = aes(x=p, y=q5, color = "5th Percentile across 50 simulations")) +
  labs(title = "Effect of Out-Degree-Preserving Perturbations on  Objective Cut Criterion",
       y = "Cut Criterion",
       x = "Probability")+
  theme(legend.position = "bottom")

tibble(p = seq(0,1,0.01), median = apply(f4_sim_results,2,median), q95 = apply(f4_sim_results, 2, function(x) quantile(x,0.95)),
       q5 = apply(f4_sim_results, 2, function(x) quantile(x,0.05))) %>%
  ggplot() + 
  geom_line(mapping = aes(x=p, y=median, color = "Median across 50 simulations")) +
  geom_line(mapping = aes(x=p, y=q95, color = "95th Percentile across 50 simulations")) +
  geom_line(mapping = aes(x=p, y=q5, color = "5th Percentile across 50 simulations")) +
  labs(title = "Effect of Edge-Preserving Perturbations on  Objective Cut Criterion",
       y = "Cut Criterion",
       x = "Probability")+
  theme(legend.position = "bottom")
