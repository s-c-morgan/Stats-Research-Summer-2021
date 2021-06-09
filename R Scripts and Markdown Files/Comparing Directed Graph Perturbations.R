library(igraph)
library(ggplot2)

set.seed(0)

rand_vec01 <- rbinom(400, 1, .25)
A <- matrix(rand_vec01, nrow = 20)
diag(A) <- 0

g <- graph_from_adjacency_matrix(A, mode = "directed")

coords <- layout_in_circle(g)

plot(g, layout = coords)

g_rewire_keepdeg <- g %>% rewire(keeping_degseq(niter = vcount(g) * 10))
g_rewire_er <- g %>% rewire(each_edge(.1))

plot(g_rewire_keepdeg, layout = coords)
plot(g_rewire_er, layout = coords)

A_rewire_keepdeg <- as.matrix(as_adjacency_matrix(g_rewire_keepdeg))
A_rewire_er <- as.matrix(as_adjacency_matrix(g_rewire_er))


tibble(deg_in = colSums(A), deg_out = rowSums(A), node = 1:20) %>%
  ggplot(mapping = aes(x=deg_in, y=deg_out, label=node)) +
  geom_text() +
  ggtitle("Original graph degree dist") +
  theme(plot.title = element_text(hjust = .5)) +
  scale_x_continuous(limits = c(0, 20), breaks = 0:10 * 2) +
  scale_y_continuous(limits = c(0, 20), breaks = 0:10 * 2)

tibble(deg_in = colSums(A_rewire_keepdeg), deg_out = rowSums(A_rewire_keepdeg), node = 1:20) %>%
  ggplot(mapping = aes(x=deg_in, y=deg_out, label=node)) +
  geom_text() +
  ggtitle("Degree-Preserving Rewired graph degree dist") +
  theme(plot.title = element_text(hjust = .5)) +
  scale_x_continuous(limits = c(0, 20), breaks = 0:10 * 2) +
  scale_y_continuous(limits = c(0, 20), breaks = 0:10 * 2)

tibble(deg_in = colSums(A_rewire_er), deg_out = rowSums(A_rewire_er), node = 1:20) %>%
  ggplot(mapping = aes(x=deg_in, y=deg_out, label=node)) +
  geom_text() +
  ggtitle("ER-Rewired graph degree dist") +
  theme(plot.title = element_text(hjust = .5)) +
  scale_x_continuous(limits = c(0, 20), breaks = 0:10 * 2) +
  scale_y_continuous(limits = c(0, 20), breaks = 0:10 * 2)