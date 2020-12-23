# Code to generate the theoretical interaction networks of figs 2 and 3

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
if(!require(DiagrammeR)) {install.packages("DiagrammeR"); library(DiagrammeR)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(reshape2)) {install.packages("reshape2"); library(reshape2)}
if(!require(DiagrammeRsvg)) {install.packages("DiagrammeRsvg"); library(DiagrammeRsvg)}
if(!require(rsvg)) {install.packages("rsvg"); library(rsvg)}
if(!require(ggsci)) {install.packages("ggsci"); library(ggsci)}

# settings ------------------------------
# whether to save plots
save_plots <- FALSE
# number of matrices
n <- 100
# connectance
connectance <- 1
# which norm to use
norm <- "l2"
# color palette
pal <- c(pal_material("red")(9)[8], pal_material("yellow")(9)[8], pal_material("blue")(9)[8])
# number of species
S <- 3
# number of fixed points to sample
fixed_points <- 100 * 2^(S-2)
# interaction strength
strength <- "low"
# type of interaction
interaction_all <- c("negative", "both", "positive")
# matrix to use
matrix_id <- 30

# plot networks for each interaction type ------------------------------
for (i in 1:length(interaction_all)) {
  # load matrices
  load(paste("data/random_matrices/", S, "sp/", n, "_random_matrices_",
             strength, "_strength_", connectance, "_connectance_", 
             interaction_all[i], "_sign.RData", sep = ""))
  # use one matrix as example
  A <- matrix_list[[matrix_id]]
  diag(A) <- NA
  # edge and nodes list
  edge_list <- melt(A)
  edge_list <- na.omit(edge_list)
  names(edge_list) <- c("from", "to", "value")
  nodes <- data.frame(id = as.integer(paste(1:S)))
  # add edge info
  edges <- edge_list %>% 
    mutate(color = "black",
           penwidth = 14 * abs(value),
           style = if_else(value > 0, 'solid', 'dashed'))
  # create network
  graph <- create_graph() %>%
    add_nodes_from_table(table = nodes) %>%
    select_nodes() %>%
    set_node_attrs_ws(node_attr = shape, value = c("oval", "square", "triangle")) %>%
    set_node_attrs_ws(node_attr = fillcolor, value = pal[i]) %>%
    set_node_attrs_ws(node_attr = color, value = "black") %>%
    set_node_attrs_ws(node_attr = fontsize, value = 18) %>%
    set_node_attrs_ws(node_attr = fontcolor, value = "black") %>%
    set_node_attrs_ws(node_attr = penwidth, value = 1) %>%
    set_node_attrs_ws(node_attr = width, value = c(0.5, 0.5, 0.6)) %>%
    set_node_attrs_ws(node_attr = height, value = c(0.5, 0.5, 0.6)) %>%
    add_edges_from_table(table = edges,
                         from_col = from,
                         to_col = to,
                         from_to_map = id_external) %>% 
    set_node_attr_to_display(attr = NULL)
  # save plot
  if (save_plots) {
    export_graph(graph = graph, 
                 file_name = paste("figs/fig2/network_", interaction_all[i], ".pdf", sep = ""))
  }
}
