# Code to generate the empirical interaction network of figs 4 and 5

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
# number of species
S <- 3

# plot network ------------------------------
# create matrix (this is system with species Ea, Pp, and Sm from Friedman et al (2017)
A <- matrix(c(-1.00, -0.82, -0.72,
              -0.87, -1.00, -0.84,
              -0.96, -0.91, -1.00), nrow = 3, byrow = TRUE)
diag(A) <- NA
# edge and nodes list
edge_list <- melt(A)
edge_list <- na.omit(edge_list)
names(edge_list) <- c("from", "to", "value")
nodes <- data.frame(id = as.integer(paste(1:S)))
# add edge info
edges <- edge_list %>% 
  mutate(color = "black",
         penwidth = 4 * abs(value),
         style = if_else(value > 0, 'solid', 'dashed'),
         tailclip = TRUE)
# create network
graph <- create_graph() %>%
  add_nodes_from_table(table = nodes) %>%
  select_nodes() %>%
  set_node_attrs_ws(node_attr = fillcolor, value = "white") %>%
  set_node_attrs_ws(node_attr = color, value = "black") %>%
  set_node_attrs_ws(node_attr = fontsize, value = 0) %>%
  set_node_attrs_ws(node_attr = fontcolor, value = "black") %>%
  set_node_attrs_ws(node_attr = penwidth, value = 0) %>%
  set_node_attrs_ws(node_attr = width, value = c(0.8, 0.8, 0.8)) %>%
  set_node_attrs_ws(node_attr = height, value = c(0.8, 0.8, 0.8)) %>%
  set_node_attrs_ws(node_attr = style, value = "invisible") %>%
  set_node_attr_to_display(attr = NULL) %>% 
  add_edges_from_table(table = edges,
                       from_col = from,
                       to_col = to,
                       from_to_map = id_external)
# save plot
if (save_plots) {
  export_graph(graph = graph, 
               file_name = "figs/fig4/Ea_Pp_Sm_network.pdf")
}
