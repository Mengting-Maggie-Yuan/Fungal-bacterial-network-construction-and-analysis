# this code calculates network topological features, node and link attributes, and generates input files for visualization using Cytoscape and Gephi.

library(igraph)
library(brainGraph)
library(bipartite)

# functions
source("/Users/maggieyuan/Documents/!AMFnetwork/GitHub/fit_power_law.R")

add_node_attribute = function(graph){
  biomarker = gsub("OTU_.*_", "", vertex_attr(graph)$name)

  OTU_name = vertex_attr(graph)$name # OTU_name[1:10]
  node_degree = centr_degree(graph)$res

  module_separation = cluster_fast_greedy(graph)
  module_membership = membership(module_separation)
  pi = part_coeff(g = graph, memb = module_membership)
  zi = within_module_deg_z_score(g = graph, memb = module_membership)
  role = rep("peripherals", gorder(graph))
  role[which(pi >= 6.2 & zi < 2.5)] = "connector"
  role[which(pi >= 6.2 & zi >= 2.5)] = "network_hub"
  role[which(pi < 6.2 & zi >= 2.5)] = "module_hub"

  graph = set_vertex_attr(graph, "biomarker", index = V(graph), biomarker)
  graph = set_vertex_attr(graph, "node_degree", index = V(graph), node_degree)
  graph = set_vertex_attr(graph, "module_membership", index = V(graph), module_membership)
  graph = set_vertex_attr(graph, "pi", index = V(graph), pi)
  graph = set_vertex_attr(graph, "zi", index = V(graph), zi)
  graph = set_vertex_attr(graph, "vertex_role", index = V(graph), role)

  return(graph)
}

add_link_attribute = function(graph, tranmatx){
  # get positive vs negative links from transition matrix
  mat_p_n = tranmatx
  mat_p_n[mat_p_n>0]<-"positive"
  mat_p_n[mat_p_n<0]<-"negative"

  # get edge attribute of negative/positive
  edge_names = unlist(strsplit(attr(E(graph), "vnames"), "|", fixed=T))
  edge_rownames = edge_names[seq(1,(length(edge_names)-1), by=2)]
  edge_colnames = edge_names[seq(2,length(edge_names), by=2)]
  edge_dir = unlist(lapply(c(1:length(edge_rownames)),
                           FUN = function(id, np_matrix, rownames, colnames){
                             dir = np_matrix[which(row.names(np_matrix) == edge_rownames[id]), which(names(np_matrix) == edge_colnames[id])]
                           }, np_matrix = mat_p_n, rownames = edge_rownames, colnames = edge_colnames))

  graph = set_edge_attr(graph, "link_sign", index=E(graph), edge_dir)

  return(graph)
}

cyto_gephi_output = function(graph_with_all_attributes){
  node_table = as.data.frame(vertex.attributes(graph_with_all_attributes))

  edge_table = as.data.frame(edge.attributes(graph_with_all_attributes))
  edge_names = unlist(strsplit(attr(E(graph_with_all_attributes), "vnames"), "|", fixed=T))
  edge_table$node1 = edge_names[seq(from = 1, to = (length(edge_names) - 1), by = 2)]
  edge_table$node2 = edge_names[seq(from = 2, to = length(edge_names), by = 2)]

  gephi_node_table = node_table
  names(gephi_node_table)[1] <- "Id"

  gephi_edge_table = edge_table
  names(gephi_edge_table)[2:3] <- c("Source", "Target")
  gephi_edge_table$Type = rep("Undirected", nrow(gephi_edge_table))

write.table(node_table, "example_data/cytoscape_node_attribute.txt", sep = "\t", row.names = F, quote = F)
write.table(edge_table, "example_data/cytoscape_edge_attribute.txt", sep = "\t", row.names = F, quote = F)
write.table(gephi_node_table, "example_data/gephi_node_attribute.csv", sep = ",", row.names = F, quote = F)
write.table(gephi_edge_table, "example_data/gephi_edge_attribute.csv", sep = ",", row.names = F, quote = F)
}


########################

setwd("/Users/maggieyuan/Documents/!AMFnetwork/GitHub/")


  # input: transition matrix
  my_tranmatx = read.table("example_data/transition-matrix_nonAMF-bacteria.txt", sep="\t", row.names=1, header=T)

  p_link = sum(my_tranmatx>0)/2 # number of positive links
  n_link = sum(my_tranmatx<0)/2 # number of negative links
  prop_p_link = p_link/(p_link+n_link) # proportion of positive links

  adj_mat=my_tranmatx
  adj_mat[adj_mat!=0]<-1 # get adjacency matrix for the network. All links are indicated by 1 in adjacency matrix (for both + and - links).
  write.table(adj_mat, "example_data/adjacency-matrix_nonAMF-bacteria.txt", sep="\t")

  # construct graph
  my_graph = graph_from_adjacency_matrix(as.matrix(adj_mat), mode="undirected", weighted = NULL, diag = FALSE, add.colnames = NULL)

  # graph properties
  conn_nodes = gorder(my_graph) # node number
  links = ecount(my_graph) # link number
  r2 = fit_power_law(my_graph) # r2 of power law fit
  avgK = mean(centr_degree(my_graph)$res) # average degree
  avgCC = transitivity(my_graph, type="average", isolates="zero") # average clustering coefficient. Zero for bipartite network (no triangular subnetwork).
  GD = mean_distance(my_graph, directed=F, unconnected=T) # geodesic distance
  gd = cluster_fast_greedy(my_graph) # greedy module separation
  modules = length(gd) # number of modules
  M = modularity(gd) # modularity
  largest_connected = round(max(component_distribution(my_graph)*conn_nodes),0) # node number in the largest connected component

  # bipartite properties from the adjacency matrix
  id_16s = which(grepl("16S",names(adj_mat)))
  id_its = which(grepl("ITS",names(adj_mat)))
  bpt_matx = adj_mat[id_16s, id_its]
  bipartite_result = networklevel(index=c("connectance", "nestedness", "NODF", "weighted connectance", "number of species", "cluster coefficient", "niche overlap", "partner diversity"), bpt_matx)

  # add graph attributes
  my_graph = add_node_attribute(graph = my_graph)
  my_graph = add_link_attribute(graph = my_graph, tranmatx = my_tranmatx)

#  output cytoscape and gephi input files for visulization
  cyto_gephi_output(graph_with_all_attributes = my_graph)
