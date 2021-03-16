rm(list=ls())
library(igraph)
library(RColorBrewer)
library(brainGraph)
library(bipartite)

# functions
source("/Users/maggieyuan/Documents/shengjing/network/functions/fit_power_law.R")

global <- function(graph) {
  
  conn_nodes = gorder(graph)
  links = ecount(graph)
  r2 = fit_power_law(graph)
  avgK_no_iso = mean(centr_degree(graph)$res)
  avgCC_no_iso = transitivity(graph, type="average", isolates="zero")
  GD = mean_distance(graph, directed=F, unconnected=T)
  gd_no_iso = cluster_fast_greedy(graph)
  modules_no_iso = length(gd_no_iso)
  M_no_iso = modularity(gd_no_iso)
  largest_connected = component_distribution(graph)
  
  values = c(conn_nodes, links, r2, avgK_no_iso, avgCC_no_iso, GD, modules_no_iso, M_no_iso)
  return(values)
}

add_node_attribute = function(graph){
  sequence = gsub("OTU_.*_", "", vertex_attr(graph)$name)
  
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
  
  graph = set_vertex_attr(graph, "sequence", index = V(graph), sequence)
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

cyto_gephi_output = function(graph_with_all_attributes, dataset_name){
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
  
#  write.table(node_table, paste("visualization/cytoscape_node_attribute_", dataset_name, ".txt", sep = ""), sep = "\t", row.names = F)
#  write.table(edge_table, paste("visualization/cytoscape_edge_attribute_", dataset_name, ".txt", sep = ""), sep = "\t", row.names = F)
#  write.table(gephi_node_table, paste("visualization/gephi_node_attribute_", dataset_name, ".csv", sep = ""), sep = ",", row.names = F)
#  write.table(gephi_edge_table, paste("visualization/gephi_edge_attribute_", dataset_name, ".csv", sep = ""), sep = ",", row.names = F)
  
  print(c(dataset_name, "output tables done."))
}


########################

setwd("/Users/maggieyuan/Documents/shengjing/network/")

# transition matrix
tranmatx_list = list.files(path="matx_tran", recursive = T)

global_title = c()

for (i in 1:length(tranmatx_list))
{
  my_tranmatx = read.table(paste("matx_tran/", tranmatx_list[i], sep=""), sep="\t", row.names=1, header=T)
  
  p_link = sum(my_tranmatx>0)/2
  n_link = sum(my_tranmatx<0)/2
  prop_p_link = p_link/(p_link+n_link)
  
  adj_mat=my_tranmatx
  adj_mat[adj_mat!=0]<-1
  
  # construct graph
  my_graph = graph_from_adjacency_matrix(as.matrix(adj_mat), mode="undirected", weighted = NULL, diag = FALSE, add.colnames = NULL)
  print("Graph constructed. Now calculating topological properties...")
  
  # global properties of the graph
  global_properties = global(graph=my_graph)
  
  # bipartite properties from the adjacency matrix
  id_16s = which(grepl("16S",names(adj_mat)))
  id_its = which(grepl("ITS",names(adj_mat)))
  bpt_matx = adj_mat[id_16s, id_its]
  if(length(dim(bpt_matx))>0){
    bipartite_result = networklevel(
      index=c("connectance", "nestedness", "NODF", "weighted connectance", "number of species", "cluster coefficient", "niche overlap", "partner diversity"), bpt_matx)
  }else{
    bipartite_result = rep("NA", 13)
  }
 
  global_title = rbind(global_title, c(global_properties, bipartite_result, p_link, n_link, prop_p_link))
  print("Global properties calculated. Now adding vertex and edge attributes...")
  
  # add graph attributes for networks without isolated nodes.
  my_graph = add_node_attribute(graph = my_graph)
  my_graph = add_link_attribute(graph = my_graph, tranmatx = my_tranmatx)
  print("Graph attributes added - without isolated nodes. Now exporting files for visualization...")
  
#  output cytoscape and gephi files
  dataset_name_clean = tranmatx_list[i]
  dataset_name_clean = gsub("/", "_", dataset_name_clean)
  cyto_gephi_output(graph_with_all_attributes = my_graph, dataset_name = dataset_name_clean)
  
  print(i)
}

colnames(global_title) <- c("size", "links", "r2_power_law", "avgK", "avgCC", "GD", "No.of_modules", "M", "connectance", "cluster coefficient", "nestedness", "NODF", "weighted connectance", "number.of.species.HL", "number.of.species.LL", "cluster.coefficient.HL", "cluster.coefficient.LL", "niche.overlap.HL", "niche.overlap.LL", "partner.diversity.HL", "partner.diversity.LL", "positive_links", "negative_links", "proportion.positive.links")
rownames(global_title) <- tranmatx_list
global_title

write.table(global_title, paste("temporary.txt", sep = ""), sep = "\t")
