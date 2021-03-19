library(igraph)
library(brainGraph)

source("/Users/maggieyuan/Documents/!AMFnetwork/GitHub/fit_power_law.R") #fit power law degree distribution

global <- function(adj){
	graph = graph_from_adjacency_matrix(as.matrix(adj), mode="undirected", weighted = NULL, diag = FALSE, add.colnames = NULL)
	size = gorder(graph)
	links = ecount(graph)
	r2 = fit_power_law(graph)
	avgK = mean(centr_degree(graph)$res)
	avgCC = transitivity(graph, type="average", isolates="zero")
	GD = mean_distance(graph, directed=F, unconnected=T)
	grdy = cluster_fast_greedy(graph)
	modules = length(grdy)
	M = modularity(grdy)

	values = c(size, links, r2, avgK, avgCC, GD, modules, M)
	names(values) = c("size", "links", "r2", "avgK", "avgCC", "GD", "modules", "M")
	return(values)
}
