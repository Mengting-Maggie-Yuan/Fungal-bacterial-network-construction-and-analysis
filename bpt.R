library(igraph)
library(brainGraph)
library(bipartite)

bpt <- function(adj){
  id_16s = which(grepl("16S",names(adj)))
  id_its = which(grepl("ITS",names(adj)))
  bpt_matx = adj[id_16s, id_its]
  
  if(length(dim(bpt_matx))>0){
    bipartite_result = networklevel(
      index=c("connectance", "nestedness", "NODF", "weighted connectance", "number of species", "cluster coefficient", "niche overlap", "partner diversity"), bpt_matx)
  }else{
    bipartite_result = rep("NA", 13)
  }

  return(bipartite_result)
}
