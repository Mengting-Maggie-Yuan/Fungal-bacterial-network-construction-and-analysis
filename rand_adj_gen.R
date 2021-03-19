# this function generates adjacency matrix for random networks based on the adjacency matrix for the empirical network.

rand_adj_gen <- function(ID, adj){
	id_16s = which(grepl("16S",names(adj)))
	id_its = which(grepl("ITS",names(adj)))
	network_L = sum(adj[id_16s, id_its])
	
	rand_adj = adj
	rand_adj[rand_adj !=0] <-0
	
	if (length(id_16s)>=length(id_its)){
		X=id_16s
		Y=id_its
		}else{
		Y=id_16s
		X=id_its			
		} #X needs to be a long edge of the rectangle
	
	rand_order = sample(c(Y, sample(Y, size=length(X)-length(Y), replace=T))) #generate random ids of the shorter edge where links should be put
	p = 1
	for (j in X){
	  rand_adj[j,rand_order[p]]<-1
	  p=p+1
	  } #put one link into each long edge, with randomized ids of the shorter edge
	
	expanded = as.numeric(as.matrix(rand_adj[X,Y])) #expand to a vector
	rand_remain_id = sample(which(expanded==0), size=(network_L-length(X))) #randomly select cells without a link
	expanded[rand_remain_id] <-1 #put the remaining links in these selected cells
	rand_adj[X,Y] = as.matrix(expanded, nrow=length(X), ncol=length(Y)) # fold back to a matrix
	
	rand_adj[Y,X]=t(rand_adj[X,Y]) #get the transpose to put in the other corner of the matrix
	if (sum(rand_adj)==sum(adj)){
		return(rand_adj)
	}else{
		return("Error")
		}
}
