# this code is to get subnetworks by removing certain nodes or keeping certain nodes.

# input: cormatx - correlation matrix
# input: cutoff - correlation cutoff
# input: keep - list of nodes to keep
# input: rmv - list of nodes to remove

# output: transition matrix,
# where positive correlation above cutoff is 1, and negative correlation above cutoff is -1.
# All correlations weaker than cutoff are 0.
# OTUs with no link to other OTUs are removed from the matrix.

setwd("/Users/maggieyuan/Documents/!AMFnetwork/GitHub/")

cor2tran <- function(cormatx, cutoff, keep, rmv){
  tranmat=cormatx
  diag(tranmat)<-0
  tranmat[abs(tranmat)<cutoff]<-0
  nod_id_tranmat = which(rowSums(tranmat)!=0)
  tranmat_noiso = tranmat[nod_id_tranmat,nod_id_tranmat]

  keep_id = which(names(tranmat_noiso) %in% keep)
  rmv_id = which(names(tranmat_noiso) %in% rmv)

  if (length(rmv)>0 & length(keep)>0){
    print("Do not accept keep and remove at the same time")
  }else{

    if (length(keep)>0 & length(keep_id)==0){
      print("Cannot find node to keep.")
    }else{
      if (length(keep_id)>0){
        nb_id = c()
        for (i in 1:length(keep_id)){
          check_col = tranmat_noiso[,keep_id[i]]
          nb_id = c(nb_id, which(check_col != 0))
        }
        keep_id = unique(c(keep_id, nb_id))
        tranmat_sub = tranmat_noiso[keep_id,keep_id]
      }

      if (length(rmv_id)>0){
        tranmat_sub = tranmat_noiso[-rmv_id, -rmv_id]
      }

      if (length(rmv)>0 & length(rmv_id)==0){
        tranmat_sub = tranmat_noiso
      }

      if (length(rmv)==0 & length(keep)==0){
        tranmat_sub = tranmat_noiso
      }

      nod_id = which(rowSums(tranmat_sub)!=0)
      tranmat_sub_noiso = tranmat_sub[nod_id,nod_id]
    }
  }
  return(tranmat_sub_noiso)
}

amf_node_list = read.table("example_data/OTU_list_AMF-in-network_22.txt")$V1

my_cor_matrix = read.table("example_data/correlation_matrix_bipartite.txt")
my_cutoff = 0.757 # the correlation cutoff was determined in MENAP (http://ieg4.rccc.ou.edu/MENA/)

tran_matrix_fun_bac = cor2tran(cormatx=my_cor_matrix, cutoff=my_cutoff, rmv=c(), keep=c()) # transition matrix for the network containing all fungal-bacterial links
write.table(tran_matrix_fun_bac, "example_data/transition-matrix_fungi-bacteria.txt", sep="\t", quote=F)

tran_matrix_nonAMF = cor2tran(cormatx=my_cor_matrix, cutoff=my_cutoff, rmv=amf_node_list, keep=c()) # transition matrix for the network containing all nonAMF-bacterial links
write.table(tran_matrix_nonAMF, "example_data/transition-matrix_nonAMF-bacteria.txt", sep="\t", quote=F)

tran_matrix_AMF = cor2tran(cormatx=my_cor_matrix, cutoff=my_cutoff, rmv=c(), keep=amf_node_list) # transition matrix for the network containing all AMF-bacterial links
write.table(tran_matrix_AMF, "example_data/transition-matrix_AMF-bacteria.txt", sep="\t", quote=F)
