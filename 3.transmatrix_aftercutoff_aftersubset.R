# this code is to get subnetworks by removing certain nodes or keeping certain nodes from a pool.
# input: cormatx - correlation matrix
# input: cutoff - correlation cutoff
# input: keep - list of nodes to keep
# input: rmv - list of nodes to remove
# input: outpath - path for output of transition matrix

rm(list=ls())
setwd("/Users/maggieyuan/Documents/shengjing/network_bacfun_maj10/")

cor2tran <- function(cormatx, cutoff, keep, rmv, outpath){
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
      
      write.table(tranmat_sub_noiso, paste(outpath, ".txt", sep=""), sep="\t")
    }
  }
  return(print(paste(outpath, "done")))
}

nwlistpaths_spearman = list.files(path="cor_matx/cormatx_spearman_bacfun/", pattern=".txt")
cut_off_table_spearman = read.table("rmt_cutoff_16S_ITS_spearman.txt", sep="\t", header=T, row.names=1)
nwlist_spearman = gsub(".txt", "", nwlistpaths_spearman)

nwlistpaths_pearson = list.files(path="cor_matx/cormatx_pearson_bacfun/", pattern=".txt")
cut_off_table_pearson = read.table("rmt_cutoff_16S_ITS_pearson.txt", sep="\t", header=T, row.names=1)
nwlist_pearson = gsub(".txt", "", nwlistpaths_pearson)

amf_node_list = as.character(as.matrix(read.table("/Users/maggieyuan/Documents/shengjing/data/OTUlist_Glom_in_network_22.txt")))

for (i in 1:18)
{
 cor_matrix_spearman = read.table(paste("cor_matx/cormatx_spearman_bacfun/", nwlistpaths_spearman[i], sep=""), sep="\t", row.names=1, header=T)
 cutoff_spearman = cut_off_table_spearman$fastest_transition[which(as.character(rownames(cut_off_table_spearman)) == nwlist_spearman[i])]
	
  cor_matrix_pearson = read.table(paste("cor_matx/cormatx_pearson_bacfun/", nwlistpaths_pearson[i], sep=""), sep="\t", row.names=1, header=T)
	cutoff_pearson = cut_off_table_pearson$fastest_transition[which(as.character(rownames(cut_off_table_pearson)) == nwlist_pearson[i])]
	
	cor2tran(cormatx=cor_matrix_spearman, cutoff=cutoff_spearman, rmv=c(), keep=c(), outpath=paste("tran_matx/tranmatx_spearman_bacfun/",nwlist_spearman[i], sep=""))
	cor2tran(cormatx=cor_matrix_pearson, cutoff=cutoff_pearson, rmv=c(), keep=c(), outpath=paste("tran_matx/tranmatx_pearson_bacfun/",nwlist_pearson[i], sep=""))
	cor2tran(cormatx=cor_matrix_pearson, cutoff=cutoff_pearson, rmv=amf_node_list, keep=c(), outpath=paste("tran_matx/tranmatx_pearson_nonamf/",nwlist_pearson[i], sep=""))
	cor2tran(cormatx=cor_matrix_pearson, cutoff=cutoff_pearson, rmv=c(), keep=amf_node_list, outpath=paste("tran_matx/tranmatx_pearson_amf/",nwlist_pearson[i], sep=""))
}
