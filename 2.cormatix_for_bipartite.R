rm(list=ls())

setwd("/Users/maggieyuan/Documents/shengjing/network_bacfun_maj10/cor_adj_matx/cormatx_spearman_full/")
cor_list = list.files()

nwlist = cor_list

for (i in 1:length(nwlist))
{
  cor_matrix = read.table(nwlist[i], sep="\t", row.names=1, header=T)
  
  id_16s = which(grepl("16S",names(cor_matrix)))
  id_its = which(grepl("ITS",names(cor_matrix)))
  
  sub_cor_matrix = cor_matrix
  sub_cor_matrix[id_16s, id_16s] <- 0
  sub_cor_matrix[id_its, id_its] <- 0
  
#  fun_cor_matrix = cor_matrix
#  fun_cor_matrix[id_16s,] <- 0
#  fun_cor_matrix[,id_16s] <- 0
    
  write.table(sub_cor_matrix, file=paste("/Users/maggieyuan/Documents/shengjing/network_bacfun_maj10/cor_adj_matx/cormatx_spearman_bacfun/", nwlist[i],sep=""), sep="\t")
#  write.table(fun_cor_matrix, file=paste("sub_cor_matx/fungi_only/", nwlist[i],sep=""), sep="\t")
  
  print(i)
}
