rm(list=ls())

code.wd="/Users/maggieyuan/Documents/shengjing/network/functions/"
source(paste(code.wd,"global.R",sep = ""))
source(paste(code.wd,"bpt.R",sep = ""))
source(paste(code.wd,"rand_adj_gen.R",sep = ""))

setwd("/Users/maggieyuan/Documents/shengjing/network/")
adj_path_list = list.files(path="matx_tran", recursive=T, full.names=T)
adj_path_list = grep("amf", adj_path_list, value=T)

adj_mat_list=list() # list containing all the transitional matrices
for (i in 1:length(adj_path_list))
{
	to_add_tran_matx = read.table(adj_path_list[i], sep="\t", header=T, row.names=1)
	to_add_adj_matx = to_add_tran_matx
	to_add_adj_matx[to_add_adj_matx !=0]  <-1
	adj_mat_list[[i]] <- to_add_adj_matx
	names(adj_mat_list)[i] <- adj_path_list[i]
}

# generate random networks
rand_adj_mat_list=list()
for (i in 1:length(adj_mat_list))
{
	times=2
	rand_adj_mat_list[((i-1)*100+1):(i*100)] = lapply(as.list(c(1:times)), FUN=rand_adj_gen, adj=adj_mat_list[[i]])
	names(rand_adj_mat_list)[((i-1)*100+1):(i*100)]=paste(names(adj_mat_list)[i], c(1:times), sep=".")
}
length(rand_adj_mat_list)

# analyze random networks
random_global = as.data.frame(lapply(rand_adj_mat_list, FUN=global))
random_bpt = as.data.frame(lapply(rand_adj_mat_list, FUN=bpt))

random_properties = as.data.frame(rbind(random_global,random_bpt))
write.table(random_properties, "network_properties/random_properties.txt", sep="\t")

# get mean and sd
random_properties2 = random_properties

means = data.frame(row.names = row.names(random_properties2))
sds = data.frame(row.names = row.names(random_properties2))
for (i in 1:(dim(random_properties2)[2]/100)){
  subdat = random_properties2[,c(((i-1)*100+1):(i*100))]
  subdat2g = matrix(as.numeric(as.matrix(subdat)),nrow=nrow(subdat))
  means = cbind(means,rowMeans(subdat2g))
  sds = cbind(sds,apply(subdat2g,1,FUN=sd))
}
names(means) = adj_path_list
names(sds) = adj_path_list

write.table(means, "network_properties/random_properties_means.txt", sep="\t")
write.table(sds, "network_properties/random_properties_sds.txt", sep="\t")
