# this code generates random bipartite networks that preserve the link and node numbers but rewire the links among nodes.
# then calculates the mean and standard deviation of network properties for multiple random networks.

code.wd="/Users/maggieyuan/Documents/!AMFnetwork/GitHub/"
source(paste(code.wd,"global.R",sep = ""))
source(paste(code.wd,"bpt.R",sep = ""))
source(paste(code.wd,"rand_adj_gen.R",sep = ""))

setwd("/Users/maggieyuan/Documents/!AMFnetwork/GitHub/")

# input: adjacency matrix
my_adj_mat = read.table("example_data/adjacency-matrix_nonAMF-bacteria.txt", sep="\t", row.names=1, header=T)

# generate random networks
times=3 # number of random networks to generate, depending on the size of the empirical network. 100 times randomization would be good for ~several hundred nodes.
rand_adj_mat_list = lapply(as.list(c(1:times)), FUN=rand_adj_gen, adj=my_adj_mat)
names(rand_adj_mat_list) = c(1:times)
length(rand_adj_mat_list) # check the number of random networks
dim(rand_adj_mat_list[[1]]) # check the dimension of the first random network

# analyze random networks
random_global = as.data.frame(lapply(rand_adj_mat_list, FUN=global))
names(random_global) = c(1:times)
random_bpt = as.data.frame(lapply(rand_adj_mat_list, FUN=bpt))
names(random_bpt) = c(1:times)
random_properties = as.data.frame(rbind(random_global,random_bpt))

# get mean and sd
means = rowMeans(random_properties)
sds = apply(random_properties, 1, FUN=sd)
