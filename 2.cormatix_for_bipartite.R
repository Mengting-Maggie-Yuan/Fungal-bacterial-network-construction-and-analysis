# This code modifies the full correlation matrix to generate a new correlation matrix for a bipartite network.

setwd("/Users/maggieyuan/Documents/!AMFnetwork/GitHub/")

cor_matrix = read.table("example_data/correlation_matrix.txt", sep="\t", header=T, row.names=1)

id_16s = which(grepl("16S",names(cor_matrix)))
id_its = which(grepl("ITS",names(cor_matrix)))
length(id_16s)
length(id_its)

sub_cor_matrix = cor_matrix
sub_cor_matrix[id_16s, id_16s] <- 0
sub_cor_matrix[id_its, id_its] <- 0

write.table(sub_cor_matrix, file="example_data/correlation_matrix_bipartite.txt", sep="\t", quote = F)
