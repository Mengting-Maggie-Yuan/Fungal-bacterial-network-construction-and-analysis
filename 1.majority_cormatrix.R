# This code calculates the full correlation matrix for constructing the network from OTU tables.

setwd("/Users/maggieyuan/Documents/!AMFnetwork/GitHub/")

otu_16s = read.table("example_data/OTU-table_16S_Season1_Week0.txt", sep="\t", head=T, row.names=1)
otu_its = read.table("example_data/OTU-table_ITS_Season1_Week0.txt", sep="\t", head=T, row.names=1)
otu_16s_its = rbind(otu_16s, otu_its)

abd = otu_16s_its
abd[is.na(abd)] <-0

keep = 10 # the frequency of occurrence to keep an OTU for correlation calculation. In this study we kept OTUs present in at least 10 of the 16 biological replicates
counts = rowSums(abd>0)
abd = abd[counts>=keep,]

## choose correlation methods from: 1) Spearman; 2) Pearson; 3) central log-ratio Pearson
# 1) Spearman
abd_cor = cor(t(abd), method="spearman")

# 2) Pearson, used in this study
abd_cor = cor(t(abd), method = "pearson")

# 3) central log-ratio Pearson
clr = apply((abd + 1), 2, function(xc){
           log(xc, 2) - sum(log(xc, 2))/length(xc)
})
abd_cor = cor(t(clr), method="pearson")
##

write.table(abd_cor, file="example_data/correlation_matrix.txt", sep="\t", quote=F)
