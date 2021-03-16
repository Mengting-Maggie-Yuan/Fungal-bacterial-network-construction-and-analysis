rm(list=ls())

maj_cut_cor_288 = function(grp, otu, keep, dataname)
{
	id = as.numeric(grp[4:19])
	abd = otu[,id]
	maj_id = apply(abd, 1, FUN=function(x){
		if (sum(x[]!=0) >= keep)
		return(1)
		else
		return(0)
	})
	maj = abd[which(maj_id[] == 1),]

	sp = cor(t(maj), method="spearman")

#	pearson = cor(t(maj), method = "pearson")

#	clr = apply((maj + 1), 2, function(xc){
#            log(xc, 2) - sum(log(xc, 2))/length(xc)
#	})
#	clr_pearson = cor(t(clr), method="pearson")
	
	write.table(sp, file=paste("cor_matx/", dataname, paste(grp[1], grp[2], grp[3], sep="_"), "_maj", as.character(keep), "_spearman.txt",sep=""), sep="\t")
#	write.table(pearson, file=paste("cor_matx/", dataname, paste(grp[1], grp[2], grp[3], sep="_"), "_maj", as.character(keep), "_pearson.txt",sep=""),sep="\t")
#	write.table(clr_pearson, file=paste("cor_matx/", dataname, paste(grp[1], grp[2], grp[3], sep="_"), "_maj", as.character(keep), "_clrpearson.txt",sep=""),sep="\t")
#	write.table(sp,file=paste("cor_matx/forJ/", dataname, paste(grp[1], grp[2], grp[3], sep="_"), "_maj", as.character(keep), "_spearman.txt",sep=""),sep="\t",row.names=F,col.names=F)
#	write.table(pearson,file=paste("cor_matx/forJ/", dataname, paste(grp[1], grp[2], grp[3], sep="_"), "_maj", as.character(keep), "_pearson.txt",sep=""),sep="\t",row.names=F, col.names =F)
#	write.table(clr_pearson,file=paste("cor_matx/forJ/", dataname, paste(grp[1], grp[2], grp[3], sep="_"), "_maj", as.character(keep), "_clrpearson.txt",sep=""),sep="\t",row.names=F, col.names =F)

}

maj_cut_cor_noplant = function(grp, otu, keep, dataname)
{
	id = as.numeric(grp[2:17])
	abd = otu[,id]
	maj_id = apply(abd, 1, FUN=function(x){
		if (sum(x[]!=0) >= keep)
		return(1)
		else
		return(0)
	})
	maj = abd[which(maj_id[] == 1),]
	
	sp = cor(t(maj), method="spearman")
#	pearson = cor(t(maj), method="pearson")
#	clr = apply((maj + 1), 2, function(xc){
#            log(xc, 2) - sum(log(xc, 2))/length(xc)
#	})
#	clr_pearson = cor(t(clr), method="pearson")

	write.table(sp, file=paste("cor_matx/", dataname, paste("S2_NP", grp[1], sep="_"), "_maj", as.character(keep), "_spearman.txt",sep=""), sep="\t")
#	write.table(pearson, file=paste("cor_matx/", dataname, paste("S2_NP", grp[1], sep="_"), "_maj", as.character(keep), "_pearson.txt",sep=""),sep="\t")
#	write.table(clr_pearson, file=paste("cor_matx/", dataname, paste("S2_NP", grp[1], sep="_"), "_maj", as.character(keep), "_clrpearson.txt",sep=""),sep="\t")
#	write.table(sp,file=paste("cor_matx_noplant/unplanted_forJ/", dataname, paste("S2_NP", grp[1], sep="_"), "_maj", as.character(keep), "_spearman.txt",sep=""),sep="\t",row.names=F,col.names=F)
#	write.table(pearson,file=paste("cor_matx_noplant/unplanted_forJ/", dataname, paste("S2_NP", grp[1], sep="_"), "_maj", as.character(keep), "_pearson.txt",sep=""),sep="\t",row.names=F, col.names =F)
#	write.table(clr_pearson,file=paste("cor_matx_noplant/unplanted_forJ/", dataname, paste("S2_NP", grp[1], sep="_"), "_maj", as.character(keep), "_clrpearson.txt",sep=""),sep="\t",row.names=F, col.names =F)
	
}

setwd("/Users/maggieyuan/Documents/shengjing/data")

otu_16s = read.table("16s_s1s2_resampled.txt", sep="\t", head=T, row.names=1)
otu_its = read.table("resampled_fungi_data.txt", sep="\t", head=T, row.names=1)
samplelist = read.table("samplelist_inorder.txt", sep="\t", head=T)

rownames(otu_16s) = paste(rownames(otu_16s), "16S", sep="_")
rownames(otu_its) = paste(rownames(otu_its), "ITS", sep="_")
otu_16s_its = rbind(otu_16s, otu_its)

# for 288 samples, balanced design
s288_id = c(which(samplelist$Plant=="P"), which(samplelist$Plant=="NP" & samplelist$week=="W0"))
otu_16s_its_s288 = otu_16s_its[, s288_id]
otu_its_s288 = otu_its[, s288_id]
otu_16s_s288 = otu_16s[, s288_id]
samplelist_s288 = samplelist[s288_id,]
samplelist_s288$ID_s288 <- c(1:nrow(samplelist_s288))

nwlist = aggregate(samplelist_s288 $ID_s288, by=list(samplelist_s288 $Season, samplelist_s288 $soil, samplelist_s288 $week), FUN=print)
nwids = apply(nwlist, 2, FUN=unlist)

setwd("/Users/maggieyuan/Documents/shengjing/16S_ITS_network")

apply(nwids, 1, otu = otu_its_s288, keep = 10, dataname = "ITS", FUN=maj_cut_cor_288)
apply(nwids, 1, otu = otu_16s_its_s288, keep = 10, dataname = "16S_ITS", FUN=maj_cut_cor_288)
apply(nwids, 1, otu = otu_16s_s288, keep = 10, dataname = "16S", FUN=maj_cut_cor_288)

# pick out unplanted
noplant_id = which(samplelist$Season=="S2" & samplelist$Plant=="NP" & samplelist$week!="W0")
otu_16s_its_noplant = otu_16s_its[,noplant_id]
otu_its_noplant = otu_its[,noplant_id]
otu_16s_noplant = otu_16s[,noplant_id]
samplelist_noplant = samplelist[noplant_id,]
samplelist_noplant$ID_noplant <- c(1:nrow(samplelist_noplant))

nwlist = aggregate(samplelist_noplant$ID_noplant, by=list(samplelist_noplant$week), FUN=print)
nwids = apply(nwlist, 2, FUN=unlist)

setwd("/Users/maggieyuan/Documents/shengjing/16S_ITS_network")

apply(nwids, 1, otu = otu_its_noplant, keep = 10, dataname = "ITS", FUN=maj_cut_cor_noplant)
apply(nwids, 1, otu = otu_16s_its_noplant, keep = 10, dataname = "16S_ITS", FUN=maj_cut_cor_noplant)
apply(nwids, 1, otu = otu_16s_noplant, keep = 10, dataname = "16S", FUN=maj_cut_cor_noplant)



