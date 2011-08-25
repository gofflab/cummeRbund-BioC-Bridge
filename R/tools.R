# TODO: Add comment
# 
# Author: lgoff
###############################################################################

JSdist<-function(mat){
	res<-matrix(0,ncol=dim(mat)[2],nrow=dim(mat)[2])
	colnames(res)<-colnames(mat)
	rownames(res)<-colnames(mat)
	for(i in 1:dim(mat)[2]){
		for(j in i:dim(mat)[2]){
			a<-mat[,i]
			b<-mat[,j]
			JSdiv<-shannon.entropy((a+b)/2)-(shannon.entropy(a)+shannon.entropy(b))*0.5
			res[i,j] = sqrt(JSdiv)
			res[j,i] = sqrt(JSdiv)
		}
	}
	as.dist(res)
}


JSdistVec<-function(p,q){
	JSdiv<-shannon.entropy((p+q)/2)-(shannon.entropy(p)+shannon.entropy(q))*0.5
	JSdist<-sqrt(JSdiv)
	JSdist
}

makeprobsvec<-function(p){
	phat<-p/sum(p)
	phat[is.na(phat)] = 0
	phat
}

shannon.entropy <- function(p) {
	if (min(p) < 0 || sum(p) <=0)
		return(Inf)
	p.norm<-p[p>0]/sum(p)
	-sum( log10(p.norm)*p.norm)
}

makeprobs<-function(a){
	colSums<-apply(a,2,sum)
	b<-t(t(a)/colSums)
	b[is.na(b)] = 0
	b
}


#THIS IS NOT MINE....I MUST REMOVE IT PRIOR TO SUBMISSION (For detailed GO analysis, check out clusterProfiler and goProfiles)
#ClusterProfiles <- function(geneClusters, onto="CC", level=3, orgPackage="org.Hs.eg.db") {
#	require(goProfiles)
#	require(plyr)
#	require(ggplot2)
#	clusterProfile <- llply(geneClusters, as.data.frame(basicProfile), onto=onto, level=level, orgPackage = orgPackage)
#	clusterProfile.df <- ldply(clusterProfile, rbind)
#	colnames(clusterProfile.df) <- c("Cluster", "Description", "GOID", "Frequency")
#	clusterProfile.df <- clusterProfile.df[clusterProfile.df$Frequency !=0,]
#	clusterProfile.df$Description <- as.character(clusterProfile.df$Description) ## un-factor
#	clusterProfile.df <- ddply(clusterProfile.df, .(Description), transform, Percent = Frequency/sum(Frequency), Total = sum(Frequency))
#	
#	x <- mdply(clusterProfile.df[, c("Description", "Total")], paste, sep=" (")
#	y <- sapply(x[,3], paste, ")", sep="")
#	clusterProfile.df$Description <- y		### label GO Description with gene counts.
#	clusterProfile.df <-  clusterProfile.df[, -6] ###drop the *Total* column##
#	mtitle <- paste(onto, "Ontology Distribution", sep = " ")
#	p <- ggplot(clusterProfile.df, aes(x = Cluster, y = Description, size = Percent))
#	p <- p + geom_point(colour="steelblue") + opts(title = mtitle) + xlab("") + ylab("")
#	p <- p + opts(axis.text.x = theme_text(colour="black", size="11", vjust = 1))
#	p <- p + opts(axis.text.y = theme_text(colour="black", size="11", hjust = 1))
#	result <- list(data=clusterProfile.df, p=p)
#	return(result)
#}



