#Add pseudo-outgroup to the resolved ingroup as input for astral 
library(phytools)

x=read.tree('ortho.tre')
for (i in 1:length(x)){
	tip<-list(edge=matrix(c(2,1),1,2),tip.label="outgroup", edge.length=1.0,Nnode=1)
	class(tip)<-"phylo"
	rt_edge=which(x[[i]]$edge[,1]==(length(x[[i]]$tip.label)+1))
	btree<-bind.tree(tip,x[[i]],where=2)
	y=multi2di(btree,random = F)
	write.tree(y,'new_ortho.trees',append=T)
}