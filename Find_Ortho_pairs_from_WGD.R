#pair orthologs derived from multiple WGD events from phylogeny trees, generate input for McScan

library(phytools)
library(base)
library(geiger)

treefile_names <- list.files(pattern = ".*inclade*")

#define function

delete_outlier<-function(phy,clade){
	nt=compute.brlen(phy, 1)
	pd=cophenetic.phylo(nt) #pairwise distance given topology only
	dis=list()
	for (i in clade){dis[[i]]=sum(pd[i,clade])}
	out=names(dis[which.max(as.integer(dis))])
	clade=clade[!clade %in% out]
	return(clade)
}

return_monophyletic_gr<-function(phy,clade){
	#remove up to three outliers in a clade
	if (is.monophyletic(phy,clade)){return(clade)}
	else{
		clade=delete_outlier(phy,clade)
		if (is.monophyletic(phy,clade)){return(clade)}	#remove 1st outlier
		else{
			clade=delete_outlier(phy,clade)
			if (is.monophyletic(phy,clade)){return(clade)}	#remove 2nd outlier
			else{
				clade=delete_outlier(phy,clade)
				return(clade)	#remove 3rd outlier, but this can still be paraphyletic
			}
		}
	}
}

pair_ortho_copies <- function(cladeA,cladeB,filename){
	for (i in cladeA){
		sapply(cladeB,function(j) write(paste(i,j),filename,append=T))
	} 
}

#loop through trees
for (i in treefile_names){ 
	tree=read.tree(i)
	taxa=tree$tip.label[which(startsWith(tree$tip.label,'Lin'))]
	taxa=return_monophyletic_gr(tree,taxa)
	if ((length(taxa)==4) | (length(taxa)==3)){
		if (length(taxa)==4){
			node=getMRCA(tree,taxa)
			daughters=which(tree$edge[,1]==node)
			clade1 = tips(tree,tree$edge[daughters[1],2])
			clade2 = tips(tree,tree$edge[daughters[2],2])
			if (length(clade1)==2){
				#((a1,a2),(b1,b2))
				pair_ortho_copies(clade1,clade2,'old_WGD.gene')
				write(paste(clade1,collapse = ' '),'new_WGD.gene',append=T)
				write(paste(clade2,collapse = ' '),'new_WGD.gene',append=T)
			}
			else {
				##(ax,(ay,(b1,b2)))
				if (length(clade1)<length(clade2)){
					sub1=clade1
					sub2=clade2
				}
				else{
					sub1=clade2
					sub2=clade1
				}
				node=getMRCA(tree,sub2)
				daughters=which(tree$edge[,1]==node)
				clade3 = tips(tree,tree$edge[daughters[1],2])
				clade4 = tips(tree,tree$edge[daughters[2],2])
				if (length(clade3)>length(clade4)){sub3=clade3}
				else{sub3=clade4}
				pair_ortho_copies(sub1,sub3,'old_WGD.gene')
				write(paste(sub3,collapse = ' '),'new_WGD.gene',append=T)
			}
		}
		else{
			#(ax,(b1,b2))
			node=getMRCA(tree,taxa)
			daughters=which(tree$edge[,1]==node)
			clade1 = tips(tree,tree$edge[daughters[1],2])
			clade2 = tips(tree,tree$edge[daughters[2],2])
			if (length(clade1)<length(clade2)){
				pair_ortho_copies(clade1,clade2,'old_WGD.gene')
				write(paste(clade2,collapse = ' '),'new_WGD.gene',append=T)
			}
			else{
				pair_ortho_copies(clade2,clade1,'old_WGD.gene')
				write(paste(clade1,collapse = ' '),'new_WGD.gene',append=T)
			}
		}
	}
}