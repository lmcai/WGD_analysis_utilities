##run under R3.3+

library(phytools)
library(geiger)
library(base)

treefiles <- list.files(pattern = "^RAx")

#define function to prune very short terminal taxa
prune_short_terminal_taxa <- function(x){
n_taxa=length(x$tip.label)

#get corresponding edge length for each tip
ee=setNames(x$edge.length[sapply(1:n_taxa,function(a,b) which(b==a),b=x$edge[,2])],x$tip.label)
#get tip label list with very small terminal branch length
short_taxa_number=which(log2(ee)<(-15))

sis=0
for (i in short_taxa_number){
	#when two sister taxa have very short branch length, drop one of them
	if (i == sis) next
	sis_nt = getSisters(x,names(ee[i]),mode='label')
	
	if ("tips" %in% attributes(sis_nt)$names){	#if sister is a tip...
		sis = which(names(ee)==sis_nt$tips) #sister is a tip
		if (sis %in% short_taxa_number){x=drop.tip(x,names(ee[i]))}
		#When the sister tip of same species have regular branch length, drop this tip
		else if (!is.na(strsplit(sis_nt$tips,'_')[[1]][1] == strsplit(names(ee[i]),'_')[[1]][1])){x=drop.tip(x,names(ee[i]))}
		#When the sister tip is different species, leave it as it is
	}
	if ("nodes" %in% attributes(sis_nt)$names){x=drop.tip(x,names(ee[i]))}
}
return(x)
}

#define function to prune long branch and its subclade
prune_long_branch <- function(x){
	#calculate trimmed mean branch length
	mean_br = mean(x$edge.length,trim=0.05)
	
	#first round set to >15 times, second time set to >10 times
	long_brs = which(x$edge.length>10*mean_br) #get list of long branches
	#prune daughter tips of long branches
	Vit = x$tip.label[which(startsWith(x$tip.label,'Vitis'))]
	tips=c()
	
	for (i in long_brs){
		daughters = tips(x,node=x$edge[i,2])
		#Vitis cannot be removed
		#if ((length(intersect(daughters,Vit)) == 0) & !(length(daughters)==1 & startsWith(daughters[1],'Podo')) & (length(daughters)<=40)){tips=c(daughters,tips)}
		if (!(length(daughters)==1 & startsWith(daughters[1],'Podo')) & (length(daughters)<=40)){tips=c(daughters,tips)}

	}
	unique(tips)
	if (length(tips)>0){x = drop.tip(x, tips)}
	return(x)
}

#define function to remove monophyletic clades that has longer root to tip length than podostemum
prune_longer_than_podo <- function(x){
	rt2tip=diag(vcv.phylo(x))
	podo_lg=0
	podo_lg=rt2tip[startsWith(names(rt2tip),'Podo')]
	if (length(podo_lg)>0){
		tips=names(rt2tip[which(rt2tip>=1.2*max(podo_lg))])
		if (length(tips)>0){x = drop.tip(x, tips)}
	}
	return(x)
}


for (i in 1:length(treefiles)){
	tree=read.tree(treefiles[i])
	
	tree=prune_short_terminal_taxa(tree)
	tree=prune_short_terminal_taxa(tree)
	#two rounds of iteration, incase of situation(A_short1,(A_long1,A_short2)). but is it too much?

	tree=prune_longer_than_podo(tree)
	#prune tips with root-to-tip length longer than podostemum

	tree = prune_long_branch(tree)
	tree = prune_long_branch(tree)
	#two rounds of long branch pruning
	
		
	newtree=strsplit(treefiles[i], '[.]')[[1]][2]
	write.tree(tree, paste(newtree,'.tre',sep=''))
}