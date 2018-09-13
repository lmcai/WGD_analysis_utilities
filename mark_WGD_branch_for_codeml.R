#read in phylogeny, generate control file for r8s, aimed at dating multiple WGDs

#load library
library(phytools)
library(base)
library(geiger)

#read in files
treefile_names <- list.files(pattern = ".*inclade*")

#define function
write_tree_section<-function(tree_str,filename){
	x=paste('	tree tree_1 = [&R] ',tree_str,sep='')
	writeLines(c("#NEXUS","begin trees;",x,"end;","\n"),filename)
}

get_align_length<-function(num){
	x=read.dna(file=paste(num,'.codon',sep=''),format='fasta')
	aln_length=length(x[1,])
	return(aln_length)
}

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

define_age_constrain <-function(phy,clade){
	if (clade=='ROOT'){
		rt_edge=which(phy$edge[,1]==(length(phy$tip.label)+1))
		taxa_1=tips(phy,node=phy$edge[rt_edge[1],2])
		taxa_2=tips(phy,node=phy$edge[rt_edge[2],2])
		rt_tips=paste(sample(taxa_1,1),sample(taxa_2,1))
		rt_node=paste('MRCA ROOT ',rt_tips,';',sep='')
		write(rt_node,ctr_filename,append=T)
	}
	
	if (clade=='euphorbioids'){
		taxa="Mani|Hev|Jatr|Endos|Ricinus|Clutia"
		member=phy$tip.label[grep(taxa,phy$tip.label)]
		member=return_monophyletic_gr(phy,member) #delete up to 3 outliers
		if (is.monophyletic(phy,member)){
			euph_node=paste('MRCA euphorbioids ',paste(member,collapse =' '),';',sep='')
			write(euph_node,ctr_filename,append=T)
		}
	}
	if (clade=='clusioids'){#stem group clusiaceae
		taxa="Clusia|Garc"
		member=phy$tip.label[grep(taxa,phy$tip.label)]
		member=return_monophyletic_gr(phy,member) #delete up to 3 outliers
		if (length(member)==1){#only one taxa
			sis=getSisters(phy,member,mode='label')
			if ('tips' %in% names(sis)){new_member=c(member,sis)}
			if ('nodes' %in% names(sis)){new_member=c(member,tips(phy,sis))}
		} else if (length(member)>1){
			node=getMRCA(phy,member)
			sis=getSisters(phy,node,mode='label')
			if ('tips' %in% names(sis)){new_member=c(member,sis)}
			if ('nodes' %in% names(sis)){new_member=c(member,tips(phy,sis))}
		} else if(length(member)==0){ #if clusia and garcinia is not in inclade, use crown group clusioids as calibration point
			taxa="Calo|Mamm|Hyper|Podo|Garc|Clusia"
			member=phy$tip.label[grep(taxa,phy$tip.label)]
			new_member=return_monophyletic_gr(phy,member) #delete up to 3 outliers
		}
		if (is.monophyletic(phy,unlist(new_member))){
			clus_node=paste('MRCA clusioids ',paste(new_member,collapse =' '),';',sep='')
			write(clus_node,ctr_filename,append=T)
		}
	}
	if (clade=='salicoids'){
		taxa="Salix|Popu|Flac|Case"
		member=phy$tip.label[grep(taxa,phy$tip.label)]
		member=return_monophyletic_gr(phy,member) #delete up to 3 outliers
		if (is.monophyletic(phy,member)){
			salx_node=paste('MRCA salicoids ',paste(member,collapse =' '),';',sep='')
			write(salx_node,ctr_filename,append=T)
		}
	}
}

add_mark2tree <-function(phy,clade){
	WGD_node=getMRCA(phy,clade)
	daughters=which(phy$edge[,1]==WGD_node)
	clade1=tips(phy,phy$edge[daughters[1],2])
	clade2=tips(phy,phy$edge[daughters[2],2])
	if (length(clade1)==1){
		new_label=paste(phy$tip.label[phy$tip.label==clade1],'xxx')
		phy$tip.label[phy$tip.label==clade1]=new_label
	}
	if (length(clade1)>1){phy=makeNodeLabel(phy,'u',nodeList = list(xxx = clade1))}
	if (length(clade2)==1){
		new_label=paste(phy$tip.label[phy$tip.label==clade2],'xxx')
		phy$tip.label[phy$tip.label==clade2]=new_label
	}
	if (length(clade2)>1){phy=makeNodeLabel(phy,'u',nodeList = list(xxx = clade2))}
	return(phy)
}

define_WGD_nodes <-function(phy,clade){
	if (clade=='MH'){
		M=phy$tip.label[which(startsWith(phy$tip.label,'Mani'))]
		H=phy$tip.label[which(startsWith(phy$tip.label,'Heve'))]
		if (length(M)+ length(H)>=2){
			MH_tips=return_monophyletic_gr(phy,c(M,H))
			M=MH_tips[which(startsWith(MH_tips,'Mani'))]
			H=MH_tips[which(startsWith(MH_tips,'Heve'))]
			if ((is.monophyletic(phy,c(M,H))) & ((length(M)>=2) | (length(H)>=2))){	#((M,H),(M,H))
				phy=add_mark2tree(phy,c(M,H))
			}
		}
		return(phy)
	}
	if (clade=='CG'){
		C=phy$tip.label[which(startsWith(phy$tip.label,'Clusia'))]
		G=phy$tip.label[which(startsWith(phy$tip.label,'Garci'))]
		if (length(C)+ length(G)>=2){
			CG_tips=return_monophyletic_gr(phy,c(C,G))
			C=CG_tips[which(startsWith(CG_tips,'Clu'))]
			G=CG_tips[which(startsWith(CG_tips,'Garc'))]
			if ((is.monophyletic(phy,c(C,G))) & ((length(C)>=2) | (length(G)>=2))){	#((C,G),(C,G))
				phy=add_mark2tree(phy,c(C,G))
			}
		}
		return(phy)
	}
	if (clade=='PS'){
		member=phy$tip.label[grep("Salix|Popu|Flac",phy$tip.label)]
		member=return_monophyletic_gr(phy,member) #delete up to 3 outliers
		if (length(member)>=2){
		node=getMRCA(phy,member)
		daughters=which(phy$edge[,1]==node)
		clade1=tips(phy,phy$edge[daughters[1],2])
		clade2=tips(phy,phy$edge[daughters[2],2])
		if (((any(startsWith(clade1,'Fla'))) & (any(startsWith(clade2,'Fla')))) | ((any(startsWith(clade1,'SalixSu_'))) & (any(startsWith(clade2,'SalixSu_')))) | ((any(startsWith(clade1,'Salix_'))) & (any(startsWith(clade2,'Salix_')))) |((any(startsWith(clade1,'Pop'))) & (any(startsWith(clade2,'Pop'))))){
			#any taxa was found duplex in two descending clades
			phy=add_mark2tree(phy,c(sample(clade1,1),sample(clade2,1)))
		}
		else if (length(c(clade1,clade2))>2){
			if (length(clade1)>length(clade2)){node=phy$edge[daughters[1],2]}
			if (length(clade1)<=length(clade2)){node=phy$edge[daughters[2],2]}
			daughters=which(phy$edge[,1]==node)
			clade1=tips(phy,phy$edge[daughters[1],2])
			clade2=tips(phy,phy$edge[daughters[2],2])
			if (((any(startsWith(clade1,'Fla'))) & (any(startsWith(clade2,'Fla')))) | ((any(startsWith(clade1,'SalixSu_'))) & (any(startsWith(clade2,'SalixSu_')))) | ((any(startsWith(clade1,'Salix_'))) & (any(startsWith(clade2,'Salix_')))) |((any(startsWith(clade1,'Pop'))) & (any(startsWith(clade2,'Pop'))))){
				phy=add_mark2tree(phy,c(sample(clade1,1),sample(clade2,1)))
			}
		}
		}
		return(phy)
	}
	if (clade=='CM'){
		member=phy$tip.label[grep("Calo|Mamm",phy$tip.label)]
		member=return_monophyletic_gr(phy,member) #delete up to 3 outliers
		if (length(member)>=2){
			node=getMRCA(phy,member)
			daughters=which(phy$edge[,1]==node)
			clade1=tips(phy,phy$edge[daughters[1],2])
			clade2=tips(phy,phy$edge[daughters[2],2])
			if (((any(startsWith(clade1,'Calo'))) & (any(startsWith(clade2,'Calo')))) | ((any(startsWith(clade1,'Mamm'))) & (any(startsWith(clade2,'Mamm'))))){
				phy=add_mark2tree(phy,c(sample(clade1,1),sample(clade2,1)))
			}
			else{
				if (length(clade1)>length(clade2)){node=tree$edge[daughters[1],2]}
				else{node=tree$edge[daughters[2],2]}
				if (length(tips(phy,node))>=2){
					daughters=which(phy$edge[,1]==node)
					clade1=tips(phy,phy$edge[daughters[1],2])
					clade2=tips(phy,phy$edge[daughters[2],2])
					if (((any(startsWith(clade1,'Calo'))) & (any(startsWith(clade2,'Calo')))) | ((any(startsWith(clade1,'Mamm'))) & (any(startsWith(clade2,'Mamm'))))){
						phy=add_mark2tree(phy,c(sample(clade1,1),sample(clade2,1)))
					}
				}
			}
		}
		return(phy)
	}
	
	
	if (clade=='Dry'){
		taxa=phy$tip.label[which(startsWith(phy$tip.label,'Dry'))]
		taxa=return_monophyletic_gr(phy,taxa) #delete up to 3 outliers
		if ((length(taxa)>=2) & (is.monophyletic(phy,taxa))){
			phy=add_mark2tree(phy,taxa)
		}
		return(phy)
	}
	if (clade=='Vio'){
		taxa=phy$tip.label[which(startsWith(phy$tip.label,'Vio'))]
		taxa=return_monophyletic_gr(phy,taxa) #delete up to 3 outliers
		if ((length(taxa)>=2) & (is.monophyletic(phy,taxa))){
			phy=add_mark2tree(phy,taxa)
		}
		return(phy)
	}
	if (clade=='Pas'){
		taxa=phy$tip.label[which(startsWith(phy$tip.label,'Pas'))]
		taxa=return_monophyletic_gr(phy,taxa) #delete up to 3 outliers
		if ((length(taxa)>=2) & (is.monophyletic(phy,taxa))){
			phy=add_mark2tree(phy,taxa)
		}
		return(phy)
	}
	if (clade=='Ery'){
		taxa=phy$tip.label[which(startsWith(phy$tip.label,'Ery'))]
		taxa=return_monophyletic_gr(phy,taxa) #delete up to 3 outliers
		if ((length(taxa)>=2) & (is.monophyletic(phy,taxa))){
			phy=add_mark2tree(phy,taxa)
		}
		return(phy)
	}
	if (clade=='Rhi'){
		taxa=phy$tip.label[which(startsWith(phy$tip.label,'Rhi'))]
		taxa=return_monophyletic_gr(phy,taxa) #delete up to 3 outliers
		if ((length(taxa)>=2) & (is.monophyletic(phy,taxa))){
			phy=add_mark2tree(phy,taxa)
		}
		return(phy)
	}
	if (clade=='Sau'){
		taxa=phy$tip.label[which(startsWith(phy$tip.label,'Sau'))]
		taxa=return_monophyletic_gr(phy,taxa) #delete up to 3 outliers
		if ((length(taxa)>=2) & (is.monophyletic(phy,taxa))){
			phy=add_mark2tree(phy,taxa)
		}
		return(phy)
	}
	if (clade=='Bis'){
		taxa=phy$tip.label[which(startsWith(phy$tip.label,'Bis'))]
		taxa=return_monophyletic_gr(phy,taxa) #delete up to 3 outliers
		if ((length(taxa)>=2) & (is.monophyletic(phy,taxa))){
			phy=add_mark2tree(phy,taxa)
		}
		return(phy)
	}
	if (clade=='Pod'){
		taxa=phy$tip.label[which(startsWith(phy$tip.label,'Podo'))]
		taxa=return_monophyletic_gr(phy,taxa) #delete up to 3 outliers
		if ((length(taxa)>=2) & (is.monophyletic(phy,taxa))){
			phy=add_mark2tree(phy,taxa)
		}
		return(phy)
	}
	if (clade=='Hyp'){
		taxa=phy$tip.label[which(startsWith(phy$tip.label,'Hyp'))]
		taxa=return_monophyletic_gr(phy,taxa) #delete up to 3 outliers
		if ((length(taxa)>=2) & (is.monophyletic(phy,taxa))){
			phy=add_mark2tree(phy,taxa)
		}
		return(phy)
	}
	if (clade=='Lin'){
		taxa=phy$tip.label[which(startsWith(phy$tip.label,'Lin'))]
		taxa=return_monophyletic_gr(phy,taxa) #delete up to 3 outliers
		if ((length(taxa)>=2) & (is.monophyletic(phy,taxa))){
			phy=add_mark2tree(phy,taxa)
		}		
		return(phy)
	}
	if (clade=='Och'){
		taxa=phy$tip.label[which(startsWith(phy$tip.label,'Och'))]
		taxa=return_monophyletic_gr(phy,taxa) #delete up to 3 outliers
		if ((length(taxa)>=2) & (is.monophyletic(phy,taxa))){
			phy=add_mark2tree(phy,taxa)
		}
		return(phy)
	}
	if (clade=='Gal'){
		taxa=phy$tip.label[which(startsWith(phy$tip.label,'Gal'))]
		taxa=return_monophyletic_gr(phy,taxa) #delete up to 3 outliers
		if ((length(taxa)>=2) & (is.monophyletic(phy,taxa))){
			phy=add_mark2tree(phy,taxa)
		}
		return(phy)
	}
	if (clade=='Bhe'){
		taxa=phy$tip.label[which(startsWith(phy$tip.label,'Bhe'))]
		taxa=return_monophyletic_gr(phy,taxa) #delete up to 3 outliers
		if ((length(taxa)>=2) & (is.monophyletic(phy,taxa))){
			phy=add_mark2tree(phy,taxa)
		}
		return(phy)
	}
}

cv_text="FIXAGE taxon=ROOT age=109;
CONSTRAIN taxon=euphorbioids min_age=107;
CONSTRAIN taxon=clusioids min_age=89;
CONSTRAIN taxon=salicoids min_age=94;
set verbose=0;
set ftol=1e-7 maxiter=2000;
DIVTIME method=pl crossv=yes cvstart=-4.5 cvinc=0.5 cvnum=19 algorithm=tn;
END;"
conduct_cv <-function(filename){
write(cv_text,filename,append=T)
}



#loop through trees
for (i in treefile_names){ 
	#tree_text=readLines(i)
	tree=read.tree(i)
	id=strsplit(i,'.inclade')
	new_filename=paste(id[[1]][1],'.marked.tre',sep='')
	
	#write tree section
	#write_tree_section(tree_text[1],ctr_filename)
	#write r8s head
	#cur_line=paste('BLFORMAT lengths=persite nsites=',get_align_length(id[[1]][1]),' ultrametric=no;',sep='')
	#d=lapply(c("BEGIN r8s;",cur_line,"COLLAPSE;\n"),function(xx){write(xx,ctr_filename,append=T)})
	
	#define age constrains #sapply?
	#for (clade in c('ROOT','euphorbioids','clusioids','salicoids')){define_age_constrain(tree,clade)}
	
	#write mark 'xxx' to WGD (foreground) branches
	for (clade in c('MH','CM','CG','PS','Dry','Vio','Pas','Ery','Rhi','Sau','Bis','Pod','Hyp','Lin','Och','Gal','Bhe')){tree=define_WGD_nodes(tree,clade)}
	write.tree(tree,new_filename)
	#conduct cv	
	#conduct_cv(ctr_filename)
}

