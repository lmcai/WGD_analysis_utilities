library(ape)
library(base)
library(geiger)

treefile_names <- list.files(pattern = ".*inclade*")
sp=c('Bergia','Bhesa','Bischofia','Calophyllum','Casearia','Chrysobalanus','Clusia','Clutia','Drypetes','Elatine','Endospermum','Erythroxylum','Flacourtia','Galearia','Galphimia','Garcinia','Hevea','Hypericum','Ixonanthes','Jatropha','Linum','Malesherbia','Mammea','Manihot','Ochna','Passiflora','Podostemum','Populus','Rhizophora','Ricinus','Rinorea','Salix_','SalixSu','Sauropus','Tristellateia','Viola')

for (i in treefile_names){
	id=strsplit(i, '[.]inclade')[[1]][1]
	spno=0
	tree=read.tree(i)
	for (j in sp){
		if (any(startsWith(tree$tip.label,j))){spno=spno+1}
	}
	write(paste(id,spno),'sp.no',append=T)
}