#extract paired na and aa sequences by entry id from all-by-all blast result
#python paired_extract.py [blast output] [aa seqfile] [na seqfile]

from Bio import SeqIO
import sys

blast_res=open(sys.argv[1],'r')
lines= blast_res.readlines()
aa_dict = SeqIO.index(sys.argv[2], "fasta")
na_dict = SeqIO.index(sys.argv[3], "fasta")
for i in range(0,len(lines)):
	id=lines[i].split()
	aa_filename=`i`+'.aa.fas'
	SeqIO.write(aa_dict[id[0]], aa_filename, "fasta")
	SeqIO.write(aa_dict[id[1]], open(aa_filename,'a'), "fasta")
	na_filename=`i`+'.na.fas'
	SeqIO.write(na_dict[id[0]], na_filename, "fasta")
	SeqIO.write(na_dict[id[1]], open(na_filename,'a'), "fasta")
