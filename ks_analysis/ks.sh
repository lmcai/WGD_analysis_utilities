########################
##Required input:
##transcriptome cds, pep sequences (header of individual fasta mush be the same) 
########################

#load environment
source new-modules.sh
module load cd-hit/4.6.4-fasrc02
module load blast/2.2.29+-fasrc01
module load Anaconda/1.9.2-fasrc01
module load trinity/2.2.0-fasrc01

########################
#reduce sequence redundancy of the transcriptome
#get longest isoform for each 'gene' and change fasta entries into numbers
/n/sw/fasrcsw/apps/Core/trinityrnaseq/2.2.0-fasrc01/util/misc/get_longest_isoform_seq_per_trinity_gene.pl trinity.nonFR.Trinity.fasta | awk '/^>/{print ">" ++i; next}{print}' > trinity_renamed.fas
#translate using transcoder
sbatch transdecoder.sh trinity_renamed.fas
#wait until it finishes

########################
#further reduce sequence redundancy of the transcriptome by removing identical sequences
cd-hit-est -T 4 -i trinity_renamed.fas.transdecoder.cds -o trinity_renamed.transdecoder.cdhit.fas -c 1 -n 10 -r 1

########################
#Internal all-by-all blast
makeblastdb -in trinity_renamed.transdecoder.cdhit.fas -dbtype nucl -out test
blastn -task dc-megablast -db test -num_threads 4 -query trinity_renamed.transdecoder.cdhit.fas -outfmt 6 -perc_identity 40 -out test.abablast.out

#get only the best non-self hit for each cds
export LC_ALL=C LC_LANG=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr test.abablast.out >test.abablast.out.srt
python filter_sorted_blast_for_top_number_hits.py test.abablast.out.srt 1
#output test.abablast.out.srt.top

#make a new subdirectory and write pairwise alignments using paired_extract.py
mkdir alignments
cd alignments
#python [path to paired_extract.py] [blast result] [pep sequence] [cds sequence]
python ../paired_extract.py ../test.abablast.out.srt.top trinity_renamed.transdecoder.pep trinity_renamed.transdecoder.cdhit.fas

#Time consuming step, alignment
#align aa using mafft, 3000 pairs in each job
sbatch ../mafft_***.sh

#wait until mafft finishes and then align na based on aa alignment using pal2nal
sbatch ../pal2nal.sh

#paml to calculate ks for all genes
cat *.codon >>paml.codon
cp ../paml.clt .
sbatch ../paml.sh
