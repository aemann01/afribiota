###############################
# Nematoda Reference tree build
###############################
# get information for each group
grep -h "Spirurida" SILVA_132_SSURef_tax_silva.fasta pr2_version_4.12.0_18S_taxo_long.fasta | sed 's/>//' | sed 's/ /\t/' | sed 's/|/\t/' > temp
grep -h "Rhabditida" SILVA_132_SSURef_tax_silva.fasta pr2_version_4.12.0_18S_taxo_long.fasta | sed 's/>//' | sed 's/ /\t/' | sed 's/|/\t/' > temp1
grep -h "Trichuris" SILVA_132_SSURef_tax_silva.fasta pr2_version_4.12.0_18S_taxo_long.fasta | sed 's/>//' | sed 's/ /\t/' | sed 's/|/\t/' > temp2
grep -h "Oxyurida" SILVA_132_SSURef_tax_silva.fasta pr2_version_4.12.0_18S_taxo_long.fasta | sed 's/>//' | sed 's/ /\t/' | sed 's/|/\t/' > temp3
grep -h "Ascari" SILVA_132_SSURef_tax_silva.fasta pr2_version_4.12.0_18S_taxo_long.fasta | sed 's/>//' | sed 's/ /\t/' | sed 's/|/\t/' > temp4
# get accession numbers
cat temp* | awk '{print $1}' > nem.ids
# get taxonomy file
cat temp* > reference_trees/annotations.nem.txt
rm temp*
# add header line in nano
# filter nematode sequences
seqtk subseq SILVA_132_SSURef_tax_silva.fasta nem.ids > all_nematoda.fa
# U to T
sed '/^[^>]/s/U/T/g' all_nematoda.fa > temp
mv temp all_nematoda.fa 
# how many entries pre-cluster?
grep ">" all_nematoda.fa -c
# 1176
# sort by length
vsearch --sortbylength all_nematoda.fa --output all_nematoda.sort
# cluster at 99% using vsearch
vsearch --cluster_fast all_nematoda.sort --centroids all_nematoda.clust --id 0.99
# how many sequences after cluster
grep ">" all_nematoda.clust -c
# 471
# clean headers
sed 's/|.*//' all_nematoda.clust > temp
mv temp all_nematoda.clust
# align with mafft
mafft --auto all_nematoda.clust > all_nematoda.align.fa
# clean up for trimal 
sed '/^[^>]/s/n/-/g' all_nematoda.align.fa > temp
mv temp all_nematoda.align.fa
# trim with trimal (as silva has degenerate bases must provide alternative matrix)
trimal -in all_nematoda.align.fa -out all_nematoda.trim.fa  -gt 0.3 -st 0.001 -matrix matrix.Degenerated_DNA 
# generate tree with raxml -- GTR and 100 bootstrap
rm *ref.tre
raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n ref.tre -s all_nematoda.trim.fa
cp RAxML_bipartitions.ref.tre nematoda.tre 
mv nematoda.tre reference_trees

################
# EPA tree build
################
cd placement_tree
# align your (unaligned) sequences to your curated reference alignment (not the trimmed one)
sina -i ../all_nematoda.align.fa --prealigned -o all_nematoda.arb
sina -i all_nem.fa -r all_nematoda.arb -o query.align.nem.fa --fs-msc 0.01 --fs-full-len=100
# concatenate query and ref
cat query.align.nem.fa ../all_nematoda.align.fa > queryPlus.align.nem.fa
#build epa tree
rm *nem.epa.tre 
raxmlHPC-PTHREADS-SSE3 -f v -G 0.2 -m GTRCAT -n nem.epa.tre -s queryPlus.align.nem.fa -t ../reference_trees/nematoda.tre -T 2
#clean up tree so you can read into figtree
sed 's/QUERY___//g' RAxML_labelledTree.nem.epa.tre | sed 's/\[I[0-9]*\]//g' > RAxML_placementTree.nem.epa.tre

#######################
#Constraint tree build
#######################
raxmlHPC-PTHREADS-SSE3 -f a -N 100 -G 0.2 -m GTRCAT -n nem.cons.tre -s queryPlus.align.nem.fa -g ../reference_trees/nematoda.tre -T 4 -x 25734 -p 25793

#now pass to phyloseq_tree.r