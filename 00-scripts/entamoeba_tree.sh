######################
# Create backbone tree
######################
# download databases
wget https://github.com/pr2database/pr2database/releases/download/v4.12.0/pr2_version_4.12.0_18S_taxo_long.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_tax_silva.fasta.gz
gzip -d *gz
# remove wordwrap on silva
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' SILVA_132_SSURef_tax_silva.fasta > temp
mv temp SILVA_132_SSURef_tax_silva.fasta
# ncbi nt query: <genus name> 18S rRNA NOT whole genome shotgun NOT mrna 
# pull all 18S Entamoeba entries from SILVA, PR2, and NCBI
grep "Entamoeba" pr2_version_4.12.0_18S_taxo_long.fasta -A 1 | sed 's/--//' > pr2_entamoeba.fa
grep "Entamoeba" SILVA_132_SSURef_tax_silva.fasta -A 1 | sed 's/--//' > silva_entamoeba.fa
# change u to t in silva
sed '/^[^>]/s/U/T/g' silva_entamoeba.fa > temp
mv temp silva_entamoeba.fa
# get ncbi, rename
mv ~/Downloads/sequence.fasta . && mv sequence.fasta ncbi_entamoeba.fa
# concatenate together
cat ncbi_entamoeba.fa silva_entamoeba.fa pr2_entamoeba.fa > all_entamoeba.fa
# manually remove two odd sequences: AAFB02001415.298.1771, X75434.1
# how many entries pre-cluster?
grep ">" all_entamoeba.fa -c
# 682
# sort by length
vsearch --sortbylength all_entamoeba.fa --output all_entamoeba.sort
# cluster at 99% using vsearch
vsearch --cluster_fast all_entamoeba.sort --centroids all_entamoeba.clust --id 0.99
# how many sequences after cluster
grep ">" all_entamoeba.clust -c
# 147
# clean headers
sed 's/|.*//' all_entamoeba.clust > temp
mv temp all_entamoeba.clust
# align with mafft
mafft --auto all_entamoeba.clust > all_entamoeba.align.fa
# clean up for trimal 
sed '/^[^>]/s/n/-/g' all_entamoeba.align.fa > temp
mv temp all_entamoeba.align.fa
# trim with trimal (as silva has degenerate bases must provide alternative matrix)
trimal -in all_entamoeba.align.fa -out all_entamoeba.trim.fa  -gt 0.3 -st 0.001 -matrix matrix.Degenerated_DNA 
# generate tree with raxml -- GTR and 100 bootstrap
rm *ref.tre
raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n ref.tre -s all_entamoeba.trim.fa
# generate taxonomy file
mkdir reference_trees
grep ">" all_entamoeba.fa | sed 's/>//' | sed 's/|/\t/' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/|/_/g' > annotations.ent.txt
mv annotations.ent.txt reference_trees/
# add header line in nano
cp RAxML_bipartitions.ref.tre entamoeba.tre 
mv entamoeba.tre reference_trees

#########################
# Generate placement tree
#########################
# BLAST unresolved reads
mkdir placement_tree
cd placement_tree
blastn -query Dal_Afribiota_18S_V4_unresolved.fasta -out Dal_Afribiota_18S_V4_unresolved.blast.out -db ~/Desktop/referenceDB/ncbi_12.5.19/nt -perc_identity 0.99 -evalue 1e-10 -max_target_seqs 500 -outfmt 6
# pull nematodes and entamoeba from megan taxonomy results
grep "Entamoeba" Dal_Afribiota_18S_V4_unresolved.taxonomy.txt | awk '{print $1}' > ent.ids
grep "Nematoda" Dal_Afribiota_18S_V4_unresolved.taxonomy.txt | awk '{print $1}' > nem.ids
# pull reads
cat ent.ids | while read line; do grep -w $line Dal_Afribiota_18S_V4_unresolved.fasta -A 1; done > ent.fa
cat nem.ids | while read line; do grep -w $line Dal_Afribiota_18S_V4_unresolved.fasta -A 1; done > nem.fa
# modify asv names and merge Dal & MI datasets
sed -i 's/ASV/Dal_ASV/' Dal_Afribiota_18S_V4_*e.fasta
sed -i 's/ASV/MI_ASV/' MI_Afribiota_18S_V4_*fasta
sed -i 's/ASV/MI_ASV/' MI_Afribiota_18S_V4_ASV_table_20_no_plants_host.txt
sed -i 's/ASV/Dal_ASV/' Dal_Afribiota_18S_V4_sequence_table.18s_R1_lwp.txt
sed -i 's/ASV/Dal_ASV/' ent.fa
sed -i 's/ASV/Dal_ASV/' nem.fa
cat MI_Afribiota_18S_V4_entamoba_ASV.fasta Dal_Afribiota_18S_V4_entamoeba_ASV.fasta ent.fa > all_ent.fa
cat MI_Afribiota_18S_V4_nematode_ASV.fasta Dal_Afribiota_18S_V4_nematode.fasta nem.fa > all_nem.fa
# blast all entamoeba reads, test for accuracy
blastn -query all_ent.fa -out all_ent.blast.out -db ~/Desktop/referenceDB/ncbi_12.5.19/nt -max_target_seqs 500 -outfmt 6

################
# EPA tree build
################
# align your (unaligned) sequences to your curated reference alignment (not the trimmed one)
sina -i ../all_entamoeba.align.fa --prealigned -o all_entamoeba.arb
sina -i all_ent.fa -r all_entamoeba.arb -o query.align.ent.fa --fs-msc 0.01 --fs-full-len=100
# concatenate query and ref
cat query.align.ent.fa ../all_entamoeba.align.fa > queryPlus.align.ent.fa
#build epa tree
rm *ent.epa.tre 
raxmlHPC-PTHREADS-SSE3 -f v -G 0.2 -m GTRCAT -n ent.epa.tre -s queryPlus.align.ent.fa -t ../reference_trees/entamoeba.tre -T 2
#clean up tree so you can read into figtree
sed 's/QUERY___//g' RAxML_labelledTree.ent.epa.tre | sed 's/\[I[0-9]*\]//g' > RAxML_placementTree.ent.epa.tre

#######################
#Constraint tree build
#######################
raxmlHPC-PTHREADS-SSE3 -f a -N 100 -G 0.2 -m GTRCAT -n ent.cons.tre -s queryPlus.align.ent.fa -g ../reference_trees/entamoeba.tre -T 4 -x 25734 -p 25793

#now pass to phyloseq_tree.r
