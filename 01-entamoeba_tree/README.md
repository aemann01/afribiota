### First create the backbone tree

Download and decompress necessary databases

```bash
wget https://github.com/pr2database/pr2database/releases/download/v4.12.0/pr2_version_4.12.0_18S_taxo_long.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_tax_silva.fasta.gz
gzip -d *gz
```

Remove wordwrap on the SILVA reference file

```bash
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' SILVA_132_SSURef_tax_silva.fasta > temp
mv temp SILVA_132_SSURef_tax_silva.fasta
```

Pull all Entamoeba 18S rRNA entries from NCBI using the following query: Entamoeba 18S rRNA NOT whole genome shotgun NOT mrna, export to file and save in your working directory (here we named it ncbi_entamoeba.fa)

Pull all Entamoeba entries from SILVA and PR2

```bash
grep "Entamoeba" pr2_version_4.12.0_18S_taxo_long.fasta -A 1 | sed 's/--//' > pr2_entamoeba.fa
grep "Entamoeba" SILVA_132_SSURef_tax_silva.fasta -A 1 | sed 's/--//' > silva_entamoeba.fa
```

Change U to T in SILVA sequences

```bash
sed '/^[^>]/s/U/T/g' silva_entamoeba.fa > temp
mv temp silva_entamoeba.fa
```

Concatenate all of your reference sequences together

```bash
cat ncbi_entamoeba.fa silva_entamoeba.fa pr2_entamoeba.fa > all_entamoeba.fa
```

Before moving on we need to remove two odd sequences from our reference datasets (confirmed in tree) -- you can manually delete these in nano (AAFB02001415.298.1771, X75434.1)

How many entries do we have before clustering?

```bash
grep ">" all_entamoeba.fa -c
682
```

Now sort by length and cluster at 99% identity

```bash
vsearch --sortbylength all_entamoeba.fa --output all_entamoeba.sort
vsearch --cluster_fast all_entamoeba.sort --centroids all_entamoeba.clust --id 0.99
```

How many sequences post clustering?

```bash
grep ">" all_entamoeba.clust -c
147
```

Clean up headers then align using mafft

```bash
sed 's/|.*//' all_entamoeba.clust > temp
mv temp all_entamoeba.clust
mafft --auto all_entamoeba.clust > all_entamoeba.align.fa
```

Clean up clustered file for trimal

```bash
sed '/^[^>]/s/n/-/g' all_entamoeba.align.fa > temp
mv temp all_entamoeba.align.fa
```

Trim with trimal (need to use alternative matrix to account for degenerate bases)

```bash
trimal -in all_entamoeba.align.fa -out all_entamoeba.trim.fa  -gt 0.3 -st 0.001 -matrix matrix.Degenerated_DNA 
```

Generate backbone tree with RAxML

```bash
rm *ref.tre
raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n ref.tre -s all_entamoeba.trim.fa
```

Generate taxonomy file for annotations

```bash
mkdir reference_trees
grep ">" all_entamoeba.fa | sed 's/>//' | sed 's/|/\t/' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/|/_/g' > annotations.ent.txt
mv annotations.ent.txt reference_trees/
```

Add header line in nano, then move reference tree to proper directory

```bash
cp RAxML_bipartitions.ref.tre entamoeba.tre 
mv entamoeba.tre reference_trees
```

### Get unassigned Entamoeba or nematodes

First we want to see if there are any reads that were unassigned that are actually Entamoeba or nematodes. You need to have a local copy of the ncbi nt database -- change your path in the command below to its location

```bash
mkdir placement_tree
cd placement_tree
blastn -query Dal_Afribiota_18S_V4_unresolved.fasta -out Dal_Afribiota_18S_V4_unresolved.blast.out -db ~/Desktop/referenceDB/ncbi_12.5.19/nt -perc_identity 0.99 -evalue 1e-10 -max_target_seqs 500 -outfmt 6
```

Get taxonomic assignments for each sequence using a LCA algorithm (here I used [MEGAN6 community edition](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004957))

Pull nematodes and Entamoeba assigned reads from your taxonomy results

```bash
grep "Entamoeba" Dal_Afribiota_18S_V4_unresolved.taxonomy.txt | awk '{print $1}' > ent.ids
grep "Nematoda" Dal_Afribiota_18S_V4_unresolved.taxonomy.txt | awk '{print $1}' > nem.ids
```

Now pull the reads by their ids

```bash
cat ent.ids | while read line; do grep -w $line Dal_Afribiota_18S_V4_unresolved.fasta -A 1; done > ent.fa
cat nem.ids | while read line; do grep -w $line Dal_Afribiota_18S_V4_unresolved.fasta -A 1; done > nem.fa
```

Now we need to modify ASV names and merge the Dal & MI datasets

```bash
sed -i 's/ASV/Dal_ASV/' Dal_Afribiota_18S_V4_*e.fasta
sed -i 's/ASV/MI_ASV/' MI_Afribiota_18S_V4_*fasta
sed -i 's/ASV/MI_ASV/' MI_Afribiota_18S_V4_ASV_table_20_no_plants_host.txt
sed -i 's/ASV/Dal_ASV/' Dal_Afribiota_18S_V4_sequence_table.18s_R1_lwp.txt
sed -i 's/ASV/Dal_ASV/' ent.fa
sed -i 's/ASV/Dal_ASV/' nem.fa
cat MI_Afribiota_18S_V4_entamoba_ASV.fasta Dal_Afribiota_18S_V4_entamoeba_ASV.fasta ent.fa > all_ent.fa
cat MI_Afribiota_18S_V4_nematode_ASV.fasta Dal_Afribiota_18S_V4_nematode.fasta nem.fa > all_nem.fa
```

Blast all Entamoeba reads to test for accuracy

```bash
blastn -query all_ent.fa -out all_ent.blast.out -db ~/Desktop/referenceDB/ncbi_12.5.19/nt -max_target_seqs 500 -outfmt 6
```

### EPA tree build

Align your (unaligned) sequences to your curated reference aligment (not the trimmed one)

```bash
cd placement_tree
sina -i ../all_entamoeba.align.fa --prealigned -o all_entamoeba.arb
sina -i all_ent.fa -r all_entamoeba.arb -o query.align.ent.fa --fs-msc 0.01 --fs-full-len=100
```

Concatenate your query sequences and reference sequences

```bash
cat query.align.ent.fa ../all_entamoeba.align.fa > queryPlus.align.ent.fa
```

Now can build your EPA tree

```bash
rm *ent.epa.tre 
raxmlHPC-PTHREADS-SSE3 -f v -G 0.2 -m GTRCAT -n ent.epa.tre -s queryPlus.align.ent.fa -t ../reference_trees/entamoeba.tre -T 2
```

Clean up tree so you can read into figtree

```bash
sed 's/QUERY___//g' RAxML_labelledTree.ent.epa.tre | sed 's/\[I[0-9]*\]//g' > RAxML_placementTree.ent.epa.tre
```

### Constraint tree build

As EPA tree does not have bootstrap information, also generate a constraint tree

```bash
raxmlHPC-PTHREADS-SSE3 -f a -N 100 -G 0.2 -m GTRCAT -n ent.cons.tre -s queryPlus.align.ent.fa -g ../reference_trees/entamoeba.tre -T 4 -x 25734 -p 25793
```

Root tree in figtree before next step

Now you can use the output of your placement tree in the script phyloseq_tree.r
