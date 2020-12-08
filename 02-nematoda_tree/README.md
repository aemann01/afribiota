### First generate the backbone tree

Get IDs for each group of interest

```bash
grep -h "Spirurida" SILVA_132_SSURef_tax_silva.fasta pr2_version_4.12.0_18S_taxo_long.fasta | sed 's/>//' | sed 's/ /\t/' | sed 's/|/\t/' > temp
grep -h "Rhabditida" SILVA_132_SSURef_tax_silva.fasta pr2_version_4.12.0_18S_taxo_long.fasta | sed 's/>//' | sed 's/ /\t/' | sed 's/|/\t/' > temp1
grep -h "Trichuris" SILVA_132_SSURef_tax_silva.fasta pr2_version_4.12.0_18S_taxo_long.fasta | sed 's/>//' | sed 's/ /\t/' | sed 's/|/\t/' > temp2
grep -h "Oxyurida" SILVA_132_SSURef_tax_silva.fasta pr2_version_4.12.0_18S_taxo_long.fasta | sed 's/>//' | sed 's/ /\t/' | sed 's/|/\t/' > temp3
grep -h "Ascari" SILVA_132_SSURef_tax_silva.fasta pr2_version_4.12.0_18S_taxo_long.fasta | sed 's/>//' | sed 's/ /\t/' | sed 's/|/\t/' > temp4
```

Get accession numbers and accompanying taxonomy file for annotations

```bash
cat temp* | awk '{print $1}' > nem.ids
cat temp* > reference_trees/annotations.nem.txt
rm temp*
```

Filter nematode sequences from the fasta files

```bash
seqtk subseq SILVA_132_SSURef_tax_silva.fasta nem.ids > all_nematoda.fa
```

Clean up the SILVA file (U to T)

```bash
sed '/^[^>]/s/U/T/g' all_nematoda.fa > temp
mv temp all_nematoda.fa 
```

How many entries pre-cluster?

```bash
grep ">" all_nematoda.fa -c
1176
```

Sort by length and then cluster at 99%

```bash
vsearch --sortbylength all_nematoda.fa --output all_nematoda.sort
vsearch --cluster_fast all_nematoda.sort --centroids all_nematoda.clust --id 0.99
```

How many sequences after cluster?

```bash
grep ">" all_nematoda.clust -c
471
```

Clean up headers

```bash
sed 's/|.*//' all_nematoda.clust > temp
mv temp all_nematoda.clust
```

Align with mafft

```bash
mafft --auto all_nematoda.clust > all_nematoda.align.fa
```

Clean up for trimal

```bash
sed '/^[^>]/s/n/-/g' all_nematoda.align.fa > temp
mv temp all_nematoda.align.fa
```

Trim with trimal

```bash
trimal -in all_nematoda.align.fa -out all_nematoda.trim.fa  -gt 0.3 -st 0.001 -matrix matrix.Degenerated_DNA 
```

Generate tree with RAxML

```bash
rm *ref.tre
raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n ref.tre -s all_nematoda.trim.fa
cp RAxML_bipartitions.ref.tre nematoda.tre 
mv nematoda.tre reference_trees
```

### EPA tree build

Align your (unaligned) sequences to your curated reference aligment (not the trimmed one)

```bash
cd placement_tree
sina -i ../all_nematoda.align.fa --prealigned -o all_nematoda.arb
sina -i all_nem.fa -r all_nematoda.arb -o query.align.nem.fa --fs-msc 0.01 --fs-full-len=100
```

Concatenate query and reference sequences

```bash
cat query.align.nem.fa ../all_nematoda.align.fa > queryPlus.align.nem.fa
```

Build the EPA tree

```bash
rm *nem.epa.tre 
raxmlHPC-PTHREADS-SSE3 -f v -G 0.2 -m GTRCAT -n nem.epa.tre -s queryPlus.align.nem.fa -t ../reference_trees/nematoda.tre -T 2
```

Clean up tree so you can read into figtree

```bash
sed 's/QUERY___//g' RAxML_labelledTree.nem.epa.tre | sed 's/\[I[0-9]*\]//g' > RAxML_placementTree.nem.epa.tre
```

### Constraint tree build

```bash
raxmlHPC-PTHREADS-SSE3 -f a -N 100 -G 0.2 -m GTRCAT -n nem.cons.tre -s queryPlus.align.nem.fa -g ../reference_trees/nematoda.tre -T 4 -x 25734 -p 25793
```

Root in figTree and fix labels

```bash
sed -i 's/_//g' RAxML_bestTree.nem.cons.root.tre
```

![nematode tree](nematoda_tree.png)

Now you can pass the tree to phyloseq_tree.r.ipynb
