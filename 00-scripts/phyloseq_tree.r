library("phyloseq");packageVersion("phyloseq")
library(ggplot2)
library(ape)


otuDal <- read.table("Dal_Afribiota_18S_V4_sequence_table.18s_R1_lwp.txt", sep="\t", header=T)
otuMI <- read.table("MI_Afribiota_18S_V4_ASV_table_20_no_plants_host.txt", sep="\t", header=T)
otu <- dplyr::full_join(otuDal, otuMI)
row.names(otu) <- otu$ASV_IDs
otu <- data.matrix(otu)
otu[is.na(otu)] <- as.numeric(0)

#read in mapping file, set as sample data type
map <- as.data.frame(read.table("afribiota_tot_finalS13_26092018_recoded.txt", sep="\t", header=T, row.names=1))
map <- sample_data(map)
map$stunted <- as.factor(map$stunted)

#entamoeba tree 
tre <- read.tree("RAxML_bestTree.ent.cons.root.tre")

otutable <- otu_table(otu, taxa_are_rows=T)
physeq <- phyloseq(otutable)
physeq <- merge_phyloseq(physeq, tre, map)

pdf("entamoeba_phyloseq_tree.pdf")
plot_tree(physeq, color="stunted", shape="pays", ladderize="left", label.tips="taxa_names", base.spacing=0.03)
dev.off()

#nematode tree 
tre <- read.tree("RAxML_bestTree.nem.cons.root.tre")

otutable <- otu_table(otu, taxa_are_rows=T)
physeq <- phyloseq(otutable)
physeq <- merge_phyloseq(physeq, tre, map)

pdf("nematode_phyloseq_tree.pdf")
plot_tree(physeq, color="stunted", shape="pays", ladderize="left", label.tips="taxa_names", base.spacing=0.03)
dev.off()









