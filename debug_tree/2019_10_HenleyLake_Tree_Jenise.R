#Scapula Beta Diversity UNIFRAC ----
#read in OTU table 
OTU2 = read.table("HLS_subsample_gg.0.05.an.shared", header=TRUE, sep="\t")
#need to use the "Group" column as the row names so that it will match our metadata
row.names(OTU2) = OTU2$Group
#we then need to remove the "label", "numOTUs", and Grou" colums as they are not OTU counts like the rest of the table 
OTU2.clean = OTU2[,-which(names(OTU2) %in% c("label", "numOtus", "Group"))]
View(OTU2.clean)

##Taxonomy of each OTU
#read in taxonomy table 
tax2 = read.table("HLS_subsample_gg.an.0.05.cons.taxonomy", header=TRUE, sep="\t")
View(tax2)
#for the taxonomy tale, we name the rows by the OTU #
row.names(tax2) = tax2$OTU
View(tax2)
#remove all OTUs that don't occur in our OTU.clean data set
tax.clean2 = tax2[row.names(tax2) %in% colnames(OTU2.clean),]
View(tax.clean2)
#we then need to separate the "taxonomy" column so that each level is in its own colun (use seprate from the tidyr package)
require("tidyr")
tax.clean2 = separate(tax.clean2, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species","Strain"), sep=";")
View(tax.clean2)
#remove the "size" and "strain" columns as well as "OTU" since these are now the row names 
tax3.clean = tax.clean2[,-which(names(tax.clean2) %in% c("Size", "Strain", "OTU"))]
View(tax3.clean)

#meta data file
meta = read.table("HenleyLake_Scapula_Meta_Data.txt", header=TRUE, sep = "\t")
#make row names the sample name
row.names(meta) = meta$Sample

#remove zero from OTU and Meta files
metareduced.nozero = meta[-c(1,20),] 
View(metareduced.nozero)
View(OTU2.clean)
OTU2.clean.nozero = OTU2.clean[-c(1,20),]
row.names(OTU2.clean.nozero)
View(OTU2.clean.nozero)

#to get orders correct and colors straight
match.order.nozero <- match(rownames(metareduced.nozero), rownames(OTU2.clean.nozero))
OTU2.clean.nozero <- OTU2.clean.nozero[match.order.nozero,]

#Begin to construct the phyloseq object
require(phyloseq)
OTU2.UF.nozero = otu_table(as.matrix(OTU2.clean.nozero), taxa_are_rows=FALSE)
tax2.UF.nozero = tax_table(as.matrix(tax3.clean))
meta.UF.nozero = sample_data(metareduced.nozero)

physeq.nozero = phyloseq(OTU2.UF.nozero, tax2.UF.nozero, meta.UF.nozero)

#Import Tree File 
treefile <- "HLS_subsample_gg.an.unique_list.0.05.rep.otu.phylip.tre"
tree <- import_mothur(mothur_tree_file = treefile)

#merge mothur imported tree file and phyloseq object
physeq.tree.nozero = merge_phyloseq(physeq.nozero,tree)

#Scapula Phylogenetic UniFrac Weighted No Zero Distance Calculation----
wUF.nozero.dist.HLS = UniFrac(physeq.tree.nozero, weighted=TRUE, normalized = TRUE)

#Scapula Phylogenetic UniFrac Unweighted No Zero Distance Calculation ----
uwUF.nozero.dist.HLS = UniFrac(physeq.tree.nozero, weighted=FALSE, normalized = TRUE)

#hclust to create tree ----
#create trees to compare (i.e weighted versus unweighted)
#create weighted unifrac tree
wUF.hclust = hclust(wUF.nozero.dist.HLS, method = "average")
plot(wUF.hclust) #issue with hclust is that it treats each sample as a leaf on the tree

#create unweighted unifrac tree
uwUF.hclust = hclust(uwUF.nozero.dist.HLS, method = "average")

#make weighted nad unweighted unifrac hclust items into dendrograms for use in the tanglegram 
wUF.dend = as.dendrogram(wUF.hclust)
uwUF.dend = as.dendrogram(uwUF.hclust)

#tanglegram in dendextend allows you to draw comparisons between two trees 
#two trees in this case are UniFrac weighted and unweighted 
require(dendextend)
tanglegram(wUF.dend, uwUF.dend)


#META.clust in RAM package----
#META.clust allows you to use a meta data file to provide grouping parameters - so we can group individual samples into ADD 
# install.packages("RAM")
require(RAM)
META.clust(meta = metareduced, type = "rectangle", group = "Collection", dist=wUF.dist, clust = "average")

## Try using META.clust grouping to create hclust object.
wUF.hclust = group.with.META.clust(wUF.nozero.dist.HLS, method = "average")
