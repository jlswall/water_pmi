## Install packages if necessary.  Then load them.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
require(phyloseq)
browseVignettes("phyloseq")
install.packages("phangorn")


## Read first table.
OTU2 = read.table("HLS_subsample_gg.0.05.an.shared", header=TRUE, sep="\t")
## need to use the "Group" column as the row names so that it will
## match our metadata
row.names(OTU2) = OTU2$Group

## we then need to remove the "label", "numOTUs", and Grou" colums as
## they are not OTU counts like the rest of the table
OTU2.clean = OTU2[,-which(names(OTU2) %in% c("label", "numOtus", "Group"))]


## Read in file with taxonomy of each OTU
## The following is the original line from Claire, but the file she
## gaves me has a different name.
tax2 = read.table("HLS_subsample_gg.an.0.05.cons.taxonomy", header=TRUE, sep="\t")
## for the taxonomy tale, we name the rows by the OTU #
row.names(tax2) = tax2$OTU

#remove all OTUs that don't occur in our OTU.clean data set
tax.clean2 = tax2[row.names(tax2) %in% colnames(OTU2.clean),]

## we then need to separate the "taxonomy" column so that each level
## is in its own colun (use seprate from the tidyr package)

require("tidyr")
tax.clean2 = separate(tax.clean2, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species","Strain"), sep=";")

##remove the "size" and "strain" columns as well as "OTU" since these are now the row names 
tax3.clean = tax.clean2[,-which(names(tax.clean2) %in% c("Size", "Strain", "OTU"))]


## meta data
meta = read.table("HenleyLake_Scapula_Season_R_microbiota.txt", header=TRUE, sep="\t")
row.names(meta) = meta$Sample


## meta$ADDCollectionDays <- factor(meta$ADDCollectionDays, levels=c("0/Baseline/0","249/Collection1/23","489/Collection2/72","727/Collection3/114", "989/Collection4/150","1242/Collection5/178","1490/Collection6/204","1738/Collection7/227", "1980/Collection8/249","2224/Collection9/269","2496/Collection10/291","2784/Collection11/310", "3103/Collection12/333","3434/Collection13/353","3666/Collection14/373","3962/Collection15/422", "4229/Collection16/477","4527/Collection17/520","4863/Collection18/556","5200/Collection19/579"), ordered = TRUE)



## meta$Sample = factor(meta$Sample, ordered = TRUE)

## JLS: I'm going to try re-ordering the rows of OTU2.clean to match
## the order of the meta data.frame.
match.order <- match(rownames(meta), rownames(OTU2.clean))
OTU2.clean <- OTU2.clean[match.order,]

require(phyloseq)
OTU2.UF = otu_table(as.matrix(OTU2.clean), taxa_are_rows=FALSE)
tax2.UF = tax_table(as.matrix(tax3.clean))
meta.UF = sample_data(meta)
sample_names(OTU2.UF)
sample_names(tax2.UF)
sample_names(meta.UF)
physeq = phyloseq(OTU2.UF, tax2.UF, meta.UF)


tree <- import_mothur(mothur_tree_file = treefile)

physeq.tree = merge_phyloseq(physeq,tree)
View(physeq)
View(tree)
View(physeq.tree)
wUF.ordu = ordinate(physeq.tree, method = "NMDS", distance = "unifrac", weighted=TRUE)
View(wUF.ordu)
par(mfrow=c(1,1))
View(meta)
plot(wUF.ordu, type="n", main="Weighted UniFrac: Henley Lake Scapula")
#warning appears: species scores not avialable 
points(wUF.ordu, pch=20, display="sites", col=c("darkblue", "darkgoldenrod1","darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", 
                                                "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", 
                                                "darkblue","royalblue4", "dodgerblue3", "steelblue1", "lightskyblue")[meta$ADDCollectionDays])
View(meta)
par(mar=c(8, 4.1, 4.1, 10), xpd=TRUE) #add extra space

legend(0.35,0.2, legend=c("0/Baseline/0","249/Collection1/23","489/Collection2/72","727/Collection3/114","989/Collection4/150","1242/Collection5/178","1490/Collection6/204","1738/Collection7/227","1980/Collection8/249","2224/Collection9/269",
                          "2496/Collection10/291","2784/Collection11/310","3103/Collection12/333","3434/Collection13/353","3666/Collection14/373","3962/Collection15/422","4229/Collection16/477","4527/Collection17/520","4863/Collection18/556","5200/Collection19/579"), 
       col=c("darkblue", "darkgoldenrod1","darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue","royalblue4", "dodgerblue3", "steelblue1", "lightskyblue")[meta$Collection],
       pch=20,
       cex = 0.75)


#same as above but with different coloring option (rainbow)
par(mfrow=c(1,1))
plot(wUF.ordu, type="n", main="Weighted UniFrac: Henley Lake")
points(wUF.ordu, pch=20, display="sites", 
       col=rainbow(20)[meta$ADDCollectionDays])
par(mar=c(6, 4.1, 4.1, 10), xpd=TRUE) #add extra space
legend(0.35,0.2, legend=c("0/Baseline/0","249/Collection1/23","489/Collection2/72","727/Collection3/114","989/Collection4/150","1242/Collection5/178","1490/Collection6/204","1738/Collection7/227","1980/Collection8/249","2224/Collection9/269",
                          "2496/Collection10/291","2784/Collection11/310","3103/Collection12/333","3434/Collection13/353","3666/Collection14/373","3962/Collection15/422","4229/Collection16/477","4527/Collection17/520","4863/Collection18/556","5200/Collection19/579"), 
       col= rainbow(20)[meta$Collection],
       pch=20,
       cex = 0.75)




####Extra in case it can help figure out why the above code isn't working
#Using ggplot this line colors them correctly, but legend is out of order
plot_ordination(physeq.tree, wUF.ordu, type = "sites", color ="DaysCollectionADD")

#These lines also worked
require(ggplot2)
meta$mylegendcolumnNew <- factor(meta$ADDCollectionDays, levels=c("0/Baseline/0","249/Collection1/23","489/Collection2/72","727/Collection3/114",
                                                                  "989/Collection4/150","1242/Collection5/178","1490/Collection6/204","1738/Collection7/227",
                                                                  "1980/Collection8/249","2224/Collection9/269","2496/Collection10/291","2784/Collection11/310",
                                                                  "3103/Collection12/333","3434/Collection13/353","3666/Collection14/373","3962/Collection15/422",
                                                                  "4229/Collection16/477","4527/Collection17/520","4863/Collection18/556","5200/Collection19/579"), 
                                 ordered = TRUE)
plot_ordination(physeq.tree, wUF.ordu, type = "sites", color="mylegendcolumnNew")+ #this keeps the order 
  labs(colour="ADD/Collection/Days")

plot_ordination(physeq.tree, wUF.ordu, type = "sites", color="mylegendcolumnNew") +
  scale_colour_manual(values = rainbow(20))+
  labs(colour="ADD/Collection/Days")
