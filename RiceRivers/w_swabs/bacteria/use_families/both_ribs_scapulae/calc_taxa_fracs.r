library("tidyverse")


# ##################################################
# Are we dealing with phlya, orders, or families?
taxalevel <- "families"

# Read in cleaned-up phyla, orders, or families taxa.
allT <- read_csv(paste0("../../", taxalevel, "_massaged.csv"))
# ##################################################


# ##################################################
# We want only rib and scapulae data.  In this case, all observations are from
# rib or scapula (no water data).  We check that here.
allT %>% select(type) %>% distinct()
# ##################################################



# ##################################################
# When we put together the "massaged" data, we only kept taxa which 
# met the 1% cutoff for more than 1 sample for the particular type (ribs or
# scapulae.)  Let's take a look at the prevalence of taxa which met the
# requirement for either ribs or scapulae but not both.

# List number of unique taxa for each group (rib or scapula).
allT %>% group_by(type) %>% summarize(uniqtaxa=n_distinct(taxon))

# Identify taxa which have sufficient numbers for both ribs and scapulae.
# There are 38 such taxa (not including the "Rare" group).
bothtaxaT <- allT %>% filter(type=="Rib") %>% distinct(taxon) %>%
    inner_join(allT %>% filter(type=="Scapula") %>% distinct(taxon))

# Take a look at the counts for taxa which have sufficient numbers for just one
# type, but not both.  After the baseline (day 0), such taxa are not present in
# large quantities (less than 10% on any one collection).
allT %>% anti_join(bothtaxaT) %>% arrange(desc(fracBySample)) %>% print(n=20)

# Need to re-classify the taxa with insufficient numbers into the "Rare
# category".  (If we just subset to include only these taxa which have
# sufficient numbers then the fractions for each sample won't sum up to one.)
newtaxon <- ifelse(allT$taxon %in% (bothtaxaT %>% pull(taxon)), allT$taxon,
  "Rare")
allT$taxon <- newtaxon
# Now, we add up the rare counts so that "Rare" appears only one for
# each sample.
allT <- allT %>%
  group_by(sampleName, taxon, type, degdays) %>%
  summarize(counts = sum(counts), fracBySample=sum(fracBySample))

## Check to make sure fractions add to 1 for each collection.
allT %>% group_by(sampleName) %>%
  summarize(sumfracBySample=sum(fracBySample)) %>% distinct(sumfracBySample)

rm(newtaxon, bothtaxaT)
# ##################################################



# ##################################################
# Save the tibble to a file for use in separate code for graphing and analysis.

write.csv(allT, file="combine_rib_scapula_massaged.csv", row.names=FALSE)
# ##################################################

