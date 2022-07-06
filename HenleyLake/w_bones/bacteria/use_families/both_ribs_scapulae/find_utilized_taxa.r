library("tidyverse")

# Read in combined rib/scapula dataset (families taxa).
allT <- read_csv("bones_combine_rib_scapula_massaged.csv")
# Get the list oif all the taxa utilized in the combined ribs/scapulae analyses.
combinedTaxa <- allT %>% filter(taxon != "Rare") %>% distinct(taxon)
rm(allT)


# #########################
# The combined model only uses taxa which meet the prevalence cutoff for both
# rib and scapula-based swabs.

# Read in list of taxa used for ribs-only models and for scapula-only models.
allribscapTaxaT <- read_csv("../../family_freq_taxa_by_type.csv")
ribonlyTaxa <-  allribscapTaxaT %>% filter(type=="Rib") %>% select(taxon)
scaponlyTaxa <-  allribscapTaxaT %>% filter(type=="Scapula") %>% select(taxon)
rm(allribscapTaxaT)

# Which taxa from the rib-based model were excluded because of this?
ribonlyTaxa %>% anti_join(combinedTaxa)
# f__p_2534_18B5

# Which taxa from the scapula-based model were excluded because of this?
scaponlyTaxa %>% anti_join(combinedTaxa)
# f__Campylobacteraceae
# f__Crenotrichaceae
# f__Desulfobulbaceae
# f__Erysipelotrichaceae
# f__Geobacteraceae
# f__JTB215
# f__Oxalobacteraceae
# f__Phormidiaceae
# f__Pirellulaceae
# f__Rhodospirillaceae
# f__Sediment_4