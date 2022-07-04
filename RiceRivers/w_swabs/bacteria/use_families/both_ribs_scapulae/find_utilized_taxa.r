library("tidyverse")

# Read in combined rib/scapula dataset (families taxa).
allT <- read_csv("combine_rib_scapula_massaged.csv")
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
# f__Bradyrhizobiaceae 
# f__Campylobacteraceae
# f__Carnobacteriaceae 
# f__Desulfobacteraceae
# f__Haliangiaceae     
# f__Neisseriaceae     
# f__Polyangiaceae     
# f__Rhodospirillaceae 
# f__Sediment_4        
# f__Thiotrichaceae    
# f__Weeksellaceae 

# Which taxa from the scapula-based model were excluded because of this?
scaponlyTaxa %>% anti_join(combinedTaxa)
# f__Chromatiaceae       
# f__Geobacteraceae      
# f__Hydrogenophilaceae  
# f__Lachnospiraceae     
# f__Marinilabiaceae     
# f__mb2424              
# f__Methylocystaceae    
# f__Mogibacteriaceae    
# f__Oxalobacteraceae    
# f__Procabacteriaceae   
# f__Rhizobiaceae        
# f__Spirochaetaceae     
# f__Syntrophaceae       
# f__Syntrophobacteraceae
# f__Tissierellaceae 