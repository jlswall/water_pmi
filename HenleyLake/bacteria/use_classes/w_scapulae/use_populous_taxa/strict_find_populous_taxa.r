library("tidyverse")


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "classes"

## Read in cleaned-up phyla, orders, or families taxa.
allT <- read_csv(paste0("../../../", taxalevel, "_massaged.csv"))
## ##################################################


## ##################################################
## Are we dealing with ribs, scapulae, or water observations?
obstype <- "Scapula"

## Filter the data to just that type.
taxaT <- allT %>% filter(type==obstype)
rm(allT)
## ##################################################


## ##################################################
## Let's try making a more stringent cutoff, so that our random forest
## model only uses taxa that are relatively prevalent.

freqCutoff <- 0.05

## Make variable denoting individual scapulae.
taxaT$scapnum <- substring(taxaT$sampleName, first=3, last=4)

## To be included a taxon must exceed 0.05 of the counts for at least
## 2 scapulae (on any date).  We make a list of those taxa below.
freqTaxaT <- taxaT %>%
  mutate(isExceed=(fracBySample>=freqCutoff)) %>%
  group_by(taxon, scapnum) %>%
  summarize(numExceed = sum(isExceed)) %>%
  filter(numExceed > 0) %>%
  group_by(taxon) %>%
  summarize(numScapulaeExceed = n()) %>% 
  filter(numScapulaeExceed > 1) %>%
  filter(taxon!="Rare") %>%
  select(taxon)


## Rename taxa that occur less than the frequency cutoff allows as
## "rare".  Then, sum all these "rare" taxa into one row for each
## sample.
commontaxaT <- taxaT
## Find rows which contain taxa which do not meet the frequency cutoff
## for that type.
indicRare <- !(commontaxaT$taxon %in% freqTaxaT$taxon)
## For these "rare" taxa, we change taxon name to rare.
commontaxaT[indicRare, "taxon"] <- "Rare"

## Now, we add up the rare counts so that "Rare" appears only one for
## each sample.
commontaxaT <- commontaxaT %>%
  group_by(sampleName, taxon, type, degdays) %>%
  summarize(counts = sum(counts), fracBySample=sum(fracBySample))


## Check that the fractions add up to 1, appropriately.
unique(
    commontaxaT %>%
    group_by(sampleName) %>%
    summarize(sumFracBySample = sum(fracBySample)) %>%
    pull(sumFracBySample)
)

## Save list of frequent taxa by type to a CSV file.
write.csv(freqTaxaT, file="class_freq_taxa_by_type.csv", row.names=F)
## ##################################################


## ##################################################
## Save the tibble to a file for use in separate code
## for graphing and analysis.

write.csv(commontaxaT, file="populous_classes.csv", row.names=FALSE)
## ##################################################



## ##################################################
## Put the data in wide format; remove days, subj, and rare taxa.

## Move back to wide format.
wideT <- taxaT %>%
  filter(taxon!="Rare") %>%
  select(degdays, sampleName, taxon, fracBySample) %>%
  spread(taxon, fracBySample) %>%
  select(-sampleName)

rm(taxaT)
## ##################################################
