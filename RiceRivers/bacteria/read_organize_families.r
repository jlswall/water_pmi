library("tidyverse")
library("readxl")
library("stringr")


## ##################################################
## Read in family-level taxonomic data from Excel sheet.
fileNm <- "orig_data_files/2019_RiceRiversCenter_Taxonomy.xlsx"

## Rows 632-636 contain summary data ("total bacteria", "total
## unclassified", etc.), so we exclude this info.  Rows 629-631 are
## empty, and row 1 is the header.  We should have 627 rows of data.
## Note that row 2 contains counts of all bacteria.
rawAllT <- read_excel(path=fileNm, sheet="TaxLevel5-Family", n_max=627)
rm(fileNm)

## The first row contains counts for all bacteria, so the totals in
## this row should equal the sums of all the rows after it.
all.equal(as.vector(as.matrix(rawAllT[1, 5:ncol(rawAllT)])),
          as.vector(as.matrix(apply(rawAllT[2:nrow(rawAllT), 5:ncol(rawAllT)], 2, sum))))
## So, we remove this first row with totals of all bacteria.
rawAllT <- rawAllT[-1,]

## For all observations, the value of "taxlevel" is 5 (denoting
## family-level), so we remove this column.  We also remove "rankID"
## and "daughterlevels", since we aren't using them.
rawAllT <- rawAllT %>% select(-taxlevel, -rankID, -daughterlevels)

## The "newtotal" column contains the totals of all the columns following.
all.equal(rawAllT$newtotal, apply(rawAllT[,4:ncol(rawAllT)], 1, sum))

## After these checks, we remove the "total" and "newtotal" columns.
rawAllT <- rawAllT %>% select(-total, -newtotal)
## ##################################################



## ##################################################
## Concentrate on the taxa counts for the individual cadavers on the
## various collection days.  Transfer from wide format to long format.
## Parse out the column names to get collection day, sample type, etc.

## Go from wide to long format.
rawIndivT <- rawAllT %>%
  gather(sampleName, counts, -taxon)


## Join with sampling information tibble to get information about
## accumulated degree days.  That tibble already has information about
## sample type (R, S, or W) and collection number (baseline,
## collection 1-19), so that we don't have to parse it out of the taxa
## spreadsheet's column names.
## NOTE: We also remove date, season, and collection variables.  This
## info is contained in the sampleName variable.
indivT <- rawIndivT %>%
  inner_join(samplingT) %>%
  select(-date, -season, -collection)
rm(rawIndivT, rawAllT, samplingT)
## ##################################################



## ##################################################
## Consider level of classification for each taxa.  (Some could not be
## classified down to the family level, but maybe the order could be
## determined.)

## Make new variable to indicate the most precise taxon which could be
## identified.
indivT$taxLvl <- ""
indivT$taxLvl[str_detect(indivT$taxon, "f__")] <- "family"
indivT$taxLvl[str_detect(indivT$taxon, "o__")] <- "order"
indivT$taxLvl[str_detect(indivT$taxon, "c__")] <- "class"
indivT$taxLvl[str_detect(indivT$taxon, "p__")] <- "phylum"
indivT$taxLvl[str_detect(indivT$taxon, "k__")] <- "kingdom"
indivT$taxLvl <- ordered(indivT$taxLvl, levels=c("family", "order", "class", "phylum", "kingdom"))


## Find percentage of counts can be classified down to the family
## level, across all types.  We see that about 88.5% of the counts can
## be attributed at the family level.
indivT %>% group_by(taxLvl) %>% summarize(sumCounts=sum(counts), percCounts=100*sum(counts)/sum(indivT$counts))

## Find percentage of counts which can be classified down to the
## family level, type by type.
indivT %>% group_by(type, taxLvl) %>% summarize(sumCountsByTypeTaxLvl=sum(counts)) %>% left_join(indivT %>% group_by(type) %>% summarize(sumCountsByType=sum(counts))) %>% mutate(perc=100*sumCountsByTypeTaxLvl/sumCountsByType)
## Percentages classified to family level: 90.5% (ribs), 88.9%
## (scapulae), 77.7% (water).


## For each sample, what percentage of counts can be classified to the
## each level?
percTaxLvlBySampleT <- indivT %>%
  group_by(degdays, type, sampleName, taxLvl) %>%
  summarize(sumTaxLvlCts=sum(counts)) %>%
  inner_join(indivT %>% group_by(sampleName) %>%
             summarize(sampleTotals=sum(counts))
             ) %>%
  mutate(percTaxLvlCts = 100*sumTaxLvlCts/sampleTotals)


## For each sample, plot percentage classified at each level
## (e.g. family, phylum, etc.) by type (rib, scapula, water) and ADD.
ggplot(percTaxLvlBySampleT) +
  geom_point(aes(x=degdays, y=percTaxLvlCts, color=taxLvl)) +
  facet_wrap(~type) +
  labs(x="Accumulated degree days", y="Percentage taxa classified each level")
## ggsave("family_perc_classif_by_add_type.pdf", width=8.5, height=6, units="in")
## dev.off()
## ##################################################



## ##################################################
## From this point forward, we consider only family-level taxa for our
## model.  So, we remove all counts of non-family-level taxa.  This is
## consistent with the anlaysis we did for Shane's and Luisa's data.
indivT <- indivT %>% filter(taxLvl=="family")


## Make a new, more readable taxon column.
## Column names with open brackets (e.g. "f__[Tissierellaceae]") cause
## problems for functions expecting traditional data frame column
## names (like randomForest function).
indivT$newtaxon <- gsub(indivT$taxon, pattern="\\[", replacement="")
indivT$newtaxon <- gsub(indivT$newtaxon, pattern="]", replacement="")
## Column names with dashes can likewise be a problem, so I replace
## dashes with underscores.
indivT$newtaxon <- gsub(indivT$newtaxon, pattern="-", replacement="_")
## Remove the column of the original, less-R-friendly taxon names from
## the tibble to avoid confusion.  Rename new column.
indivT <- indivT %>% select(-taxon)
indivT <- rename(indivT, taxon = newtaxon)

## Count the number of unique family-level taxa observed for the
## various types (rib, scapula, water).
indivT %>% filter(taxLvl=="family") %>% filter(counts>0) %>% group_by(type) %>% distinct(taxon) %>% summarize(n=n())
##   type        n
##   <chr>   <int>
## 1 Rib       154
## 2 Scapula   210
## 3 Water     216


## For use in graphs and in calculating percentages later, we need
## total counts for each sample.
ctBySampleT <- indivT %>%
  group_by(sampleName) %>%
  summarize(totals=sum(counts))
## ##################################################



## ##################################################
## Some taxa don't occur frequently.  It's hard to make a hard cutoff
## for what constitutes "frequently".  There are 245 family-level taxa
## represented, and a lot of them appear in less than 1% of
## samples.

## I'm going to set the frequency cutoff at 1% (0.01).  This means
## that in order to meet the cutoff, a specific family-level taxa must
## make up at least 1% of the total family-level counts for that
## sample.
freqCutoff <- 0.01

## Get list of maximum taxa percentages for various types (rib,
## scapula, water) sorted in descending order:
## indivT %>%
##   left_join(ctBySampleT) %>%
##   mutate(fracBySubjDay = counts/totals) %>%
##   group_by(type, taxon) %>%
##   summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
##   ## filter(maxFracBySubjDay >= freqCutoff) %>%
##   arrange(type, desc(maxFracBySubjDay)) %>%
##   print(n = Inf)


## Count how many times each taxa beats the cutoff (freqCutoff) within
## each type.  Then, keep only those which meet cutoff for more than 1
## sample.  Save those taxa names for each type.
freqTaxaByTypeT <- indivT %>%
  left_join(ctBySampleT) %>%
  mutate(fracBySubjDay = counts/totals,
         isExceed=(fracBySubjDay>=freqCutoff)) %>%
  group_by(type, taxon) %>%
  summarize(numExceed = sum(isExceed)) %>%
  filter(numExceed > 1) %>%
  ## arrange(type, desc(numExceed)) %>%
  select(type, taxon)

## For each type, how many taxa meet the cutoff?
freqTaxaByTypeT %>% group_by(type) %>% summarize(n=n())
## Save list of frequent taxa by type to a CSV file.
write.csv(freqTaxaByTypeT, file="family_freq_taxa_by_type.csv", row.names=F)


## Rename taxa that occur less than the frequency cutoff allows as
## "rare".  Then, sum all these "rare" taxa into one row for each
## sample.
commontaxaT <- indivT
## Find rows which contain taxa which do not meet the frequency cutoff
## for that type.
indicRare <- !( with(commontaxaT, paste(taxon, type, sep=":")) %in% with(freqTaxaByTypeT, paste(taxon, type, sep=":")) )
## For these "rare" taxa, we change taxon name to rare.
commontaxaT[indicRare, "taxon"] <- "Rare"
## Now, we add up the rare counts so that "Rare" appears only one for
## each sample.
commontaxaT <- commontaxaT %>%
  group_by(sampleName, taxon, type, degdays) %>%
  summarize(counts = sum(counts))


## Use the table of total counts by subj/day to find the fraction
## represented by each taxa for each subj/day.
commontaxaT <- commontaxaT %>%
  left_join(ctBySampleT) %>%
  mutate(fracBySample=counts/totals) %>%
  select(-totals)


## Check that the fractions add up to 1, appropriately.
unique(
    commontaxaT %>%
    group_by(sampleName) %>%
    summarize(sumFracBySample = sum(fracBySample)) %>%
    pull(sumFracBySample)
)


## Remove the list of taxa names that satisfied the frequency cutoff.
rm(indicRare, freqCutoff, freqTaxaByTypeT)
## ##################################################



## ##################################################
## Save the tibble to a file for use in separate code
## for graphing and analysis.

write.csv(commontaxaT, file="families_massaged.csv", row.names=FALSE)
## ##################################################
