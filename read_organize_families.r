library("tidyverse")
library("readxl")
library("stringr")


## ##################################################
## Read in family-level taxonomic data from Excel sheet.
fileNm <- "orig_data/HenleyLake_Taxonomy.xlsx"

## Rows 523-527 contain summary data ("total bacteria", "total
## unclassified", etc.), so we exclude this info.  Row 522 is empty,
## and row 1 is the header.  We should have 520 rows of data.
rawAllT <- read_excel(path=fileNm, sheet="TaxLevel 5-Family", n_max=520)
rm(fileNm)

## For all observations, the value of "taxlevel" is 5 (denoting
## family-level), so we remove this column.  We also remove "rankID"
## and "daughterlevels", since we aren't using them.
rawAllT <- rawAllT %>% select(-taxlevel, -rankID, -daughterlevels)

## We also check the "total" column.  For most rows, the "total" is
## the sum of all the subsequent columns.  However, for 107 of the
## rows, the "total" is higher than the sum of the subsequent rows,
## due to some control samples, etc. which were purposefully omitted.
## In these cases, the "total" number should match the number for that
## taxon given in the "Original" tab.  (I've spotchecked several of
## these cases, and that seems to be true.)
sumCols <- apply(rawAllT[,-c(1, 2)], 1, sum)
sumMinusTotal <- sumCols - (rawAllT %>% pull(total)) 
percDiff <- 100 * sumMinusTotal/(rawAllT %>% pull(total))
whichDiffer <- which(abs(sumMinusTotal) >= 1)
## There are 107 rows (out of 520) for which the totals in the "total"
## column (column "E") do not match the sum of the all the
## observations in that row.
sum(abs(sumMinusTotal) > 0)
## Some of these differences are quite big percentage differences, but
## in all cases the calculated total is less than the total given in
## column "E".
summary(percDiff)
## After these checks, we can remove the "total" column".
rawAllT <- rawAllT %>% select(-total)

rm(sumCols, sumMinusTotal, percDiff, whichDiffer)
## ##################################################




## ##################################################
## Read information about collections, sample types, dates, etc. from
## spreadsheet in "HenleyLake_SampleInformation.xlsx".

## The first 3 rows of the spreadsheet are all part of the header, so
## skip those.
fileNm <- "orig_data/HenleyLake_SampleInformation.xlsx"
samplingT <- read_excel(fileNm, skip=3,
                        col_names=c("date", "degdays", "season",
                                    "location", "type", "collection",
                                    "extractMethod", "sampleName"))

## Switch date to "Date" class.  Default is POSIX class, which is for
## datetimes.  Here, we only have the dates, not a time of day.
samplingT$date <- as.Date(samplingT$date)

## All of these samples were taken at the same location (Henley Lake)
## and using the sampe extraction method, so "location" and
## "extractMethod" columns are not necessary.
samplingT <- samplingT %>% select(-location, -extractMethod)


## The type can be either scapula, rib, or water.  Simplify by using
## just S, R, or W.
samplingT$type <- substring(samplingT$type, first=1, last=1)


## The "collection" column contains one of "Baseline", "Collection 1",
## ..., "Collection 19".  We save space by refering to these as "B",
## "1", ..., "19". 
tmp <- samplingT$collection
tmp <- str_remove(tmp, pattern="Collection ")
tmp[tmp=="Baseline"] <- "B"
samplingT$collection <- tmp
## To keep the chronological order (and avoid ABC order), we specify
## this as an ordered factor.
samplingT$collection <- ordered(samplingT$collection, levels=c("B", as.character(1:19)))


## Some of the sample names have spaces in them.  That is not true of
## the corresonding sample names in the taxonomy file.  So, we take
## the spaces out for consistency.
## Take out spaces and store elsewhere.
tmp <- str_replace_all(samplingT$sampleName, " ", "")
## Check to see which sample names were changed, and make sure all looks ok.
cbind(samplingT$sampleName, tmp)[samplingT$sampleName!=tmp,]
## Replace sample name column with adjusted names.
samplingT$sampleName <- tmp

## #######################
## Ideally, we would have 5 rib samples, 5 scapula samples, and 1
## water sample for each collection day.  That would 11 samples in
## total for each of 20 collection days.  In reality, all samples are
## not available for every collection day.

## Count the number of such samples for each collection day.
samplingT %>% group_by(collection, type) %>% summarize(nObs=n())

## Identify collection days on which samples are missing.
samplingT %>%
  group_by(collection, type) %>%
  summarize(nObs=n()) %>%
  filter( ((type=="R" | type=="S") & nObs<5) | (type=="W" & nObs<1) )
## This table shows which collections have missing scapula or rib
## samples.  Water samples are not missing on any of the collection
## days.
##    collection type   nObs
##    <ord>      <chr> <int>
##  1 B          S         2
##  2 3          R         4
##  3 4          R         4
##  4 6          R         4
##  5 9          R         3
##  6 10         R         4
##  7 11         R         3
##  8 12         R         4
##  9 13         R         4
## 10 14         R         3
## 11 15         R         3
## 12 16         R         4
## 13 17         R         3
## 14 19         R         4
## #######################

rm(fileNm, tmp)
## ##################################################



## ##################################################
## Concentrate on the taxa counts for the individual cadavers on the
## various collection days.  Transfer from wide format to long format.
## Parse out the column names to get collection day, sample type, etc.

## #######################
## Go from wide to long format.
rawIndivT <- rawAllT %>%
  gather(sampleName, counts, -taxon)


## Join with sampling information tibble to get information about
## accumulated degree days.  That tibble already has information about
## sample type (R, S, or W) and collection number (baseline,
## collection 1-19), so that we don't have to parse it out of the taxa
## spreadsheet's column names.
indivT <- rawIndivT %>% inner_join(samplingT)
rm(rawIndivT, rawAllT, samplingT)
## #######################
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
## level.  We see that about 88.5% of the counts can be attributed at
## the family level.
indivT %>% group_by(taxLvl) %>% summarize(sumCounts=sum(counts), percCounts=100*sum(counts)/sum(indivT$counts))


## For each collection day, what percentage of counts can be
## classified to the each level?
percTaxLvlByDayT <- indivT %>%
  group_by(degdays, collection, taxLvl) %>%
  summarize(sumTaxLvlCts=sum(counts)) %>%
  inner_join(indivT %>% group_by(collection) %>%
             summarize(dayTotals=sum(counts))
             ) %>%
  mutate(percTaxLvlCts = 100*sumTaxLvlCts/dayTotals)


## Plot percentage classified at each level (e.g. family, phylum,
## etc.) by ADD.
ggplot(percTaxLvlByDayT) +
  geom_point(aes(x=degdays, y=percTaxLvlCts, color=taxLvl)) +
  labs(x="Accumulated degree days", y="Percentage taxa classified each level")
## ##################################################




## ##################################################
## For use in graphs and in calculating percentages later, we need
## total counts (over all taxa, unclassified taxa excluded) for each
## sample (sample is individual rib or scapula on a particular
## collection time).

## ########## WORKING HERE!
## May want to remove all taxa that couldn't be classified at family
## level from these counts, or we may want to try it both ways.

## Total taxa counts by sample.
ctBySampleT <- indivT %>%
  group_by(degdays, type, sampleName) %>%
  summarize(totals=sum(counts))
## ##################################################




## #######################
## We need to separate out the baseline samples, which were taken
## before the cadavers were ever placed in the water.
baselineT <- allT %>% filter(collection=="Baseline") 
## Now, remove these from the rest of the tibble (which includes data
## from each collection day).
indivT <- indivT %>% filter(collection!="Baseline")

rm(isBaseline)
## #######################






## ############## OLD STUFF FOLLOWS ##############



## ##################################################
## Make other adjustments to the dataset so that it's easier to use.

## Remove the counts associated with unclassifed taxa,
## "Incertae_Sedis", and "uncultured" (*_unclassified, "_uncultured",
## "Incertae_Sedis").  Also, include accum. degree days in the tibble.
indivT <- indivT %>%
  filter(!grepl("_unclassified|_uncultured|Incertae_Sedis", taxon)) %>%
  left_join(timeDF, by="days")

## Remove the counts associated with order "Mammalia" (which is
## probably the pig's DNA), and with order "Aves" (birds, not sure why
## that's in there).
indivT <- indivT %>%
  filter(!grepl("Mammalia|Aves", taxon))
## ##################################################




## ##################################################
## For use in graphs and in calculating percentages later, we need
## total counts (over all taxa, unclassified taxa excluded) by:
##   Each pig and each day 
##   Each day (all pigs combined)

## Total taxa counts by day and subject (each pig separately).
ctBySubjDayT <- indivT %>%
  group_by(days, degdays, subj) %>%
  summarize(totals=sum(counts))

## Total taxa counts by day (all pigs combined).
ctByDayT <- indivT %>%
  group_by(days, degdays) %>%
  summarize(totals = sum(counts))
## ##################################################




## ##################################################
## Some taxa don't occur frequently.  It's hard to make a hard cutoff
## for what constitutes "frequently".  There are 164 (classified) taxa
## in the dataset, and a lot of them appear in less than 0.1% of
## samples.

## I'm going to set the cutoff at 1% (0.01).  This means that in order
## to be included in the dataset, a specific taxa must make up at
## least 1% of the total counts on at least 1 day for at least 1
## cadaver.
freqCutoff <- 0.01

## Get list of maximum taxa percentages sorted in descending order:
data.frame(indivT %>%
  left_join(ctBySubjDayT) %>%
  mutate(fracBySubjDay = counts/totals) %>%
  group_by(taxon) %>%
  summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
  arrange(desc(maxFracBySubjDay))
)


## Save the taxa names (in a tibble) which satisfy the frequency
## cutoff.
freqTaxaT <- indivT %>%
  left_join(ctBySubjDayT) %>%
  mutate(fracBySubjDay = counts/totals) %>%
  group_by(taxon) %>%
  summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
  filter(maxFracBySubjDay >= freqCutoff) %>%
  arrange(desc(maxFracBySubjDay)) %>%
  select(taxon)


## Rename taxa that occur less than the frequency cutoff allows as
## "rare".  Then, sum all these "rare" taxa into one row.
commontaxaT <- indivT
commontaxaT[!(commontaxaT$taxon %in% freqTaxaT$taxon), "taxon"] <- "Rare"
commontaxaT <- commontaxaT %>%
  group_by(days, degdays, subj, taxon) %>%
  summarize(counts = sum(counts))

## Remove the list of taxa names that satisfied the frequence cutoff.
rm(freqTaxaT)
## ##################################################




## ##################################################
## Add percentages by subj/day to the commontaxaT table.

## Use the table of total counts by subj/day to find the fraction
## represented by each taxa for each subj/day.
commontaxaT <- commontaxaT %>%
  left_join(ctBySubjDayT) %>%
  mutate(fracBySubjDay=counts/totals) %>%
  select(-totals)


## Check that the fractions add up to 1, appropriately.
unique(
    unlist(commontaxaT %>%
           group_by(days, subj) %>%
           summarize(sumFracBySubjDay = sum(fracBySubjDay)) %>%
           ungroup() %>%
           select(sumFracBySubjDay))
)
## ##################################################




## ##################################################
## Save the tibble to a file for use in separate code
## for graphing and analysis.

## I have to use the base R write.csv() routine, because write_csv
## will write out scientific notation, which read_csv() doesn't read
## in properly.
write.csv(commontaxaT, file="families_massaged.csv", row.names=FALSE)
## ##################################################
