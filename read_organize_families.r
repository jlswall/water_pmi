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


## Some of the sample names have spaces in them.  That is not true of
## the corresonding sample names in the taxonomy file.  So, we take
## the spaces out for consistency.
## Take out spaces and store elsewhere.
tmp <- str_replace_all(samplingT$sampleName, " ", "")
## Check to see which sample names were changed, and make sure all looks ok.
cbind(samplingT$sampleName, tmp)[samplingT$sampleName!=tmp,]
## Replace sample name column with adjusted names.
samplingT$sampleName <- tmp

rm(fileNm, tmp)
## ##################################################



## ##################################################
## Concentrate on the taxa counts for the individual cadavers on the
## various collection days.  Transfer from wide format to long format.
## Parse out the column names to get collection day, sample type, etc.

## #######################
## Go from wide format to long format. 

## Go from wide to long format.
indivT <- rawAllT %>%
  gather(origColNm, counts, -taxon)
## Each column started with "HL" (abbrev. for location), so we "HL".
indivT$reviseNm <- str_remove(indivT$origColNm, pattern="HL")
## Each column ends with "CS", so we remove "CS". 
indivT$reviseNm <- str_remove(indivT$reviseNm, pattern="CS")
## #######################

## #######################
## We need to separate out the baseline samples, which were taken
## before the cadavers were ever placed in the water.  These are the
## columns with names ending with "B".  They do include the pattern
## "C#" (number of collection day).
isBaseline <- str_detect(indivT$reviseNm, pattern="B$")
baselineT <- indivT[isBaseline,]
## Now, remove these from the rest of the tibble (which includes data
## from each collection day).
indivT <- indivT[!isBaseline,]

rm(isBaseline)
## #######################


## #######################
## Determine collection day.

## For the non-baseline data, each column name contained the pattern
## "C#", where the "#" denotes the number of the collection day.
indivT <- separate(indivT, reviseNm, sep="C", into=c("siteindiv", "day"), convert=T)
## #######################


## #######################
## Figure out site of each sample.

## The first letter in the "siteindiv" column denotes the site of the
## sample, either water ("W"), rib ("R"), or scapula ("S").

## Separate water samples from submerged cadaver samples.  There is
## only one water sample for each collection day.
waterT <- indivT %>% filter(siteindiv=="W") %>% rename(site=siteindiv)

## After separating the water samples, separate the "siteindiv" column
## into a site column ("R" or "S") and a subj column (number of the
## cadaver).
mainT <- indivT %>%
  filter(siteindiv!="W") %>%
  separate(siteindiv, sep=1, into=c("site", "subj"), convert=T)

rm(indivT)
## #######################



## ######### WORKING HERE!!!!!


## #######################
## Put individual counts and daily sums into different tables.

## Identify column names starting with "A". Save these as the counts
## for individual pigs on the various days.
namesA <- colnames(rawAllT)[substring(first=1, last=1, colnames(rawAllT))=="A"]
wideIndivT <- rawAllT[,c("taxon", namesA)]

## Identify column names starting with "T".  The values in these
## columns are the sums of the individual pigs at each time point.
## The names of these columns contain the number of days since death
## and the accumulated degree days.  The number of days since death
## immediately follows the "T", and the number of accumulated degree
## days follows the "_".
namesT <- colnames(rawAllT)[substring(first=1, last=1, colnames(rawAllT))=="T"]
## Separate the days from the accumulated degree days.
timeDF <- separate(data.frame(x=substring(namesT, first=2),
                              stringsAsFactors=F),
                   x, sep=" - ", into=c("days", "degdays"), convert=T)
## Note: number of days and accum. degree days are strongly correlated.
with(timeDF, cor(degdays, days))


## Extract the columns with the taxa names and the sums across pigs
## for each time point and taxon.
wideSumsT <- rawAllT[,c("taxon", namesT)]


## The last few rows are special cases, so I exclude them from the
## tables of counts for the taxa.  They are: "Eukaryota",
## '"Unclassified"', "% Unclassified", "Total Classified"
wideIndivT <- wideIndivT %>% filter(!(taxon %in% c('Eukaryota', '"Unclassified"', "% Unclassified", "Total Classified")))
wideSumsT <- wideSumsT %>% filter(!(taxon %in% c('Eukaryota', '"Unclassified"', "% Unclassified", "Total Classified")))


rm(namesA, namesT)
## #######################


## #######################
## Go from wide format to long format.  Check the columns and rows
## containing totals.

## Go from wide to long format, separating column names into subject,
## day, and swab number.
indivT <- wideIndivT %>%
  gather(indiv_time, counts, -taxon) %>%
  separate(indiv_time, sep="T", into=c("subj", "days"), convert=T) %>%
  separate(days, sep="S", into=c("days", "swab"), convert=T)


## ######
## Make sure the totals for each day and taxa (across subjects) match
## those I get from doing the sums.
mysumsT <- indivT %>% group_by(taxon, days) %>% summarize(counts=sum(counts))
## Adjust the times so that they match up with the first few
## characters of the original column names for comparison.
mysumsT[,"compareDays"] <- paste0("T", as.vector(as.matrix(mysumsT[,"days"])), " - ")
mysumsT <- mysumsT %>% select(-days) %>% spread(compareDays, counts)
## Put column names of mysumsT into same order as those for
## wideSums. It's sufficient to just check that the "T##" part
## matches.
match.order <- match(substring(colnames(wideSumsT), first=1, last=4), substring(colnames(mysumsT), first=1, last=4))
mysumsT <- mysumsT[,match.order]
## Put rows of mysumsT into same order as those for wideSums.
match.order <- match(wideSumsT$taxon, mysumsT$taxon)
mysumsT <- mysumsT[match.order,]
unique(as.vector(as.matrix(wideSumsT[,-1]) - as.matrix(mysumsT[,-1])))

rm(match.order, wideSumsT, mysumsT)
## ######
## #######################
## ##################################################




## ##################################################
## Find the total percentage of counts which are unclassified.

## About 32.76% are unclassified.
pull(indivT %>%
     filter(grepl("_unclassified|_uncultured|Incertae_Sedis", taxon)) %>%
     summarize(total_uncl=sum(counts)), "total_uncl") / sum(indivT[,"counts"])
## Unclassified: 0.3276325
## ##################################################




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
