library("tidyverse")
library("readxl")
library("stringr")


## ##################################################
## Read in family-level taxonomic data from Excel sheet.
fileNm <- "orig_data_files/Rice_River_18S_Taxonomy.xlsx"

## Rows 178-182 contain summary data ("Eukaryota", "Total
## Unclassified", etc.), so we exclude this info.  Rows 176-177 are
## empty, and row 1 is the header.  We should have 174 rows of data.
## Note that row 2 contains counts of all eukaryotes.
rawAllT <- read_excel(path=fileNm, sheet="Tax level 5- Family", n_max=174)
rm(fileNm)

## The first row contains counts for all eukaryotes, so the totals in
## this row should equal the sums of all the rows after it.
all.equal(as.vector(as.matrix(rawAllT[1, 6:ncol(rawAllT)])),
          as.vector(as.matrix(apply(rawAllT[2:nrow(rawAllT), 6:ncol(rawAllT)], 2, sum))))
## So, we remove this first row with totals of all eukaryotes.
rawAllT <- rawAllT[-1,]

## For all observations, the value of "taxlevel" is 5 (denoting
## family-level), so we remove this column.  We also remove "rankID"
## and "daughterlevels", since we aren't using them.
rawAllT <- rawAllT %>% select(-taxlevel, -rankID, -daughterlevels)

## The `New total` column contains the totals of all the columns following.
all.equal(rawAllT %>% pull(`New total`), apply(rawAllT[,3:ncol(rawAllT)], 1, sum))

## After these checks, we remove the "newtotal" column.
rawAllT <- rawAllT %>% select(-`New total`)
## ##################################################



## ##################################################
## Dealing with unusual columns.

## Column "RRWBCS" (denoted EV in the spreadsheet) doesn't have a
## match in the sampling information file.)  I think this column
## contains baseline (day 0, when ribs and scapulae were placed in the
## river) eukaryote counts from the water.  We have no other baseline
## data, so we remove this column.
rawAllT <- rawAllT %>% select(-RRWBCS)

## We have some column names in which the ending "S" was omitted.  We
## correct this so that the columns can be matched with the sample
## names from the sample information file. Column names affected are:
## RRR2C16C, RRR3C16C, RRR4C16C, RRR5C16C, RRS2C17C
tmp <- colnames(rawAllT)
tmp[tmp=="RRR2C16C"] <- "RRR2C16CS"
tmp[tmp=="RRR3C16C"] <- "RRR3C16CS"
tmp[tmp=="RRR4C16C"] <- "RRR4C16CS"
tmp[tmp=="RRR5C16C"] <- "RRR5C16CS"
tmp[tmp=="RRS2C17C"] <- "RRS2C17CS"
colnames(rawAllT) <- tmp
## ##################################################



## ##################################################
## Concentrate on the taxa counts for the individual cadavers on the
## various collection days.  Transfer from wide format to long format.

## Go from wide to long format.
rawIndivT <- rawAllT %>%
  gather(sampleName, counts, -taxon)
## ##################################################


## ##################################################
## Join with sampling information to get information about accumulated
## degree days.  That tibble already has information about sample type
## (R, S, or W) and collection number, so that we don't have to
## parse it out of the taxa spreadsheet's column names.

## Read in sampling information for Henley Lake eukaryotes.  Note that
## the baseline data isn't included here, though it was for the
## bacterial data.
samplingT <- read_csv("sampling_info.csv")

## NOTE: We also remove date, season, and collection variables.  This
## info is contained in the sampleName variable.
indivT <- rawIndivT %>%
  left_join(samplingT) %>%
  select(-date, -season, -collection)
rm(rawIndivT, rawAllT, samplingT)
## ##################################################



## ##################################################
## Make a new, more readable taxon column.  Count number of unique
## taxa observed for the various types (rib, scapula, water).

## Column names with dashes and periods can be a problem, so I replace
## them with underscores.  This affects taxa such as "MAST-2_fa" and
## "P34.45_fa".
indivT$newtaxon <- indivT$taxon
indivT$newtaxon <- gsub(indivT$newtaxon, pattern="-", replacement="_")
indivT$newtaxon <- gsub(indivT$newtaxon, pattern=".", replacement="_", fixed=TRUE)
## Remove the column of the original, less-R-friendly taxon names from
## the tibble to avoid confusion.  Rename new column.
indivT <- indivT %>% select(-taxon)
indivT <- rename(indivT, taxon = newtaxon)


## Count the number of unique family-level taxa observed for the
## various types (rib, scapula, water).
indivT %>% filter(counts>0) %>% group_by(type) %>% distinct(taxon) %>% summarize(n=n())
##   type        n
##   <chr>   <int>
## 1 Rib       116
## 2 Scapula   146
## 3 Water      74
## ##################################################




## ##################################################
## The eukaryote data includes many taxa that couldn't be classified
## down to the family level.  So, we consider the percentage of each
## sample that is made up of each taxon (no matter the classification
## level).  Then, all taxa that do not meet a certain criteria will be
## grouped together into a "Rare" category.  The taxa which do meet
## the criteria will be considered in the random forest model.

## To calculate percentages, we need total counts for each sample.
## We find that the total for each sample is the same (5048).
ctBySampleT <- indivT %>%
  group_by(sampleName) %>%
  summarize(totals=sum(counts))


## #########################
## Some taxa don't occur frequently.  It's hard to make a hard cutoff
## for what constitutes "frequently", but many of the taxa make up
## less than 1% of the counts for any sample.

## I'm going to set the frequency cutoff at 1% (0.01).  This means
## that in order to meet the cutoff, a specific family-level taxa must
## make up at least 1% of the total family-level counts for that
## sample.
freqCutoff <- 0.01

## ## Get list of maximum taxa percentages for various types (rib,
## ## scapula, water, mud) sorted in descending order:
## indivT %>%
##   left_join(ctBySampleT) %>%
##   mutate(fracBySampleName = counts/totals) %>%
##   group_by(type, taxon) %>%
##   summarize(maxFracBySampleName = max(fracBySampleName)) %>%
##   ## filter(maxFracBySampleName >= freqCutoff) %>%
##   arrange(type, desc(maxFracBySampleName)) %>%
##   print(n = Inf)


## Count how many times each taxa beats the cutoff (freqCutoff) within
## each type.  Then, keep only those which meet cutoff for more than 1
## sample.  Save those taxa names for each type.
freqTaxaByTypeT <- indivT %>%
  left_join(ctBySampleT) %>%
  mutate(fracBySampleName = counts/totals,
         isExceed=(fracBySampleName>=freqCutoff)) %>%
  group_by(type, taxon) %>%
  summarize(numExceed = sum(isExceed)) %>%
  filter(numExceed > 1) %>%
  ## arrange(type, desc(numExceed)) %>%
  select(type, taxon)


## For each type, how many taxa meet the cutoff?
freqTaxaByTypeT %>% group_by(type) %>% summarize(n=n())
##   type        n
##   <chr>   <int>
## 1 Rib        16
## 2 Scapula    27
## 3 Water      27

## Save list of frequent taxa by type to a CSV file.
write.csv(freqTaxaByTypeT, file="family_freq_taxa_by_type.csv", row.names=F)
## #########################
## ##################################################



## ##################################################
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
