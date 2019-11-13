library("tidyverse")
library("readxl")
library("stringr")


## ##################################################
## Read in family-level taxonomic data from Excel sheet.
fileNm <- "orig_data_files/HenleyLake_Taxonomy.xlsx"

## Rows 181-185 contain summary data ("total bacteria", "total
## unclassified", etc.), so we exclude this info.  Row 180 is empty,
## and row 1 is the header.  We should have 178 rows of data.
rawAllT <- read_excel(path=fileNm, sheet="TaxLevel 3-Class", n_max=179)
rm(fileNm)

## For all observations, the value of "taxlevel" is 3 (denoting
## class-level), so we remove this column.  We also remove "rankID"
## and "daughterlevels", since we aren't using them.
rawAllT <- rawAllT %>% select(-taxlevel, -rankID, -daughterlevels)

## We also check the "total" column.  For most rows, the "total" is
## the sum of all the subsequent columns.  However, for 35 of the
## rows, the "total" is higher than the sum of the subsequent rows,
## due to some control samples, etc. which were purposefully omitted.
## In these cases, the "total" number should match the number for that
## taxon given in the "Original" tab.  (I've spotchecked several of
## these cases, and that seems to be true.)
sumCols <- apply(rawAllT[,-c(1, 2)], 1, sum)
sumMinusTotal <- sumCols - (rawAllT %>% pull(total)) 
percDiff <- 100 * sumMinusTotal/(rawAllT %>% pull(total))
whichDiffer <- which(abs(sumMinusTotal) >= 1)
## There are 35 rows (out of 178) for which the totals in the "total"
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
## When we organized the family taxonomic info, we also read in and
## organized the sample information.  Now, we can use the CSV file we
## made of the sample information.  (This info originially came from
## file "HenleyLake_SampleInformation.xlsx").
samplingT <- read_csv(file="family_sampling_info.csv")

## To keep the chronological order (and avoid ABC order), we specify
## this as an ordered factor.
samplingT$collection <- ordered(samplingT$collection, levels=c("B", as.character(1:19)))
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
indivT$taxLvl[str_detect(indivT$taxon, "c__")] <- "class"
indivT$taxLvl[str_detect(indivT$taxon, "p__")] <- "phylum"
indivT$taxLvl[str_detect(indivT$taxon, "k__")] <- "kingdom"
indivT$taxLvl <- ordered(indivT$taxLvl, levels=c("class", "phylum", "kingdom"))


## Find percentage of counts can be classified down to the class
## level, across all types.  We see that about 93.3% of the counts can
## be attributed at the class level.
indivT %>% group_by(taxLvl) %>% summarize(sumCounts=sum(counts), percCounts=100*sum(counts)/sum(indivT$counts))

## Find percentage of counts which can be classified down to the
## class level, type by type.
indivT %>% group_by(type, taxLvl) %>% summarize(sumCountsByTypeTaxLvl=sum(counts)) %>% left_join(indivT %>% group_by(type) %>% summarize(sumCountsByType=sum(counts))) %>% mutate(perc=100*sumCountsByTypeTaxLvl/sumCountsByType)
## Percentages classified to class level: 94.4% (ribs), 93.5%
## (scapulae), 87.2% (water).


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
## ggsave("class_perc_classif_by_add_type.pdf", width=8.5, height=6, units="in")
## dev.off()
## ##################################################



## ##################################################
## From this point forward, we consider only class-level taxa for our
## model.  So, we remove all counts of non-class-level taxa.  This is
## consistent with the anlaysis we did for Shane's and Luisa's data.
indivT <- indivT %>% filter(taxLvl=="class")


## Make a new, more readable taxon column.
## Column names with open brackets (e.g., "c__[Saprospirae]") cause
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

## Count the number of unique class-level taxa observed for the
## various types (rib, scapula, water).
indivT %>% filter(taxLvl=="class") %>% filter(counts>0) %>% group_by(type) %>% distinct(taxon) %>% summarize(n=n())
##   type        n
##   <chr>   <int>
## 1 Rib       95
## 2 Scapula   124
## 3 Water     127


## For use in graphs and in calculating percentages later, we need
## total counts for each sample.
ctBySampleT <- indivT %>%
  group_by(sampleName) %>%
  summarize(totals=sum(counts))
## ##################################################



## ##################################################
## Some taxa don't occur frequently.  It's hard to make a hard cutoff
## for what constitutes "frequently".  There are 145 family-level taxa
## listed (some with 0 counts), and a lot of them appear in less than
## 1% of samples.

## I'm going to set the frequency cutoff at 1% (0.01).  This means
## that in order to meet the cutoff, a specific family-level taxa must
## make up at least 1% of the total family-level counts for that
## sample.
freqCutoff <- 0.01

## ## Get list of maximum taxa percentages for various types (rib,
## ## scapula, water) sorted in descending order:
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
write.csv(freqTaxaByTypeT, file="class_freq_taxa_by_type.csv", row.names=F)


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

write.csv(commontaxaT, file="classes_massaged.csv", row.names=FALSE)
## ##################################################
