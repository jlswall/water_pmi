library("tidyverse")
library("readxl")
library("stringr")


# ##################################################
# Read in family-level taxonomic data from Excel sheet.
fileNm <- "orig_data_files/HVJ_HL_Taxonomy.xlsx"

# The first 454 rows contain the bone swabs.  The rows after than contain
# Claire's bone samples, which we dealt with in a previous study.  Row 1
# contains the header info.  We go up to column "BJ".  The columns after this
# contain summary tables.
rawAllT <- read_excel(path=fileNm, sheet="Family", range="A1:BJ454")
rm(fileNm)

# For all observations, the value of "taxlevel" is 5 (denoting family-level), so
# we remove this column.  We also remove "rankID", since we're not using it.
# The total column looks like it was supposed to get the sum for each bacterial
# family over the various collections, but the range of the sum appears to have
# been extended past the raw data into the summary table after column BJ.
rawAllT <- rawAllT %>% select(-taxlevel, -rankID, -total)

# Go from wide to long format.
rawIndivT <- rawAllT %>%
  gather(sampleName, counts, -taxon)
# ##################################################




# ##################################################
# Read information about collections, sample types, ADD, etc. from
# sheet "ADD_Collections".

fileNm <- "orig_data_files/HVJ_HL_Taxonomy.xlsx"
rawsamplingT <- read_excel(fileNm, sheet="ADD_Collections")

# We only need the actual ADD and the codes for the rib ane scapula samples
# (e.g., "R1B", "R1C2", "S1C17")
rawsamplingT <- rawsamplingT %>% select(`Actual ADD`, `Rib Swab Samples`,
              `Scapulae Swab Samples`)
# Separate the `Rib Swab Samples` column, which has 3 codes listed for each
# collection day, so that each code is listed in its own column.
rawsamplingT <- rawsamplingT %>% separate(`Rib Swab Samples`,
              into=c("R1", "R2", "R3")) 
# Do the same thing for the scapula samples.  However, the first and last
# collection days are missing one sample each.  This means we'll have two
# missing values and a warning message.
rawsamplingT <- rawsamplingT %>% separate(`Scapulae Swab Samples`,
              into=c("S1", "S2", "S3")) 

# Now, we go from wide to long format.  There are 2 rows with missing values,
# from the missing scapula swabs on the first and last days, and we drop those.
samplingT <- rawsamplingT %>% gather(RiOrSc, sampleName, -`Actual ADD`) %>%
              select(-RiOrSc) %>%
              arrange(`Actual ADD`) %>%
              rename(degdays=`Actual ADD`) %>%
              drop_na(sampleName)

# Make "type" variable to make it easy to tell which came from ribs and which
# came from scapulae.  This is also consistent with how I set up the dataset for
# the analysis with bones (rather than swabs).
samplingT$type <- ifelse(substring(samplingT$sampleName, first=1, last=1)=="S",
  "Scapula", "Rib") 

rm(rawsamplingT)
# ##################################################




# ##################################################
# Join the sampling information with ADD to the taxa counts.

indivT <- rawIndivT %>%
  inner_join(samplingT)
rm(rawIndivT, rawAllT, samplingT)
# ##################################################




# ##################################################
# Consider level of classification for each taxa.  (Some could not be
# classified down to the family level, but maybe the order could be
# determined.)

# Make new variable to indicate the most precise taxon which could be
# identified.
indivT$taxLvl <- ""
indivT$taxLvl[str_detect(indivT$taxon, "f__")] <- "family"
indivT$taxLvl[str_detect(indivT$taxon, "o__")] <- "order"
indivT$taxLvl[str_detect(indivT$taxon, "c__")] <- "class"
indivT$taxLvl[str_detect(indivT$taxon, "p__")] <- "phylum"
indivT$taxLvl[str_detect(indivT$taxon, "k__")] <- "kingdom"
indivT$taxLvl <- ordered(indivT$taxLvl,
      levels=c("family", "order", "class", "phylum", "kingdom"))


# Find percentage of counts can be classified down to the family
# level, across all types.  We see that about 85.1% of the counts can
# be attributed at the family level.
indivT %>%
  group_by(taxLvl) %>%
  summarize(sumCounts=sum(counts),
            percCounts=100*sum(counts)/sum(indivT$counts))

# Find percentage of counts which can be classified down to the
# family level, type by type.
indivT %>%
  group_by(taxLvl, type) %>%
  summarize(sumCountsByTypeTaxLvl=sum(counts)) %>%
  left_join(indivT %>%
            group_by(type) %>%
            summarize(sumCountsByType=sum(counts))
            ) %>%
  mutate(perc=100*sumCountsByTypeTaxLvl/sumCountsByType)
# Percentages classified to family level: 86.8% (ribs), 83.2%
# (scapulae)


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
ggsave("swabs_family_perc_classif_by_add_type.pdf", width=8.5, height=6,
  units="in")
dev.off()
## ##################################################




# ##################################################
# From this point forward, we consider only family-level taxa for our
# model.  So, we remove all counts of non-family-level taxa.  This is
# consistent with the anlaysis we did for Shane's and Luisa's data.
indivT <- indivT %>% filter(taxLvl=="family")


# Make a new, more readable taxon column.
# Column names with open brackets (e.g. "f__[Tissierellaceae]") cause
# problems for functions expecting traditional data frame column
# names (like randomForest function).
indivT$newtaxon <- gsub(indivT$taxon, pattern="\\[", replacement="")
indivT$newtaxon <- gsub(indivT$newtaxon, pattern="]", replacement="")
# Column names with dashes can likewise be a problem, so I replace
# dashes with underscores.
indivT$newtaxon <- gsub(indivT$newtaxon, pattern="-", replacement="_")
# Remove the column of the original, less-R-friendly taxon names from
# the tibble to avoid confusion.  Rename new column.
indivT <- indivT %>% select(-taxon)
indivT <- rename(indivT, taxon = newtaxon)

# Count the number of unique family-level taxa observed for the
# various types (rib, scapula, water).
indivT %>%
  filter(taxLvl=="family") %>%
  filter(counts>0) %>%
  group_by(type) %>%
  distinct(taxon) %>%
  summarize(n=n())
#  type        n
#  <chr>   <int>
# Rib       172
# Scapula   191


# For use in graphs and in calculating percentages later, we need
# total counts for each sample.
ctBySampleT <- indivT %>%
  group_by(sampleName) %>%
  summarize(totals=sum(counts))
# ##################################################




# ##################################################
# Some taxa don't occur frequently.  It's hard to make a hard cutoff
# for what constitutes "frequently".  There are 224 family-level taxa
# represented, and a lot of them appear in less than 1% of
# samples.

# I'm going to set the frequency cutoff at 1% (0.01).  This means
# that in order to meet the cutoff, a specific family-level taxa must
# make up at least 1% of the total family-level counts for that
# sample.
freqCutoff <- 0.01

# Get list of maximum taxa percentages for various types (rib,
# scapula, water) sorted in descending order:
# indivT %>%
#   left_join(ctBySampleT) %>%
#   mutate(fracBySubjDay = counts/totals) %>%
#   group_by(type, taxon) %>%
#   summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
#   ## filter(maxFracBySubjDay >= freqCutoff) %>%
#   arrange(type, desc(maxFracBySubjDay)) %>%
#   print(n = Inf)


# Count how many times each taxa beats the cutoff (freqCutoff) within
# each type.  Then, keep only those which meet cutoff for more than 1
# sample.  Save those taxa names for each type.
freqTaxaByTypeT <- indivT %>%
  left_join(ctBySampleT) %>%
  mutate(fracBySubjDay = counts/totals,
        isExceed=(fracBySubjDay>=freqCutoff)
        ) %>%
  group_by(type, taxon) %>%
  summarize(numExceed = sum(isExceed)) %>%
  filter(numExceed > 1) %>%
  ## arrange(type, desc(numExceed)) %>%
  select(type, taxon)

# For each type, how many taxa meet the cutoff?
freqTaxaByTypeT %>% group_by(type) %>% summarize(n=n())
#   type        n
#   <chr>   <int>
# 1 Rib        21
# 2 Scapula    37
# Save list of frequent taxa by type to a CSV file.
write.csv(freqTaxaByTypeT, file="family_freq_taxa_by_type.csv", row.names=F)


# Rename taxa that occur less than the frequency cutoff allows as
# "rare".  Then, sum all these "rare" taxa into one row for each
# sample.
commontaxaT <- indivT
# Find rows which contain taxa which do not meet the frequency cutoff
# for that type.
indicRare <- !( with(commontaxaT, paste(taxon, type, sep=":")) %in%
  with(freqTaxaByTypeT, paste(taxon, type, sep=":")) )
# For these "rare" taxa, we change taxon name to rare.
commontaxaT[indicRare, "taxon"] <- "Rare"
# Now, we add up the rare counts so that "Rare" appears only one for
# each sample.
commontaxaT <- commontaxaT %>%
  group_by(sampleName, taxon, type, degdays) %>%
  summarize(counts = sum(counts))


# Use the table of total counts by subj/day to find the fraction
# represented by each taxa for each subj/day.
commontaxaT <- commontaxaT %>%
  left_join(ctBySampleT) %>%
  mutate(fracBySample=counts/totals) %>%
  select(-totals)


# Check that the fractions add up to 1, appropriately.
unique(
    commontaxaT %>%
    group_by(sampleName) %>%
    summarize(sumFracBySample = sum(fracBySample)) %>%
    pull(sumFracBySample)
)


# Remove the list of taxa names that satisfied the frequency cutoff.
rm(indicRare, freqCutoff, freqTaxaByTypeT)
# ##################################################



# ##################################################
# Save the tibble to a file for use in separate code
# for graphing and analysis.

write.csv(commontaxaT, file="families_massaged.csv", row.names=FALSE)
# ##################################################
