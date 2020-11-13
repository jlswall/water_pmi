library("tidyverse")
library("readxl")
library("stringr")


## ##################################################
## Read information about collections, sample types, dates, etc. from
## spreadsheet in "HenleyLake18S_SampleInformation.xlsx".

## The header takes up the first 3 lines.  We skip reading these
## provide our own column names.  We should end up with 272 rows.
fileNm <- "orig_data_files/HL_RR_Combined_Data_Information.xlsx"
samplingT <- read_excel(fileNm, skip=1,
                        col_names=c("date", "degdays", "season",
                                    "location", "type", "collection",
                                    "extractMethod", "sampleName"))
rm(fileNm)

## Switch date to "Date" class.  Default is POSIX class, which is for
## datetimes.  Here, we only have the dates, not a time of day.
samplingT$date <- as.Date(samplingT$date)
## ##################################################


## ##################################################
## Remove unnecessary columns.

## All of these samples were taken using the same extraction method, so the
## "extractMethod" column is not necessary.
samplingT <- samplingT %>% select(-extractMethod)
## ##################################################


## ##################################################
## The "collection" column contains "Collection 1", ..., "Collection
## 19".  We save space by refering to these as "1", ..., "19"

tmp <- samplingT$collection
tmp <- str_remove(tmp, pattern="Collection ")
## We have some rows where collection is misspelled.
tmp <- str_remove(tmp, pattern="Colelction ")
samplingT$collection <- tmp

rm(tmp)
## ##################################################


## ##################################################
## The "sampleName" column contains a few entries which include extra
## space.  We remove those spaces.

tmp <- samplingT$sampleName
tmp <- str_remove_all(tmp, pattern=" ")
samplingT$sampleName <- tmp

rm(tmp)
## ##################################################


## ##################################################
## Ideally, we would have 5 rib samples, 5 scapula samples, and 1
## water sample, for each collection day at each location.  That would mean 11
## samples in total for each collection.  In reality, all samples are not
## available for every collection day.

## If we had the full number of samples for each collection, the
## counts would look like this table.
## WORKING HERE!
fullObsT <- as_tibble(expand.grid(location=c("Henley Lake", "James River"),
              type=c("Rib", "Scapula", "Water"),
              collection=c("Baseline", 1:19))
              )
fullObsT$type <- as.character(fullObsT$type)
fullObsT$collection <- as.character(fullObsT$collection)
fullObsT$fullObsn <- 5
fullObsT$fullObsn[fullObsT$type=="Water"] <- 1


## Identify count how many samples we have for rib, scapula, water for
## each collection day.  In the next step, we join this with fullObsT
## and compare to see which samples are missing.
samplingT %>%
  group_by(collection, type) %>%
  summarize(nObs=n())


fullObsT %>%
  left_join(samplingT %>%
            group_by(collection, type) %>%
            summarize(nObs=n())
            ) %>%
  filter((nObs < fullObsn) | is.na(nObs)) %>%
  select(-fullObsn) %>%
  print(n=Inf)
## This table shows which collections have missing scapula or rib
## samples. "NA" means that there are 0 observations for that combination.
##    type    collection  nObs
##    <chr>        <dbl> <int>
##  1 Rib              1     3
##  2 Scapula          1     4
##  3 Rib              3     4
##  4 Rib              4     3
##  5 Scapula          4     4
##  6 Rib              5     4
##  7 Scapula          5     4
##  8 Rib              6     3
##  9 Scapula          6     4
## 10 Rib              7     4
## 11 Rib              8     4
## 12 Scapula          8     3
## 13 Rib              9     2
## 14 Rib             10     4
## 15 Scapula         10     4
## 16 Rib             11     2
## 17 Scapula         11     4
## 18 Rib             12     4
## 19 Scapula         12     3
## 20 Rib             13     3
## 21 Scapula         13     4
## 22 Rib             14     3
## 23 Rib             15     3
## 24 Rib             16     4
## 25 Water           16    NA
## 26 Rib             17     3
## 27 Rib             18     4
## 28 Rib             19     3
## ##################################################


## ##################################################
## Write out the sample information into a CSV file.
write.csv(samplingT, file="sampling_info.csv", row.names=F)
## ##################################################
