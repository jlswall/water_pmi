library("tidyverse")
library("readxl")
library("stringr")


## ##################################################
## Read information about collections, sample types, dates, etc. from
## spreadsheet in "HenleyLake18S_SampleInformation.xlsx".

## The header takes up the first 3 lines.  We skip reading these
## provide our own column names.  We should end up with 272 rows.
fileNm <- "orig_data_files/HenleyLake18S_SampleInformation.xlsx"
samplingT <- read_excel(fileNm, skip=3, n_max=167,
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

## All of these samples were taken at the same location ("Henley
## Lake") and using the same extraction method, so the "location" and
## "extractMethod" columns are not necessary.
samplingT <- samplingT %>% select(-location, -extractMethod)
## ##################################################


## ##################################################
## The "collection" column contains "Collection 1", ..., "Collection
## 19".  We save space by refering to these as "1", ..., "19"

tmp <- samplingT$collection
tmp <- str_remove(tmp, pattern="Collection ")
samplingT$collection <- as.numeric(tmp)

rm(tmp)
## ##################################################


## ##################################################
## Ideally, we would have 5 rib samples, 5 scapula samples, and 1
## water sample, for each collection day.  That would mean 11 samples
## in total for each collection.  In reality, all samples are not
## available for every collection day.

## Count the number of such samples for each collection day.
samplingT %>% group_by(collection) %>% summarize(nObs=n())
## Count the number of such samples for each collection day and type.
samplingT %>% group_by(collection, type) %>% summarize(nObs=n())


## If we had the full number of samples for each collection, the
## counts would like this table.
fullObsT <- as_tibble(expand.grid(type=c("Rib", "Scapula", "Water"), collection=1:19))
fullObsT$type <- as.character(fullObsT$type)
fullObsT$fullObsn <- 5
fullObsT$fullObsn[fullObsT$type=="Water"] <- 1


## Identify count how many samples we have for rib, scapula, water for
## each collection day.
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

## WORKING HERE.  Need to remove NA for water on collection 16.  Also,
## need to replace table below with the table for this situation.


## This table shows which collections have missing scapula or rib
## samples.  Water and mud samples are not missing on any of the collection
## days.
##    collection type     nObs
##    <ord>      <chr>   <int>
##  1 5          Rib         2
##  2 6          Rib         1
##  3 9          Rib         3
##  4 12         Rib         3
##  5 18         Scapula     4
##  6 20         Rib         2
##  7 21         Rib         3
##  8 22         Rib         4
##  9 23         Rib         4
## 10 24         Rib         1
## ##################################################


## ##################################################
## Write out the sample information into a CSV file.
write.csv(samplingT, file="sampling_info.csv", row.names=F)
## ##################################################

