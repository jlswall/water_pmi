library("tidyverse")
library("readxl")
library("stringr")


## ##################################################
## Read information about collections, ADD, dates, etc. from
## spreadsheet in "Rice_RiversCenter_18S_SampleInformation.xlsx".

## Skip the header in row 1 and use our own column names.  We should
## end up with 152 rows of data.
fileNm <- "orig_data_files/Rice_RiversCenter_18S_SampleInformation.xlsx"
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

## All of these samples were taken at the same location ("James
## River") and using the same extraction method, so the "location" and
## "extractMethod" columns are not necessary.
samplingT <- samplingT %>% select(-location, -extractMethod)
## ##################################################


## ##################################################
## The "collection" column contains "Collection 1", ..., "Collection
## 19".  We save space by refering to these as "1", ..., "19"

tmp <- samplingT$collection
tmp <- str_remove(tmp, pattern="Collection ")
tmp <- str_remove(tmp, pattern="Colelction ")
samplingT$collection <- tmp

rm(tmp)
## ##################################################


## ##################################################
## Make capitalization consistent in type column.  In 2 cases, "Rib"
## was written as "RIb".
samplingT$type[samplingT$type=="RIb"] <- "Rib"
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
## counts would look like this table.
fullObsT <- as_tibble(expand.grid(type=c("Rib", "Scapula", "Water"), collection=1:19))
fullObsT$type <- as.character(fullObsT$type)
fullObsT$fullObsn <- 5
fullObsT$fullObsn[fullObsT$type=="Water"] <- 1
fullObsT$collection <- as.character(fullObsT$collection)


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
##    <chr>   <chr>      <int>
##  1 Rib     1              1
##  2 Scapula 1              4
##  3 Water   1             NA
##  4 Rib     2              1
##  5 Water   2             NA
##  6 Rib     3              4
##  7 Water   3             NA
##  8 Scapula 4              4
##  9 Water   4             NA
## 10 Rib     5              2
## 11 Rib     6              1
## 12 Water   6             NA
## 13 Rib     7             NA
## 14 Water   7             NA
## 15 Water   8             NA
## 16 Rib     9              2
## 17 Water   9             NA
## 18 Water   10            NA
## 19 Rib     11             4
## 20 Rib     12             3
## 21 Water   12            NA
## 22 Rib     13             4
## 23 Rib     14             4
## 24 Scapula 14             4
## 25 Water   14            NA
## 26 Rib     15             4
## 27 Water   15            NA
## 28 Water   16            NA
## 29 Rib     17            NA
## 30 Scapula 17             1
## 31 Rib     18             4
## 32 Scapula 18             4
## ##################################################


## ##################################################
## Write out the sample information into a CSV file.
write.csv(samplingT, file="sampling_info.csv", row.names=F)
## ##################################################
