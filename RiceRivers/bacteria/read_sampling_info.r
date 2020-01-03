library("tidyverse")
library("readxl")
library("stringr")


## ##################################################
## Read information about collections, sample types, dates, etc. from
## spreadsheet in "HenleyLake_SampleInformation.xlsx".

## Skip reading the header row, because we're providing our own column
## names.  Rows 274-275 in the spreadsheet contain NAs for control
## filter and swab, so we omit those rows.  We should end up with 272
## rows.
fileNm <- "orig_data_files/2019_RiceRiversCenter_SampleSheet.xlsx"
samplingT <- read_excel(fileNm, skip=1, n_max=272,
                        col_names=c("date", "degdays", "season",
                                    "location", "type", "collection",
                                    "extractMethod", "sampleName"))
rm(fileNm)

## Switch date to "Date" class.  Default is POSIX class, which is for
## datetimes.  Here, we only have the dates, not a time of day.
samplingT$date <- as.Date(samplingT$date)
## ##################################################


## ##################################################
## All of these samples were taken at the same location ("James
## River") and using the same extraction method, so the "location" and
## "extractMethod" columns are not necessary.
samplingT <- samplingT %>% select(-location, -extractMethod)


## The "collection" column contains one of "Baseline", "Collection 1",
## ..., "Collection 24", "Collection 24 Disturbed".  Note that the
## word "collection" is misspelled for collection 14 and that we may
## or may not have space between the word "collection" and the number
## (e.g. "Collection7", "Collection15").  We save space by refering to
## these as "B", "1", ..., "24", "24D".
tmp <- samplingT$collection
tmp <- str_remove(tmp, pattern="Collection ")
tmp <- str_remove(tmp, pattern="Collection")
tmp <- str_remove(tmp, pattern="Colelction ")
tmp[tmp=="Baseline"] <- "B"
tmp[tmp=="24 Disturbed"] <- "24Dist"
samplingT$collection <- tmp
## To keep the chronological order (and avoid ABC order), we specify
## this as an ordered factor.
samplingT$collection <- ordered(samplingT$collection, levels=c("B", as.character(1:24), "24Dist"))
rm(tmp)
## ##################################################


## ##################################################
## Ideally, we would have 5 rib samples, 5 scapula samples, 1 water
## sample, and 1 mud sample for each collection day.  That would mean
## 12 samples in total for each collection.  In reality, all samples
## are not available for every collection day.

## Count the number of such samples for each collection day.
samplingT %>% group_by(collection) %>% summarize(nObs=n())
## Count the number of such samples for each collection day and type.
samplingT %>% group_by(collection, type) %>% summarize(nObs=n())


## Identify collection days on which samples are missing.
samplingT %>%
  group_by(collection, type) %>%
  summarize(nObs=n()) %>%
  filter( ((type=="Rib" | type=="Scapula") & nObs<5) | (type=="Water" & nObs<1) | (type=="Mud" & nObs<1) )
## This table shows which collections have missing scapula or rib
## samples.  Water and mud samples are not missing on any of the collection
## days.
##   collection type     nObs
##   <ord>      <chr>   <int>
## 1 6          Rib         1
## 2 9          Rib         3
## 3 12         Rib         3
## 4 18         Scapula     4
## 5 20         Rib         2
## 6 21         Rib         3
## 7 22         Rib         4
## 8 23         Rib         4
## 9 24         Rib         1
## ##################################################


## ##################################################
## Write out the sample information into a CSV file.
write.csv(samplingT, file="sampling_info.csv", row.names=F)
## ##################################################

