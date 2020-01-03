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

## Switch date to "Date" class.  Default is POSIX class, which is for
## datetimes.  Here, we only have the dates, not a time of day.
samplingT$date <- as.Date(samplingT$date)

## All of these samples were taken at the same location ("James
## River") and using the same extraction method, so the "location" and
## "extractMethod" columns are not necessary.
samplingT <- samplingT %>% select(-location, -extractMethod)


## ##### WORKING HERE!


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
  filter( ((type=="Rib" | type=="Scapula") & nObs<5) | (type=="Water" & nObs<1) )
## This table shows which collections have missing scapula or rib
## samples.  Water samples are not missing on any of the collection
## days.
##    collection type     nObs
##    <ord>      <chr>   <int>
##  1 B          Scapula     2
##  2 3          Rib         4
##  3 4          Rib         4
##  4 6          Rib         4
##  5 9          Rib         3
##  6 10         Rib         4
##  7 11         Rib         3
##  8 12         Rib         4
##  9 13         Rib         4
## 10 14         Rib         3
## 11 15         Rib         3
## 12 16         Rib         4
## 13 17         Rib         3
## 14 19         Rib         4
## #######################

## Write out the sample information into a CSV file.
write.csv(samplingT, file="sampling_info.csv", row.names=F)

rm(fileNm, tmp)
## ##################################################

