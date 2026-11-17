
### This file searches through the NDNS files and pulls out the individual blood weights and strata for each participant.

### Note 16 March 2026 
### Why no PSU/cluster information?
### 

library(data.table)

getIMD <- function(){

  filenames <- dir("c:/data/UKDA-6533-tab/tab/",pattern = "indiv")
  
  ## This is a check to test that the filenames appear as I expect them.
  stopifnot(identical(filenames,c("ndns_rp_yr1-4a_indiv_uk.tab", 
                                  "ndns_rp_yr5-6a_indiv.tab",
                                  "ndns_rp_yr7-8a_indiv.tab",
                                  "ndns_rp_yr9-11a_indiv_20211020.tab",
                                  "ndns_yr1-3indiva_weights.tab",
                                  "ndns_yr12-15a_indiv.tab")))
  
  tabs <- lapply(filenames, \(f){
    fread(file = sprintf("c:/data/UKDA-6533-tab/tab/%s",f))
  }) |> setNames(filenames)
  
  
IMD_dat <-  lapply(tabs[c(1:4)] , \(t) t[ , IMD5 := fcoalesce(EIMD_2015_quintile,
                                             WIMD_2014_quintile,
                                             SIMD_2016_quintile,
                                             NIMD_2017_quintile) ][, .(seriali,IMD5) ]  ) |> rbindlist() |> rbind(
  
  tabs[[6]][ , IMD5 := fcoalesce(qeimd,qwimd,qnimd,qsimd)][, .(seriali=SerialI,IMD5) ] )

fwrite(x = IMD_dat,"ndns_imd_dat.csv")  
  
}
getIMD()

tabs <- list()
dir("c:/data/UKDA-6533-tab/tab/",pattern="12")

filenames <- dir("c:/data/UKDA-6533-tab/tab/",pattern = "indiv")

## This is a check to test that the filenames appear as I expect them.
stopifnot(identical(filenames,c("ndns_rp_yr1-4a_indiv_uk.tab", 
                                "ndns_rp_yr5-6a_indiv.tab",
                                "ndns_rp_yr7-8a_indiv.tab",
                                "ndns_rp_yr9-11a_indiv_20211020.tab",
                                "ndns_yr1-3indiva_weights.tab",
                                "ndns_yr12-15a_indiv.tab")))

## Load the datasets
tabs <- lapply(filenames, \(f){
  fread(file = sprintf("c:/data/UKDA-6533-tab/tab/%s",f))
}) |> setNames(filenames)

## Checking what information we have on participant age
#tabs[[1]][ , .SD , .SDcols=grep(x = names(tabs[[1]]),"age") , ]
#tabs[[1]] |> names() |> grep(pattern="dob")

## Identify which variables we need to pull
weightvars <- lapply(names(tabs), \(t) {  grep(pattern="wtb|strata2|astrata2|Area|area",names(tabs[[t]]),value=TRUE) }) |> setNames(filenames)

## Now stack the files with only the variables we want
weightsDat <- rbindlist(list(
tabs[[1]][ , .(seriali,age,.SD) , .SDcols = weightvars[[1]]],
tabs[[2]][ , .(seriali,age,.SD) , .SDcols = weightvars[[2]]],
tabs[[3]][ , .(seriali,age,.SD) , .SDcols = weightvars[[3]]],
tabs[[4]][ , .(seriali,age=AgeR,.SD) , .SDcols = weightvars[[4]]],
tabs[[6]][ , .(seriali=SerialI,age=AgeR,.SD) , .SDcols = weightvars[[6]]]),fill=TRUE,idcol="file")

## Check there's only one weight per person
stopifnot(all(
weightsDat[ , as.numeric(!is.na(.SD.wtb_UKY1234))+
              as.numeric(!is.na(.SD.wtb_Y56))+
              as.numeric(!is.na(.SD.wtb_Y78))+
              as.numeric(!is.na(.SD.wtb_Y911))+
              as.numeric(!is.na(.SD.wtb_Y1215))]<=1
))

## Check there's only one stratum per person
stopifnot(all(
weightsDat[ , as.numeric(!is.na(.SD.astrata2))+
              as.numeric(!is.na(.SD.strata2))] <=1
))

## Make weight and strata variables
## (use strata2 now, previously I used strata1)
weightsDat[ !is.na(.SD.wtb_UKY1234), weight := .SD.wtb_UKY1234 ]
weightsDat[ !is.na(.SD.wtb_Y56), weight := .SD.wtb_Y56 ]
weightsDat[ !is.na(.SD.wtb_Y78), weight := .SD.wtb_Y78 ]
weightsDat[ !is.na(.SD.wtb_Y911), weight := .SD.wtb_Y911 ]
weightsDat[ !is.na(.SD.wtb_Y1215), weight := .SD.wtb_Y1215 ]

weightsDat[ !is.na(.SD.strata2) , strata:=.SD.strata2]
weightsDat[ !is.na(.SD.astrata2) , strata:=.SD.astrata2]

## Make Area:
weightsDat[ !is.na(.SD.Area) , PSU := .SD.Area]
weightsDat[ !is.na(.SD.area) , PSU := .SD.area]


## A bit of checking
weightsDat[weight>0,weight] |> hist(breaks = 100)
weightsDat[weight>0,] |> nrow()
weightsDat[weight>0,mean(weight),by=file]

## Remove the over 18s and those with no blood weight
weightsDat <- weightsDat[ age <= 18 & !is.na(weight) & weight>0 , ]

## Now merge the weights as per NDNS documentation
## This code block explicitly follows the instructions for re-weighting in the NDNS documentation here:
## https://doc.ukdataservice.ac.uk/doc/6533/mrdoc/pdf/6533_ndns_y12_15_user_guide.pdf (page 22/23)
weightsDat[ ,   sumweight := sum(weight), by=file]
weightsDat[ , .(sumweight = sum(weight)), by=file]
weightsDat[ , totalsumweight := sum(weight) ]
weightsDat[ , .(totalsumweight = sum(weight)) ]
fileweights <- c(4,2,2,3,4)
weightsDat[ , finalweight := weight * (totalsumweight / sumweight) * (fileweights[file] / sum(fileweights))  ]
weightsDat[ , finalweight := finalweight / mean(finalweight)  ]


## Write out the weights/strata file to load later
weightsDat[ , .(seriali , finalweight, strata,PSU)] |> fwrite(file="childBloodWeights.csv")

#  Notes from the meeting with the Cambridge Epidemiology Unit team
## How do GOOS cohort compare to the obese population
#  Use ethnic groups: Indian (Asian?), White, Black
#  Recode NDNS similarly
## How does GOOS exclude mixed population
###### Look over time, and by region?
###### Reassure yourself that it doesn't make too much difference,
###### Also ethnicity and sex.
###### Make a template for the other biomarkers?
#### Percentile chart can be done using GAMLSS
##  I think I'd want to see the distributions in both cohorts.

## PSUs are unique to files, so this is OK.
weightsDat[, length(unique(file)), by=PSU]$V1 |> table()

