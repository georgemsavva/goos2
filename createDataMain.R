rm(list=ls())

library(readxl)      # read excel data
library(data.table)  # data wrangling
library(ggplot2)     # plotting
library(gamlss)      # centile modelling
library(gamlss.cens) # modelling CRP with a censored distribution
library(survival)    # also to help with the censored data
library(splines)     # basis splines (not used?)
library(ggbeeswarm)  # beeswarm plots
library(survey)      # complex survey methods
library(gtsummary)   # nice tabulations
library(flextable)   # nice tabulations
library(officer)     # making powerpoints
library(patchwork)   # making multipanel graphs
library(childsds)    # calculate bmi_sds
require(emmeans)     # model mean estimation
library(modelsummary)# tabulate regression models
library(gt)          # tabulation
library(Hmisc)       # survey estimation
library(stringr)     # string functions
library(ggridges)    # ridge plots
library(hrbrthemes)  # plot theming
library(viridis)     # plot colours

## Loads functions to make the datasets we need
source("makeDatasets.R")

## Loads functions to produce (most of) the tables and figures for the paper
#source("list of tables.r")

## We need this for censored quantile modelling
# I don't think we use this ultimately.
#gen.cens("ST3","left")

## vitamin E is defined by (γ- and α-tocopherol)
## total carotenoids is lutein, α-cryptoxanthin, β-cryptoxanthin, 
## lycopene, α-carotene, and β-carotene

listOfOutcomes <- fread("listofoutcome.csv")
units <- c("CRP"="mg/L",
           "PFerritin"="μg/L", # The NDNS report is inconsistent but I think this is right
           "holoTC"="pmol/L",
           "SVitB12"="pmol/L",
           "RCFolate"="nmol/L",
           "v25OHD_Std"="nmol/L",
           "Retinol"="μmol/L",
           "TotalCarotenoids"="μmol/L",
           "vitE"="μmol/L",
           "alphatocopherol"="μmol/L",
           "gammatocopherol"="μmol/L",
           "pivka"="AU/L",
           "vitK"="µg/L",
           "VitC"="µmol/L",
           "VitB1"="ETKAC",
           "VitB2"="EGRAC",
           "VitB6PLP"="nmol/L",
           "VitB6PA"="nmol/L",
           "Zn"="µmol/L",
           "Se"="µmol/L",
           "gammaVsAlpha"="μmol/μmol",
           "gammaVsChol"="μmol/mmol",
           "alphaVsChol"="μmol/mmol",
           "ACRY"="μmol/L",
           "ACAR"="μmol/L",
           "BCRY"="μmol/L",
           "BCAR"="μmol/L",
           "LYCO"="μmol/L",
           "LEUT"="μmol/L",
           "FerritinCRP"="μg/mg",
           "SerumFolate"="μg/L",
           "RetinolVsChol"="μmol/mmol",
           "CarotVsChol"="μmol/mmol",
           "cholesterol"="mmol"
)

# THIS IS THE ORDER THAT EMMA WANTS THE OUTPUTS IN:
# Vitamin B1
# Vitamin B2
# Vitamin B6 PLP nmol/L
# Vitamin B6 PA nmol/L
# RC Folate nmol/L
# vitamin B12 pmol/L
# Holo-TC pmol/L
# CRP mg/L
# Ferritin μg/L
# Vitamin C µmol/L
# Retinol μmol/L
# 25-OH-D nmol/L
# vitamin E μmol/L
# a-Tocopherol μmol/L
# γ-Tocopherol:α-tocopherol
# γ-Tocopherol μmol/L
# Total Carotenoids μmol/L
# Selenium µmol/L
# Zinc µmol/L

## SHORT FOR FIGURES
order_of_outputs <- c("VitB1","VitB2","VitB6PLP","VitB6PA",
                      "RCFolate","SerumFolate",
                      "SVitB12","holoTC",
                      "CRP","PFerritin","FerritinCRP","VitC",
                      "Retinol","RetinolVsChol",
                      "v25OHD_Std",
                      "alphatocopherol","gammaVsAlpha","gammatocopherol",
                      "TotalCarotenoids","BCRY","BCAR","LYCO","LEUT","Se","Zn")

order_of_outputs <- c("VitB1","VitB2","VitB6PLP","VitB6PA",
                      "RCFolate","SerumFolate",
                      "SVitB12","holoTC",
                      "PFerritin","FerritinCRP","VitC",
                      "Retinol",
                      "v25OHD_Std",
                      "alphatocopherol","gammatocopherol",
                      "gammaVsChol","alphaVsChol",
                      "TotalCarotenoids","CarotVsChol",
                      "BCRY","BCAR","LYCO","LEUT","Se","Zn")


## LONG FOR TABLES / FOR COMPLETENESS
order_of_outputs2 <- c("VitB1","VitB2","VitB6PLP","VitB6PA",
                       "RCFolate","SerumFolate",
                       "SVitB12","holoTC",
                       "CRP","PFerritin","FerritinCRP","VitC",
                       "Retinol","RetinolVsChol",
                       "v25OHD_Std",
                       "alphatocopherol","gammaVsAlpha","gammatocopherol",
                       "gammaVsChol","alphaVsChol",
                       "TotalCarotenoids","CarotVsChol","BCRY","BCAR","LYCO","LEUT","Se","Zn")

### Remove two outcomes that are below the limit of detection almost everywhere.
listOfOutcomes <- listOfOutcomes[ varName %notin% c("ACRY","ACAR"), ]

#### Create all datasets
ndns <-         makeNDNSdata()
GOOSdata <-     makeGOOSdata()
NNUHdata <-     makeNNUHdata()
NNUHhealthy <-  makeNorfolkHealthyData2()

combinedData <- makeCombinedData(GOOSdata,ndns,NNUHdata,NNUHhealthy)

## Make design for SVY analysis
combined_design <- svydesign(ids=~PSU,
                             weights=~finalweight, 
                             strata=~strata,
                             data=combinedData
)
options(survey.lonely.psu="adjust")

#### Make the small datasets for GAMLSS model building (not used subsequently)
datasets <- makeIndividualOutcomeDatasets(listOfOutcomes,ndns,NNUHhealthy)

## Drop pivka and vitamin K at this stage, they are not used as the assay is thought to be unreliable
listOfOutcomes <- listOfOutcomes[varName %notin% c("pivka","vitK")]

#### Identify which outcomes are present in each dataset (helpful for loops and outputs)
presentInGOOS <- intersect(listOfOutcomes$varName,names(GOOSdata))
presentInNDNS <- intersect(listOfOutcomes$varName,names(ndns))
presentInNNUH <- intersect(listOfOutcomes$varName,names(NNUHdata))
presentInNNUHhealthy <- intersect(listOfOutcomes$varName,names(NNUHhealthy))

#### Make some nice labelling vectors
outcomevector <- listOfOutcomes$varName |> setNames(listOfOutcomes$labels)
labelsvector <- listOfOutcomes$labels |> setNames(listOfOutcomes$varName)
labelsvector2 <- paste(labelsvector,units[listOfOutcomes$varName]) |> 
  setNames(listOfOutcomes$varName)
labelsvector3 <- paste0(labelsvector," (",units[listOfOutcomes$varName],")") |> 
  setNames(listOfOutcomes$varName)

labelsvector4 <- c(labelsvector,
                   "anyDeficiency"="At least one deficiency",
                   "anyDeficiency2"="At least one deficiency (GOOS vars only)",
                   "multipleDeficiency"="More than one deficiency",
                   "multipleDeficiency2"="More than one deficiency (GOOS vars only)")

labelsvector3b <- ifelse(nchar(labelsvector3)>30, gsub(" ","\n",labelsvector3),labelsvector3) |> 
  gsub("l\nC","l C",x=_)


groupLabelsVector <- c("NDNS Normal-weight"="NDNS <85th centile",
                       "NDNS Over-weight"="NDNS 85th-95th centile",
                       "NDNS Obese"="NDNS >95th centile",
                       "GOOS"="GOOS")

#  scale_x_discrete(labels=c("NDNS <85th centile","NDNS 85th-95th centile","NDNS >95th centile","GOOS"))

