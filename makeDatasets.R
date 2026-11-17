## makeDatasets

stars<-function(p.value) ifelse(p.value<0.05,
                                ifelse(p.value<0.01,
                                       ifelse(p.value<0.001,"***","**"),
                                       "*"),
                                "")
which.nonnum <- function(x) {
  badNum <- is.na(suppressWarnings(as.numeric(as.character(x))))
  x[which(badNum & !is.na(x))]
}

formatp <- function(p) ifelse(p<0.001,"p<0.001***",sprintf("p=%0.3f%s",p,stars(p)))


makeNDNSdata <- function(){
  
##### THERE ARE NO NON-NUMERICS IN THE NUMERICAL DATA HERE.
  
## Load weights information
weightsDat <- fread("dat/childBloodWeights.csv")

## Load NDNS dat
ndns_raw <- read_excel("dat/NDNS_Child_micronutrientdata.xlsx")
setDT(ndns_raw)

## Load energy intake data
energyKcal= fread("dat/energyKcal.csv")
energyKcal=energyKcal[ seriali %in% ndns_raw$seriali]

ndns_raw <- ndns_raw[  energyKcal , on="seriali"]

## Create vitamin E variable
ndns_raw[ , vitE := ATC + GTC...42] # There are two GTC variables in the sheet but they are identical

ndns_raw[ , alphatocopherol := ATC] 
ndns_raw[ , gammatocopherol := `GTC...42`]
ndns_raw[ , gammaVsAlpha := gammatocopherol / alphatocopherol] 
ndns_raw[ , gammaVsChol := G_chol ]
ndns_raw[ , alphaVsChol := A_chol ]
ndns_raw[ , RetinolVsChol := Retinol / Chol]
ndns_raw[ , TotalCarotenoids:=TotalCarotenoids...43] # two Total Carotenoids variables but they are identical
ndns_raw[ , CarotVsChol := TotalCarotenoids / Chol]

ndns_raw[ , LEUT := TotalCarotenoids - (ACRY+BCRY+ACAR+BCAR+LYCO) ]

ndns_raw[ , FerritinCRP := PFerritin / CRP]

## Check the weights are all present
setdiff( ndns_raw$seriali , weightsDat$seriali )

## Cleaning the ethnicity codes
ndns_raw[ ethgr5 == "Asian or Asian British" , ethgr5 := "Asian"]
ndns_raw[ ethgr5 == "Asian or asian British" , ethgr5 := "Asian"]
ndns_raw[ ethgr5 == "Black or Black British" , ethgr5 := "Black"]
ndns_raw[ ethgr5 == "Mixed ethnic group" ,     ethgr5 := "Mixed"]
ndns_raw[ ethgr5 == "Any other group" ,     ethgr5 := "Other"]

# fix age
ndns_raw[ , bmi_sds := childsds::sds(bmival,age,Sex,
                                     ref = childsds::ukwho.ref,
                                     male = "Male",
                                     female="Female",
                                     item = "bmi")]

ndns_raw |> ggplot(aes(x=bmicat218,y=bmi_sds)) + geom_boxplot() + 
  coord_flip()
ndns_raw[ , cholesterol := Chol]
sdcols = c("seriali", "SurveyYear", "Sex", "age", "ethgr5", 
           "vitE_intake"="VitaminEmg",
           "vitA_intake"="VitaminAretinolequivalentsµg",
           "vitB12_intake"="VitaminB12µg",
           "TotalCarot_intake"="Totalcaroteneµg",
           "eqv3", "region", "bmival","bmi_sds", "RCFolate","Energykcal",
           "bmicat218","cholesterol",
           intersect(listOfOutcomes$varName,names(ndns_raw)))

ndns <- ndns_raw[ ,.SD, .SDcols = sdcols]
ndns[ , bmi:=bmival]
#ndns_raw[ , cor(VitaminDµg,v25OHD_Std,use="pairwise")]

## Check whether any NDNS cases are without weights
#### We actually know there are four rows without blood data..
#stopifnot(length(setdiff( ndns$seriali , weightsDat$seriali ))==0)

## Load the weights into the nutrient data
ndns <- weightsDat[ ndns , on="seriali"]
ndns <- ndns[!is.na(finalweight)]

ndns
}


makeIndividualOutcomeDatasets <- function(listOfOutcomes,ndns,NNUHhealthy){
  
  datasets <- list()
  for(i in c(1:nrow(listOfOutcomes))){
    if(listOfOutcomes$varName[i] %notin% c("pivka","vitK")){
      ndns[ , target := get(listOfOutcomes$varName[i])]
      ndns[ , target := ifelse(target==0,
                               min(target[target!=0],na.rm=TRUE),
                               target)]
      ndns[ , logtarget := log(target)]
      ## Create a 'healthy' dataset to build the model
      ndnsHealthy <- ndns[ bmicat218 == "Normal-weight", 
                           .(age,target,logtarget,finalweight,Sex,ethgr5,region,SurveyYear)]
      
      datasets[[i]] <- copy(na.omit(ndnsHealthy))
    }
    
  }
  datasets <- setNames(datasets,listOfOutcomes$varName )
  ### CRP is highly censored so use a survival analysis
  datasets[["CRP"]][ , event := as.numeric(target>1)]
  datasets[["CRP"]][ , targetCen := ifelse(event==1,target,1)]
  datasets[["CRP"]][ , logtargetCen := log(targetCen)]
  
  ### These two are only in NNUH, so we don't have an NDNS dataset to work with.
#  datasets[[9]] <- NNUHhealthy[ !is.na(pivka), .(age,target=pivka,logtarget=log(pivka),
#                                    finalweight=1,Sex=sex,ethgr5,region="Norfolk",SurveyYear=14)]
#  datasets[[10]] <- NNUHhealthy[ !is.na(vitK), .(age,target=vitK,logtarget=log(vitK),
#                                                 finalweight=1,Sex=sex,ethgr5,region="Norfolk",SurveyYear=14)]
  #datasets[[""]] <- "empty"  # Don't look at PIVKA/Vitamin K any more.
  #datasets[[10]] <- "empty"

  datasets
  
}


## File including the functions to load and encode the datasets.

## If less than 5% censored low then I'm not going to worry about it.

makeGOOSdata <- function(){
  ## Load the GOOS data
  GOOSvitaminK <- read_excel(path="dat/VitK1GOOSresults.xlsx") |> 
    setNames(c("subjectID","vitK"))
  setDT(GOOSvitaminK)
  GOOSvitaminK[ , vitKCensored := vitK %like% "<" | vitK == "ND"]
  GOOSvitaminK[ , vitK := as.numeric(gsub("<","",vitK))]
  
  ## I'm replacing the value of vitamin K with 0.1 if its ND, and calling it censored.  The interpretation in the modelling is "less than 0.1".
  ## Since this is less than the deficiency threshold anyway I don't think it matters.
  GOOSvitaminK[ is.na(vitKCensored), vitK == 0.1]

  GOOSimd <- read_excel(path="dat/GOOSimd.xlsx",sheet=1,range="A1:D237")
  setDT(GOOSimd)
  GOOSimd[ , imd2 := as.numeric(`Index of Multiple Deprivation (IMD) Decile`)]
  GOOSimd[ , inc2 := as.numeric(`Income Decile`)]
  GOOSimd[ , subjectID := `Subject ID`]
  
  GOOSbiometric <- read_excel(path = "dat/GOOS biometric data.xlsx") |> as.data.table()
  GOOSbiometric <- GOOSbiometric[ , .(subjectID=`Subject ID`, bmi=`BMI (kg/m2)`, bmi_sds=`BMI SDS`)]
  GOOSdata <- read_excel("dat/GOOS data summay-original.xlsx" , 
                         sheet="Sheet1") |> 
    setNames(c("subjectID","sex","ethgr5","age","v25OHD_Std","Retinol",
               "gammatocopherol","alphatocopherol","lutein",
               "alphacryptoxanthin",
               "betacryptoxanthin","lycopene",
               "alphacarotene","betacarotene",
               "TotalCarotenoids","PFerritin","CRP",
               "SVitB12","holoTC"))
  
  setDT(GOOSdata)
  
  GOOSdata[ , gammatocopherol := as.numeric(gsub("<","",gammatocopherol)) ] # five non-numerics
  GOOSdata[ , alphatocopherol := as.numeric(alphatocopherol) ]
  GOOSdata[ , vitE := alphatocopherol + gammatocopherol]
  GOOSdata[ , gammaVsAlpha := gammatocopherol/alphatocopherol]
  
  
  
  GOOSdata[ , CRP := as.numeric(gsub("<","",CRP))]  #  one non-numeric
  GOOSdata[ , SVitB12 := as.numeric(SVitB12)]
  GOOSdata[ , PFerritin := as.numeric(PFerritin)]
  GOOSdata[ , v25OHD_Std := as.numeric(v25OHD_Std)]
  GOOSdata[ , holoTC := as.numeric(gsub(">","",holoTC))]  #  (23 of these are 'high' since we care mostly about deficiency happy to ignore)
  GOOSdata[ , Retinol := as.numeric(Retinol)]
  GOOSdata[ , TotalCarotenoids := as.numeric(`TotalCarotenoids`)]
  
  GOOSdata[ , ACRY := as.numeric(gsub("<","",alphacryptoxanthin))]
  GOOSdata[ , BCRY := as.numeric(gsub("<","",betacryptoxanthin))]
  GOOSdata[ , ACAR := as.numeric(gsub("<","",alphacarotene))]
  GOOSdata[ , BCAR := as.numeric(gsub("<","",betacarotene))]
  GOOSdata[ , LEUT := as.numeric(gsub("<","",lutein))]
  GOOSdata[ , LYCO := as.numeric(gsub("<","",lycopene))]
  
  GOOSdata[ , FerritinCRP := PFerritin / CRP]
  
  GOOSdata <- GOOSbiometric[ GOOSdata , on="subjectID" ]
  GOOSdata <- GOOSvitaminK[ GOOSdata , on="subjectID"]
  GOOSdata <- GOOSimd[ , .(subjectID,imd2,inc2)][GOOSdata, on="subjectID"]
  
  ## Clean ethnicity in GOOS
  GOOSdata[ ethgr5 == "Caucasian" , ethgr5 := "White"]
  GOOSdata[ ethgr5 == "Indo-origin" , ethgr5 := "Asian"]
  GOOSdata[ ethgr5 == "Afro-caribbean" , ethgr5 := "Black"]
  
  GOOSdata$ethgr5 |> table()
  
  GOOSdata
}


makeNNUHdata <- function(){

## Load the NNUH data

NNUHdata_raw <- read_excel(path = "dat/Copy of NNUH Obese CYP Anonomised micronutrient data.xlsx",range="A1:AB170")

## Emma sent me this file later, I'm just extracting the extra columns from it.
NNUHdata_CholCRPDep <- read_excel(path = "dat/Copy of Anon micronutrient data with CRP and Chol for George.xlsx",range="A1:D170") |> as.data.table()


setDT(NNUHdata_raw)

NNUHdata_raw[ `Vit A`=="1.41 (06/07/22)"  , `Vit A`:="1.41"]
NNUHdata_raw[ `Vit E (13-24) umol/L)`=="22.9 (06/07/2022)", `Vit E (13-24) umol/L)`:="22.9"]

## All non-numerics dealt with apart from folate.  holoTC and K are just replaced with the top value
NNUHdata <-  NNUHdata_raw[ , .(subjectID=paste0("NNUH_",`...1`),
                               age=`Age at blood test`,
                               sex=Sex,
                               Ethnicity,
                               bmi=as.numeric(BMI),
                               bmi_sds=`BMI SDS`,
                               vitK=`Vit K1 (0.15-1.55)`,
                               holoTC=as.numeric(gsub(">","",`Holotranscobalamin (>34)`)),
                               # Holo TC TOP VALUES CENSORED
                               vitE=as.numeric(`Vit E (13-24) umol/L)`),
                               SerumFolate=as.numeric(gsub(x=`Folate (3.1--20)`,"<","")),
                               Retinol=as.numeric(`Vit A`),
                               PFerritin=as.numeric(Ferritin),
                               Zn=as.numeric(`Zinc (9.8-19) umol/L`), ## Happy that all the missing Zn are really missing
                               pivka=as.numeric(`PIVKA II (17.36-50.9)`),
                               v25OHD_Std=as.numeric(`25 OHD Vit D`))][!is.na(subjectID)]

NNUHdata <- NNUHdata_CholCRPDep[ , 
                                 .(subjectID=paste0("NNUH_",`...1`),
                                   CRP=as.numeric(gsub("<","",CRP)),
                                   cholesterol,imd10=`Deprivation Decile (-1 is unknown)`)] |> 
  _[NNUHdata,on="subjectID"]

#NNUHdata[ , alphaVsChol := vitE/ ]
NNUHdata[ , FerritinCRP := PFerritin / CRP]
NNUHdata[ , RetinolVsChol := Retinol / cholesterol]

NNUHdata[ , vitKCensored := vitK %like% "<"]
NNUHdata[ , vitK := as.numeric(gsub("<","",vitK))]

NNUHdata[,ethgr5:="Unknown"]

NNUHdata[Ethnicity=="White-British",ethgr5:="White"]
NNUHdata[Ethnicity=="white",ethgr5:="White"]
NNUHdata[Ethnicity=="White-British - British",ethgr5:="White"]
NNUHdata[Ethnicity=="white british",ethgr5:="White"]
NNUHdata[Ethnicity=="Lithuanian",ethgr5:="White"]
NNUHdata[Ethnicity=="white other",ethgr5:="White"]
NNUHdata[Ethnicity=="White other",ethgr5:="White"]
NNUHdata[Ethnicity=="?portugese",ethgr5:="White"]
NNUHdata[Ethnicity=="Turkish",ethgr5:="White"]

NNUHdata[Ethnicity=="asian",ethgr5:="Asian"]
NNUHdata[Ethnicity=="Asian",ethgr5:="Asian"]
NNUHdata[Ethnicity=="Asian -Filipino",ethgr5:="Asian"]
NNUHdata[Ethnicity=="asian",ethgr5:="Asian"]

NNUHdata[Ethnicity=="Mixed white and blcak caribbean",ethgr5:="Mixed"]
NNUHdata[Ethnicity=="Mixed white and black african",ethgr5:="Mixed"]
NNUHdata[Ethnicity=="Portugese + Black",ethgr5:="Mixed"]

NNUHdata[Ethnicity=="African",ethgr5:="Black"]
NNUHdata[Ethnicity=="Black",ethgr5:="Black"]
#NNUHdata[, table(Ethnicity,ethgr5,useNA = "ifany")]

NNUHdata[ imd10<0 , imd10 := NA]

imd2 <- read_excel("dat/Deprivation scores CEW BMI 3SDS.xlsx")

NNUHdata$imd2 = imd2$`Index of Multiple Deprivation (IMD) Decile`

NNUHdata

}


makeNorfolkHealthyData2 <- function(){
  logNumber <- function(s,x) print(sprintf("NNUH Vit D: %s %d",s,x))
  healthyDat <- read_excel(path = "dat/Vitamin K and Vitamin D data combined.xlsx")
  logNumber("Load dataset",nrow(healthyDat))
  healthyDatIMD <- read_excel(path = "dat/Multiple depiravation data anon Vitamin D NNUH cohort_also includes missing patients.xlsx")
  setDT(healthyDat)
  setDT(healthyDatIMD)
  
  healthyDatIMD <- 
    healthyDatIMD[ , .(subjectID=`0`,
                       imd2=`Index of Multiple Deprivation Decile`,
                       bmi_sds2=as.numeric(`BMI SDS`),
                       age2=Age)]
  
  #hist(healthyDatIMD$age2,breaks=100)
  
  ## Non-numerics dealt with.  Vitmain K censoring indicator, none for vitmain D or PIVKA.
  healthyDat <- healthyDat[ , .(subjectID=`0`, 
                                #pivka=as.numeric(`PIVKA-II`), 
                                #vitK = `Serum\r\nVitamin K`,
                                #svk=`Serum\r\nVitamin K`,
                                v25OHD_Std=
                                  as.numeric(`Total 25(OH)D nmol/L`),
                                Ethnicity,
                                SexNumeric=`Gender \r\n(1=male)`,
                                age=Age,
                                #height=`Height (m)`,
                                #weight=`Weight (kg)`, 
                                bmi=BMI)]
  
  healthyDat <- healthyDatIMD[healthyDat, on="subjectID"]
  
  logNumber("Merge IMD",nrow(healthyDat))
  
  healthyDat[ age<0,  age:= NA] # added this because there was a negative age and I have no other information about the ages.
  
  #healthyDat[ , vitKcensored := vitK %like% "<"]
  #healthyDat[ , vitK := as.numeric(gsub("<","",vitK))]

  healthyDat[ , ethgr5:="Unknown"]
  
  healthyDat[ Ethnicity=="White British", ethgr5:="White"]
  healthyDat[ Ethnicity=="White-Any other background", ethgr5:="White"]
  healthyDat[ Ethnicity=="White-any other background", ethgr5:="White"]
  
  healthyDat[ Ethnicity=="Asian-any other background", ethgr5:="Asian"]
  healthyDat[ Ethnicity=="Asian or Asian British Indian", ethgr5:="Asian"]
  healthyDat[ Ethnicity=="Other Ethnic group-Chinese", ethgr5:="Asian"]
  
  healthyDat[ Ethnicity=="Black-Any other background", ethgr5:="Black"]
  healthyDat[ Ethnicity=="Black or Black British Caribbean", ethgr5:="Black"]
  healthyDat[ Ethnicity=="Black or Black British African", ethgr5:="Black"]
  
  healthyDat[ Ethnicity %like% "Mixed" ,  ethgr5:="Mixed"]
  healthyDat[ Ethnicity %like% "Mied" ,  ethgr5:="Mixed"]
  
  healthyDat[ SexNumeric==1 , sex := "Male"]
  healthyDat[ SexNumeric==2 , sex := "Female"]
  
  healthyDat[ is.na(ethgr5) , ethgr5:="Unknown"]
  healthyDat[ , subjectID:=paste0("NNUH_healthy_",subjectID)]

  healthyDat[ , bmi_sds := childsds::sds(bmi,age,sex,ref = childsds::ukwho.ref,
                                             male = "Male",female="Female",item = "bmi")]
  healthyDat[ , bmi_sds3 := childsds::sds(bmi,floor(age),sex,ref = childsds::ukwho.ref,
                                         male = "Male",female="Female",item = "bmi")]
  
  # Exclude where there is no blood, or no overlap in age with the clinical cases
  # Also exclude where the BMI is greated than 1.04 SDs
  
  logNumber("Just before exclusions",nrow(healthyDat))
  
  
  healthyDat <- healthyDat[ age>=2 ]
  
  logNumber("Exclude younger than 2 years",nrow(healthyDat))
  
  healthyDat <- healthyDat[ bmi_sds < 1.04  ]
  
  logNumber("Exclude BMI SDS greater than 1.04",nrow(healthyDat))
  
  healthyDat <- healthyDat[  !is.na(v25OHD_Std) ]
  
  logNumber("Vitamin D missing",nrow(healthyDat))
  
  
  copy(healthyDat[, .(subjectID,v25OHD_Std,ethgr5,age,sex,bmi,bmi_sds,imd2)])
  
}


## This function picks which is the descriptive to use for each statistic
weightedDescriptive <- function(x,w,digits=2) {
  res=NA
  if(all(is.na(x)) | length(x)==0) return("-")
  w <- w[!is.na(x)]
  x <- x[!is.na(x)]
  
  if(length(unique(x))<=2){
    #print(length(x))
    #print(length(w))
    #print(wtd.mean(x,weights =w))
    #print(table(x,useNA = "ifany"))
    res=sprintf("%d/%d (%d%%)",sum(x==1,na.rm = TRUE),sum(!is.na(x),na.rm = TRUE),round(100*wtd.mean(x,weights =w)))
    return(res)
  } else{
    wq=round(wtd.quantile(x,weights = w,probs = c(0.25,0.5,0.75),na.rm = FALSE),digits = digits)
    res=sprintf("%2g (%2g-%2g)",wq[2],wq[1],wq[3])  }
  res
}



makeCombinedData <- function(GOOSdata ,ndns, NNUHdata, NNUHhealthy){
    combinedData <-  rbindlist( list(GOOS=GOOSdata , NDNS=ndns, NNUH=NNUHdata, NNUHhealthy=NNUHhealthy ), 
                                idcol = "Source", fill=TRUE)
  
    combinedData[ Source=="GOOS" , Group := "GOOS"]
    combinedData[ Source=="NNUH" , Group := "NNUH"]
    combinedData[ Source=="NNUHhealthy" , Group := "NNUHhealthy"]
    combinedData[ Source=="NDNS" , Group := paste("NDNS",bmicat218)]
    combinedData[ Source=="GOOS" | Source=="NNUH" | Source=="NNUHhealthy" , finalweight := 1] # I'm not completely sure how to incorporate weights for the GOOS variables to ensure a fair comparison, but this seems reasonable.
    combinedData[ Source=="GOOS" , weight_goos := 1] #
    combinedData[ Source=="GOOS" , strata := 1] # 1 is the unique name for the GOOS stratum
    combinedData[ Source=="NNUH" , strata := 2] # 2 is the unique name for the NNUH stratum
    combinedData[ Group=="NNUHhealthy" , strata := 3] # 3 is the unique name for the NNUH healthy cohort
    
    combinedData[ is.na(subjectID) , subjectID:=seriali ]
    
    ## Subject IDs are unique to groups.
    combinedData[ , length(unique(Group)) , by=(subjectID) ]$V1 |> table()
    combinedData[ , PSU := ifelse(Source=="NDNS",PSU,subjectID)]
    
    combinedData[ ethgr5 %like% "Qn" , ethgr5:=NA] # Set ethnicity to missing in this case
    combinedData[ , ethgr5:= factor(ethgr5,levels=c("White","Asian","Black","Mixed","Other"))]
    combinedData[ ,Group := factor(Group,levels=c("NDNS Normal-weight", "NDNS Over-weight","NDNS Obese","GOOS","NNUH","NNUHhealthy")) ]
    combinedData[ !is.na(Sex), sex := Sex]
    combinedData[ is.na(Sex),  Sex := sex]
    
    ### Dropping one participant from Norfolk healthy group with missing age:
    ### Now - doesn't seem to do anything?
    # nrow(combinedData)
    # combinedData <- combinedData[!is.na(age)]
    # nrow(combinedData)
    
    ### Add the threshold for each outcome
    ### This is done BY REFERENCE as its a data.table
    
    
    # We should add half a year to the NDNS data because we are comparing 
    # the whole year to exact year in other sources
    
    combinedData[ Source=="NDNS", age := age + 0.5]
    
    combinedData[ , Source2 := ifelse(Source=="NNUHhealthy","NNUH",Source) ]
    combinedData[ , Source2 := factor(Source2 , levels=c("NDNS","GOOS","NNUH")) ]
    combinedData[ , Group2 := gsub("NDNS ","",Group) |> 
                    gsub("NNUHhealthy","Healthy weight",x=_) |>
                    gsub("NNUH","CYP\nObese",x=_) |>
                    factor(levels=c("Normal-weight","Over-weight",
                                    "Obese","GOOS",
                                    "Healthy weight",
                                    "CYP\nObese")) ]
    
    ## Convert the Serum folate unit to ng/Ml
    combinedData[Source %like% "NDNS" , SerumFolate := SerumFolate/2.26]
    
    
    ## Make an approximate IMD tertile group 
    combinedData[ , imd3 := cut(imd2,c(0,4,7,10))]
    combinedData[ , imd3 := cut(imd2,c(0,4,7,10))]
    
    combinedData[ , imd5 := floor((imd2+1)/2)]
    IMD_dat <- fread("dat/ndns_imd_dat.csv")[, subjectID:=as.character(seriali)][, seriali:=NULL]
    
    
    combinedData <- IMD_dat[ combinedData , on="subjectID"]
    combinedData[ Source %like% "NDNS", imd5 := IMD5]
    
    combinedData[ , imd5F := factor(imd5, levels=c(1,2,3,4,5),
                                    labels = c("1-2","3-4","5-6","7-8","9-10"))]
    
    combinedData[ , imd5C := as.numeric(imd5F)]
    combinedData[ , imd5O := ordered(imd5F)]
    
    
    
    combinedData[ Source %like% "NNUH" & imd3=="(0,4]"  , eqv3 := "Lowest Tertile" ]
    combinedData[ Source %like% "NNUH" & imd3=="(4,7]"  , eqv3 := "Middle Tertile"  ]
    combinedData[ Source %like% "NNUH" & imd3=="(7,10]" , eqv3 := "Highest Tertile"  ]
    
    combinedData[ Source %like% "GOOS" & imd3=="(0,4]"  , eqv3 := "Lowest Tertile" ]
    combinedData[ Source %like% "GOOS" & imd3=="(4,7]"  , eqv3 := "Middle Tertile"  ]
    combinedData[ Source %like% "GOOS" & imd3=="(7,10]" , eqv3 := "Highest Tertile"  ]
    
    combinedData[, eqv3F := factor(eqv3, levels=c("Highest Tertile","Middle Tertile","Lowest Tertile"))]
    combinedData
 
       
}
    
    
addPredictionsAndDeficiencies <- function(combinedData){
   
    # everything here should happen by reference
    # so we don't need to return the dataset
  
    ## This makes the centiles and the deficiency indicator
    ## for each outcome and places them within combineddat.
    
    ############################################################
    #### This makes a set of predicted values for each     #####
    #### outcome for each individual in the dataset.       #####
    ############################################################
    
    ## Doing this again to use 'combinedData' as the reference dataset.
    #predictions_All <-lapply(listOfOutcomes$varName,\(i){
    #  print(i)
    #  predictAll(mods[[i]], type="response" , 
    #             newdata = data.frame(age=combinedData$age))
    #})
    
    for(i in seq_along(outcomevector)){
      
      # FOR THE CENSORED OUTCOMES I AM NOT USING THE CENSORED DISTRIBUTIONS RIGHT NOW 
      # pFN <- pST3
      # if(listOfOutcomes$varName[i]%in% c("CRP","vitK")) pFN <- pST3lc
      # qFN <- qST3
      # if(listOfOutcomes$varName[i]%in% c("CRP","vitK")) qFN <- qST3lc
      
      v=outcomevector[i]
      print(v)
      combinedData[[paste0(v,"_deficient")]] <- 
        combinedData[[v]] < combinedData[[paste0(v,"_deficient_threshold")]]
      combinedData[[paste0(v,"_possible_deficient")]] <- 
        (combinedData[[v]] < combinedData[[paste0(v,"_possible_deficient_threshold")]] )&
        (combinedData[[v]] >= combinedData[[paste0(v,"_deficient_threshold")]] )
      combinedData[[paste0(v,"_elevated")]] <- 
        combinedData[[v]] > combinedData[[paste0(v,"_elevated_threshold")]]
      
      combinedData[[paste0(v,"_natural")]] <- combinedData[[v]] 
      
      #preddat <- copy(cbind(combinedData,as.data.table(predictions_All[[i]])))
      #preddat$logtarget <- log(preddat[[listOfOutcomes$varName[i]]])
      #preddat <- preddat[!is.na(logtarget)]
      
      #preddat[ , centile2 := pFN(q=logtarget, mu=mu, sigma=sigma,nu=nu,tau=tau)]
      #preddat[ , prediction := qFN(p = 0.5, mu=mu, sigma=sigma,nu=nu,tau=tau)]
      
      #combinedData[ !is.na(combinedData[[v]]) , 
      #              eval(paste0(v,"_centile")) := preddat$centile2 ]
      
      #combinedData[ !is.na(combinedData[[v]]) , 
      #              eval(paste0(v,"_prediction")) := preddat$prediction ]
    }
    
    
    ## These two are special, the only outcomes with a 'possible elevated' category:
    combinedData[ , VitB1_possible_elevated := VitB1 > VitB1_possible_elevated_threshold &
                    VitB1 <= VitB1_elevated_threshold   ]
    #As of June 2026 we don't use the VitB2 possible elevated threshold
    #combinedData[ , VitB2_possible_elevated := VitB2 > VitB2_possible_elevated_threshold  ]
    
    
    sumDefCols <- c(
      
      "PFerritin_deficient", 
      "v25OHD_Std_deficient",
      "SVitB12_deficient",
      "Retinol_deficient",
      "TotalCarotenoids_deficient",
      "vitE_deficient",
      "SerumFolate_deficient",
      "VitC_deficient",
      "VitB1_elevated",
      "VitB2_elevated",
      "VitB6PLP_deficient",
      "Zn_deficient"
      
    )
    
    sumDefCols
    
    ### ONLY INCLUDE DEFICIENCIES MEASURED IN GOOD AND NDNS.
    ### Request for B12 to be included instead of Holo TC for the GOOS/NDNS comparison
    sumDefCols2 <- c(
      "SVitB12_deficient",
      "PFerritin_deficient", 
      "v25OHD_Std_deficient",
      "Retinol_deficient",
      "TotalCarotenoids_deficient",
      "vitE_deficient"
      )
    
    # Sum how many deficiencies from the list are present in each participant
    deficientCount0 <- combinedData[ , .SD , .SDcols = sumDefCols] |> 
      lapply(as.numeric) |> 
      as.data.table() |>
      rowSums(na.rm = TRUE) # 
    combinedData[ , deficientCount := deficientCount0]
    combinedData[ , anyDeficiency := deficientCount>0]
    combinedData[ , multipleDeficiency := deficientCount>1]
    
    # Sum how many deficiencies from the list are present in each participant
    deficientCount20 <- combinedData[ , .SD , .SDcols = sumDefCols2] |> 
      lapply(as.numeric) |> 
      as.data.table() |>
      rowSums(na.rm = TRUE) # 
    combinedData[ , deficientCount2 := deficientCount20]
    combinedData[ Source=="NNUH", deficientCount2:=NA]  # This variable not defined for NNUH
    combinedData[ Source=="NNUHhealthy", deficientCount2:=NA]  # This variable not defined for NNUH
    
    combinedData[ , anyDeficiency2 := deficientCount2>0]
    combinedData[ , multipleDeficiency2 := deficientCount2>1]
    
    
    denominatorCount <- combinedData[ , .SD , .SDcols = sumDefCols] |> 
      lapply(as.numeric) |> 
      as.data.table() |> 
      lapply(\(v)  !is.na(v)) |> 
      as.data.table() |> 
      rowSums(na.rm = TRUE)
    combinedData[ , denominatorCount := denominatorCount]
    
    # table(combinedData$denominatorCount, combinedData$Group)
    # table(combinedData$deficientCount, combinedData$Group)
    # table(combinedData$deficientCount2, combinedData$Group)
    # 
    # table(combinedData$anyDeficiency, combinedData$Group)
    # table(combinedData$multipleDeficiency, combinedData$Group)
    # table(combinedData$anyDeficiency2, combinedData$Group)
    # table(combinedData$multipleDeficiency2, combinedData$Group)
    
    ## There are seven deficiencies measured in NNUH.
    ## Ferritin, Vitamin D, Holo-TC, Retinol, Vitamin E, Zinc and Serum Folate
    
    ## Note which deficiencies are included in each cohort.
    ## It is possible to impute the missing values in a lot of cases..
    ## We are treating 'missing' deficiencies as absent.
    ## Ideally we would attempt this..
    
    ### Global assignments here because I need these labels later..
    nnuhSumVars <<- 
      labelsvector[sapply(presentInNNUH , \(w) length(grep(w , sumDefCols))>0) |> which() |> names()] |> 
      paste(collapse=", ")
    GOOSSumVars <<- 
      labelsvector[sapply(presentInGOOS , \(w) length(grep(w , sumDefCols))>0) |> which() |> names()]|> 
      paste(collapse=", ")
    ndnsSumVars <<- 
      labelsvector[sapply(presentInNDNS , \(w) length(grep(w , sumDefCols))>0) |> which() |> names()]|> setdiff("Selenium") |>
      paste(collapse=", ")
  
    
    
    GOOSSumVars2 <<- 
      labelsvector[sapply(presentInGOOS , \(w) length(grep(w , sumDefCols2))>0) |> which() |> names()]|>   paste(collapse=", ")
    ndnsSumVars2 <<- 
      labelsvector[sapply(presentInNDNS , \(w) length(grep(w , sumDefCols2))>0) |> which() |> names()]|> paste(collapse=", ")
    
    combinedData[ , ethgr2:=ifelse(ethgr5=="White","White","Non-white")]
    
    
    combinedData

}

####### FUNCTIONS TO MAKE GRAPHS

# This function adds pairwise p-values, only allows three groups (NDNS only)
makeObesityGraph <- function(i,modslist){
  emgrid <- emmeans(modslist[[i]],
                    ~Group,
                    type="response",
                    data=modslist[[i]]$survey.design$variables)
  plotdat <- emgrid |> contrast("pairwise",adjust="none",infer=c(TRUE,TRUE)) |> 
    as.data.frame() |> as.data.table() |>  _[ , c("group1","group2") := tstrsplit(contrast, " / ")][]
  plotdat[, group1 := gsub(pattern = "(","",group1,fixed = TRUE)]
  plotdat[, group1 := gsub(pattern = ")","",group1,fixed = TRUE)]
  plotdat[, group2 := gsub(pattern = "(","",group2,fixed = TRUE)]
  plotdat[, group2 := gsub(pattern = ")","",group2,fixed = TRUE)]
  plotdat2 <- emgrid  |> as.data.frame() |> as.data.table() 
  height=max(plotdat2$upper.CL) - min(plotdat2$lower.CL)
  plotdat$ypos <- max(plotdat2$upper.CL) + height*c(0.2,0.6,0.4)
  plotdat[ , p.value := ifelse(p.value<0.001,"p<0.001",sprintf("p=%0.3f",p.value))]
  ggplot(plotdat2) + 
    aes(x=Group,y=response,ymin=lower.CL,ymax=upper.CL) + 
    geom_errorbar(width=0.2) + geom_point() + 
    ggpubr::stat_pvalue_manual(data=plotdat , y.position = "ypos",label = "p.value") + 
    theme_bw() + 
    labs(x="BMI SDS Group",y=labelsvector2[i])
} 

## This function is designed to include GOOS as well where needed, 
## and only use treatment contrasts
makeObesityGraph2 <- function(i,modslist){
  emgrid <- emmeans(modslist[[i]],
                    ~Group,
                    type="response",
                    at=list(age=12,sex="Male",ethgr5="White"),
                    data=modslist[[i]]$survey.design$variables)
  plotdat <- emgrid |> contrast("trt.vs.ctrl",adjust="none",infer=c(TRUE,TRUE)) |> 
    as.data.frame() |> as.data.table() |>  _[ , c("group1","group2") := tstrsplit(contrast, " / ")][]
  plotdat[, group1 := gsub(pattern = "(","",group1,fixed = TRUE)]
  plotdat[, group1 := gsub(pattern = ")","",group1,fixed = TRUE)]
  plotdat[, group2 := gsub(pattern = "(","",group2,fixed = TRUE)]
  plotdat[, group2 := gsub(pattern = ")","",group2,fixed = TRUE)]
  plotdat2 <- emgrid  |> as.data.frame() |> as.data.table() 
  height=max(plotdat2$upper.CL) - min(plotdat2$lower.CL)
  plotdat$ypos <- max(plotdat2$upper.CL) + height*(1:nrow(plotdat))/10
  plotdat[ , p.value := ifelse(p.value<0.001,"p<0.001",sprintf("p=%0.3f",p.value))]
  ggplot(plotdat2) + 
    aes(x=Group,y=response,ymin=lower.CL,ymax=upper.CL) + 
    geom_errorbar(width=0.2) + geom_point() + 
    ggpubr::stat_pvalue_manual(data=plotdat , y.position = "ypos",label = "p.value",size=3.5) + 
    theme_bw() + 
    labs(x="BMI SDS Group",y=labelsvector3b[i])
} 






makeObesityTable2 <- function(i,modslist){
  emgrid <- emmeans(modslist[[i]],
                    ~Group,
                    type="response",
                    at=list(age=12,sex="Male",ethgr5="White"),
                    data=modslist[[i]]$survey.design$variables)
  # plotdat <- emgrid |> contrast("trt.vs.ctrl",adjust="none",infer=c(TRUE,TRUE)) |> 
  #   as.data.frame() |> as.data.table() |>  _[ , c("group1","group2") := tstrsplit(contrast, " / ")][]
  # plotdat[, group1 := gsub(pattern = "(","",group1,fixed = TRUE)]
  # plotdat[, group1 := gsub(pattern = ")","",group1,fixed = TRUE)]
  # plotdat[, group2 := gsub(pattern = "(","",group2,fixed = TRUE)]
  # plotdat[, group2 := gsub(pattern = ")","",group2,fixed = TRUE)]
  plotdat2 <- emgrid  |> as.data.frame() |> as.data.table() 
  #plotdat[ , p.value := ifelse(p.value<0.001,"p<0.001",sprintf("p=%0.3f",p.value))]
  plotdat2
} 



makeIncomeGraph <- function(i,modslist){
  emgrid <- emmeans(modslist[[i]],
                    ~eqv3F,
                    type="response",
                    data=modslist[[i]]$survey.design$variables)
  plotdat <- emgrid |> contrast("pairwise",adjust="none",infer=c(TRUE,TRUE)) |> 
    as.data.frame() |> as.data.table() |>  _[ , c("group1","group2") := tstrsplit(contrast, " / ")][]
  plotdat2 <- emgrid  |> as.data.frame() |> as.data.table() 
  height=max(plotdat2$upper.CL) - min(plotdat2$lower.CL)
  plotdat$ypos <- max(plotdat2$upper.CL) + height*c(0.2,0.6,0.4)
  plotdat[ , p.value := ifelse(p.value<0.001,"p<0.001",sprintf("p=%0.3f",p.value))]
  ggplot(plotdat2) + 
    aes(x=eqv3F,y=response,ymin=lower.CL,ymax=upper.CL) + 
    geom_errorbar(width=0.2) + geom_point() + 
    ggpubr::stat_pvalue_manual(data=plotdat , y.position = "ypos",label = "p.value") + 
    theme_bw() + 
    labs(x="Income tertile",y=labelsvector3b[i])
} 

makeIMDGraph <- function(i,modslist){
  emgrid <- emmeans(modslist[[i]],
                    ~imd5F,
                    type="response",
                    data=modslist[[i]]$survey.design$variables)
  plotdat <- emgrid |> contrast("trt.vs.ctrl",adjust="none",infer=c(TRUE,TRUE)) |> 
    as.data.frame() |> as.data.table() |>  
    _[ , c("group1","group2") := tstrsplit(contrast, " / ")] |>
    _[ , group1 := gsub(pattern = "\\(","",group1)] |>
    _[ , group1 := gsub(pattern = "\\)","",group1)]|>
    _[ , group2 := gsub(pattern = "\\(","",group2)]|>
    _[ , group2 := gsub(pattern = "\\)","",group2)][]
  
  plotdat2 <- emgrid  |> as.data.frame() |> as.data.table() 
  height=max(plotdat2$upper.CL) - min(plotdat2$lower.CL)
  
  #plotdat$ypos <- max(plotdat2$upper.CL) + height*c(0.2,0.3,0.4,.5)
  plotdat$ypos <- max(plotdat2$upper.CL) + height*c(0.2,0.1,0.2,0.1)
  plotdat[ , p.value := ifelse(p.value<0.001,"p<0.001",sprintf("p=%0.3f",p.value))]
  print(plotdat)
  ggplot(plotdat2) + 
    aes(x=imd5F,y=response,ymin=lower.CL,ymax=upper.CL) + 
    geom_errorbar(width=0.2) + geom_point() + 
    #ggpubr::stat_pvalue_manual(data=plotdat , y.position = "ypos",label = "p.value") + 
    geom_text(aes(y=ypos,x=group1,label=p.value),inherit.aes = FALSE,data=plotdat,size = 10/.pt)+
    theme_bw() + 
    labs(x="IMD quintile (lower = more deprived)",y=labelsvector3b[i])
} 

makeEthGraph <- function(i,modslist){
  emgrid <- emmeans(modslist[[i]],
                    ~ethgr5,
                    type="response",
                    data=modslist[[i]]$survey.design$variables)
  plotdat <- emgrid |> contrast("trt.vs.ctrl",adjust="none",infer=c(TRUE,TRUE)) |> 
    as.data.frame() |> as.data.table() |>  _[ , c("group1","group2") := tstrsplit(contrast, " / ")][]
  plotdat2 <- emgrid  |> as.data.frame() |> as.data.table() 
  height=max(plotdat2$upper.CL) - min(plotdat2$lower.CL)
  plotdat$ypos <- max(plotdat2$upper.CL) + height*c(0.2,0.4,0.6,0.8)
  plotdat[ , p.value := ifelse(p.value<0.001,"p<0.001",sprintf("p=%0.3f",p.value))]
  ggplot(plotdat2) + 
    aes(x=ethgr5,y=response,ymin=lower.CL,ymax=upper.CL) + 
    geom_errorbar(width=0.2) + geom_point() + 
    ggpubr::stat_pvalue_manual(data=plotdat , y.position = "ypos",label = "p.value") + 
    theme_bw() + 
    labs(x="Ethnicity",y=labelsvector3b[i])
} 


## make a standardisation weight (by reference)
makeStandardisationWeight <- function(combinedData){
  combinedData[, ageband := cut(age,breaks = c(0, 5, 10, 15, 20, Inf),right = FALSE)]
  combinedData[, std_stratum := interaction(ageband, sex, drop = TRUE)]
  combinedData[, base_wt := finalweight]
  target_dist <- combinedData[Group%like%"NDNS" &
                                !is.na(std_stratum) & !is.na(Group),.(target_n = sum(base_wt, na.rm = TRUE)),by = std_stratum] |>
    _[,  target_prop := target_n / sum(target_n)] |>
    _[, .(std_stratum, target_prop)]
  group_dist <- combinedData[!is.na(std_stratum) & !is.na(Group),
                             .(group_n = sum(base_wt, na.rm = TRUE)),  by = .(Group, std_stratum)] |>
    _[, group_prop := group_n / sum(group_n), by = Group]
  std_weights <- group_dist[target_dist, on = "std_stratum"] |> 
    _[, std_factor := target_prop / group_prop] |> 
    _[,  .(Group, std_stratum, std_factor) ]
 
   #combinedData <- std_weights[combinedData, on = .(Group, std_stratum)]
  
  combinedData[ # this is a target join, which makes the file work by reference
    std_weights,
    on = .(Group, std_stratum),
    std_factor := i.std_factor
  ]
  
  combinedData[, analysis_wt := base_wt * std_factor]
}
