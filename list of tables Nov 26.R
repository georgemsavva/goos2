

### List of tables for the paper along with code:


########################################################
## Table 2: Make the biomarkers descriptive table    ###
########################################################
# Adjusted means of the biomarkers (perhaps replace with medians, and then use the figure to show the means? )

bioMarkerDescriptiveTable <- function(dir="june26figs",filename="simpleMedians.docx"){

  metricsTab2 <- combinedData[ , .SD , 
                              .SDcols = c("Group","finalweight","strata",order_of_outputs2)] |> 
    melt( id.vars = c("finalweight","strata","Group") ) |> 
    _[ , weightedDescriptive(value,w=finalweight,digits = 2) , by=.(Group, variable)] |> 
    dcast(variable ~ Group, value.var="V1")  
  
  ## Replace NAs with "-"
  metricsTab2 <- metricsTab2[ !apply((metricsTab2[, -1, with = FALSE]), 1, \(r) all(r == "-")) ] 
  metricsTab2[ , variable := factor(variable,levels=order_of_outputs2)]
  metricsTab2 <- metricsTab2[ order(variable)]
  metricsTab2[ , variable := labelsvector3[as.character(variable)]]
  
  metricsTab2 |> 
    gt() |> 
    tab_header("Weighted median and IQR for each outcome in each group") |>
    gtsave(file.path(dir,filename))

}

bioMarkerDescriptiveTableWeighted <- function(dir="june26figs",filename="weightedMedians.docx"){
  
  metricsTab2 <- combinedData[ , .SD , 
                               .SDcols = c("Group","analysis_wt","strata",order_of_outputs2)] |> 
    melt( id.vars = c("analysis_wt","strata","Group") ) |> 
    _[ , weightedDescriptive(value,w=analysis_wt,digits = 2) , by=.(Group, variable)] |> 
    dcast(variable ~ Group, value.var="V1")  
  
  ## Replace NAs with "-"
  metricsTab2 <- metricsTab2[ !apply((metricsTab2[, -1, with = FALSE]), 1, \(r) all(r == "-")) ] 
  metricsTab2[ , variable := factor(variable,levels=order_of_outputs2)]
  metricsTab2 <- metricsTab2[ order(variable)]
  metricsTab2[ , variable := labelsvector3[as.character(variable)]]
  
  metricsTab2 |> 
    gt() |> 
    tab_header("Weighted median and IQR for each outcome in each group") |>
    gtsave(file.path(dir,filename))
  
}


#####################################################
# Table 1: Participant demographic characteristics  #
#####################################################

makeDemographicTable <- function(dir="oct26figs" , filename="demographics.docx"){

  combinedData[,weight2 := finalweight / sum(finalweight) * .N , by=Group]
  
  sdc = c("Age"="age",
          "Sex"="Sex",
          "BMI_SDS"="bmi_sds",
          "Income"="eqv3F",
          "Deprivation"="imd5F",
          "Ethnicity"="ethgr5",
          "weight2"="weight2",
          "strata"="strata",
          "PSU"="PSU",
          "Group"="Group")
  svy <- combinedData[ , .SD , .SDcols = sdc] |> 
    setNames(names(sdc))  |>
    svydesign(ids=~PSU,weights=~weight2, strata=~strata,data=_)|> 
    tbl_svysummary(by=Group,
                   #missing = "no",
                   include= -c(strata,weight2,PSU),
                   statistic = list(all_categorical() ~
                                      "{n_unweighted} ({p}%)"),
                   missing_stat = "{N_miss_unweighted}") |> 
    as.data.table() 
  
  
  
  # |>
  #   setNames(names(metricsTab))
  #svy[ GOOS == "0 (NA%)", GOOS := "-"]
  
  svy |> gt() |> gtsave(file.path(dir,filename))
  
}
  
#######################################################################################
# Table 3: Deficiency / Excess (prevalence) with p-values (and possibly odds ratios)  #
#######################################################################################

prevalenceOfDeficiencyTable <- function(dir="june26figs" , filename="prevalence_of_deficiency.docx"){
  
  svy2 <- combinedData[ , .SD , .SDcols = 
                          c(deficiencyoutcomes,"finalweight","strata","PSU","Group")] |> 
    #setNames(names(sdc))  |>
    svydesign(ids=~PSU,weights=~finalweight, strata=~strata,data=_) |> 
    tbl_svysummary(by=Group,
                   missing = "no",
                   include= -c(strata,finalweight,PSU),
                   statistic = list(all_categorical() ~
                                      "{n_unweighted} ({p}%)")) |> 
    as.data.table() |>
    setNames(names(metricsTab))
  #svy[ GOOS == "0 (NA%)", GOOS := "-"]
  
  svy3 <- svy2[ , marker := stringr::word(variable,1,sep=fixed("_"))] |>
    _[ marker=="v25OHD", marker := "v25OHD_Std"] |> 
    _[ marker %in% order_of_outputs] |> 
    _[ , marker := factor(marker, levels=order_of_outputs)] |> 
    _[ variable %like% "deficient" , threshold := "Deficient"] |>
    _[ variable %like% "possible_deficient" , threshold := "Insufficient"] |>
    _[ variable %like% "elevated" , threshold := "Elevated"] |>
    _[ variable %like% "elevated" & marker %in% c("VitB1","VitB2") , threshold := "Deficient"] |>
    _[ variable %like% "possible_elevated" & marker %in% c("VitB1","VitB2") , threshold := "Insufficient"] |>
    _[order(marker,threshold)] |>
    _[ marker %notin% presentInGOOS , GOOS := ""] |>
    _[ marker %notin% presentInNNUH , NNUH := ""] |>
    _[ marker %notin% presentInNNUHhealthy , NNUHhealthy := ""] |>
    _[ , .(marker=labelsvector[as.character(marker)],
           threshold,`NDNS Normal-weight`,`NDNS Over-weight`,`NDNS Obese`,`GOOS`,`NNUH`,`NNUHhealthy`)] 
      
    svy3 |> gt() |> tab_header("Prevalence of deficiency (N, weighted %) across groups") |>
    gtsave(filename = file.path(dir,filename))

}
  
########################################
# Table S1:  Sample sizes per outcome  #
########################################

sampleSizePerOutcomeTable <- function(dir="june26figs",filename="listOfOutcomes.docx"){
  
  combinedData[ , .SD , .SDcols = c("subjectID","Group",order_of_outputs2,"cholesterol")] |> 
    melt(id.vars = c("subjectID","Group"),variable.name = "varName") |> 
    _[, varName := factor(varName,levels=c(order_of_outputs2,"vitE","cholesterol"))] |> 
    _[!is.na(value)] |> 
    dcast(varName ~ Group , fun.aggregate = length) |>
    _[, unit:=units[as.character(varName)]] |>
    _[, outcome:=c(labelsvector,"cholesterol"="Cholesterol")[as.character(varName)]] |>
    _[,varName:=NULL] |> setcolorder(c("outcome","unit")) |>
    gt() |>
    tab_header("List of outcomes with number of non-missing entries per dataset") |>
    gtsave(file.path(dir,filename))
}
  

#### Supplemental figure - ridge plots of distributions (replaces box plots)

makeRidgePlot <- function(dir="june26figs",filename="ridgePlots.png",...){
  combinedData[ !(Group%like%"NNUH"), .SD , .SDcols = c(presentInNDNS, "Group", "analysis_wt","subjectID")] |>   melt(id.vars = c("subjectID","analysis_wt","Group"),variable.name = "outcome",value.name = "value") |>
    _[ , outcome2 := factor(outcome,levels=order_of_outputs)] |>
    _[ !is.na(outcome2)] |>
    ggplot() + 
    aes(x=value,y=droplevels(Group),weight = analysis_wt,fill=Group) + 
    geom_density_ridges(scale=1.5,quantile_lines = TRUE,alpha=0.6,quantiles = c(0.5)) + theme_bw() +
    scale_fill_viridis_d() +
    facet_wrap(~outcome2, scale = "free_x",labeller=labeller(outcome2=labelsvector3b,nrow=5),
               strip.position="bottom")+
    scale_x_log10() + 
    labs(y=NULL,x=NULL)+
    theme(legend.position ="none", 
          strip.background = element_blank(),
          strip.placement="outside", 
          strip.text.x.bottom = element_text(
            margin = margin(t = -1, r = 0, b = 0, l = 0)
          ))+
    scale_y_discrete(labels=groupLabelsVector)
  ggsave(filename=file.path(dir,filename),...)
}

  
# Table S2:  Assays

# Table S3:  Assays 

# Table S4:  Thresholds

# Table S5:  Unadjusted medians (replace Table 2?)

# Table S6:  Medians and centiles (very similar to S5)

# Table S7:  Daily dietary nutrient intakes??
  
  
  #####################################################################       
  ###### CALCULATE THE REGRESSION TABLE FOR THE BINARY OUTCOMES.  #####
  #####################################################################       


makeDeficiencyTable <- function(dir="june26figs",filename="propresults3_SHORT.docx",analysisweight=FALSE){
  
  deficiencyoutcomes <- grep("(elevated|deficient|Deficiency2)$", names(combinedData),value=TRUE)
  outcomevector <- listOfOutcomes$varName |> setNames(listOfOutcomes$labels)
  outcomevector2 <- sort(c(paste0(outcomevector,"_natural")|>
                             setNames(paste(names(outcomevector),"Natural")),
                           paste0(outcomevector,"_deficient")|>
                             setNames(paste(names(outcomevector),"Deficient")) ,
                           paste0(outcomevector,"_possible_deficient")|>
                             setNames(paste(names(outcomevector),"Possible Deficient")) ,
                           paste0(outcomevector,"_elevated")|>
                             setNames(paste(names(outcomevector),"Elevated")) ,
                           paste0(c("VitB1"),"_possible_elevated")|>
                             setNames(paste(c("Vitamin B1"),"Possible Elevated")) ,
                           paste0(outcomevector,"_centile")  |>
                             setNames(paste(names(outcomevector),"Centiles"))))
  
  ## Identify the rows with at least one positive (binary outcome) in each group
  rowsToInclude <- combinedData[ , .SD, .SDcols=c("Group",deficiencyoutcomes)] |> 
    melt(idcols="Group") |> 
    _[!is.na(value)][ , sum(value) , by=.(Group,variable)] |> 
    _[ , min(V1, na.rm=TRUE), by=variable][V1>0][, as.character(variable)]
  rowsToInclude <- setdiff(rowsToInclude,"Zn_elevated")
  
  ## Identify the rows with at least one positive binary outcome overall.
  rowsToInclude2 <- combinedData[ , .SD, .SDcols=c("Group",deficiencyoutcomes)] |> 
    melt(idcols="Group") |> 
    _[!is.na(value)][ , sum(value) , by=.(Group,variable)] |> 
    _[ , max(V1), by=variable][V1>0][, as.character(variable)]
  rowsToInclude2 <- setdiff(rowsToInclude2,"Zn_elevated")
  rowsToInclude2 <- setdiff(rowsToInclude2,"vitE_elevated")
  
  if(analysisweight) combinedData[ , useweight := analysis_wt]
  else combinedData[ , useweight := finalweight]
  
  metricsTab <- combinedData[ , .SD , 
                              .SDcols = c("Group","useweight","strata",outcomevector2,
                                          "anyDeficiency","multipleDeficiency",
                                          "anyDeficiency2","multipleDeficiency2")] |> 
    melt( id.vars = c("useweight","strata","Group") ) |> 
    _[ , weightedDescriptive(value,w=useweight) , by=.(Group, variable)] |> 
    dcast(variable ~ Group, value.var="V1")  
  
  ## Replace NAs with "-"
  metricsTab <- metricsTab[ !apply((metricsTab[, -1, with = FALSE]), 1, \(r) all(r == "-")) ] 
  

  
  propresults3 <- lapply( propmods[rowsToInclude], \(m){
    #m=propmods[[1]]
    emmeans::emmeans(m, ~Group , data=m$survey.design$variables,at=list(age=12,sex="Male",ethgr2="White")) |> 
      contrast("revpairwise",adjust="none") |> 
      as.data.frame() |> 
      subset(contrast %in% c(
        "(NDNS Over-weight) - (NDNS Normal-weight)", 
        "NDNS Obese - (NDNS Normal-weight)",
        "GOOS - (NDNS Normal-weight)",
        #"NNUHhealthy - NNUH"
        "NNUH - (NDNS Normal-weight)"
      ))
  }) |> rbindlist(idcol="outcome") |>
    _[ , string := formatp(p.value)][] |>
    dcast( outcome ~ contrast , value.var="string",fill="-" ) |>
    _[ , .(variable=outcome,
           OWp=`(NDNS Over-weight) - (NDNS Normal-weight)`,
           OBp=`NDNS Obese - (NDNS Normal-weight)`,
           GOp=`GOOS - (NDNS Normal-weight)`,
           #NNUHp=`NNUHhealthy - NNUH`
           NNUHp=`NNUH - (NDNS Normal-weight)`
    )]  |>
    _[metricsTab,on="variable"] |>
    _[, .(outcome=variable, 
          `NDNS Normal-weight`, 
          `NDNS Over-weight`=paste(`NDNS Over-weight`,ifelse(is.na(OWp),"",OWp)),
          `NDNS Obese`=paste(`NDNS Obese`,ifelse(is.na(OBp),"",OBp)),
          `GOOS`=paste(`GOOS`,ifelse(is.na(GOp),"",GOp)),
          `NNUH`=paste(`NNUH`,ifelse(is.na(NNUHp),"",NNUHp)),
          `NNUH Healthy weight`=NNUHhealthy
          #`NNUH`=paste(`NNUH`,ifelse(is.na(NNUHp),"",NNUHp),ifelse(is.na(NNUHp2),"",NNUHp2))
    )] |>
    _[outcome %like% "deficient|elevated|Deficiency" ] |>
    
    _[outcome %in% c(rowsToInclude2,"anyDeficiency","multipleDeficiency")] |>  # uncomment this to get the smaller table
    _[ , outcome := sub(x=outcome,replacement = "v25OHDStd",pattern = "v25OHD_Std",fixed = TRUE)] |>
    _[ , c("var","outcome"):=tstrsplit(sub(x=outcome,replacement = "__",pattern = "_",fixed = TRUE),split="__",fixed=TRUE)] |>
    _[ , var := sub(x=var,replacement = "v25OHD_Std",pattern = "v25OHDStd",fixed = TRUE)] |>
    _[ var=="VitB1" & outcome=="elevated", outcome:="deficient"]  |> #  BECAUSE THE MEASURE IS BACKWARDS
    _[ var=="VitB1" & outcome=="possible_elevated", outcome:="possible_deficient"]  |> #  BECAUSE THE MEASURE IS BACKWARDS
    _[ var=="VitB2" & outcome=="elevated", outcome:="deficient"]  |> #  BECAUSE THE MEASURE IS BACKWARDS
    _[ var=="VitB2" & outcome=="possible_elevated", outcome:="possible_deficient"]  |> #  BECAUSE THE MEASURE IS BACKWARDS
    _[ ,var:=factor(var,levels=c(order_of_outputs,"vitE","CRP","anyDeficiency","anyDeficiency2","multipleDeficiency","multipleDeficiency2"))] |>
    _[ ,outcome:=factor(outcome,levels=c("possible_deficient","deficient","elevated"))] |>
    _[ order(var,outcome)] |>
    _[,var:=labelsvector4[as.character(var)]] |>
    setcolorder(c("var","outcome")) |>
    _[]
  
  
  # gt(propresults3) |> cols_align(align = "left") |>
  #   tab_header("Deficiency/elevated statistics by group","With p-values corresponding to survey weighted comparison with NDNS healthy group / NNUH healthy group") |>
  #   gtsave(filename = "graphs/propresults3.html")
  # 
  
  ###### WRITE EVERYTHING OUT.
  
  # gt(propresults3) |>cols_align(align = "left")|>
  #   tab_header("Deficiency/elevated statistics by group","Prevalences are weighted proportions, p-values are from regression models adjusted for age, sex and ethnicity") |>
  #   gtsave(filename = "finalFigures/propresults3_ALL.docx")

  
  
footnote1 <-  sprintf("The deficiencies included in the 'at least one' or 'more than one' deficiency outcomes vary across cohorts, and are: \n For NDNS: %s.  \n For GOOS: %s.  \n For NNUH: %s.  For these rows, missing values are counted as not deficient." , ndnsSumVars, GOOSSumVars , nnuhSumVars) 
  
footnote2 <-  sprintf("For the GOOS only row, the variables are:%s", GOOSSumVars2) 

  gt(propresults3) |>cols_align(align = "left")|>

    tab_header("Deficiency/elevated statistics by group","Prevalences are weighted proportions, p-values are from regression models adjusted for age, sex and ethnicity") |> 
    tab_footnote(footnote1) |>
    tab_footnote(footnote2) |>
    gtsave(filename = file.path(dir,filename))
  
  }

propmodsVitmainDspecial_deficient |> summary()
  

marginalMeansTable <- function(dir="june26figs",filename="marginalMeansTable.docx"){
  
  tabs = lapply(order_of_outputs2, makeObesityTable2, naturalMods2)
  tabs = setNames(tabs,order_of_outputs2)
  rbindlist(tabs,idcol = "outcome") |>
    _[ , GroupLab := groupLabelsVector[Group]] |>
    _[ , outcome2 := factor(labelsvector3[outcome],levels=labelsvector3[order_of_outputs2])] |>
    _[ , outstr := sprintf("%0.2f (%0.2f, %0.2f)",response,lower.CL,upper.CL)] |>
    dcast(outcome2 ~ GroupLab , value.var = "outstr" ) |>
    _[ , .SD , .SDcols=c("outcome2",groupLabelsVector)] |>
    gt()|>cols_align(align = "left")|>    
    sub_missing(missing_text = "-") |>
    tab_header("Marginal means, adjusted for age, sex, and ethnicity") |> 
    gtsave(filename = file.path(dir,filename))
  }


