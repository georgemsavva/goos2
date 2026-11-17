
###################################################################
#### MODELLING THE PREVALENCE OF DEFICIENCY BETWEEN GROUPS.
###################################################################


### We need to collapse the ethncity groups into two for the proportion models, 
### because of the sparsity of some of the outcomes leading to nonsensical values

combined_design <- svydesign(ids=~PSU,weights=~finalweight,strata=~strata,data=copy(combinedData))

deficiencyoutcomes <- grep("(elevated|deficient|Deficiency2)$", names(combinedData),value=TRUE)

## Remove vitamin K 
deficiencyoutcomes <- deficiencyoutcomes[!deficiencyoutcomes%like%"vitK"]
deficiencyoutcomes <- deficiencyoutcomes[!deficiencyoutcomes%like%"pivka"]

## Use only the *defined* deficiency outcomes.
deficiencyoutcomes <- sapply(deficiencyoutcomes, \(s) !all(is.na(combinedData[[s]]))) |> 
  which() |> names()

## Identify the rows with at least one positive binary outcome overall.
rowsToInclude2 <- combinedData[ , .SD, .SDcols=c("Group",deficiencyoutcomes)] |> 
  melt(idcols="Group") |> 
  _[!is.na(value)][ , sum(value) , by=.(Group,variable)] |> 
  _[ , max(V1), by=variable][V1>0][, as.character(variable)]

### We estimate all the outcomes, but we don't use the models where there is not enough data
### So we can ignore warnings here
### But you can't use Zn elevated, because that only exists for NNUH
propmods <- lapply(deficiencyoutcomes, function(y) {
  tryCatch({
    ## Using linear age here because the data is very sparse
    formula <- as.formula(paste(y, "~ Group + age + sex + ethgr2"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design, cc & Source %in% c("NDNS","GOOS","NNUH"))
    svyglm(formula, design =  design_cc, family = quasibinomial())
  }, error = function(e) NULL )
}) |> setNames(deficiencyoutcomes)

propmods <- propmods[!sapply(propmods,is.null)]

propmodsVitmainDspecial_deficient <- 
  glm(v25OHD_Std_deficient ~ Group + age + sex + ethgr2 , family="quasibinomial",data=combinedData[Source %like% "NNUH"])
propmodsVitmainDspecial_possible_deficient <- 
  glm(v25OHD_Std_possible_deficient ~ Group + age + sex + ethgr2 , family="quasibinomial",data=combinedData[Source %like% "NNUH"])

propmodsVitmainDspecial_deficient |> summary()
propmodsVitmainDspecial_possible_deficient |> summary()

# Extract the results from the models.  This creates a table of odds ratios related to the reference group.
propresults <- propmods |> 
  lapply(broom::tidy,exponentiate=TRUE) |> 
  rbindlist(idcol="outcome") |> 
  _[ , string := sprintf("%0.2f (%0.2f) p=%0.3f %s",estimate,std.error,p.value,stars(p.value))][] |>
  dcast( outcome ~ term , value.var="string" ) |>
  
  _[ , marker := stringr::word(outcome,1,sep=fixed("_"))] |>
  _[ marker=="v25OHD", marker := "v25OHD_Std"] |> 
  _[ marker %in% order_of_outputs, marker := factor(marker, levels=order_of_outputs)] |> 
  _[ outcome %like% "deficient" , threshold := "Deficient"] |>
  _[ outcome %like% "possible_deficient" , threshold := "Insufficient"] |>
  _[ outcome %like% "elevated" , threshold := "Elevated"] |>
  _[ outcome %like% "elevated" & marker %in% c("VitB1","VitB2") , threshold := "Deficient"] |>
  _[ outcome %like% "possible_elevated" & marker %in% c("VitB1","VitB2") , threshold := "Insufficient"] |>
  _[order(marker,threshold)] |>
  _[ , .(marker,threshold,`GroupNDNS Over-weight`,`GroupNDNS Obese`,`GroupGOOS`,`GroupNNUH`)]

propresultsCI <- propmods |> 
  lapply(broom::tidy,exponentiate=TRUE, conf.int=TRUE) |> 
  rbindlist(idcol="outcome") |> 
  _[ , string := sprintf("%0.2f (%0.2f-%0.2f) p=%0.3f %s",estimate,conf.low,conf.high,p.value,stars(p.value))][] |>
  dcast( outcome ~ term , value.var="string" ) |>
  
  _[ , marker := stringr::word(outcome,1,sep=fixed("_"))] |>
  _[ marker=="v25OHD", marker := "v25OHD_Std"] |> 
  _[ marker %in% order_of_outputs, marker := factor(marker, levels=order_of_outputs)] |> 
  _[ outcome %like% "deficient" , threshold := "Deficient"] |>
  _[ outcome %like% "possible_deficient" , threshold := "Insufficient"] |>
  _[ outcome %like% "elevated" , threshold := "Elevated"] |>
  _[ outcome %like% "elevated" & marker %in% c("VitB1","VitB2") , threshold := "Deficient"] |>
  _[ outcome %like% "possible_elevated" & marker %in% c("VitB1","VitB2") , threshold := "Insufficient"] |>
  _[order(marker,threshold)] |>
  _[ , .(labelsvector[marker],threshold,`GroupNDNS Over-weight`,`GroupNDNS Obese`,`GroupGOOS`,`GroupNNUH`)]

 propresultsCI |> 
   gt() |>
   tab_header("Odds ratios, confidence intervals and p-values comparing overweight groups with NDNS healthy") |> 
   gtsave("finalFigures/odds_ratios_ci_p_values.docx")

# propresults |> 
#   gt() |>
#   tab_header("Odds ratios, standard errors and p-values comparing overweight groups with NDNS healthy") |> 
#   gtsave("finalFigures/odds_ratios_p_values.docx")

#propresults |> fwrite(file="propresults.csv")

propresults


## For the abstract
propmods$anyDeficiency2 |> emmeans::emmeans( trt.vs.ctrl ~ Group, 
                                             adjust="none", 
                                             type="response", infer=c(TRUE,TRUE),
                                             data=propmods$anyDeficiency2$survey.design$variables,
                                             at=list(age=12,sex="Male",ethgr2="White"))


