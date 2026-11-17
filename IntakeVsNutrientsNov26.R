


### This file supposes that we have the combinedData file already made and in memory?
### So it should be called as an analysis file after the main file has been run.

#### We load in the food intake data separately in this file, then create the 

if(!exists("combinedData")){
  source("U://work/qib projects/mithen - food and health programme/emmawebb/goos/main_analysis_file_June_22.R")
}


format_cor <- function(x, y, method = "pearson", digits = 3) {
  # Run the correlation test
  test_result <- cor.test(x, y, method = method)
  
  # Extract components
  rho <- round(as.numeric(test_result$estimate), digits)
  p_val <- test_result$p.value
  
  # Format p-value to avoid looking like "0" if it's highly significant
  if (p_val < 0.001) {
    p_str <- "< 0.001"
  } else {
    p_str <- paste0("= ", round(p_val, digits))
  }
  
  # Return the formatted string: rho (p)
  return(paste0(rho, " (p ", p_str, ")"))
}


## Now load the person level data:
personlevel <- dir("c:/data/UKDA-6533-tab/tab/",full.names = TRUE) |> 
  grep(pattern="personleveldiet", value=TRUE)
personleveldat <- lapply(personlevel, fread) 

namesOfNutrients <- c("seriali",
  "Energykcal",
  "Retinolµg","Totalcaroteneµg","Alphacaroteneµg",
  "Betacaroteneµg","Betacryptoxanthinµg","VitaminAretinolequivalentsµg",
  "VitaminDµg","VitaminEmg","Thiaminmg",
  "Riboflavinmg","Niacinequivalentmg","VitaminB6mg",
  "VitaminB12µg","Folateµg","Pantothenicacidmg",
  "Biotinµg","VitaminCmg", "Ironmg",
  "Haemironmg","Nonhaemironmg","Zincmg","Seleniumµg")

labelsOfNutrients <- c("Energykcal"="Energy (kcal)",
                       "Retinolµg"="Retinol (µg)",
                       "Totalcaroteneµg"="Total carotene (µg)",
                       "Alphacaroteneµg"="Alpha-carotene (µg)",
                       "Betacaroteneµg"="Beta-carotene (µg)",
                       "Betacryptoxanthinµg"="Beta-cryptoxanthin (µg)",
                       "VitaminAretinolequivalentsµg"="Vitamin A retinol equivalents (µg)",
                       "VitaminDµg"="Vitamin D (µg)",
                       "VitaminEmg"="Vitamin E (mg)", # Fixed unit typo from original snippet
                       "Thiaminmg"="Thiamin (mg)",
                       "Riboflavinmg"="Riboflavin (mg)",
                       "Niacinequivalentmg"="Niacin equivalent (mg)",
                       "VitaminB6mg"="Vitamin B6 (mg)",
                       "VitaminB12µg"="Vitamin B12 (µg)",
                       "Folateµg"="Folate (µg)",
                       "Pantothenicacidmg"="Pantothenic acid (mg)",
                       "Biotinµg"="Biotin (µg)",
                       "VitaminCmg"="Vitamin C (mg)",
                       "Ironmg"="Iron (mg)",
                       "Haemironmg"="Haem iron (mg)",
                       "Nonhaemironmg"="Non-haem iron (mg)",
                       "Zincmg"="Zinc (mg)",
                       "Seleniumµg"="Selenium (µg)")


namesOfNutrients5 <- c("seriali",paste0(namesOfNutrients[-1],"_UI"))
nutrientList <- setdiff(namesOfNutrients,"seriali")

personLevelNutrients <- rbindlist(
  list(
    yr14 = personleveldat[[1]][ , .SD , .SDcols = namesOfNutrients],
    yr56 = personleveldat[[2]][ , .SD , .SDcols = namesOfNutrients],
    yr78 = personleveldat[[3]][ , .SD , .SDcols = namesOfNutrients],
    yr911 = personleveldat[[4]][ , .SD , .SDcols = namesOfNutrients],
    yr1215 = personleveldat[[5]][ , .SD , .SDcols = namesOfNutrients5] |> setnames(old = namesOfNutrients5 , new = namesOfNutrients)
    ),
  idcol = "wave"
)


ndns_nutrients <- personLevelNutrients[combinedData , on="seriali"]
ndns_nutrients <- ndns_nutrients[Source=="NDNS"]
#ndns_nutrients[ , .(wave,seriali,i.Energykcal,Energykcal)]
ndns_nutrients[ , EnergykcalResid := resid(lm(Energykcal ~ Sex+age+I(age^2)))]
ndns_nutrients[ , bmicat218 := factor(bmicat218,c("Normal-weight","Over-weight","Obese"))]


### For each nutrient, we need:

## Distribution of nutrient intake over weight (bmi_sds) and age.
## Adjusted for energy intake?

### Aim:
##
## Understand how reported micronutrient status is related to serum levels of corresponding micronutrient.
## Understand how intake is related to BMI_SDS, with and without correction for weight?

ggplot(ndns_nutrients , aes(x=Group, y=Energykcal)) + geom_boxplot(fill="red") + facet_wrap(~ageband, ncol=1) + 
  coord_flip()+ theme_bw() + ggtitle("Energy intake by age group and weight status") + 
  scale_x_discrete(labels=groupLabelsVector)

ggplot(ndns_nutrients , aes(x=Group, y=EnergykcalResid)) + geom_boxplot(fill="red") + facet_wrap(~ageband, ncol=1) +  coord_flip() + 
  theme_bw() + ggtitle("Residualised energy intake by age group and weight status","(Residualised by age and sex)")+ 
  scale_x_discrete(labels=groupLabelsVector)

## There is no apparent effect of BMI SDS group on reported energy intake.
lm(data=ndns_nutrients,Energykcal ~ Sex+age+I(age^2)+bmicat218+wave) |> summary()


## Energy intake goes up with age, but not with BMI SDS

ggplot(ndns_nutrients , aes(x=age, y=Energykcal)) + 
  geom_point(position = position_jitter(width=0.1)) + 
  stat_smooth() + 
  facet_wrap(~Sex) + 
  theme_bw() + 
  ggtitle("Reported energy intake by age and sex")

ggplot(ndns_nutrients , aes(x=bmi_sds, y=EnergykcalResid)) + 
  geom_point() + 
  theme_bw() + 
  stat_smooth(method="lm") + 
  ggtitle("Residualised reported energy intake vs BMI SDS","Residualised for age, age^2 and sex")


## This is the absolute values of nutrients by BMI_SDS
ndns_nutrients[ , .SD , .SDcols = c("seriali","age","Sex","bmi_sds","bmicat218",nutrientList)] |> 
  melt( id.vars = c("seriali","age","Sex","bmi_sds","bmicat218","Energykcal")  )  |>
  ggplot() + 
  aes(x=bmi_sds , y=value) + 
  geom_point(alpha=0.5) + 
  facet_wrap(~variable, scale="free_y",labeller=labeller(variable=labelsOfNutrients)) + 
  scale_y_log10() + 
  theme_bw() + 
  stat_smooth() + 
  labs(y="Intake",x="BMI SDS")
  ggpubr::stat_cor() + 
  ggtitle("Reported micronutrient intake by BMI SDS (crude correlation)")

## Residualised nutrient values per energy intake
ndns_nutrients[ , .SD , .SDcols = c("seriali","age","Sex","bmi_sds","bmicat218",nutrientList)] |> 
  melt( id.vars = c("seriali","age","Sex","bmi_sds","bmicat218","Energykcal")  )  |>
  subset(value>0) |>
  ggplot() + 
  aes(x=bmi_sds , y=resid(lm(log(value)~Energykcal))) + 
  facet_wrap(~variable, scale="free_y",labeller=labeller(variable=labelsOfNutrients)) + 
  geom_point() + 
  labs(y="log-intake (residualised)",x="BMI SDS")+
  theme_bw() + 
  stat_smooth(method="gam") + 
  ggpubr::stat_cor()+
  ggtitle("Reported micronutrient intake (residualised to Energy intake) by BMI SDS")

## Figure S5 - normalised intake values per energy intake
ndns_nutrients[ , .SD , .SDcols = c("seriali","age","Sex","bmi_sds","bmicat218",nutrientList)] |> 
  melt( id.vars = c("seriali","age","Sex","bmi_sds","bmicat218","Energykcal")  )  |>
  subset(value>0) |>
  ggplot() + 
  aes(x=bmi_sds , y= 1000*value/Energykcal) + 
  geom_point(col="darkgrey") + 
  facet_wrap(~variable, scale="free_y",labeller=labeller(variable=labelsOfNutrients)) + 
  labs(y="Intake per 1000 kcal",x="BMI SDS")+
  scale_y_log10() +
  theme_bw() + 
  stat_smooth(method="gam") + 
  ggpubr::stat_cor(method="spearman")+
  ggtitle("Reported micronutrient intake (normalised to Energy intake) by BMI SDS","Spearman correlation (p-value)")
ggsave("ReportedMicronutrientDensityBySDS_S5.png",width=14,height=12,dpi="retina")


### Models for the relationship between BMI SDS and normalised nutrient intake
combined_design_nutrients <- svydesign(ids=~PSU,weights=~finalweight, 
                                       strata=~strata,data=copy(ndns_nutrients))
## no evidence for an association between energy intake and bmi
svyglm( formula = Energykcal ~ bmi_sds + age + sex + ethgr5 , design = combined_design_nutrients) |> summary()

svyglm( formula = log(VitaminB12µg/Energykcal) ~ bmi_sds + age + sex + ethgr5 , design = combined_design_nutrients) |> summary()

s="VitaminB12µg"
nutrient_bmi_mods <- lapply(nutrientList[-1] , \(s){
  svyglm( formula = as.formula(sprintf("log(%s/Energykcal) ~ bmi_sds + age + sex + ethgr5",s)) , 
          design = svydesign(ids=~PSU,weights=~finalweight, 
                             strata=~strata,data=copy(ndns_nutrients[is.finite(get(s)) & get(s)>0])))
}) |> setNames(nutrientList[-1])

Nutrient_BMI_Model_Table <- nutrient_bmi_mods |> lapply(broom::tidy,conf.int=TRUE, exp=TRUE) |> 
  rbindlist(idcol="Outcome") |> 
  subset(term %like% "bmi")
Nutrient_BMI_Model_Table[ , q.value := p.adjust(p.value,method="fdr")]
Nutrient_BMI_Model_Table[ , resultsString := sprintf("%0.2f (%0.2f - %0.2f) p=%0.3f ; q=%0.3f" , estimate,conf.low,conf.high,p.value,q.value ) ]
Nutrient_BMI_Model_Table <- Nutrient_BMI_Model_Table[ , .(labelsOfNutrients[Outcome],resultsString)]
gt(Nutrient_BMI_Model_Table) |> gtsave("oct26figs/nutrient_bmi_table.docx")


ndns_nutrients[ , .SD , .SDcols = c("seriali","age","Sex","bmi_sds","Group",nutrientList)] |> 
  melt( id.vars = c("seriali","age","Sex","bmi_sds","Group","Energykcal")  )  |>
  subset(value>0) |>
  ggplot() + 
  aes(x=Group , y= 1000*value/Energykcal) + 
  geom_boxplot() + 
  facet_wrap(~variable, scale="free_y") + 
  scale_y_log10() +
  scale_x_discrete(labels=groupLabelsVector)+
  theme_bw() + 
  stat_smooth() + 
  ggtitle("Reported micronutrient intake (normalised to Energy intake) by BMI SDS Group (boxplot)") + 
  theme(axis.text.x = element_text(angle=30,hjust=1))
ggsave("ReportedMicronutrientDensityByGroupsBoxplot.png",width=14,height=12,dpi="retina")

naturalMods3[[1]] |> performance::check_model()

### Correlations are miniscule between BMI SDS and reported nutrient intake.  We have no way of 
ndns_nutrients[ , .SD , .SDcols = c("seriali","age","Sex","bmi_sds","Group",nutrientList)] |> 
  melt( id.vars = c("seriali","age","Sex","bmi_sds","Group","Energykcal")  )  |>
  subset(value>0) |>
  ggplot() + 
  aes(x=Group , y= 1000*value/Energykcal) + 
  stat_summary() + 
  facet_wrap(~variable, scale="free_y") + 
  scale_y_log10() +
  scale_x_discrete(labels=groupLabelsVector)+
  theme_bw() + 
  stat_smooth() + 
  ggtitle("Reported micronutrient intake (normalised to Energy intake) by BMI SDS Group (boxplot)") + 
  theme(axis.text.x = element_text(angle=30,hjust=1))
ggsave("oct26figs/ReportedMicronutrientDensityByGroupsMeanSE.png",width=14,height=12,dpi="retina")




ndns_nutrients[ , .SD , .SDcols = c("seriali","age","Sex","bmi_sds","bmicat218",nutrientList)] |> 
  melt( id.vars = c("seriali","age","Sex","bmi_sds","bmicat218","Energykcal")  )  |>
  _[  variable != "Energykcal"] |> 
  _[ , .(corSimple=cor(bmi_sds,value,method="spearman",use="pairwise"),
           cor=format_cor(bmi_sds,value,method="spearman"),
         cor_density=format_cor(bmi_sds,value/Energykcal,method="spearman")) , 
         by=variable][order(corSimple)][, corSimple := NULL] |>
  _[,variable := labelsOfNutrients[as.character(variable)]][]

ndns_nutrients[ , .SD , .SDcols = c("seriali","age","Sex","bmi_sds","bmicat218",nutrientList)] |> 
  melt( id.vars = c("seriali","age","Sex","bmi_sds","bmicat218","Energykcal")  )  |>
  _[  variable != "Energykcal"] |> 
  _[ , .(pval=cor.test(bmi_sds,value,method="spearman",use="pairwise")$p.value),
         by=variable][ , padjust:=p.adjust(pval,method="fdr")][]
ndns_nutrients[ , .SD , .SDcols = c("seriali","age","Sex","bmi_sds","bmicat218",nutrientList)] |> 
  melt( id.vars = c("seriali","age","Sex","bmi_sds","bmicat218","Energykcal")  )  |>
  _[  variable != "Energykcal"] |> 
  _[ , .(pval=cor.test(bmi_sds,value/Energykcal,method="spearman",use="pairwise")$p.value),
     by=variable][ , padjust:=p.adjust(pval,method="fdr")][padjust<0.05]




## Standardised intakes by BMI SDS
ndns_nutrients[ , .SD , .SDcols = c("seriali","age","Sex","bmi_sds","bmicat218",nutrientList)] |> 
  melt( id.vars = c("seriali","age","Sex","bmi_sds","bmicat218","Energykcal")  )  |>
  subset(value>0) |>
  ggplot() + 
  aes(x=bmicat218 , y= 1000*value/Energykcal) + 
  stat_summary() + 
  facet_wrap(~variable, scale="free_y") + 
  scale_y_log10() +
  theme_bw() + 
  stat_smooth() + 
  ggtitle("Reported micronutrient intake (normalised to Energy intake) by BMI SDS Group (mean, se)")

stat_test <- ndns_nutrients[ , .SD , .SDcols = c("seriali","age","Sex","bmi_sds","Group",nutrientList)] |> 
  melt( id.vars = c("seriali","age","Sex","bmi_sds","Group","Energykcal") ,variable.name="Nutrient" )  |>
  subset(value>0) |>
  dplyr::group_by(Nutrient) |>
  rstatix::wilcox_test(value ~ Group) |>
  rstatix::adjust_pvalue(method = "bonferroni") |>
  # CRUCIAL: Tell rstatix that scales vary across the groups
  rstatix::add_xy_position(x = "Group", fun = "max", scales = "free")

stat_test$groups <- NULL

library(rstatix)

## Intakes by BMI SDS group
ndns_nutrients[ , .SD , .SDcols = c("seriali","age","Sex","bmi_sds","Group",nutrientList)] |> 
  melt( id.vars = c("seriali","age","Sex","bmi_sds","Group","Energykcal")  )  |>
  subset(value>0) |>
  ggplot() + 
  aes(x=Group , y= value) + 
  stat_summary() +
  facet_wrap(~variable, scale="free_y") + 
  #scale_y_log10() +
  theme_bw() + 
  ggpubr::stat_pvalue_manual(data=stat_test, label = "p.adj", tip.length = 0.01) +
  ggtitle("Reported micronutrient intake by BMI SDS Group (mean, se)") + 
  theme(axis.text.x=element_text(angle=30,hjust=1))


### Now check the relationship between nutrient intakes and serum nutrients

pairList <- fread("pairedNutrients.txt",fill = TRUE,sep=",",header = TRUE)
pairListLong <- melt(pairList,id.vars = "NDNSfoodintake") |> subset(value!="")

mapply(FUN=\(X,Y) with(ndns_nutrients,c(X,Y,cor(get(X),get(Y),use = "pairwise"))), 
       X=pairListLong$NDNSfoodintake , Y=pairListLong$value ) |>
  t()

corGraphs <- mapply(FUN=\(X,Y) {
  ggplot(ndns_nutrients[`VitaminB12µg`/Energykcal>(1e-4) & `VitaminB12µg`/Energykcal<(.1)]) + 
    aes(x=get(X),y=get(Y)) +
    labs(x=X,y=Y) +  
    theme_bw() + 
    geom_point(col="grey") + 
    ggpubr::stat_cor(method="spearman") +
    stat_smooth(method="gam")
  }, 
       X=pairListLong$NDNSfoodintake , Y=pairListLong$value )

corGraphs |> wrap_plots() & scale_y_log10() & scale_x_log10() 


getListOfIntakeVars <- function(outcome){
  pairListLong[value==outcome,NDNSfoodintake] |> paste(collapse="+")
}

getListOfIntakeVars("Se")
getListOfIntakeVars("FerritinCRP")

corGraphsE <- mapply(FUN=\(X,Y) {
  ggplot(ndns_nutrients[(`VitaminB12µg`/Energykcal)>(1e-4) & (`VitaminB12µg`/Energykcal)<(.1)]) + 
    aes(x=1000*get(X)/Energykcal,y=get(Y)) +
    labs(x=labelsOfNutrients[X],y=labelsvector3b[Y]) +  
    theme_bw() + 
    geom_point(col="grey") + 
    ggpubr::stat_cor(method="spearman") +
    stat_smooth(method="lm")
}, 
X=pairListLong$NDNSfoodintake , Y=pairListLong$value )

corGraphsE |> wrap_plots()  + 
  plot_annotation("Nutrient density (per 1000kcal) vs corresponding biochemical markers",
                  subtitle = "Spearman's correlation (p-value); remember Vit B1 and Vit B2 are measured in reverse")& 
  scale_y_log10() & 
  scale_x_log10()

ggsave("oct26figs/IntakeDensityVsBiochemistry.png",width=14,height=12,dpi="retina")

### Next up, can you explain any of the adiposity effect by intake?

## For each model:

allMods <- mapply(
  FUN=\(X,Y) {
    form1 <- as.formula( sprintf("log(%s) ~ Group + age + Sex " , Y  ))
    form2 <- as.formula( sprintf("log(%s) ~ Group + log(%s)  + age + Sex " , Y,X  ))
    form3 <- as.formula( sprintf("log(%s) ~ Group + I(log(1000*%s/Energykcal)) + age + Sex " , Y,X  ))
    
    tempdat <- ndns_nutrients[get(Y)>0 & is.finite(get(Y)) & get(X)>0 & is.finite(get(X))]
    
    Mod1 <- lm(data=tempdat , form1 ) |> summary()
    Mod2 <- lm(data=tempdat , form2 ) |> summary()
    Mod3 <- lm(data=tempdat , form3 ) |> summary()

    list(Mod1,Mod2,Mod3) |> setNames(c(Y,paste(Y,"incl. intake"),paste(Y,"(incl. density)")))
    
    },
  X=pairListLong$NDNSfoodintake , 
  Y=pairListLong$value ,
  SIMPLIFY = FALSE
) |> setNames(sprintf("%s vs %s",pairListLong$NDNSfoodintake,pairListLong$value))
allMods[[1]]
names(allMods)


for(i in seq_along(pairListLong$value)){
  print(i)
  file_path <- paste0("intake_model_outputs/", pairListLong$NDNSfoodintake[i],pairListLong$value[i], "_intake_mod.html")
  modelsummary::modelsummary(allMods[[i]][1:3],stars = TRUE,
                             drop=c("age","SexMale","(Intercept)"),output = file_path)
}



## This is the main block that makes the obesity model in NDNS only
listOfNutrientMods <- pairListLong$value |> unique() |> setdiff("LEUT")
## Need to remake combineddesign since I added new variables.
combined_design_nutrients <- svydesign(ids=~PSU,weights=~finalweight, 
                             strata=~strata,data=copy(ndns_nutrients))
naturalMods_intake <- lapply(listOfNutrientMods, function(y) {
  tryCatch({
    print(y)
    formula <- as.formula(paste("log(" , y, ") ~ Group + ns(age,df=4) + sex + ethgr5 +" , getListOfIntakeVars(y)))
    print(formula)
    vars <- all.vars(formula)
    vars_dt <- combined_design_nutrients$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design_nutrients , cc &  Source %in% c("NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(listOfNutrientMods)

obesityGraphs_intakeAdj <- lapply(listOfNutrientMods , makeObesityGraph2,naturalMods_intake) |> 
  setNames(listOfNutrientMods)
#obesityGraphs_intakeAdj <- wrap_plots(obesityGraphs_intakeAdj)
meansOverObeseGraph_nutrientAdj <- wrap_plots(obesityGraphs_intakeAdj[intersect(order_of_outputs,listOfNutrientMods)],ncol = 5) + 
  plot_layout(axis_titles = "collect") & 
  scale_x_discrete(labels=c("NDNS <85th centile","NDNS 85th-95th centile","NDNS >95th centile","GOOS")) & 
  scale_y_continuous(expand = expansion(0.1)) &
  theme(axis.text.x=element_text(angle=20,hjust=1)) &
  theme(plot.margin = margin(l = 20,t=10,r=10,b = 10, unit = "pt"))
ggsave(plot = meansOverObeseGraph_nutrientAdj,
       filename = file.path("oct26figs","obeseGraph3_main_nutrientAdj.png"),
       width=14,height=14,dpi="retina")

## This is the list of models just with NDNS 
#naturalMods3

## This is the list of models just with NDNS and the intake variables adjusted for.
## (The sensitivity analysis - how much effect is explained by intake)
#naturalMods_intake

#outcome = listOfNutrientMods[[1]]
tableOfEffects <- lapply(listOfNutrientMods,\(l){
  list(natural3 = naturalMods3[[l]] |> broom::tidy(exponentiate=TRUE, conf.int=TRUE) ,
  withintake = naturalMods_intake[[l]] |> broom::tidy(exponentiate=TRUE, conf.int=TRUE)) |>
    rbindlist(idcol = "type") |> subset(term %like% "Group")
}) |> setNames(listOfNutrientMods) |> rbindlist(idcol="outcome")

tableOfEffects[ , stringresults := sprintf("%0.2f (%0.2f) p=%0.3f",estimate,std.error,p.value)]
tableOfEffects[ , .(stringresults,outcome=labelsvector[outcome],type,term)] |> 
  dcast(outcome ~ term + type, value.var = "stringresults") |> 
  gt() |> gtsave(filename = "oct26figs/intake_moderation_table.docx")



tableOfEffects2 <- lapply(order_of_outputs2,\(l){
  naturalMods2[[l]] |> broom::tidy(exponentiate=TRUE, conf.int=TRUE)
  }) |> setNames(order_of_outputs2) |> 
  rbindlist(idcol = "outcome") |> subset(term %like% "Group")

tableOfEffects2[ , q.value := p.adjust(p.value,"fdr")]
tableOfEffects2[ , stringresults := sprintf("%0.2f (%0.2f) p=%0.3f; q=%0.3f",estimate,std.error,p.value,q.value)]
tableOfEffects2[ , stringresults2 := sprintf("%0.2f (%0.2f-%0.2f) p=%0.3f; q=%0.3f",estimate,conf.low,conf.high,p.value,q.value)]

tableOfEffects2[ , stringresults := gsub("=0.000","<0.001",stringresults)]
tableOfEffects2[ , stringresults2 := gsub("=0.000","<0.001",stringresults2)]


# tableOfEffects2[ , .(stringresults,outcome,term)] |> dcast(outcome ~ term , value.var = "stringresults") |>
#   _[, .(labelsvector[outcome],`GroupNDNS Over-weight`,`GroupNDNS Obese`,`GroupGOOS`,
#         p.value=sprintf("%0.3f",bmi_effectsNDNSonlySPLINE[outcome]),
#         q.value=sprintf("%0.3f",p.adjust(bmi_effectsNDNSonlySPLINE[outcome],"fdr"))
#    )]

tableOfEffects2[ , .(stringresults,outcome,term)] |> dcast(outcome ~ term , value.var = "stringresults") |>
  _[, .(labelsvector[outcome],`GroupNDNS Over-weight`,`GroupNDNS Obese`,`GroupGOOS`,
        `trend p.value`=sprintf("%0.3f",bmi_effectsNDNSonlySPLINE[outcome]),
        `trend q.value`=sprintf("%0.3f",p.adjust(bmi_effectsNDNSonlySPLINE[outcome],"fdr"))
  )] |> gt() |> gtsave(filename="oct26figs/fdr_table.docx")

tableOfEffects2[ , .(stringresults2,outcome,term)] |> dcast(outcome ~ term , value.var = "stringresults2") |>
  _[, .(labelsvector[outcome],`GroupNDNS Over-weight`,`GroupNDNS Obese`,`GroupGOOS`,
        `trend p.value`=sprintf("%0.3f",bmi_effectsNDNSonlySPLINE[outcome]),
        `trend q.value`=sprintf("%0.3f",p.adjust(bmi_effectsNDNSonlySPLINE[outcome],"fdr"))
  )] |> gt() |> gtsave(filename="oct26figs/fdr_table_ci.docx")


# bmi_effectsNDNSonlySPLINE
# 
# 
# contMods2NDNSonlySPLINE |> lapply(broom::tidy) |> rbindlist(outcome)
# 
# 
# 
