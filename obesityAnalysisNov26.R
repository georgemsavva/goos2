
## Restate the svy design at the start of this file
combined_design <- svydesign(ids=~PSU,weights=~finalweight, 
                             strata=~strata,data=copy(combinedData))

###################################################################################
################### MAKE THE MAIN SVY GLM REGRESSION MODELS #######################
###################################################################################

#### natural mods 2 doesn't include survey year as a covariate because its not in GOOS.
#### the GOOS samples were collected years, but analysed recently, so it's not clear what the 'year' should be.
#### So we could assume that the survey year is various different things.
nsknots <- quantile(combinedData$age,c(0.25,.5,.75))


### These models include BMI as a continuous variable

contMods2 <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ bmi_sds + ns(age,df=4) + sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("GOOS","NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)
bmi_effects <- contMods2 |> lapply(broom::tidy) |> rbindlist(idcol = "outcome") |> subset(term=="bmi_sds")
bmi_effects[ , qval:=p.adjust(p.value, method="BH")][order(p.value)]


### We only want to run this in NDNS since the GOOS group and BMI is confounded, and GOOS may not be as comparable as we'd like.

contMods2NDNSonly <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ bmi_sds + ns(age,df=4) + sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("NDNS") & abs(bmi_sds)<4)
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)
bmi_effectsNDNSonly <- contMods2NDNSonly |> lapply(broom::tidy) |> rbindlist(idcol = "outcome") |> subset(term=="bmi_sds")
bmi_effectsNDNSonly[ , qval:=p.adjust(p.value, method="BH")][order(p.value)]

combinedData[Source %in% c("NDNS") & abs(bmi_sds)<4, .SD, .SDcols = c("bmi_sds","Group",presentInNDNS)] |> 
  melt(id.vars=c("bmi_sds","Group"), variable.name = "outcome") |> 
  ggplot(aes(x=bmi_sds,y=value)) + facet_wrap(~outcome,scale="free_y") + geom_point() + 
  scale_y_log10() + 
  stat_smooth(method="gam",se=FALSE) + 
  stat_smooth(method="lm",se = FALSE,col="red")

##
## We can use these sets of model to describe the variance explained by the bmi term.
##

contMods2NDNSonlySPLINE <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ ns(bmi_sds,df=4) + ns(age,df=4) + sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)
contMods2NDNSonlyNOBMI <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ ns(age,df=4) + sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("NDNS") & !is.na(bmi_sds))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)

####
#### Estimating partial R^2 associated with BMI SDS for each of the biomarkers
####

wrss <- function(m) sum(m$prior.weights * residuals(m, type = "response")^2)
partial_R2 <- (wrss(contMods2NDNSonlyNOBMI[[1]]) - 
               wrss(contMods2NDNSonlySPLINE[[1]])) / wrss(contMods2NDNSonlyNOBMI[[1]])
sapply(presentInNDNS , \(s) (wrss(contMods2NDNSonlyNOBMI[[s]]) - wrss(contMods2NDNSonlySPLINE[[s]])) / wrss(contMods2NDNSonlyNOBMI[[s]])) |>
  cbind()


#### 
#### Calculate the p-values and the q-values for the effect of "BMI" as a predictor of the biomarker (spline)
####

bmi_effectsNDNSonlySPLINE <- contMods2NDNSonlySPLINE |> 
  lapply(regTermTest,~ns(bmi_sds,df=4)) |> 
  lapply(\(x) x[["p"]][1,1]) |> unlist() |> sort()
cbind(
  pval=bmi_effectsNDNSonlySPLINE,
  qval=p.adjust(bmi_effectsNDNSonlySPLINE, method="fdr")) |> 
  as.data.frame() |> 
  data.table(keep.rownames = "outcome") |> gt()


### These are the main regression models linking BMI group to biomarkers in NDNS and GOOS
### Adjusted for ethnicity, age (spline) and sex.

naturalMods2 <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ Group + ns(age,df=4) + sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("GOOS","NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)

linearMods2 <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste(y, " ~ Group + ns(age,df=4) + sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("GOOS","NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)

### Now the model without adjustment for ethnicity

naturalMods2_sens_noethnicity <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ Group + ns(age,df=4) + sex"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("GOOS","NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)

### Now the model adjusted for income tertile (with NDNS and GOOS)
### Evidently I have make eqv3F into a factor which has mappings of the SES variables across all of the different groups.

naturalMods_sens_income <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ Group + ns(age,df=4) + eqv3F +
                                sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("GOOS","NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)

naturalMods_sens_IMD <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ Group + ns(age,df=4) + imd5F +
                                sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("GOOS","NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)

### Models with income but not with ethnicity

naturalMods_sens_income_no_ethnicity <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ Group + ns(age,df=4) + eqv3F +
                                sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("GOOS","NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)

##
###
#### Sensitivity analysis, assume GOOS is survey year 2,8,14?
###  In these cases, apply a spline for the survey year effect
##

### Remove possibly null models.
naturalMods2 <- naturalMods2[!sapply(naturalMods2,is.null)] # no models are null
### Output the diagnostic plots - observe no real problems except heavy tails everywhere.
for(i in 1:length(naturalMods2)){
  png(sprintf("diagnostics/check_model_%s.png",names(naturalMods2)[i]),width=1000,height=1000)
  par(mfrow=c(2,2))
  plot(naturalMods2[[i]])
  dev.off()
}

##
## These are NDNS only models.
##

naturalMods3 <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ Group + ns(age,df=4) + sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc &  Source %in% c("NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)



### This is the main output for the graphical comparison of average concentrations over obesity
### naturalMods2


######## THE CODE BELOW IS FOR THE GRAPHS WITHOUT THE FDR CORRECTION.

######## THE FDR CORRECTED TABLES ARE IN THE FILE WITH THE NUTRIENT ADJUSTMENT ANALYSIS (IntakeVsNutrients.R)
######## THE FDR CORRECTED GRAPHS ARE IN MakeLargeModelsFigure.R
######## BECAUSE THERE IS WHERE I SET UP THE GRAPHS WITH THE Q-VALUE ADJUSTMENTS

if(FALSE){

obesityGraphs3 <- lapply(presentInNDNS , makeObesityGraph2,naturalMods2) |> 
  setNames(presentInNDNS)
meansOverObeseGraph3 <- wrap_plots(obesityGraphs3[order_of_outputs],ncol = 5) + 
  plot_layout(axis_titles = "collect") & 
  scale_x_discrete(labels=c("NDNS <85th centile","NDNS 85th-95th centile","NDNS >95th centile","GOOS")) & 
  scale_y_continuous(expand = expansion(0.1)) &
  theme(axis.text.x=element_text(angle=20,hjust=1)) &
  theme(plot.margin = margin(l = 20,t=10,r=10,b = 10, unit = "pt"))
ggsave(plot = meansOverObeseGraph3,
       filename = file.path("june26figs","obeseGraph3_main.png"),
       width=14,height=14,dpi="retina")

#################################################
# SENSITIVITY ANALYSIS GRAPHS IF YOU WANT THEM  #
#################################################

### Income adjustment main sensitivity analysis

obesityGraphs3_sens_inc <- lapply(presentInNDNS , makeObesityGraph2,naturalMods_sens_income) |> 
  setNames(presentInNDNS)
meansOverObeseGraph3_sens_income <- wrap_plots(obesityGraphs3_sens_inc[order_of_outputs],ncol = 5) + 
  plot_layout(axis_titles = "collect") & 
  scale_x_discrete(labels=c("NDNS <85th centile","NDNS 85th-95th centile","NDNS >95th centile","GOOS")) & 
  scale_y_continuous(expand = expansion(0.1)) &
  theme(axis.text.x=element_text(angle=20,hjust=1)) &
  theme(plot.margin = margin(l = 20,t=10,r=10,b = 10, unit = "pt"))
ggsave(plot = meansOverObeseGraph3_sens_income,
       filename = file.path("june26figs","obeseGraph3_sens_inc_.png"),
       width=14,height=14,dpi="retina")

### Income adjustment main sensitivity analysis - no ethnicity

obesityGraphs3_sens_inc_no_ethnicity <- lapply(presentInNDNS , makeObesityGraph2,naturalMods_sens_income_no_ethnicity) |> 
  setNames(presentInNDNS)
meansOverObeseGraph3_sens_income_no_ethnicity <- 
  wrap_plots(obesityGraphs3_sens_inc_no_ethnicity[order_of_outputs],ncol = 5) + 
  plot_layout(axis_titles = "collect") & 
  scale_x_discrete(labels=c("NDNS <85th centile","NDNS 85th-95th centile","NDNS >95th centile","GOOS")) & 
  scale_y_continuous(expand = expansion(0.1)) &
  theme(axis.text.x=element_text(angle=20,hjust=1)) &
  theme(plot.margin = margin(l = 20,t=10,r=10,b = 10, unit = "pt"))
ggsave(plot = meansOverObeseGraph3_sens_income_no_ethnicity,
       filename = file.path("june26figs","obeseGraph3_sens_inc_no_ethnicity.png"),
       width=14,height=14,dpi="retina")

### Income adjustment main sensitivity analysis - no ethnicity, no income

obesityGraphs3_sens_no_ethnicity <- 
  lapply(presentInNDNS , makeObesityGraph2,naturalMods2_sens_noethnicity) |> 
  setNames(presentInNDNS)
meansOverObeseGraph3_sens_no_ethnicity <- 
  wrap_plots(obesityGraphs3_sens_no_ethnicity[order_of_outputs],ncol = 5) + 
  plot_layout(axis_titles = "collect") & 
  scale_x_discrete(labels=c("NDNS <85th centile","NDNS 85th-95th centile","NDNS >95th centile","GOOS")) & 
  scale_y_continuous(expand = expansion(0.1)) &
  theme(axis.text.x=element_text(angle=20,hjust=1)) &
  theme(plot.margin = margin(l = 20,t=10,r=10,b = 10, unit = "pt"))
ggsave(plot = meansOverObeseGraph3_sens_no_ethnicity,
       filename = file.path("june26figs","obeseGraph3_sens_no_ethnicity.png"),
       width=14,height=14,dpi="retina")


combinedData[ , SurveyYear2:=ifelse(Source=="GOOS",2,SurveyYear)]
combined_design <- svydesign(ids=~PSU,weights=~finalweight, 
                             strata=~strata,data=copy(combinedData))

naturalMods_sens1 <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ Group + ns(age,df=4) + ns(SurveyYear2,df=4) +
                                sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("GOOS","NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)


## Pretend GOOS was recruited/analysed at survey year 14 

combinedData[ , SurveyYear2:=ifelse(Source=="GOOS",14,SurveyYear)]
combined_design <- svydesign(ids=~PSU,weights=~finalweight, 
                             strata=~strata,data=copy(combinedData))
naturalMods_sens_goos_14 <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ Group + ns(age,df=4) + ns(SurveyYear2,df=4) +
                                sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("GOOS","NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)



## Pretend GOOS was recruited/analysed at survey year 2 

combinedData[ , SurveyYear2:=ifelse(Source=="GOOS",2,SurveyYear)]
combined_design <- svydesign(ids=~PSU,weights=~finalweight, 
                             strata=~strata,data=copy(combinedData))
naturalMods_sens_goos_2 <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ Group + ns(age,df=4) + ns(SurveyYear2,df=4) +
                                sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("GOOS","NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)


## Pretend GOOS was recruited/analysed at survey year 8 

combinedData[ , SurveyYear2:=ifelse(Source=="GOOS",8,SurveyYear)]
combined_design <- svydesign(ids=~PSU,weights=~finalweight, 
                             strata=~strata,data=copy(combinedData))
naturalMods_sens_goos_8 <- lapply(presentInNDNS, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ Group + ns(age,df=4) + ns(SurveyYear2,df=4) +
                                sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("GOOS","NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(presentInNDNS)


### Suppose GOOS was collected in year 8

obesityGraphs3_sens_goos_8 <- 
  lapply(presentInNDNS , makeObesityGraph2,naturalMods_sens_goos_8) |> 
  setNames(presentInNDNS)
meansOverObeseGraph3_sens_goos_8 <- 
  wrap_plots(obesityGraphs3_sens_goos_8[order_of_outputs],ncol = 5) + 
  plot_layout(axis_titles = "collect") & 
  scale_x_discrete(labels=c("NDNS <85th centile","NDNS 85th-95th centile","NDNS >95th centile","GOOS")) & 
  scale_y_continuous(expand = expansion(0.1)) &
  theme(axis.text.x=element_text(angle=20,hjust=1)) &
  theme(plot.margin = margin(l = 20,t=10,r=10,b = 10, unit = "pt"))
ggsave(plot = meansOverObeseGraph3_sens_goos_8,
       filename = file.path("june26figs","obeseGraph3_sens_goos_8.png"),
       width=14,height=14,dpi="retina")

### Suppose GOOS was collected in year 14

obesityGraphs3_sens_goos_14 <- 
  lapply(presentInNDNS , makeObesityGraph2,naturalMods_sens_goos_14) |> 
  setNames(presentInNDNS)
meansOverObeseGraph3_sens_goos_14 <- 
  wrap_plots(obesityGraphs3_sens_goos_8[order_of_outputs],ncol = 5) + 
  plot_layout(axis_titles = "collect") & 
  scale_x_discrete(labels=c("NDNS <85th centile","NDNS 85th-95th centile","NDNS >95th centile","GOOS")) & 
  scale_y_continuous(expand = expansion(0.1)) &
  theme(axis.text.x=element_text(angle=20,hjust=1)) &
  theme(plot.margin = margin(l = 20,t=10,r=10,b = 10, unit = "pt"))
ggsave(plot = meansOverObeseGraph3_sens_goos_14,
       filename = file.path("june26figs","obeseGraph3_sens_goos_14.png"),
       width=14,height=14,dpi="retina")


#### Suppose GOOS was collected in year 2

obesityGraphs3_sens_goos_2 <- 
  lapply(presentInNDNS , makeObesityGraph2,naturalMods_sens_goos_2) |> 
  setNames(presentInNDNS)
meansOverObeseGraph3_sens_goos_2 <- 
  wrap_plots(obesityGraphs3_sens_goos_2[order_of_outputs],ncol = 5) + 
  plot_layout(axis_titles = "collect") & 
  scale_x_discrete(labels=c("NDNS <85th centile","NDNS 85th-95th centile","NDNS >95th centile","GOOS")) & 
  scale_y_continuous(expand = expansion(0.1)) &
  theme(axis.text.x=element_text(angle=20,hjust=1)) &
  theme(plot.margin = margin(l = 20,t=10,r=10,b = 10, unit = "pt"))
ggsave(plot = meansOverObeseGraph3_sens_goos_2,
       filename = file.path("june26figs","obeseGraph3_sens_goos_2.png"),
       width=14,height=14,dpi="retina")



}


MakeRegressionDiagnostics <- FALSE
if(MakeRegressionDiagnostics==TRUE){

## Regression diagnostics


naturalMods2_diag <- lapply(presentInNDNS, function(y) {
  tryCatch({
    #y=presentInNDNS[1]
    print(y)
    formula <- as.formula(paste("log(" , y, ") ~ Group + ns(age,df=4) + sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("GOOS","NDNS"))
    m=svyglm(formula, design =  design_cc, family = gaussian())
    design_cc_2 <- update(design_cc,fit=fitted(m),r=residuals(m,type="response"))
    
    png(filename=sprintf("oct26figs/diags/reg_diag_%s.png",y),width=800,height=400)
    par(mfrow=c(1,2))
    print("Got here")
    svyplot(r ~ fit, design=design_cc_2,main=sprintf("%s (fitted vs residual)",labelsvector[y]),ylab="Residual",xlab="Fitted")
    abline(h=0, lty=2)
    lines(svysmooth(r ~ fit, design = design_cc_2,bandwidth = .1), col = "red",lwd=2)
    
    svyqqmath(~r,design=design_cc_2,main=sprintf("%s (residual QQ plot)",labelsvector[y]))
    dev.off()    
    
    
  }, error = function(e) print(e))
}) 


linearMods2_diag <- lapply(presentInNDNS, function(y) {
  tryCatch({
    #y=presentInNDNS[1]
    print(y)
    formula <- as.formula(paste( y, " ~ Group + ns(age,df=4) + sex + ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design , cc & Source %in% c("GOOS","NDNS"))
    m=svyglm(formula, design =  design_cc, family = gaussian())
    design_cc_2 <- update(design_cc,fit=fitted(m),r=residuals(m,type="response"))
    
    png(filename=sprintf("oct26figs/diags/linear_diag_%s.png",y),width=800,height=400)
    par(mfrow=c(1,2))
    print("Got here")
    svyplot(r ~ fit, design=design_cc_2,main=sprintf("%s (fitted vs residual)",labelsvector[y]),ylab="Residual",xlab="Fitted")
    abline(h=0, lty=2)
    lines(svysmooth(r ~ fit, design = design_cc_2), col = "red",lwd=2)
    
    svyqqmath(~r,design=design_cc_2,main=sprintf("%s (residual QQ plot)",labelsvector[y]))
    dev.off()    
    
    
  }, error = function(e) print(e))
}) 

vars <- all.vars(VitC ~ Group + ns(age,df=4) + sex + ethgr5)
vars_dt <- combined_design$variables[, ..vars]
cc <- stats::complete.cases(vars_dt)
design_cc <- subset(combined_design , cc & Source %in% c("GOOS","NDNS"))

vitCnaturalScaleModel <- svyglm( VitC ~ Group + ns(age,df=4) + sex + ethgr5 , design=design_cc, family=gaussian())
emmeans(vitCnaturalScaleModel , trt.vs.ctrl~ Group, data=vitCnaturalScaleModel$data)
emmeans(naturalMods2$VitC , trt.vs.ctrl~ Group, data=naturalMods2$VitC$data, type="response")


}