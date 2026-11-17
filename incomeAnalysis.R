
#####  THIS FILE PERFORMS ANALYSIS TO DESCRIBE THE EFFECTS OF INCOME AND ETHNICITY ON OUTCOMES.

#####  WE USE THE NDNS DATA ONLY FOR MANY ANALYSES BECAUSE ITS EASIER TO UNDERSTAND THE SAMPLING AND ITS RELATIONSHIP WITH THE UK POPULATION.

##### INCOME ANALYSIS IS ONLY DONE IN THE NDNS GROUPS
##### IMD ANALYSIS CAN BE DONE ACROSS ALL GROUPS

combined_design <- svydesign(ids=~PSU,
                             weights=~finalweight, 
                             strata=~strata,
                             data=combinedData)
options(survey.lonely.psu="adjust")

ggplot(combinedData[!is.na(imd5F)], aes(x=imd5F,y=bmi_sds)) + 
  geom_violin() + 
  coord_flip()+ 
  theme_bw() + ggridges::geom_ridgeline_gradient()

# 1. Helper function for weighted median in base R
# From Gemini - I'm not using this for anything important!
wtd_median <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cum_w <- cumsum(w) / sum(w)
  x[which(cum_w >= 0.5)[1]]
}

# 2. Compute exact weighted medians per IMD group using data.table
medians_dt <- combinedData[Source=="NDNS" & !is.na(imd5F), .(
  wtd_med = wtd_median(bmi_sds, finalweight)
), by = imd5F]

ggplot(combinedData[!is.na(imd5F) & Source=="NDNS"], 
       aes(y=imd5F,x=bmi_sds,fill=imd5F)) + 
  theme_bw() + 
  ggridges::geom_density_ridges(alpha=1,
                                aes(weight=finalweight),
                                #quantile_lines=TRUE,
                                quantiles=2,
                                )  + 
  geom_vline(xintercept=c(-1,0,1,2),lty="dotted")+
  #geom_point(data=medians_dt,aes(x=wtd_med,y=imd5F),inherit.aes=FALSE,shape="|",size=8) + 
  coord_cartesian(xlim=c(-3,4.2))+
  labs(x="BMI SDS (UK90 Reference)",y="IMD Quintile (1–2 = Most Deprived)" )+
  scale_fill_brewer(palette = "YlOrRd", direction = -1)

# ggplot(combinedData[!is.na(imd5F) & Source=="NDNS"]) + 
#   aes(x=imd5F,fill=cut(bmi_sds,c(-Inf,-1,0,1,2,Inf)),y=1) + geom_col(position="fill")

combinedData[Source=="NDNS",sd(bmi_sds),.(SurveyYear,imd5F)] |> 
  ggplot(aes(x=SurveyYear,y=V1,col=imd5F)) + 
  geom_point() + 
  geom_line()
combinedData[Source=="NDNS",sd(bmi_sds, na.rm = TRUE),.(imd5F)] 

wtd_median <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cum_w <- cumsum(w) / sum(w)
  x[which.max(cum_w >= 0.5)]
}

diag_dt <- combinedData[Source=="NDNS"][!is.na(bmi_sds) & !is.na(imd5F), .(
  n_raw       = .N,
  eff_n       = (sum(finalweight)^2) / sum(finalweight^2),
  unwtd_med   = median(bmi_sds),
  wtd_med     = wtd_median(bmi_sds, finalweight),
  wtd_obese_p = sum(finalweight[bmi_sds >= 2]) / sum(finalweight) * 100
), by = imd5F][order(imd5F)]

print(diag_dt)

tails_dt <- combinedData[Source=="NDNS"][!is.na(bmi_sds) & !is.na(imd5F) & !is.na(finalweight), .(
  n_raw          = .N,
  n_underwtd     = sum(bmi_sds <= -2.0),
  underweight_p  = sum(finalweight[bmi_sds <= -2.0]) / sum(finalweight) * 100,
  low_bmi_p      = sum(finalweight[bmi_sds <= -1.33]) / sum(finalweight) * 100,
  healthy_p      = sum(finalweight[bmi_sds > -1.33 & bmi_sds < 1.33]) / sum(finalweight) * 100,
  overweight_p   = sum(finalweight[bmi_sds >= 1.33]) / sum(finalweight) * 100,
  obese_p        = sum(finalweight[bmi_sds >= 2.0]) / sum(finalweight) * 100
), by = imd5F][order(imd5F)]

## Calculate models
eqv_mods_year <- lapply(listOfOutcomes$v, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ 
                                eqv3F + 
                                factor(SurveyYear) + 
                                bs(age) + 
                                sex + 
                                ethgr5 + 
                                Group"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design, cc)
    design_cc <- subset(design_cc, Source == "NDNS")
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(listOfOutcomes$v)

## Calculate models for ethnicity
eth_mods <- lapply(listOfOutcomes$v, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~
                                eqv3F +
                                factor(SurveyYear) +
                                bs(age) +
                                sex +
                                ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design, cc)
    design_cc <- subset(design_cc, Source == "NDNS")
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(listOfOutcomes$v)

## Calculate models for ethnicity
eth_mods2 <- lapply(listOfOutcomes$v, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~
                                eqv3F +
                                factor(SurveyYear) +
                                bs(age) +
                                sex +
                                ethgr5 + 
                                Group"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design, cc)
    design_cc <- subset(design_cc, Source %in% c("NDNS"))
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(listOfOutcomes$v)

## Calculate models excluding BMI group as a covariate
eqv_mods_year_noGroup <- lapply(listOfOutcomes$v, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ 
                                eqv3F + 
                                factor(SurveyYear) + 
                                bs(age) + 
                                sex + 
                                ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design, cc)
    design_cc <- subset(design_cc, Source == "NDNS")
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) NULL)
}) |> setNames(listOfOutcomes$v)



## Calculate models excluding BMI group as a covariate
imd_mods_year_noGroup <- lapply(listOfOutcomes$v, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ 
                                imd5F + 
                                factor(SurveyYear) + 
                                bs(age) + 
                                sex + 
                                ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design, cc)
    design_cc <- subset(design_cc, Source == "NDNS")
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) e)
}) |> setNames(listOfOutcomes$v)

imd_mods_continuous <- lapply(listOfOutcomes$v, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ 
                                imd5C + 
                                factor(SurveyYear) + 
                                bs(age) + 
                                sex + 
                                ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design, cc)
    design_cc <- subset(design_cc, Source == "NDNS")
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) e)
}) |> setNames(listOfOutcomes$v)

imd_mods_ordered <- lapply(listOfOutcomes$v, function(y) {
  tryCatch({
    formula <- as.formula(paste("log(" , y, ") ~ 
                                imd5O + 
                                factor(SurveyYear) + 
                                bs(age) + 
                                sex + 
                                ethgr5"))
    vars <- all.vars(formula)
    vars_dt <- combined_design$variables[, ..vars]
    cc <- stats::complete.cases(vars_dt)
    design_cc <- subset(combined_design, cc)
    design_cc <- subset(design_cc, Source == "NDNS")
    svyglm(formula, design =  design_cc, family = gaussian())
  }, error = function(e) e)
}) |> setNames(listOfOutcomes$v)

imdresults <- imd_mods_continuous |> lapply(broom::tidy, conf.int=TRUE) |> 
  rbindlist(idcol="Outcome") |> 
  _[term=="imd5C"]
setkey(imdresults,Outcome)
imdresults[order_of_outputs]

imdOresults <- imd_mods_ordered |> lapply(broom::tidy, conf.int=TRUE) |> 
  rbindlist(idcol="Outcome") |> 
  _[term%like%"imd5O"]
setkey(imdOresults,Outcome)
imdOresults[order_of_outputs] |> _[term%like%"L" &p.value<0.05]


### Export models to spreadsheet
eqv_mods_year[!sapply(eqv_mods_year,is.null)] |> 
  modelsummary(
    estimate  = "{estimate} {p.value} {stars}",
    statistic = "({conf.low}, {conf.high})",
    fmt = fmt_decimal(digits = 2, pdigits = 3),
    exponentiate = TRUE,
    coef_omit = "bs(age)",
    output = "supplementary/eqvmodels.xlsx"
  )

eqv_mods_year_noGroup[!sapply(eqv_mods_year_noGroup,is.null)] |> 
  modelsummary(
    estimate  = "{estimate} {p.value} {stars}",
    statistic = "({conf.low}, {conf.high})",
    fmt = fmt_decimal(digits = 2, pdigits = 3),
    exponentiate = TRUE,
    coef_omit = "bs(age)",
    output = "graphs/eqvmodels_not_adjusted_for_group.xlsx"
  )

notNullGraphs <- seq_along(eqv_mods_year)[!sapply(eqv_mods_year,is.null)]

## Make income means graphs
incomeGraphs <- lapply(notNullGraphs , makeIncomeGraph,modslist=eqv_mods_year) |> 
  setNames(listOfOutcomes$varName[notNullGraphs])
meansOverEqvGraph <- wrap_plots(incomeGraphs[order_of_outputs],ncol = 5) + 
  plot_layout(axes="collect",guides="collect") & 
  scale_x_discrete(labels=c("High","Mid","Low")) & 
  scale_y_continuous(expand = expansion(0.1))

## Make income means graphs with no adjustment for obesity
incomeGraphs2 <- lapply(notNullGraphs , makeIncomeGraph,modslist=eqv_mods_year_noGroup) |> 
  setNames(listOfOutcomes$varName[notNullGraphs])
meansOverEqvGraph2 <- wrap_plots(incomeGraphs2[order_of_outputs],ncol = 5) + 
  plot_layout(axes="collect",guides="collect") & 
  scale_x_discrete(labels=c("High","Mid","Low")) & 
  scale_y_continuous(expand = expansion(0.1))


IMDGraphs2 <- lapply(seq_along(imd_mods_year_noGroup) , 
                     makeIMDGraph,modslist=imd_mods_year_noGroup) |> 
  setNames(listOfOutcomes$varName)
meansOverIMDGraph2 <- wrap_plots(IMDGraphs2[order_of_outputs],ncol = 5) + 
  plot_layout(axes="collect",guides="collect") 

meansOverIMDGraph2

## Make ethnicity graph
#plot_annotation(title="NDNS only, adjusting for income but not weight")
notNullEthGraphs <- seq_along(eth_mods)[!sapply(eth_mods,is.null)]
ethGraphs <- lapply(notNullEthGraphs, makeEthGraph,modslist=eth_mods)|> 
  setNames(listOfOutcomes$varName[notNullEthGraphs])
meansOverEth <- wrap_plots(ethGraphs[order_of_outputs],ncol = 5) + 
  plot_layout(axes="collect",guides="collect") & 
  theme(axis.text.x=element_text(angle=20,hjust=1)) &
  scale_y_continuous(expand = expansion(0.1))
ggsave(plot = meansOverEth,
       filename = "finalFigures/ethGraph.png",
       width=12,height=12,dpi="retina")

#plot_annotation(title="NDNS only, adjusting for weight not income and not year (GOOS doesn't have income or year)")
notNullEthGraphs2 <- seq_along(eth_mods2)[!sapply(eth_mods2,is.null)]
ethGraphs2 <- lapply(notNullEthGraphs2, makeEthGraph,modslist=eth_mods2)|> 
  setNames(listOfOutcomes$varName[notNullEthGraphs2])
meansOverEth2 <- wrap_plots(ethGraphs2[order_of_outputs],ncol = 5) + 
  plot_layout(axes="collect",guides="collect") & 
  theme(axis.text.x=element_text(angle=20,hjust=1)) &
  scale_y_continuous(expand = expansion(0.1)) 

ggsave(plot = meansOverEth2,
       filename = "finalFigures/ethGraph2.png",
       width=12,height=12,dpi="retina")

vitDeth <- combinedData[Source %in% c("NDNS","GOOS") & !is.na(ethgr5)] |> 
  ggplot(aes(x=ethgr5,y=v25OHD_Std,fill=ethgr5)) + ggbeeswarm::geom_beeswarm(cex=0.3) + geom_boxplot(alpha=0.5) + 
  facet_wrap(~Source)  + theme_bw()
ggsave(plot = vitDeth, filename="finalFigures/vitDeth.png", height=6,width=10,dpi="retina")

ggsave(plot = meansOverEqvGraph,
       filename = file.path("june26figs","income_graph_no_bmi_group_adjust.png"),
       width=12,height=12,dpi="retina")

ggsave(plot = meansOverEqvGraph2,
       filename = file.path("june26figs","income_graph_with_bmi_group_adjust.png"),
       width=12,height=12,dpi="retina")

ggsave(plot = meansOverIMDGraph2,
       filename = "finalFigures/IMD_Graph2.png",
       width=12,height=12,dpi="retina")

######################################################
### Similar analysis within the NNUH cohort only.  ###
######################################################

## We wouldn't necessarily expect the relationship between imd10 and nutrition to exist within the 
## very obese cohorts, because of the selection bias.  
## If we condition on being very overweight we remove the link between IMD and nutrition.
## Having said that, I haven't adjusted for anything in the simple graphical association..
NNUHdata |> melt(id.vars = c("subjectID","imd10","age","sex","Ethnicity","bmi","bmi_sds","ethgr5","imd2")) |>
  remove_missing() |> 
  ggplot(aes(x=factor(imd10),y=value)) + 
  geom_boxplot() + 
  facet_wrap(~variable,scale="free_y") + 
  scale_y_log10()

NNUHdata |> melt(id.vars = c("subjectID","imd10","age","sex","Ethnicity","bmi",
                             "bmi_sds","ethgr5","imd2","vitKCensored")) |>
  remove_missing() |> 
  ggplot(aes(x=imd10,y=value)) + 
  geom_point() + 
  facet_wrap(~variable,scale="free_y") + 
  scale_y_log10() + 
  stat_smooth(method="lm") + 
  ggpubr::stat_cor()

NNUHdataModels <- NNUHdata |> melt(id.vars = c("subjectID","imd10","age","sex","Ethnicity","bmi","bmi_sds","ethgr5","imd2")) |>
  remove_missing() |> 
  split(by="variable") |> 
  lapply(\(d){
    lm( data=d , value ~ imd10 + age + sex + ethgr5 )
  })

## Yep, no evidence for an age effect in that dataset.
NNUHdataModels |> modelsummary::modelsummary(stars=TRUE)