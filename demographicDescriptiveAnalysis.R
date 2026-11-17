
combinedData[, eqv3F := factor(eqv3, levels=c("Highest Tertile","Middle Tertile","Lowest Tertile"))]
combinedData |> ggplot(aes(y=SVitB12, x=eqv3F))+ 
  geom_boxplot() + theme_bw()
combinedData |> ggplot(aes(y=SVitB12, x=eqv3F))+ stat_summary()

listOfOutcomes$varName

table(ndns$eqv3)

combinedData[ Source=="NDNS"& !is.na(eqv3F),.SD, 
              .SDcols = c("ethgr5","eqv3F","Group","age","sex",presentInNDNS)] |> 
  melt(measure.vars = presentInNDNS) |>
  ggplot(aes(x=eqv3F, y=value)) + 
  facet_wrap(~variable,scales="free_x") + 
  geom_boxplot() + theme_bw() + scale_y_log10() + 
  ggpubr::stat_compare_means() + coord_flip()

combinedData[ Source=="NDNS"& !is.na(eqv3F),.SD, 
              .SDcols = c("ethgr5","eqv3F","Group","age","sex",presentInNDNS)] |> 
  melt(measure.vars = presentInNDNS) |>
  ggplot(aes(x=Group, y=value)) + 
  facet_wrap(~variable,scales="free_x") + 
  geom_boxplot() + theme_bw() + scale_y_log10() + 
  ggpubr::stat_compare_means() + coord_flip()

combinedData[ Source=="NDNS"& !is.na(eqv3F),.SD, 
              .SDcols = c("ethgr5","eqv3F","Group","age","sex",presentInNDNS)] |> 
  melt(measure.vars = presentInNDNS) |>
  ggplot(aes(x=Group, y=value)) + 
  facet_wrap(~variable,scales="free_x") + 
  geom_boxplot() + theme_bw() + scale_y_log10() + 
  ggpubr::stat_compare_means() + coord_flip()

combinedData[ Source=="NDNS"& !is.na(eqv3F),.SD, 
              .SDcols = c("ethgr5","eqv3F","Group","age","sex",presentInNDNS)] |> 
  melt(measure.vars = presentInNDNS) |>
  ggplot(aes(y=Group, x=value)) + 
  facet_wrap(~variable,scales="free_x") + 
  stat_summary() + theme_bw() + scale_x_log10() + 
  ggpubr::stat_compare_means() 

combinedData[ Source=="NDNS"& !is.na(eqv3F),.SD, 
              .SDcols = c("ethgr5","eqv3F","Group","age","sex",presentInNDNS)] |> 
  melt(measure.vars = presentInNDNS) |>
  ggplot(aes(y=ethgr5, x=value)) + 
  facet_wrap(~variable,scales="free_x") + 
  stat_summary() + theme_bw() + scale_x_log10() + 
  ggpubr::stat_compare_means() 


combinedData$imd10 |> table()
combinedData[ !is.na(imd10),.SD, 
      .SDcols = c("ethgr5","imd10","Group","age","sex",presentInNNUH)] |>
  melt(measure.vars = c(presentInNNUH,"age")) |>
  ggplot(aes(x=imd10, y=value)) + 
  facet_wrap(~variable,scales="free_y") + 
  stat_smooth(method="lm") + theme_bw() + 
  scale_y_log10() + geom_point() + 
  ggpubr::stat_cor()


#### Obesity explains only a small proportion of the effect of income on carotenoids, ferritin

list(
lm(data=combinedData[Source=="NDNS"], log(TotalCarotenoids) ~ Group + age + sex) ,
lm(data=combinedData[Source=="NDNS"], log(TotalCarotenoids) ~ eqv3F + age + sex) ,
lm(data=combinedData[Source=="NDNS"], log(TotalCarotenoids) ~ Group + eqv3F + age + sex)
) |> modelsummary(stars=TRUE)

list(
  lm(data=combinedData[Source=="NDNS"], log(FerritinCRP) ~ Group + age + sex) ,
  lm(data=combinedData[Source=="NDNS"], log(FerritinCRP) ~ eqv3F + age + sex) ,
  lm(data=combinedData[Source=="NDNS"], log(FerritinCRP) ~ Group + eqv3F + age + sex)
) |> modelsummary(stars = TRUE)

list(
  lm(data=combinedData[Source=="NDNS"], log(v25OHD_Std) ~ Group + age + sex) ,
  lm(data=combinedData[Source=="NDNS"], log(v25OHD_Std) ~ eqv3F + age + sex) ,
  lm(data=combinedData[Source=="NDNS"], log(v25OHD_Std) ~ Group + eqv3F + age + sex)
) |> modelsummary(stars = TRUE)


lm(data=combinedData[Source=="NDNS"], 
   log(v25OHD_Std) ~ Group + eqv3F + age + sex) |> performance::check_model()
lm(data=combinedData[Source=="NDNS"], 
   (v25OHD_Std) ~ Group + eqv3F + ethgr5 + age + sex) |> performance::check_model()

#### No evidence for an interaction here..
modinteract <- lm(data=combinedData[Source=="NDNS"], 
   log(TotalCarotenoids) ~ Group * eqv3F + age + sex)
modinteract |> car::Anova(type=3)

emmeans::emmeans(modinteract , ~eqv3F|Group) |> emmeans::contrast("pairwise")

deficient_names <- names(combinedData) |> grep(pattern="deficient$", value=TRUE) |> 
  setdiff(names(combinedData) |> grep(pattern="possible_deficient$", value=TRUE))

names(combinedData) |> grep(pattern="deficient$", value=TRUE)
names(combinedData) |> grep(pattern="elevated$", value=TRUE)

## (SF (µg/l) and CRP (mg/l))
## Vitamins A, B1, B6, B12, C, D, E, FerritinVsCRP, Carotenoids, Zinc, Selenium 

deficiencyCols <- 
  c("Retinol_deficient",
  "VitB1_elevated",
  "VitB6PA_deficient", 
  "SVitB12_deficient", 
  "VitC", 
  "v25OHD_Std_deficient", 
  "vitE_deficient", 
  "FerritinCRP_deficient", 
  "TotalCarotenoids_deficient")

combinedData[ Source=="NDNS", (.SD) , .SDcols = deficiencyCols]

## To make the total number of deficiencies, combine:
#### 

deficient_names <- names(combinedData) |> grep(pattern="deficient$", value=TRUE)

combinedData[ , totalDeficient := rowSums(.SD, na.rm = TRUE) , .SDcols = deficient_names ]
combinedData[ , totalMeasured := rowSums(!is.na(.SD), na.rm = TRUE) , .SDcols = deficient_names ]
combinedData[ , .(totalDeficient,totalMeasured)]

combinedData[ , table(Group,totalMeasured)]
combinedData[ ]

combinedData[ Source=="NDNS", table(as.character(Group),totalDeficient)] |> 
  prop.table(1)

combinedData[ Source=="NDNS", .SD, .SDcols = c("Group", deficient_names)] |>
  _[ , lapply(.SD , \(v) table(v,Group)) , .SDcols=deficient_names]

#### What is the correlation between all of the deficiencies?
#### What is the multi-deficiency outcome..





############################################
### PCA ANALYSIS THAT WE ARE NOT DOING  ####
############################################

combinedData[Source=="GOOS", .SD , .SDcols=presentInGOOS] |> pca()
combinedData[Source=="GOOS", .SD , .SDcols=presentInGOOS] |> princomp()
combinedData[Source=="GOOS", .SD , .SDcols=presentInGOOS] |> remove_missing() |> princomp()
combinedData[Source=="GOOS", .SD , .SDcols=presentInGOOS] |> remove_missing() |> princomp() |> screeplot()
combinedData[Source=="GOOS", .SD , .SDcols=presentInGOOS] |> remove_missing() |> prcomp() |> screeplot()
combinedData[Source=="GOOS", .SD , .SDcols=presentInGOOS] |> remove_missing() |> prcomp(scale. = TRUE) |> screeplot()
combinedData[Source=="GOOS", .SD , .SDcols=presentInGOOS] |> remove_missing() |> prcomp(scale. = TRUE,center = TRUE) |> screeplot()
combinedData[Source=="GOOS", .SD , .SDcols=presentInGOOS] |> remove_missing() |> prcomp(scale. = TRUE) |> screeplot()
library(psych)
pca1 <- combinedData[Source=="GOOS", .SD , .SDcols=presentInGOOS] |> remove_missing() |> principal() |> screeplot()
pca1 <- combinedData[Source=="GOOS", .SD , .SDcols=presentInGOOS] |> remove_missing() |> principal()
pca1
pca1 <- combinedData[Source=="GOOS", .SD , .SDcols=presentInGOOS] |> remove_missing() |> principal(nfactors = 5)
pca1
pca1 <- combinedData[Source=="GOOS", .SD ] |> remove_missing() |> principal(nfactors = 5)
presentInGOOS
presentInGOOS |> dput()
presentInGOOS2 <- c("SVitB12", "v25OHD_Std", "holoTC", "Retinol",
                    "TotalCarotenoids", "CRP", "vitK", "alphatocopherol",
                    "gammatocopherol", "FerritinCRP")
pca1 <- combinedData[Source=="GOOS", .SD , .SDcols = presentInGOOS] |> remove_missing() |> principal(nfactors = 5)
pca1 <- combinedData[Source=="GOOS", .SD , .SDcols = presentInGOOS2] |> remove_missing() |> principal(nfactors = 5)
pca1
presentInGOOS2 <- c("SVitB12", "v25OHD_Std", "Retinol",
                    "TotalCarotenoids", "CRP", "vitK", "alphatocopherol",
                    "gammatocopherol", "FerritinCRP")
pca1 <- combinedData[Source=="GOOS", .SD , .SDcols = presentInGOOS2] |> remove_missing() |> principal(nfactors = 5)
pca1
pca1 <- combinedData[Source=="GOOS", .SD , .SDcols = presentInGOOS2] |> corrplot::corrplot()
pca1 <- combinedData[Source=="GOOS", .SD , .SDcols = presentInGOOS2] |> cor() |> corrplot::corrplot()
pca1 <- combinedData[Source=="GOOS", .SD , .SDcols = presentInGOOS2] |> cor(use = "everything") |> corrplot::corrplot()
pca1 <- combinedData[Source=="GOOS", .SD , .SDcols = presentInGOOS2] |> cor(use = "pairwise") |> corrplot::corrplot()
pca1 <- combinedData[Source=="GOOS", .SD , .SDcols = presentInGOOS2] |> cor(use = "pairwise", method="sprearman") |> corrplot::corrplot()
pca1 <- combinedData[Source=="GOOS", .SD , .SDcols = presentInGOOS2] |> cor(use = "pairwise", method="spearman") |> corrplot::corrplot()
pca1 <- combinedData[Source=="NDNS", .SD , .SDcols = presentInGOOS2] |> cor(use = "pairwise", method="spearman") |> corrplot::corrplot()
pca1 <- combinedData[Source=="GOOS", .SD , .SDcols = presentInNDNS] |> cor(use = "pairwise", method="spearman") |> corrplot::corrplot()
pca1 <- combinedData[Source=="NDNS", .SD , .SDcols = presentInNDNS] |> cor(use = "pairwise", method="spearman") |> corrplot::corrplot()