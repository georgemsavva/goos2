

mods <- list()



mods[["SVitB12"]] <- gamlss( data = datasets[["SVitB12"]] , 
        logtarget ~ pb(age) , 
        sigma.formula = ~ 1 , 
        nu.formula = ~1 , 
        tau.formula = ~1,
        family="ST3C" ,
        weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)
centiles(mods[[1]], xvar = datasets[[1]]$age, 
         cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["SVitB12"]])

mods[["PFerritin"]] <- gamlss( data = datasets[["PFerritin"]] , 
                     logtarget ~ pb(age) , 
                     sigma.formula = ~ 1 , 
                     nu.formula = ~1 , 
                     tau.formula = ~1,
                     family="ST3C" ,
                     weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try=FALSE)
centiles(mods[["PFerritin"]], xvar = datasets[["PFerritin"]]$age, 
         cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[[2]])
summary(mods[[2]])


mods[["v25OHD_Std"]] <- gamlss( data = datasets[["v25OHD_Std"]] , 
                     logtarget ~ pb(age) , 
                     sigma.formula = ~ 1 , 
                     nu.formula = ~1 , 
                     tau.formula = ~1,
                     family="ST3C" ,
                     weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try=FALSE)
centiles(mods[["v25OHD_Std"]], xvar = datasets[["v25OHD_Std"]]$age, 
         cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["v25OHD_Std"]])
summary(mods[["v25OHD_Std"]])



mods[["holoTC"]] <- gamlss( data = datasets[["holoTC"]] , 
                     logtarget ~ pb(age) , 
                     sigma.formula = ~ 1 , 
                     nu.formula = ~1 , 
                     tau.formula = ~1,
                     family="ST3C" ,
                     weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try=FALSE)
centiles(mods[["holoTC"]], xvar = datasets[["holoTC"]]$age, 
         cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["holoTC"]])
summary(mods[["holoTC"]])


mods[["Retinol"]] <- gamlss( data = datasets[["Retinol"]] , 
                     logtarget ~ pb(age) , 
                     sigma.formula = ~ 1 , 
                     nu.formula = ~1 , 
                     tau.formula = ~1,
                     family="ST3C" ,
                     weights = finalweight)|>stepGAICAll.A(scope= ~age,mu.try=FALSE)
centiles(mods[["Retinol"]], xvar = datasets[["Retinol"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["Retinol"]])
summary(mods[["Retinol"]])


mods[["TotalCarotenoids"]] <- gamlss( data = datasets[["TotalCarotenoids"]] , 
                     logtarget ~ pb(age,control = pb.control(method = "GAIC", k = log(1297))) , 
                     sigma.formula = ~ 1 , 
                     nu.formula = ~1 , 
                     tau.formula = ~1,
                     family="ST3C" ,
                     weights = finalweight)|>stepGAICAll.A(scope= ~age,mu.try=FALSE)
centiles(mods[["TotalCarotenoids"]], xvar = datasets[["TotalCarotenoids"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["TotalCarotenoids"]])
summary(mods[["TotalCarotenoids"]])

##### CRP needs a separate modelling procedure because its highly censored distribution
gen.cens("ST3",type="left")
gen.cens("BCPE",type="left")
datasets[["CRP"]][ , event := as.numeric(target>1)]
datasets[["CRP"]][ , targetCen := ifelse(event==1,target,1)]
datasets[["CRP"]][ , logtargetCen := log(targetCen)]

mods[["CRP"]] <- gamlss( data = datasets[["CRP"]] , 
                     Surv( time=logtargetCen , type="left" , 
                           event=event) ~ pb(age) ,
                     sigma.formula = ~ pb(age) , 
                     nu.formula = ~1 , 
                     tau.formula = ~1,
                     family="ST3lc", ,
                     weights = finalweight,
                     control = gamlss.control(n.cyc = 200, trace = TRUE))


mods[["vitE"]] <- gamlss( data = datasets[["vitE"]] , 
                     logtarget ~ pb(age) , 
                     sigma.formula = ~ 1 , 
                     nu.formula = ~1 , 
                     tau.formula = ~1,
                     family="ST3C" ,
                     weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)
centiles(mods[["vitE"]], xvar = datasets[["vitE"]]$age, 
         cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["vitE"]])



mods[["pivka"]] <- gamlss( data = datasets[["pivka"]] , 
                     logtarget ~ pb(age) , 
                     sigma.formula = ~ 1 , 
                     nu.formula = ~1 , 
                     tau.formula = ~1,
                     family="ST3C" ,
                     weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)
centiles(mods[["pivka"]], xvar = datasets[["pivka"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["pivka"]])

datasets[["vitK"]][ , event := as.numeric(target>.1)]
datasets[["vitK"]][ , targetCen := ifelse(event==1,target,.1)]
datasets[["vitK"]][ , logtargetCen := log(targetCen)]
mods[["vitK"]] <- gamlss( data = datasets[["vitK"]] , 
                      Surv( time=logtargetCen , type="left" , 
                            event=event) ~ pb(age) ,
                     sigma.formula = ~ 1 , 
                     nu.formula = ~1 , 
                     tau.formula = ~1,
                     family="ST3lc" 
                    )|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)
centiles(mods[[10]], xvar = datasets[[10]]$age, cent = c(3, 10, 25, 50, 75, 90, 97))
wp(mods[[10]])


mods[["RCFolate"]] <- gamlss( data = datasets[["RCFolate"]] , 
                     logtarget ~ pb(age) , 
                     sigma.formula = ~ 1 , 
                     nu.formula = ~1 , 
                     tau.formula = ~1,
                     family="ST3C" ,
                     weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)
centiles(mods[["RCFolate"]], xvar = datasets[["RCFolate"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["RCFolate"]])

mods[["alphatocopherol"]] <- gamlss( data = datasets[["alphatocopherol"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)
centiles(mods[["alphatocopherol"]], xvar = datasets[["alphatocopherol"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["alphatocopherol"]])

n=nrow(datasets[["gammatocopherol"]])
mods[["gammatocopherol"]] <- gamlss( data = datasets[["gammatocopherol"]] , 
                      logtarget ~ pb(age,control = pb.control(method = "GAIC", k = log(1303))) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)
centiles(mods[["gammatocopherol"]], xvar = datasets[["gammatocopherol"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["gammatocopherol"]])

mods[["VitC"]] <- gamlss( data = datasets[["VitC"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["VitC"]], xvar = datasets[["VitC"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["VitC"]])

mods[["VitB1"]] <- gamlss( data = datasets[["VitB1"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["VitB1"]], xvar = datasets[["VitB1"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["VitB1"]])

mods[["VitB2"]] <- gamlss( data = datasets[["VitB2"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["VitB2"]], xvar = datasets[["VitB2"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["VitB2"]])

mods[["VitB6PLP"]] <- gamlss( data = datasets[["VitB6PLP"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["VitB6PLP"]], xvar = datasets[["VitB6PLP"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["VitB6PLP"]])

mods[["VitB6PA"]] <- gamlss( data = datasets[["VitB6PA"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["VitB6PA"]], xvar = datasets[["VitB6PA"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["VitB6PA"]])

mods[["Zn"]] <- gamlss( data = datasets[["Zn"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["Zn"]], xvar = datasets[["Zn"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["Zn"]])

mods[["Se"]] <- gamlss( data = datasets[["Se"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["Se"]], xvar = datasets[["Se"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["Se"]])

mods[["gammaVsAlpha"]] <- gamlss( data = datasets[["gammaVsAlpha"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["gammaVsAlpha"]], xvar = datasets[["gammaVsAlpha"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["gammaVsAlpha"]])

mods[["gammaVsChol"]] <- gamlss( data = datasets[["gammaVsChol"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["gammaVsChol"]], xvar = datasets[["gammaVsChol"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["gammaVsChol"]])

mods[["alphaVsChol"]] <- gamlss( data = datasets[["alphaVsChol"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["alphaVsChol"]], xvar = datasets[["alphaVsChol"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["alphaVsChol"]])

mods[["BCRY"]] <- gamlss( data = datasets[["BCRY"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["BCRY"]], xvar = datasets[["BCRY"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["BCRY"]])

mods[["BCAR"]] <- gamlss( data = datasets[["BCAR"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["BCAR"]], xvar = datasets[["BCAR"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["BCAR"]])

mods[["LYCO"]] <- gamlss( data = datasets[["LYCO"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["LYCO"]], xvar = datasets[["LYCO"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["LYCO"]])

mods[["LEUT"]] <- gamlss( data = datasets[["LEUT"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["LEUT"]], xvar = datasets[["LEUT"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["LEUT"]])

mods[["FerritinCRP"]] <- gamlss( data = datasets[["FerritinCRP"]] , 
                      logtarget ~ pb(age) , 
                      sigma.formula = ~ 1 , 
                      nu.formula = ~1 , 
                      tau.formula = ~1,
                      family="ST3C" ,
                      weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["FerritinCRP"]], xvar = datasets[["FerritinCRP"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["FerritinCRP"]])


mods[["SerumFolate"]] <- gamlss( data = datasets[["SerumFolate"]] , 
                                 logtarget ~ pb(age) , 
                                 sigma.formula = ~ 1 , 
                                 nu.formula = ~1 , 
                                 tau.formula = ~1,
                                 family="ST3C" ,
                                 weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["SerumFolate"]], xvar = datasets[["SerumFolate"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["FerritinCRP"]])

mods[["RetinolVsChol"]] <- gamlss( data = datasets[["RetinolVsChol"]] , 
                                 logtarget ~ pb(age) , 
                                 sigma.formula = ~ 1 , 
                                 nu.formula = ~1 , 
                                 tau.formula = ~1,
                                 family="ST3C" ,
                                 weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["RetinolVsChol"]], xvar = datasets[["RetinolVsChol"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["RetinolVsChol"]])


mods[["CarotVsChol"]] <- gamlss( data = datasets[["CarotVsChol"]] , 
                                   logtarget ~ pb(age) , 
                                   sigma.formula = ~ 1 , 
                                   nu.formula = ~1 , 
                                   tau.formula = ~1,
                                   family="ST3C" ,
                                   weights = finalweight)|>
  stepGAICAll.A(scope= ~age,mu.try = FALSE)

centiles(mods[["CarotVsChol"]], xvar = datasets[["CarotVsChol"]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
wp(mods[["CarotVsChol"]])


# mods[[29]] <- gamlss( data = datasets[[29]] , 
#                       logtarget ~ pb(age) , 
#                       sigma.formula = ~ 1 , 
#                       nu.formula = ~1 , 
#                       tau.formula = ~1,
#                       family="ST3C" ,
#                       weights = finalweight)|>
#   stepGAICAll.A(scope= ~age,mu.try = FALSE)
# 
# centiles(mods[[29]], xvar = datasets[[29]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
# wp(mods[[29]])
# 
# mods[[30]] <- gamlss( data = datasets[[30]] , 
#                       logtarget ~ pb(age) , 
#                       sigma.formula = ~ 1 , 
#                       nu.formula = ~1 , 
#                       tau.formula = ~1,
#                       family="ST3C" ,
#                       weights = finalweight)|>
#   stepGAICAll.A(scope= ~age,mu.try = FALSE)
# 
# centiles(mods[[30]], xvar = datasets[[30]]$age, cent = c(3, 10, 25, 50, 75, 90, 97),legend = FALSE)
# wp(mods[[30]])



saveRDS(mods,file = "mods.rds")

png("graphs/diagnostics1.png",width = 800,height=1200)
par(mfrow=c(5,2))
for(i in 1:5){
  wp(mods[[i]])
  centiles(mods[[i]])
}
dev.off()

png("graphs/diagnostics2.png",width = 800,height=1200)
par(mfrow=c(5,2))
for(i in 6:10){
  wp(mods[[i]])
  centiles(mods[[i]])
}
dev.off()

png("graphs/diagnostics3.png",width = 800,height=1200)
par(mfrow=c(5,2))
for(i in 11:15){
  wp(mods[[i]])
  centiles(mods[[i]])
}
dev.off()

png("graphs/diagnostics4.png",width = 800,height=1200)
par(mfrow=c(6,2))
for(i in 16:21){
  wp(mods[[i]])
  centiles(mods[[i]])
}
dev.off()


png("graphs/diagnostics5.png",width = 800,height=1200)
par(mfrow=c(4,2))
for(i in 22:25){
  wp(mods[[i]])
  centiles(mods[[i]])
}
dev.off()

png("graphs/diagnostics6.png",width = 800,height=1200)
par(mfrow=c(5,2))
for(i in 26:30){
  wp(mods[[i]])
  centiles(mods[[i]])
}
dev.off()



