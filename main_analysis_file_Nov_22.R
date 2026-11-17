
#### THIS IS THE MAIN FILE.. RUN THIS ONE.

#### If you 'source' this file it will generate the datasets then stop with an 'error' caused by the stop() function.  This is expected!  Then you can source whichever of the analysis files you need to.


source("createDataMain.R")

#################################################
## Make GAMLSS models or load existing models  ##
#################################################

## RUN THIS IF YOU WANT TO REBUILD THE CENTILE MODELS 
## LEAVE WITH COMMENT IF YOU WANT JUST TO LOAD THE EXISTING MODELS.
# source("buildModels.R")
# Do we use the centile models for anything?
# mods <- readRDS("mods.rds")

## Now add the binary outcomes and centiles based on the thresholds and mods
source("add thresholds_june26.R")

add_threshold(combinedData) #by reference
makeStandardisationWeight(combinedData) #by reference (adds the 'analysis weight')

## As of Nov 26 no longer makes predictions, because we don't use that information.
combinedData <- addPredictionsAndDeficiencies(combinedData) # not by reference?

###################################################################################
## Uncomment this to stop once all the data is loaded, 
###################################################################################

 stop("Created file, now you can run whatever graphs/analyses you need")

###################################################################
#### Files and functions to produce outputs
###################################################################


##### Table 1 : Participant demographic characteristics
makeDemographicTable()

source("prevalenceAnalysis.R") # Make propmods for the prevalence analysis.
##### Table 2 : Weighted prevalence of micronutrient deficiencies, with models
makeDeficiencyTable()  # - relies on 'propmods' existing

##### 
propmodsVitmainDspecial_deficient |> summary()
propmodsVitmainDspecial_possible_deficient |> summary()


##### Table 3 : Sample weighted medians of biomarkers across groups (no modelling)
bioMarkerDescriptiveTable()

##### Table 3 (alternative) : Population weighted medians of biomarkers across groups (no modelling)
bioMarkerDescriptiveTableWeighted()

##### Figure 1 : NDNS normal centiles plus GOOS distributions
source("predictionsAndPlots.R")

##### Figure Supplement (ridgelines) : 
##### Replaces supp figure 2 : Boxplots of centiles, per biomarker
makeRidgePlot(width=9.5,height=10)

##########################################################
##### MAIN OBESITY ANALYSIS              #################
##### Figure 2 : Marginal means across NSNS/GOOS groups  #
##########################################################

source("obesityAnalysis.R") # produces the figure 2 graphs and associated sensitivity analysis

###########################################################
##### Consider the impact of Household Income (tertiles) ##
##### Figure 3 : Marginal means across income groups     ##
###########################################################

source("incomeAnalysis.R")  # produces the figure 3 graphs and associated sensitivity analysis

##### Supp table 3 : Sample sizes
sampleSizePerOutcomeTable()

##### Supp table 5 : Raw medians and centiles relative to NDNS healthy
## Better to use this for model means
marginalMeansTable()


# Session information for reproducibility
sink("june26figs/sessionInfo.txt")
sessionInfo()
sink(NULL)

# citations
citation("emmeans")
citation("survey")
