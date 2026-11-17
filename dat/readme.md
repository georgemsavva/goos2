
# How does the GOOS analysis file work

# Create the datasets

Before we run any analysis we need to create a dataset.  This is done using the latest main_analysis_file_[date].R

This in turn calls various files to load and code datasets, load various functions, establish thresholds for deficiency, including

createDataMain.R in turn calls makedatasets.R
add_threshold.R


Then, to get specific analyses

* obesityanalysis[date].R for the regression models across graphics across BMI groups

* IntakeVsNutrientsNov26.R does the adjustment for intake variables, and an explanation of how much the effects of obesity are explained by intake (NDNS only as only NDSD has that information).