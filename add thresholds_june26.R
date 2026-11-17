

## 23 June 2026 - This file updated with the definitions from table S4 in the paper.
## Add person-specific deficiency etc threshold to the combined data file:

add_threshold <- function(combinedData){
  # We don't need to return anything, because a data.table is modified in place by reference
 
    ## my UNIT here IS ug/l
    combinedData[ CRP<=5 & age<5, PFerritin_deficient_threshold := 12 ] #  table S4
    combinedData[ CRP<=5 & age>=5, PFerritin_deficient_threshold := 15 ]  #  table S4
    combinedData[ CRP>5 & age<5, PFerritin_deficient_threshold := 30 ] #  table S4
    combinedData[ CRP>5 & age>=5, PFerritin_deficient_threshold := 70 ] #  table S4

    combinedData[ CRP>5 & age>=5   , PFerritin_elevated_threshold := 500 ]  #  table S4
    
    combinedData[ CRP<=5 & age>=5 & sex=="Male"  , PFerritin_elevated_threshold := 200 ]  #  table S4
    combinedData[ CRP<=5 & age>=5 & sex=="Female", PFerritin_elevated_threshold := 150 ] #  table S4
    
    combinedData[ , PFerritin_possible_deficient_threshold := NA ]
    
    
    
  # Vitamin B12: 
  
    combinedData[ age<=16, SVitB12_deficient_threshold := 150  ]  # table S4
    combinedData[ age>16, SVitB12_deficient_threshold := 133  ] # table S4
    
    combinedData[ age>16, SVitB12_possible_deficient_threshold := 258  ] # table S4
    combinedData[ , SVitB12_elevated_threshold := 675  ] # table S4
    
  # Holo TC: 
    
    combinedData[ age<=16, holoTC_deficient_threshold := 32  ]
    combinedData[ age>16, holoTC_deficient_threshold := 25  ]
    ## Don't use this because it excludes most of the cohort.
    ## There are 571 over 16s in the whole dataset, so maybe stick with it.
    combinedData[ age>16, holoTC_possible_deficient_threshold := 70  ] #  table S4
    #combinedData[ age>16, holoTC_possible_deficient_threshold := NA  ]
    combinedData[ age>2 & age<11, holoTC_elevated_threshold := 206  ]  #  table S4
    combinedData[ age>=11, holoTC_elevated_threshold := 178  ] #  table S4
    
  # Total 25-OH VD: 
 
    combinedData[ , v25OHD_Std_deficient_threshold := 30  ] # table s4
    combinedData[ , v25OHD_Std_possible_deficient_threshold := 50  ] # table s4
    combinedData[ , v25OHD_Std_elevated_threshold := 125  ] # table s4
    
  # Retinol：
    
    combinedData[ , Retinol_deficient_threshold := 0.35  ]
    combinedData[ , Retinol_possible_deficient_threshold := 0.7  ]
    
    combinedData[ , Retinol_elevated_threshold := 3.4  ]
    combinedData[ age<13, Retinol_elevated_threshold := 2.8  ]
    combinedData[ age<7, Retinol_elevated_threshold := 2.3  ]
    
    #combinedData[ , Retinol_elevated_threshold := NA] # note from Albert - no good ref for elevated retinol.
    
    combinedData[ , RetinolVsChol_deficient_threshold := NA  ]
    combinedData[ , RetinolVsChol_possible_deficient_threshold := NA  ]
    combinedData[ , RetinolVsChol_elevated_threshold := NA  ]
    
  ## Vitamin E
    
    combinedData[ , vitE_deficient_threshold := 13  ] # Confirmed with Albert by Emma in March 25 meeting
    combinedData[ , vitE_possible_deficient_threshold := NA  ]
    combinedData[ , vitE_elevated_threshold := NA ] # FROM Albert - no reference

    
    combinedData[ , alphatocopherol_deficient_threshold := 11.6  ] # table s4
    combinedData[ , alphatocopherol_possible_deficient_threshold := NA  ]
    combinedData[ , alphatocopherol_elevated_threshold := NA ] # table s4
    
    combinedData[ , gammatocopherol_deficient_threshold := NA  ] 
    combinedData[ , gammatocopherol_possible_deficient_threshold := NA  ]
    combinedData[ , gammatocopherol_elevated_threshold := NA ] 
    
    combinedData[ , gammaVsAlpha_deficient_threshold := NA ]  # from Emma's table
    combinedData[ , gammaVsAlpha_possible_deficient_threshold := NA ]
    combinedData[ , gammaVsAlpha_elevated_threshold := NA ]
    
  # Vitamin K1
    
    combinedData[ , vitK_deficient_threshold := NA  ] # modify to 0.15 from 0.1 following meeting with Emma
    combinedData[ , vitK_possible_deficient_threshold := NA  ]
    combinedData[ , vitK_elevated_threshold := NA  ]
    
  # CRP:

    combinedData[ , CRP_deficient_threshold := NA  ]
    combinedData[ , CRP_possible_deficient_threshold := NA  ]
    combinedData[ , CRP_elevated_threshold := 5  ]
    
  # PIVKA
    
    combinedData[ , pivka_deficient_threshold := NA  ]
    combinedData[ , pivka_possible_deficient_threshold := NA  ]
    combinedData[ , pivka_elevated_threshold := 50.9  ]
    
  # Folate threshold
    
    combinedData[ , RCFolate_deficient_threshold := 305  ] # NDNS publication defined.
    combinedData[ , RCFolate_possible_deficient_threshold := NA  ]
    combinedData[ , RCFolate_elevated_threshold := NA  ]

    combinedData[ , SerumFolate_deficient_threshold := 3  ] # table S4
    combinedData[ , SerumFolate_possible_deficient_threshold := 5.9  ] # table S4
    combinedData[ , SerumFolate_elevated_threshold := 20  ] # table S4
    
  # Total Carotenoids
  # 1000 nmol/L = 1 µmol/L
    
    combinedData[ , TotalCarotenoids_deficient_threshold := 1 ] # s4
    combinedData[ , TotalCarotenoids_possible_deficient_threshold := NA  ]
    combinedData[ , TotalCarotenoids_elevated_threshold := NA ]
    
  # Vitamin C
      
    combinedData[ , VitC_deficient_threshold := 11 ] # table S4
    combinedData[ , VitC_possible_deficient_threshold := 23 ] # table S4
    combinedData[ , VitC_elevated_threshold := NA ]
    
  ## B1 and B2 measures are "high=bad", so we use "elevated" not "deficient".
  
    combinedData[ , VitB1_deficient_threshold := NA  ] # from Emma's table S4
    combinedData[ , VitB1_possible_deficient_threshold := NA ]
    combinedData[ , VitB1_elevated_threshold := 1.25 ]
    combinedData[ , VitB1_possible_elevated_threshold := 1.15 ]
    
    combinedData[ , VitB2_deficient_threshold := NA ] # from Emma's table S4
    combinedData[ , VitB2_possible_deficient_threshold := NA ]
    combinedData[ , VitB2_elevated_threshold := 1.3 ]      ### NOTE FROM EMMA 
    combinedData[ , VitB2_possible_elevated_threshold := NA ]
  
  ## B6
      
    combinedData[ , VitB6PLP_deficient_threshold := 20 ] # table S4
    combinedData[ , VitB6PLP_possible_deficient_threshold := NA ]
    combinedData[ , VitB6PLP_elevated_threshold := NA ]
    
    combinedData[ , VitB6PA_deficient_threshold := NA ]
    combinedData[ , VitB6PA_possible_deficient_threshold := NA ]
    combinedData[ , VitB6PA_elevated_threshold := NA ]
  
  ## Se
    
    combinedData[ , Se_deficient_threshold := NA ]  # from Albert - no ref
    combinedData[ , Se_possible_deficient_threshold := NA ]
    combinedData[ , Se_elevated_threshold := NA ]
  
  ## Zn
    
    combinedData[ , Zn_deficient_threshold := 9.9 ] # from Emma table S4
    combinedData[ , Zn_possible_deficient_threshold := NA_real_ ]
    combinedData[ , Zn_elevated_threshold := NA_real_ ]
   
  ## Everything else should be set to NA.
     
    combinedData[ , gammaVsChol_deficient_threshold := NA ]  
    combinedData[ , gammaVsChol_possible_deficient_threshold := NA ]
    combinedData[ , gammaVsChol_elevated_threshold := NA ]
    
    combinedData[ , alphaVsChol_deficient_threshold := NA ]  
    combinedData[ , alphaVsChol_possible_deficient_threshold := NA ]
    combinedData[ , alphaVsChol_elevated_threshold := NA ]
    
    combinedData[ , ACRY_deficient_threshold := NA ]  
    combinedData[ , ACRY_possible_deficient_threshold := NA ]
    combinedData[ , ACRY_elevated_threshold := NA ]
    
    combinedData[ , BCRY_deficient_threshold := NA ]  
    combinedData[ , BCRY_possible_deficient_threshold := NA ]
    combinedData[ , BCRY_elevated_threshold := NA ]
    
    combinedData[ , ACAR_deficient_threshold := NA ]  
    combinedData[ , ACAR_possible_deficient_threshold := NA ]
    combinedData[ , ACAR_elevated_threshold := NA ]
    
    combinedData[ , BCAR_deficient_threshold := NA ]  
    combinedData[ , BCAR_possible_deficient_threshold := NA ]
    combinedData[ , BCAR_elevated_threshold := NA ]
    
    combinedData[ , LYCO_deficient_threshold := NA ]  
    combinedData[ , LYCO_possible_deficient_threshold := NA ]
    combinedData[ , LYCO_elevated_threshold := NA ]
    
    combinedData[ , LEUT_deficient_threshold := NA ]  
    combinedData[ , LEUT_possible_deficient_threshold := NA ]
    combinedData[ , LEUT_elevated_threshold := NA ]

    combinedData[ , FerritinCRP_deficient_threshold := NA ]  
    combinedData[ , FerritinCRP_possible_deficient_threshold := NA ]
    combinedData[ , FerritinCRP_elevated_threshold := NA ]
    
    combinedData[ , CarotVsChol_deficient_threshold := NA ]  
    combinedData[ , CarotVsChol_possible_deficient_threshold := NA ]
    combinedData[ , CarotVsChol_elevated_threshold := NA ]
   
}

