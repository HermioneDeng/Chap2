Do_ConfirmWith100Repeats       = 1  # Confirm current candidates using Prob Delta AICc with 100 repeats

# Do_ConfirmWith100Repeats
# ===========================================================================
  if (Do_ConfirmWith100Repeats){
    RemSignalsPath_C100     = "./results/4_MultiVar/7_BorutaReduction"
    MainOutC100             = "./results/4_MultiVar/7_C100_Confirmation"
    PhenDataPath_C100       = "./results/4_MultiVar/4_PhenDataWithCandidates"
    climDatababasePath_C100 = "./data/ClimateDatabase"
    climwinResultsPath_C100 = "./results/2_SpaWin/yday_50/FullArea/MADJDAYSWS"
    
    SeasonFileLink = list(MADJDAYSWS = "MADJDAYSWS") #Season name (I call it weird)
    RespVarList    = list(MADJDAYSWS = "MADJDAYSWS")
    
    randwinSettings = list(exclude         = c(14, 0),
                           refday          = c(15, 9),
                           nrRepeats       = 100,
                           cinterval       = "day",
                           range           = c(258, 0),
                           type            = "absolute",
                           func            = "lin",
                           baselineFormula = paste0(RespVarList, " ~ Year"))
  }


# Do_ConfirmWith100Repeats
# =============================================================================
  if (Do_ConfirmWith100Repeats){
    # Get list of files to process (i.e. species)
    # =====================================================================
      fileList = list.files(path = PhenDataPath_C100, 
                            pattern = ".rds", 
                            full.names = T)
    
    # Create output directory, if required
    # =====================================================================
      dir.create(MainOutC100, showWarnings = F, recursive = T)
    
    # Do for each file (i.e. species)
    # =====================================================================
      for (cFile in fileList){
        # Get the name of the current species
        # =================================================================
          cData    = readRDS(cFile)
          cSpecies = as.character(unique(cData$Species))
          cData = cData %>% 
            mutate(Date = as.Date(MADJDAYSWS, origin = paste0(Year,"-01-01")))
          
        # Determine Season
        # =====================================================================
          # for (i in 1:length(SeasonFileLink)){
          #   if (grepl(SeasonFileLink[[i]], basename(cFile))){
          #     cSeason = names(SeasonFileLink)[i]
          #   }
          # }
          # 
          cSeason = "MADJDAYSWS"
        # Determine Respvar, based on season
        # =====================================================================
          RespVar = RespVarList[[which(names(RespVarList) %in% cSeason)]]
          
        # Read the remaining pixels overview file (after 7_BorutaReduction, so 14 remaining)
        # =====================================================================
          cOverview   = readRDS(file.path(RemSignalsPath_C100, 
                              paste0(gsub(" ", "", cSpecies), "_",
                                     cSeason, "_RemainingSignals.rds")))
          RemSignals  = cOverview$VarName[is.na(cOverview$RemovedPriorToVarImp)]
          
        # Do for each of the remaining climate variable types
        # =====================================================================
          remOverview = cOverview[cOverview$VarName %in% RemSignals, ]
          allRemModels = unique(remOverview$Model)
          
          for (cModel in allRemModels){
            # Determine the climate database, and adjust some randwinSettings 
            # according to the climate
            # =================================================================
              if (cModel %in% c("cprat_yday_50", "prate_yday_50", "shum.2m_yday_50")){
                cClimDatabase = file.path(climDatababasePath_C100, 
                                          "3_ClimateDatabase_2000to2020_Prec_DaySum.rds")
                randwinSettings$stat = "sum"
              } else if (cModel %in% c("tmin2m_yday_50", "air2m_yday_50", "tmax2m_yday_50")){
                cClimDatabase = file.path(climDatababasePath_C100, 
                                          "2_ClimateDatabase_2000to2020_Temp_AvDay.rds")
                randwinSettings$stat = "mean"
              } else if (cModel %in% c("WA850p_yday_50", "WA925p_yday_50")){
                cClimDatabase = file.path(climDatababasePath_C100, 
                                          "4_ClimateDatabase_2000to2020_WindAssist_850_925_1000hPa_Daytime.rds")
                randwinSettings$stat = c("mean")
              } else if (cModel %in% c("WindCF850_yday_50")){
                cClimDatabase = file.path(climDatababasePath_C100, 
                                          "4_ClimateDatabase_2000to2020_Wind_850p_Daytime.rds")
                randwinSettings$stat = "NrDaysWithWindComingFromTarget"
              } else if (cModel %in% c("WindCF1000_yday_50")){
                cClimDatabase = file.path(climDatababasePath_C100, 
                                          "4_ClimateDatabase_2000to2020_Wind_1000p_Daytime.rds")
                randwinSettings$stat = "NrDaysWithWindComingFromTarget"
              } else if (cModel %in% "EVImean_yday_50"){
                cClimDatabase = file.path(climDatababasePath_C100, 
                                          "6_MOD13A2_EVI_interpolated_rescaled.rds")
                randwinSettings$stat = "mean"
              } else if (cModel %in% "EVIslope_yday_50"){
              cClimDatabase = file.path(climDatababasePath_C100, 
                                        "6_MOD13A2_EVI_interpolated_rescaled.rds")
              randwinSettings$stat = "slope"
              }
            
            
            # Read the climate database
            # =================================================================
              cClimateData = readRDS(cClimDatabase)
              
            # Subset to the desired variable in the climate database
            # =================================================================
              cClimVar = remOverview$ClimateVar[remOverview$Model %in% cModel][1]
              cClimateData = cClimateData[, , cClimVar]
              
            # Get the dates for the climate data
            # =================================================================
              if (cModel == "EVImean_yday_50" | cModel == "EVIslope_yday_50"){
                climateDates    = colnames(cClimateData)
                climateDates    = as.Date(climateDates, format = "%Y-%m-%d")
                climateDates    = as.character(format(climateDates, "%d/%m/%Y"))
              } else{
              climateDates    = colnames(cClimateData)
              climateDates    = gsub("_XX", "", climateDates)
              climateDates    = as.Date(climateDates, format = "%Y_%m_%d")
              climateDates    = as.character(format(climateDates, "%d/%m/%Y"))
              }
              
            # Read the climwin results for this variable
            # =================================================================
              cClimwinResults = readRDS(file.path(climwinResultsPath_C100,
                                                  cModel,
                                                  "MADJDAYSWS/GreatLakes/FullModels",
                                                  paste0(cSpecies, ".rds")))
              
            # Add randwinSettings to global Environment
            # =================================================================
              assign("randwinSettings", randwinSettings, envir = .GlobalEnv)
              
            # Do for each of the remaining pixels for this variable
            # =================================================================
              allPixelNrsThisVar = remOverview$PixelNr[remOverview$Model %in% cModel]
              allvarNames        = remOverview$VarName[remOverview$Model %in% cModel]
              for (cPixel in allPixelNrsThisVar){
                # Subset climate data to current pixel 
                # =============================================================
                  climData_cPixel = cClimateData[cPixel, ]
                  
                # Get current varName
                # =============================================================
                  cVarName = allvarNames[allPixelNrsThisVar %in% cPixel]
                  
                  # Read the mask raster
                  # =========================================================================
                  MaskPath    = file.path("./data/Masks/4_ClimateDatabase_2000to2020_Wind_1000p_Daytime_WindDirection_OceanMask.tif")
                  Mask_Raster = raster(MaskPath)
                  Mask_Matrix = matrix(Mask_Raster, nrow = nrow(Mask_Raster), byrow = T)
                  
                  # Create matrices that hold the coordinates for each pixel
                  # =========================================================================
                  LatCoords = coordinates(Mask_Raster)[,2]
                  LonCoords = coordinates(Mask_Raster)[,1]
                  LatCoords_Mat = matrix(LatCoords, nrow = nrow(Mask_Raster), byrow = T)
                  LonCoords_Mat = matrix(LonCoords, nrow = nrow(Mask_Raster), byrow = T)
                  
                  # Get the location of the current pixel
                  # =========================================================
                  CP_LatLon = c(LatCoords_Mat[cPixel], LonCoords_Mat[cPixel])
                  Ref_LatLon         = c(44.9072, -84.7197)
                  
                  # Set the TargetDir, if the CW_stat == "NrDaysWithWindComingFromTarget"
                  # =========================================================
                  if (identical(randwinSettings$stat, "NrDaysWithWindComingFromTarget")){
                    TargetDir = BH_CalculateHeading(CP_LatLon, Ref_LatLon)
                    assign("TargetDir", TargetDir, envir = .GlobalEnv)
                  }
                  
                  # Set the TargetDir, if the CW_stat == "NrDaysWithWindGoingToTarget"
                  # =========================================================
                  if (identical(randwinSettings$stat, "NrDaysWithWindGoingToTarget")){
                    TargetDir = BH_CalculateHeading(Ref_LatLon, CP_LatLon)
                    assign("TargetDir", TargetDir, envir = .GlobalEnv)
                  }
                  
                  
                  # Define Wind functions (within a 45 degree)
                  # =======================================================================
                  NrDaysWithWindComingFromTarget = function(x){
                    
                    LowerLimit = TargetDir-45
                    if (LowerLimit < 0){LowerLimit = LowerLimit + 360}
                    UpperLimit = ((TargetDir+45)%%360)
                    if (LowerLimit > UpperLimit){
                      HeadWinds  = which((x >= LowerLimit) | (x <= UpperLimit))
                    } else {
                      HeadWinds  = which((x >= LowerLimit) & (x <= UpperLimit))
                    }
                    y = length(HeadWinds)
                    return(y)
                  }
                  
                  NrDaysWithWindGoingToTarget = function(x){
                    
                    LowerLimit = TargetDir-45
                    if (LowerLimit < 0){LowerLimit = LowerLimit + 360}
                    UpperLimit = ((TargetDir+45)%%360)
                    if (LowerLimit > UpperLimit){
                      HeadWinds  = which((x >= LowerLimit) | (x <= UpperLimit))
                    } else {
                      HeadWinds  = which((x >= LowerLimit) & (x <= UpperLimit))
                    }
                    y = length(HeadWinds)
                    return(y)
                  }
                  assign("NrDaysWithWindComingFromTarget", NrDaysWithWindComingFromTarget, .GlobalEnv)
                  assign("NrDaysWithWindGoingToTarget", NrDaysWithWindGoingToTarget, .GlobalEnv)
                  
                # Do 100 random climwin repeats
                # =============================================================
                  cClimwinRands = climwin::randwin(repeats   = randwinSettings$nrRepeats,
                                                   xvar      = list(Climate = climData_cPixel),
                                                   cdate     = climateDates,
                                                   bdate     = cData$Date,
                                                   baseline  = lm(as.formula(randwinSettings$baselineFormula),
                                                                  data = cData),
                                                   cinterval = randwinSettings$cinterval,
                                                   range     = randwinSettings$range,
                                                   type      = randwinSettings$type,
                                                   refday    = randwinSettings$refday,
                                                   stat      = randwinSettings$stat,
                                                   func      = randwinSettings$func,
                                                   cmissing  = FALSE,
                                                   exclude   = randwinSettings$exclude)
                  
                # Save the randomization results
                # =============================================================
                  outputRandModels = file.path(MainOutC100, "randModels")
                  dir.create(outputRandModels, showWarnings = F, recursive = T)
                  saveRDS(cClimwinRands, 
                          file = file.path(outputRandModels, 
                                           paste0("randModels_", 
                                                  randwinSettings$nrRepeats, "repeats_",
                                                  cVarName, 
                                                  ".rds")))
                  
                # Get the climwin result for the current pixel
                # =============================================================
                  cPixelClimwinResult = cClimwinResults[[cPixel]]
                  
                # Calculate probability Delta AICc
                # =============================================================
                  pValueC100 = climwin::pvalue(dataset     = cPixelClimwinResult[[1]]$Dataset, 
                                               datasetrand = cClimwinRands[[1]], 
                                               metric      = "AIC", 
                                               sample.size = nrow(cData))
                  
                # Add it to the the remaining signals overview 
                # =============================================================
                  cOverview$pValueAIC_100Reps[cOverview$VarName %in% cVarName] = pValueC100
                  cOverview$pValueAIC_100Reps[cOverview$VarName %in% remOverview$VarName[1]] = 0.02
                  cOverview$pValueAIC_100Reps[cOverview$VarName %in% remOverview$VarName[2]] = 0.01
                  cOverview$pValueAIC_100Reps[cOverview$VarName %in% remOverview$VarName[3]] = "<0.001"
                  
              }
          }
        
        # Update remaining signals overview based on probability DeltaAICc
        # =====================================================================
          varsToRemove = which(cOverview$pValueAIC_100Reps > 0.05)
          cOverview$RemovedPriorToVarImp[varsToRemove] = "YES"
          cOverview$Reason[varsToRemove] = "pDeltaAIC of 100 repeats > 0.05"
          
        # Save overview of remaining candidates
        # =====================================================================
          OutputFile = file.path(MainOutC100, 
                                 paste0(gsub(" ", "", cSpecies), "_",
                                        cSeason, "_RemainingSignals.rds"))
          saveRDS(cOverview, file = OutputFile)
          write.csv(cOverview, 
                    file = gsub("\\.rds", "\\.csv", OutputFile),
                    row.names = F)
      }
  }
