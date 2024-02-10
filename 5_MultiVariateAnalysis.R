# =============================================================================
# Author:           Birgen Haest, birgen.haest@protonmail.com
# Last changed:     February 10, 2020
# =============================================================================
# Functionality:
#   Script to perform multivariate analysis on the output of Spatclimwin
# =============================================================================
# =============================================================================
# Clear memory
# =============================================================================
  env = parent.frame()
  rm(list = ls(all.names = TRUE, env = env), envir = env)
  
# Package requirements and files to source
# =============================================================================
  ReqPackages = c("lme4", "arm", "meifly", "mgcv", "MuMIn", "binhf",
                  "ggplot2", "ggrepel", "ggpmisc", "R.utils", "rgdal", 
                  "tidyr", "plyr", "dplyr", "BHGeometryTools", "relaimpo")
  invisible(lapply(ReqPackages, require, character.only = T, quietly = T))
  source("./functions/BH_AD_Do_GatherResults_Spatial.R")
  source("./functions/BH_CreatePixelsToAddFiles.R")
  source("./functions/BH_AddCorrectPixelNrToCandidatePixels.R")
  source("./functions/BH_AddCandidatesToPhenData.R")
  source("./functions/BH_AD_Do_VarImpAll.R")  
  source("./functions/BH_ProgressReport_DoubleForLoop.R")

# Map the current working directory to a drive to avoid issues with long paths
# =============================================================================
  # if (!("X:" %in% names(System$getMappedDrivesOnWindows()))){
  #   currWD = getwd()
  #   setwd(System$mapDriveOnWindows("X:"))
  # }
  
# =============================================================================
# Global Input Settings
# =============================================================================
  MainInputDir  = "./results/3_SpaWin_Adj"
  MainOutputDir = "./results/4_MultiVar"
  
  RespVar       = "MADJDAYSWS"
  # RespVar       = "WinterPop"
  # PredictorVariables = c("SpatSST_North_SST_pixel24_Open58_Close37_mean")
  
# =============================================================================
# Option Settings
# =============================================================================
  Do_JoinModelRasters            = 0 #have to run this to get the data structure
  Do_AICOverview                 = 0 #one-by-one
  Do_GetCandidates               = 0 #get the valley
  Do_CreatePixelsToAddFiles      = 0 #use the shapefile to get timing & location - in csv
  Do_AddCandidatesToPhenDatabase = 0   
  Do_DeltaAICcIntOnly            = 0 #phenology ~ climate + int (auto mark the invalid candidates)
  Do_AnalyseCollinearity         = 0 #correlate the remaining variables 
  Do_ReduceWithBoruta            = 0 #get the last 14 winners (for the accuracy for the variable importance analysis)
  Do_VarImpAll                   = 1 #ranking, relative importance (game theory version is fun to look at)
  Do_GetFinalSignals             = 0
  
# =============================================================================
# Option-dependent Input Settings
# =============================================================================
  # Do_JoinModelRasters
  # ===========================================================================
    if(Do_JoinModelRasters){
      Settings_JoinRas  = list(
                 InputFolder  = NA_character_,
                 OutputFolder = file.path(MainOutputDir, "1_JoinResults"),
                 TypesToJoin = c("AICc", 
                                 "R-squares", 
                                 "Slopes",
                                 "_Open.tif",
                                 "_Close.tif",
                                 "tiffs"),
                 ClimVarList  = c("air2m", #name before the model
                                  "tmax2m",
                                  "tmin2m",
                                  "shum.2m",  
                                  "cprat",
                                  "prate",
                                  "WindCF850",
                                  "WindGT850",
                                  "WindCF925",
                                  "WindGT925" , 
                                  "WindCF1000",
                                  "WindGT1000", 
                                  "WA850p",
                                  "WA925p",
                                  "WA1000p",
                                  "aquatemp",
                                  "EVImean",
                                  "EVIslope"), 
                 AllModels    = list(air2m.Models     = c("air.2m_yday50_Climwin"), #name of the folder
                                     tmax2m.Models    = c("air.max_yday50_Climwin"),
                                     tmin2m.Models    = c("air.min_yday50_Climwin"),
                                     shum2m.Models  = c("shum.2m_yday50_Climwin"),
                                     cprat.Models = c("cprat.sfc_yday50_Climwin"),
                                     prate.Models = c("prate.sfc_yday50_Climwin"),
                                     WindCF850.Models = c("WindComingFromTarget_850p_yday50_Climwin"),
                                     WindGT850.Models = c("WindGoingToTarget_850p_yday50_Climwin"),
                                     WindCF925.Models = c("WindComingFromTarget_925p_yday50_Climwin"),
                                     WindGT925.Models = c("WindGoingToTarget_925p_yday50_Climwin"),
                                     WindCF1000.Models     = c("WindComingFromTarget_1000p_yday50_Climwin"),
                                     WindGT1000.Models     = c("WindGoingToTarget_1000p_yday50_Climwin"),
                                     WA850.Models     = c("WindAssist_850p_yday50_Climwin"),
                                     WA925.Models     = c("WindAssist_925p_yday50_Climwin"),
                                     WA1000.Models     = c("WindAssist_1000p_yday50_Climwin"),
                                     aquatemp.Models     = c("AquaTemp_yday50_Climwin"),
                                     EVImean.Models     = c("EVI_mean_yday50_Climwin"),
                                     EVIslope.Models     = c("EVI_slope_yday50_Climwin"))
                 )
    }
  
  # Do_AICOverview
  # ===========================================================================
    if(Do_AICOverview){
      Settings_AIC = list(Path_Main        = file.path(MainOutputDir, "1_JoinResults"),
                          Path_Output      = file.path(MainOutputDir, "2_AICOverview"), 
                          DirsToInclude    = c("AICc", 
                                               "R-squares", 
                                               "Slopes",
                                               "WOpen",
                                               "WClose"),
                          RefDate           = "2019-09-15",  # SummerPop
                          # RefDate           = "2019-02-01",  # WinterPop
                          ShapeWithLocation = "./data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp", 
                          ShapeLayerColumn  = c("ne_10m_admin_0_countries", 
                                                "SOVEREIGNT"), 
                          Do_CI             = FALSE,
                          CI_SubsetToIgnore = NA_character_)
    }
  
  # Do_GetCandidates
  # ===========================================================================
    if(Do_GetCandidates){
      Settings_Candid = list(Path_Mosaics      = file.path(MainOutputDir, "1_JoinResults"),
                             Path_Output       = file.path(MainOutputDir, "3_CandidatePixels"), 
                             DirsToInclude     = c("AICc"), 
                             Path_AICcOverview = file.path(MainOutputDir, "2_AICOverview"))
    }
 
  # Do_CreatePixelsToAddFiles
  # ===========================================================================
    if(Do_CreatePixelsToAddFiles){
      Settings_PixelsAdd = list(Path_CandPixels  = file.path(MainOutputDir, "3_CandidatePixels"),
                                OutputFile       = file.path(MainOutputDir, "3_CandidatePixels"), 
                                TargetCoords     = "44.9072, -84.7197")   # Latlon of KAPX station location (mid of GreatLakes)
    }
  
  # Do_AddCandidatesToPhenDatabase
  # ===========================================================================
    if (Do_AddCandidatesToPhenDatabase){
      InputFileBase_AddCand  = "./data/PhenologyDataset.rds"
      OutputFilebase_AddCand = file.path("./results/4_MultiVar/4_PhenDataWithCandidates", 
                                         paste0("3_PhenologyDataset_Updated.rds"))
      ModelClimDataLink      = list(wVar   = c("air2m_yday_50",
                                               "tmax2m_yday_50",
                                               "tmin2m_yday_50",
                                               "shum.2m_yday_50",  
                                               "cprat_yday_50",
                                               "prate_yday_50",
                                               "WindCF850_yday_50",
                                               "WindGT850_yday_50",
                                               "WindCF925_yday_50",
                                               "WindGT925_yday_50" , 
                                               "WindCF1000_yday_50",
                                               "WindGT1000_yday_50", 
                                               "WA850p_yday_50",
                                               "WA925p_yday_50",
                                               "WA1000p_yday_50",
                                               "aquatemp_yday_50",
                                               "EVImean_yday_50",
                                               "EVIslope_yday_50"),
                                    ClimData = c("./data/ClimateDatabase/2_ClimateDatabase_2000to2020_Temp_AvDay.rds",
                                                 "./data/ClimateDatabase/2_ClimateDatabase_2000to2020_Temp_AvDay.rds",
                                                 "./data/ClimateDatabase/2_ClimateDatabase_2000to2020_Temp_AvDay.rds",
                                                 "./data/ClimateDatabase/3_ClimateDatabase_2000to2020_Prec_DaySum.rds",
                                                 "./data/ClimateDatabase/3_ClimateDatabase_2000to2020_Prec_DaySum.rds",
                                                 "./data/ClimateDatabase/3_ClimateDatabase_2000to2020_Prec_DaySum.rds",
                                                 "./data/ClimateDatabase/4_ClimateDatabase_2000to2020_Wind_850p_Daytime.rds",
                                                 "./data/ClimateDatabase/4_ClimateDatabase_2000to2020_Wind_850p_Daytime.rds",
                                                 "./data/ClimateDatabase/4_ClimateDatabase_2000to2020_Wind_925p_Daytime.rds",
                                                 "./data/ClimateDatabase/4_ClimateDatabase_2000to2020_Wind_925p_Daytime.rds",
                                                 "./data/ClimateDatabase/4_ClimateDatabase_2000to2020_Wind_1000p_Daytime.rds",
                                                 "./data/ClimateDatabase/4_ClimateDatabase_2000to2020_Wind_1000p_Daytime.rds",
                                                 "./data/ClimateDatabase/4_ClimateDatabase_2000to2020_WindAssist_850_925_1000hPa_Daytime.rds",
                                                 "./data/ClimateDatabase/4_ClimateDatabase_2000to2020_WindAssist_850_925_1000hPa_Daytime.rds",
                                                 "./data/ClimateDatabase/4_ClimateDatabase_2000to2020_WindAssist_850_925_1000hPa_Daytime.rds",
                                                 "./data/ClimateDatabase/5_MYD21A2_GreatLakesTemp_interpolated_rescaled.rds",
                                                 "./data/ClimateDatabase/6_MOD13A2_EVI_interpolated_rescaled.rds",
                                                 "./data/ClimateDatabase/6_MOD13A2_EVI_interpolated_rescaled.rds"))
      SelectionFile          = file.path("./results/4_MultiVar/3_CandidatePixels")
      
      # MultiRefDates          = list(SummerPop = c(1, 9),     # Name of list variable should be your seasons
      #                               WinterPop = c(1, 2))
      RefDate                = c(15,9)      # c(day, month), e.g. July 1 -> c(1,7)    # Summer population
      # RefDate                = c(1,2)      # c(day, month), e.g. July 1 -> c(1,7)    # Winter population
      # NewVarNamePrefix       = "Temp_"    
      SpeciesList            = c("yday_50")
    }
  
  # Do_DeltaAICcIntOnly
  # ===========================================================================
    if (Do_DeltaAICcIntOnly){
      PhenData_Path         = "./results/4_MultiVar/4_PhenDataWithCandidates"
      OutputFolderPath_DIO  = "./results/4_MultiVar/5_DeltaAICcIntOnly"
      CandPixelOverviewPath = "./results/4_MultiVar/3_CandidatePixels/REPLACERAND/4_CandidatePixels_ToAddToPhen.csv"
    }
  
  # Do_AnalyseCollinearity
  # ===========================================================================
    if (Do_AnalyseCollinearity){
      PhenData_Path        = "./results/4_MultiVar/4_PhenDataWithCandidates"
      RemainSignals_Path   = "./results/4_MultiVar/5_DeltaAICcIntOnly"
      OutputFolderPath_Col = "./results/4_MultiVar/6_CollinearityAnalysis"
    }
  
  # Do_ReduceWithBoruta
  # ===========================================================================
    if (Do_ReduceWithBoruta){
      RemSignalsPath_Bor = "./results/4_MultiVar/6_CollinearityAnalysis"
      MainOutBor         = "./results/4_MultiVar/7_BorutaReduction"
      PhenDataPath_Bor   = "./results/4_MultiVar/4_PhenDataWithCandidates"
      MaxSignalsOut      = 14
      BorutaAll_Settings = list(InputDatabasePath_AD     = NULL,
                                OutputFolder             = NULL,
                                # ResponseVariables        = RespVar,
                                ResponseVariables        = NULL,
                                PredictorVariables       = NULL,
                                SubsetVarCond_One        = list(TRUE, "Species", ""),
                                SubsetVarCond_Two        = list(FALSE, "", ""),
                                PerSpecies               = FALSE, 
                                Do_Standardise_PredVars  = TRUE,
                                Do_AIC_VarImp            = FALSE,
                                AIC_VarImp_Settings      = list(Do_CompareAllPerm       = TRUE,
                                                                VarsToInvert            = "",
                                                                VarsToLog10             = "",
                                                                Interactions            = "", 
                                                                Do_Standardise_RespVar  = FALSE,
                                                                Do_MakePositive         = FALSE, 
                                                                LimitModels             = list(1, "NrVariables", 4),
                                                                MakeAddedVariablePlots  = FALSE,
                                                                MakePredictionPlot      = FALSE,
                                                                MakeDiagnosticsPlot     = FALSE,
                                                                LeaveOneOutVarImp       = FALSE, 
                                                                DoDredgeBest            = FALSE),  
                                Do_RandomForest          = FALSE, 
                                RF_VarImp_Settings       = NULL, 
                                Do_LinRelaImpo           = FALSE,
                                RelaImpo_VarImp_Settings = list(RI_type     = c("lmg"), 
                                                                BootstrapCI = FALSE), 
                                Do_Boruta                = TRUE,
                                Boruta_Settings          = list(MaxRuns           = 1000,
                                                                doTraceSetting    = 1),
                                Do_LASSO                 = FALSE,
                                LASSO_Settings           = NULL,
                                Gather_VarImp            = FALSE)
    }
  
  # Do_VarImpAll
  # ===========================================================================
    if (Do_VarImpAll){
      RemSignalsPath = "./results/4_MultiVar/7_C100_Confirmation/"
      MainOutVarImp  = "./results/4_MultiVar/8_VarImp"
      PhenDataPath   = "./results/4_MultiVar/4_PhenDataWithCandidates"
      VarImpAll_Settings  = list(InputDatabasePath_AD     = NULL,
                                 OutputFolder             = NULL,
                                 # ResponseVariables        = RespVar,
                                 ResponseVariables        = NULL,
                                 PredictorVariables       = NULL,
                                 SubsetVarCond_One        = list(FALSE, "Species", ""),
                                 SubsetVarCond_Two        = list(FALSE, "", ""),
                                 PerSpecies               = FALSE, 
                                 Do_Standardise_PredVars  = TRUE,
                                 Do_AIC_VarImp            = TRUE, #METHOD 1
                                 AIC_VarImp_Settings      = list(Do_CompareAllPerm       = TRUE,
                                                                 VarsToInvert            = "",
                                                                 VarsToLog10             = "",
                                                                 Interactions            = "", 
                                                                 Do_Standardise_RespVar  = FALSE,
                                                                 Do_MakePositive         = FALSE, 
                                                                 LimitModels             = list(1, "NrVariables", 3),
                                                                 MakeAddedVariablePlots  = FALSE,
                                                                 MakePredictionPlot      = FALSE,
                                                                 MakeDiagnosticsPlot     = FALSE,
                                                                 LeaveOneOutVarImp       = FALSE, 
                                                                 DoDredgeBest            = FALSE),  
                                 Do_RandomForest          = FALSE, 
                                 RF_VarImp_Settings       = NULL, 
                                 Do_LinRelaImpo           = TRUE, #METHOD 2
                                 RelaImpo_VarImp_Settings = list(RI_type     = c("lmg"), 
                                                                 BootstrapCI = FALSE), 
                                 Do_Boruta                = TRUE, #METHOD 3
                                 Boruta_Settings          = list(MaxRuns           = 1000,
                                                                 doTraceSetting    = 1),
                                 Do_LASSO                 = FALSE,
                                 LASSO_Settings           = NULL,
                                 Gather_VarImp            = TRUE)
    }
  
  # Do_GetFinalSignals
  # ===========================================================================
    if (Do_GetFinalSignals){
      PhenData_Path        = "./results/4_MultiVar/4_PhenDataWithCandidates"
      RemainSignals_Path   = "./results/4_MultiVar/7_BorutaReduction"
      VarImp_Path          = "./results/4_MultiVar/8_VarImp"
      ShapeFile_Path       = "./results/4_MultiVar/3_CandidatePixels/RANDTOCHANGE/2_CandidatePixels_Filtered"
      OutputDir_FinalSig   = "./results/5_FinalSignals"
      
      # NrMostImp = 4
      NrMostImp = 3   # This is dependent on the multivariate analysis, how much you set to include maximum in each model
    }
 
# =============================================================================
# =============================================================================
# START PROCESSING
# =============================================================================
# Get list of the seasons
# =============================================================================
  allSeasons = list.files(MainInputDir, include.dirs = T)

# =============================================================================
# Do_JoinModelRasters, if requested
# =============================================================================
  if(Do_JoinModelRasters){
    # Do for each of the seasons
    # =========================================================================
      for (cSeason in allSeasons){
        Settings_JoinRas$InputFolder  = file.path(MainInputDir, cSeason)
        Settings_JoinRas$OutputFolder = file.path(MainOutputDir, "1_JoinResults", cSeason)
        BH_AD_Do_GatherResults_Spatial(GatherType      = "JoinModelRasters", 
                                       Gather_Settings = Settings_JoinRas)
      }
  }

# Do_AICOverview
# =============================================================================
  if(Do_AICOverview){
    # Do for each of the seasons
    # =========================================================================
      for (cSeason in allSeasons){
        Settings_AIC$Path_Main   = file.path(MainOutputDir, "1_JoinResults", cSeason)
        Settings_AIC$Path_Output = file.path(MainOutputDir, "2_AICOverview", cSeason)
        BH_AD_Do_GatherResults_Spatial(GatherType      = "AICOverview", 
                                       Gather_Settings = Settings_AIC)
      }
  }

# Do_GetCandidates
# =============================================================================
  if(Do_GetCandidates){
    # Do for each of the seasons
    # =========================================================================
      for (cSeason in allSeasons){
        Settings_Candid$Path_Mosaics      = file.path(MainOutputDir, "1_JoinResults", cSeason)
        Settings_Candid$Path_AICcOverview = file.path(MainOutputDir, "2_AICOverview", cSeason)
        Settings_Candid$Path_Output       = file.path(MainOutputDir, "3_CandidatePixels", cSeason)
        BH_AD_Do_GatherResults_Spatial(GatherType      = "SelectCandidatePixels", 
                                       Gather_Settings = Settings_Candid)
      }
  }

# Do_CreatePixelsToAddFiles
# =============================================================================
  if(Do_CreatePixelsToAddFiles){
    # Do for each of the seasons
    # =========================================================================
      for (cSeason in allSeasons){
        Settings_PixelsAdd$Path_CandPixels = file.path(MainOutputDir, "3_CandidatePixels", cSeason, "2_CandidatePixels_Filtered")
        Settings_PixelsAdd$OutputFile      = file.path(MainOutputDir, "3_CandidatePixels", cSeason, "4_CandidatePixels_ToAddToPhen.csv")
        BH_CreatePixelsToAddFiles(Settings_PixelsAdd)
      }
  }

# Do_AddCandidatesToPhenDatabase
# =============================================================================
  if (Do_AddCandidatesToPhenDatabase){
    # Do for each of the seasons
    # =========================================================================
      for (cSeason in allSeasons){
        SelectionFile          = file.path("./results/4_MultiVar/3_CandidatePixels", 
                                           cSeason, "4_CandidatePixels_ToAddToPhen.csv")
        NewVarNamePrefix       = cSeason
        
        RefDate = RefDate #MultiRefDates[[cSeason]]
        
        BH_AddCandidatesToPhenData(InputFileBase_AddCand, 
                                   OutputFilebase_AddCand,
                                   ModelClimDataLink, 
                                   SelectionFile, 
                                   RefDate, 
                                   NewVarNamePrefix, 
                                   cSeason, 
                                   allSeasons)
        }
  }

# Do_DeltaAICcIntOnly
# =============================================================================
  if (Do_DeltaAICcIntOnly){
    # Get list of files to process (i.e. seasons)
    # =======================================================================
      fileList = list.files(path = PhenData_Path, 
                            pattern = ".rds", 
                            full.names = T)
      
    # Do for each file (i.e. species)
    # =======================================================================
      for (cFile in fileList){
        # Check from which season the file is
        # ===================================================================
          found = 0
          for (i in 1:length(allSeasons)){
            tmp = grepl(allSeasons[i], cFile)
            if ((tmp == TRUE) & (found == 0)){
              cSeason = allSeasons[i]
              found == 1
            } else if ((tmp == TRUE) & (found == 1)){
              error("Something is wrong with finding the season..")
            }
          }
          RespVar = cSeason
        
        # Read the file 
        # ===================================================================
          cData = readRDS(cFile)
          cData = cData[with(cData, order(Year)), ]
          
        # Subset to variables pertaining to the current season only
        # ===================================================================
          RandResInd = grep(cSeason, colnames(cData))
          cData = cData[, c(1:(RandResInd[1]-1), grep(cSeason, colnames(cData)))]
          
        # Subset to variables pertaining to the one Species only
        # ===================================================================
          cData = cData %>% filter(Species == "yday_50")
          
        # Rename response variable in data to "RespVar"
        # ===================================================================
          colnames(cData)[colnames(cData) %in% RespVar] <- "RespVar"
        
        # Get AICc for intercept only model
        # ===================================================================
          BaseModel = lm("RespVar ~ 1", data = cData)
          # BaseModel = MASS::glm.nb("RespVar ~ 1", data = cData)
          BaseAICc  = AICc(BaseModel)  
          
        # Set names of Predictorvariables
        # ===================================================================
          PredVars = colnames(cData)[grep(cSeason, colnames(cData))]
          
        # Create dataframe to hold all aiccs
        # ===================================================================
          AllAICc = data.frame(VarName         = c("Intercept Only", PredVars), 
                               AICc            = NA_real_, 
                               DiffWithIntOnly = NA_real_)
          AllAICc$AICc[1] = BaseAICc
          
        # Get all AICc 
        # ===================================================================
          for (cVar in PredVars){
            tmpData = cData
            colnames(tmpData)[colnames(tmpData) %in% cVar] <- "cVar"
            cModel  = lm("RespVar ~ cVar", data = tmpData)
            #cModel  = MASS::glm.nb("RespVar ~ cVar", data = tmpData)
            AllAICc$AICc[AllAICc$Var %in% cVar] = AICc(cModel)
          }  
          
        # Calculate Delta AICcs
        # ===================================================================
          AllAICc$DiffWithIntOnly = AllAICc$AICc - AllAICc$AICc[1]
        
        # Get the name of the current species
        # ===================================================================
          cSpecies = as.character(unique(cData$Species))
          
        # Save to file
        # ===================================================================
          dir.create(OutputFolderPath_DIO, showWarnings = F, recursive = T)
          OutputFile = file.path(OutputFolderPath_DIO, 
                                 paste0(gsub(" ", "", cSpecies), "_",
                                        cSeason, "_DeltaAICcWithIntOnly.rds"))
          saveRDS(AllAICc, file = OutputFile)
          write.csv(AllAICc, file = gsub("\\.rds", "\\.csv", OutputFile), 
                    row.names = F)
          
        # Create overview of remaining candidates
        # ===================================================================
          cCandPixFile      = gsub("REPLACERAND", cSeason, CandPixelOverviewPath)
          cOverview         = read.csv(cCandPixFile)
          cOverview         = cOverview[cOverview$Species %in% cSpecies, ]
          cOverview$VarName = gsub(" ", "", paste0(cSeason, "_", cOverview$Model, "_",
                                            # cOverview$GridLoc, "_pixel", cOverview$PixelNr,
                                            "pixel", cOverview$PixelNr,
                                            "_Open", cOverview$WOpen, 
                                            "_Close", cOverview$WClose, 
                                            "_", cOverview$OperationOverTime))
          
        # Flag the ones removed because of Delta AICc with IntOnly > -2
        # ===================================================================
          RemainingSignals = as.character(AllAICc$Var[AllAICc$DiffWithIntOnly <= -2])
          SignalsToRemove  = as.character(AllAICc$Var[AllAICc$DiffWithIntOnly > -2])
          cOverview = merge(cOverview, AllAICc[, c("VarName", "DiffWithIntOnly")])
          cOverview$RemovedPriorToVarImp = NA_character_
          cOverview[cOverview$VarName %in% SignalsToRemove, "RemovedPriorToVarImp"] = "YES"
          cOverview$Reason = NA_character_
          cOverview[cOverview$VarName %in% SignalsToRemove, "Reason"] = "Delta AICc compared to Intercept-Only model > -2"
          
        # Save overview of remaining candidates
        # ===================================================================
          OutputFile = file.path(OutputFolderPath_DIO, 
                                 paste0(gsub(" ", "", cSpecies), "_",
                                        cSeason, "_RemainingSignals.rds"))
          saveRDS(cOverview, file = OutputFile)
          write.csv(cOverview, file = gsub("\\.rds", "\\.csv", OutputFile),
                    row.names = F)
      }
  }

# Do_AnalyseCollinearity
# =========================================================================
  if (Do_AnalyseCollinearity){
    # Get list of files to process (i.e. seasons)
    # =====================================================================
      fileList = list.files(path = PhenData_Path, 
                            pattern = ".rds", 
                            full.names = T)
      
    # Create output directory, if required
    # =====================================================================
      dir.create(OutputFolderPath_Col, showWarnings = F, recursive = T)
      
    # Do for each file (i.e. species)
    # =====================================================================
      for (cFile in fileList){
        # Check from which season the file is
        # ===================================================================
          found = 0
          for (i in 1:length(allSeasons)){
            tmp = grepl(allSeasons[i], cFile)
            if ((tmp == TRUE) & (found == 0)){
              cSeason = allSeasons[i]
              found == 1
            } else if ((tmp == TRUE) & (found == 1)){
              error("Something is wrong with finding the season..")
            }
          }
        
        # Read the phenology file 
        # =================================================================
          cData = readRDS(cFile)
          cData = cData[with(cData, order(Year)), ]
          
        # Get the name of the current species
        # =================================================================
          cSpecies = as.character(unique(cData$Species))
          
        # Read the remaining signals overview file
        # =================================================================
          cRemSignalsFile = file.path(RemainSignals_Path, 
                                      paste0(gsub(" ", "", cSpecies), "_",
                                             cSeason, "_RemainingSignals.rds"))
          cOverview = readRDS(cRemSignalsFile)
          
        # Determine current amount of remaining signals
        # =================================================================
          RemSignals   = cOverview$VarName[is.na(cOverview$RemovedPriorToVarImp)]
          NrRemSignals = length(RemSignals)
          
        # CASE: no remaining signals, just copy the overview file
        # =================================================================
          if (NrRemSignals == 0){
            # Nothing is done here, because the saving is done after the 
            # if statement
        # CASE: Signals remaining, check for collinearity
        # =================================================================
          } else {
            # Subset phenology data to remaining variables
            # =============================================================
              cData = cData[, RemSignals]
              
            # Sort remaining signals according to DeltaAICc compared to IntOnly model
            # =============================================================
              RemSignals = RemSignals[order(cOverview$DiffWithIntOnly[cOverview$VarName %in% RemSignals])]
              
            # Make sure all variables are numeric, and sort columns (similar to above)  
            # =============================================================
              cDataNumeric = as.data.frame(apply(cData, 2, "as.numeric"))
              cDataNumeric = cDataNumeric[, RemSignals]
              
            # Calculate Perason correlation coefficients
            # =============================================================
              CurrentCorrCoef = cor(cDataNumeric, 
                                    method = "pearson", 
                                    use = "complete.obs")
              
            # Check on which variables to exclude from analysis based on collinearity
            # =============================================================
              SignalsToExclude = NULL
              for (cVarLoc in 1:(NrRemSignals-1)){
                cVar           = RemSignals[cVarLoc]
                # Check whether the variable has been excluded already,
                # and if so, skip
                # =========================================================
                  if (cVar %in% SignalsToExclude$VarName){
                    next()
                  }
                LocsToCheck    = (cVarLoc+1):NrRemSignals
                SignalsToCheck = RemSignals[LocsToCheck]   
                PosCorrLoc     = CurrentCorrCoef[cVarLoc, LocsToCheck] >= 0.7
                NegCorrLoc     = CurrentCorrCoef[cVarLoc, LocsToCheck] <= -0.7
                if (length(SignalsToCheck[as.vector(PosCorrLoc | NegCorrLoc)]) > 0){
                  tmpSignalstoExclude = data.frame(VarName = SignalsToCheck[as.vector(PosCorrLoc | NegCorrLoc)], 
                                                   ColWith = cVar)
                  SignalsToExclude    = rbind(SignalsToExclude, tmpSignalstoExclude)
                }
              }
            
            # Remove duplicates 
            # =============================================================
              if (!is.null(SignalsToExclude)){
                SignalsToExclude = SignalsToExclude[!duplicated(SignalsToExclude$VarName), ]  
              }
              
            # Update the overview of candidate signals
            # =============================================================
              if (!is.null(SignalsToExclude)){
                cOverview[cOverview$VarName %in% SignalsToExclude$VarName, "RemovedPriorToVarImp"] = "YES"
                for (cLoc in 1:nrow(SignalsToExclude)){
                  cOverview[cOverview$VarName %in% SignalsToExclude$VarName[cLoc], "Reason"] = paste0("Highly collinear with ", 
                                                                                                      SignalsToExclude$ColWith[cLoc])
                }
              }
              
            # Save collinearity values to file
            # =============================================================
              OutputFile = file.path(OutputFolderPath_Col, 
                                     paste0(gsub(" ", "", cSpecies), "_",
                                            cSeason, "_Collinearity.rds"))
              saveRDS(CurrentCorrCoef, file = OutputFile)
              write.csv(CurrentCorrCoef, 
                        file = gsub("\\.rds", "\\.csv", OutputFile), 
                        row.names = F)
          }
          
        # Save overview of remaining candidates
        # =================================================================
          OutputFile = file.path(OutputFolderPath_Col, 
                                 paste0(gsub(" ", "", cSpecies), "_",
                                        cSeason, "_RemainingSignals.rds"))
          saveRDS(cOverview, file = OutputFile)
          write.csv(cOverview, 
                    file = gsub("\\.rds", "\\.csv", OutputFile),
                    row.names = F)
      }
  }

# Do_ReduceWithBoruta, if requested
# =========================================================================
  if (Do_ReduceWithBoruta){
    # Get list of files to process (i.e. species)
    # =====================================================================
      fileList = list.files(path = PhenDataPath_Bor, 
                            pattern = ".rds", 
                            full.names = T)
      
    # Create output directory, if required
    # =====================================================================
      dir.create(MainOutBor, showWarnings = F, recursive = T)
      
    # Do for each file (i.e. species)
    # =====================================================================
      for (cFile in fileList){
        # Check from which season the file is
        # ===================================================================
          found = 0
          for (i in 1:length(allSeasons)){
            tmp = grepl(allSeasons[i], cFile)
            if ((tmp == TRUE) & (found == 0)){
              cSeason = allSeasons[i]
              found == 1
            } else if ((tmp == TRUE) & (found == 1)){
              error("Something is wrong with finding the season..")
            }
          }
        
        # Get the name of the current species
        # =================================================================
          cData    = readRDS(cFile)
          cSpecies = as.character(unique(cData$Species))
          
        # Read the remaining pixels overview file
        # =================================================================
          cOverview = readRDS(file.path(RemSignalsPath_Bor, 
                              paste0(gsub(" ", "", cSpecies), "_",
                                     cSeason, "_RemainingSignals.rds")))
          RemSignals = cOverview$VarName[is.na(cOverview$RemovedPriorToVarImp)]
          
        # CASE: No remaining signals <= MaxSignalsOut -> Just copy the file with the 
        #       remaining signals
        # =================================================================
          if (length(RemSignals) <= MaxSignalsOut){
            # Do NOTHING (copying is done later)
        # CASE: Remaining signals -> Do Boruta reduction
        # =================================================================
          } else {
            # Update Boruta settings
            # =============================================================
              BorutaAll_Settings$InputDatabasePath_AD   = cFile
              BorutaAll_Settings$PredictorVariables     = RemSignals
              BorutaAll_Settings$SubsetVarCond_One[[3]] = cSpecies
              BorutaAll_Settings$OutputFolder           = file.path(MainOutBor,
                                                                    cSpecies,
                                                                    cSeason)
              BorutaAll_Settings$ResponseVariables      = cSeason

            # Do Boruta
            # =============================================================
              BH_AD_Do_VarImpAll(BorutaAll_Settings)
              
            # Read Boruta output file
            # =============================================================
              borFile = file.path(MainOutBor, cSpecies, cSeason, 
                                  "VarImp", "Boruta", "borutaDecision.csv")
              borResults = read.csv(borFile)
              
            # Get locations of those that were rejected
            # =============================================================
              indRejected = which(borResults$decision %in% "Rejected")
              
            # Update the overview of candidate signals
            # =============================================================
              if (length(indRejected) != 0){
                VarsToRemove = borResults$Variable[indRejected]
                cOverview[cOverview$VarName %in% VarsToRemove, "RemovedPriorToVarImp"] = "YES"
                cOverview[cOverview$VarName %in% VarsToRemove, "Reason"] = paste0("Boruta reduction of more than 14 remaining candidates")
              }
              
            # Check if number of signals has been reduced to MaxSignalsOut,
            # if not reduce to MaxSignalsOut based on boruta varimp
            # =============================================================
              cRemSignals = cOverview$VarName[is.na(cOverview$RemovedPriorToVarImp)]
              cBorRestuls = borResults[-indRejected, ]
              if (length(cRemSignals) >= MaxSignalsOut){
                VarsToRemove = cBorRestuls$Variable[which((max(rank(cBorRestuls$meanImp)) + 1 - rank(cBorRestuls$meanImp)) > MaxSignalsOut)]
                cOverview[cOverview$VarName %in% VarsToRemove, "RemovedPriorToVarImp"] = "YES"
                cOverview[cOverview$VarName %in% VarsToRemove, "Reason"] = paste0("Boruta reduction of more than 14 remaining candidates")
              }  
          } 
           
        # Save overview of remaining candidates
        # =================================================================
          OutputFile = file.path(MainOutBor, 
                                 paste0(gsub(" ", "", cSpecies), "_",
                                        cSeason, "_RemainingSignals.rds"))
          saveRDS(cOverview, file = OutputFile)
          write.csv(cOverview, 
                    file = gsub("\\.rds", "\\.csv", OutputFile),
                    row.names = F)
      }
  }

# Do_VarImpAll, if requested
# =========================================================================
  if (Do_VarImpAll){
    # Get list of files to process (i.e. species)
    # =====================================================================
      fileList = list.files(path = PhenDataPath, 
                            pattern = ".rds", 
                            full.names = T)
      
    # Create output directory, if required
    # =====================================================================
      dir.create(MainOutVarImp, showWarnings = F, recursive = T)
      
    # Do for each file (i.e. species)
    # =====================================================================
      for (cFile in fileList){
        # Check from which season the file is
        # ===================================================================
          found = 0
          for (i in 1:length(allSeasons)){
            tmp = grepl(allSeasons[i], cFile)
            if ((tmp == TRUE) & (found == 0)){
              cSeason = allSeasons[i]
              found == 1
            } else if ((tmp == TRUE) & (found == 1)){
              error("Something is wrong with finding the season..")
            }
          }
        
        # Get the name of the current species
        # =================================================================
          cData    = readRDS(cFile)
          cSpecies = as.character(unique(cData$Species))
          
        # Read the remaining pixels overview file
        # =================================================================
          # cOverview = readRDS(file.path(RemSignalsPath, 
          #                     paste0(gsub(" ", "", cSpecies), "_",
          #                            cSeason, "_RemainingSignals.rds")))
          cOverview = read.csv(file.path(RemSignalsPath, 
                                        paste0(gsub(" ", "", cSpecies), "_",
                                               cSeason, "_RemainingSignals.csv")))
          RemSignals = cOverview$VarName[is.na(cOverview$RemovedPriorToVarImp)]
          
        # CASE: No or only 1 remaining signal(s) -> Don't do anything
        # =================================================================
          if (length(RemSignals) <= 1){
            # Do NOTHING
        # CASE: Remaining signals -> Analyse variable importance    
        # =================================================================
          } else {
            # Update Varimp settings
            # =============================================================
              VarImpAll_Settings$InputDatabasePath_AD   = cFile
              VarImpAll_Settings$PredictorVariables     = RemSignals
              VarImpAll_Settings$SubsetVarCond_One[[3]] = cSpecies
              VarImpAll_Settings$OutputFolder           = file.path(MainOutVarImp,
                                                                    cSpecies,
                                                                    cSeason)
              VarImpAll_Settings$ResponseVariables      = cSeason

            # Do VarImp
            # =============================================================
              BH_AD_Do_VarImpAll(VarImpAll_Settings)
          }
      }
  }

# Do_GetFinalSignals
# ===========================================================================
  if (Do_GetFinalSignals){
    # Get list of files to process (i.e. species)
    # =====================================================================
      fileList = list.files(path = PhenData_Path, 
                            pattern = ".rds", 
                            full.names = T)
      
    # Create output directory, if required
    # =====================================================================
      dir.create(OutputDir_FinalSig, showWarnings = F, recursive = T)
      
    # Do for each file (i.e. species)
    # =====================================================================
      for (cFile in fileList){
        # Check from which season the file is
        # ===================================================================
          found = 0
          for (i in 1:length(allSeasons)){
            tmp = grepl(allSeasons[i], cFile)
            if ((tmp == TRUE) & (found == 0)){
              cSeason = allSeasons[i]
              found == 1
            } else if ((tmp == TRUE) & (found == 1)){
              error("Something is wrong with finding the season..")
            }
          }
          RespVar = cSeason
          
        # Get the name of the current species
        # =================================================================
          cData    = readRDS(cFile)
          cSpecies = as.character(unique(cData$Species))
          
        # Read the remaining pixels overview file
        # =================================================================
          filePath_Overview = file.path(RemainSignals_Path, 
                                        paste0(gsub(" ", "", cSpecies), "_",
                                               cSeason, "_RemainingSignals.rds"))
          cOverview = readRDS(filePath_Overview)
          RemSignals = cOverview[is.na(cOverview$RemovedPriorToVarImp), ]
          
        # CASE: No remaining signals -> Don't do anything
        # =================================================================
          if (nrow(RemSignals) == 0){
            # Do NOTHING
        # CASE: Remaining signals
        # =================================================================
          } else {
            # SUBCASE 1: # remaining signals == 1
            # =============================================================
              if (nrow(RemSignals) == 1){
                # Do not update cOverview
            # SUBCASE 2: # remaining signals > 1: Add varimp to cOverview
            # =============================================================  
              } else {
                # Read the variable importance overview file
                # =============================================================
                  VarImp_File = file.path(VarImp_Path, cSpecies, cSeason, 
                                          "VarImp", "Overview_Varimp.csv")
                  cVarImp = read.csv(VarImp_File)
                  colnames(cVarImp)[1] <- "VarName"
                  
                # Merge the varimp with the overview 
                # =============================================================
                  cOverview  = merge(cOverview, cVarImp, all.x = TRUE)
                  RemSignals = merge(RemSignals, cVarImp)
              }
              
            # Save the (updated) overview file 
            # =============================================================
              cOutputFile = file.path(OutputDir_FinalSig, "1_OverviewFiles",
                                      paste0(gsub(" ", "", cSpecies), "_",
                                               cSeason, "_FinalOverview.rds"))
              dir.create(dirname(cOutputFile), showWarnings = F, recursive = T)
              saveRDS(cOverview, cOutputFile)
              write.csv(cOverview, gsub("\\.rds", "\\.csv", cOutputFile), 
                        row.names = F)
              
            # Update the shapefile path to the current randomisation
            # =============================================================
              cShapeFile_Path = gsub("RANDTOCHANGE", cSeason, ShapeFile_Path)
              
            # Do for each of the remaining signals
            # =============================================================
              for (cSignalNr in 1:nrow(RemSignals)){
                # Get the weather variable and model name
                # =========================================================
                  cVarType = RemSignals[cSignalNr, "ClimateVar"]
                  cModel   = RemSignals[cSignalNr, "Model"]
                  bname = strsplit(cModel, "_", fixed = T)[[1]][1]
                  
                # Determine the filename of the shapefile
                # =========================================================
                  if (bname %in% c("air2m", "tmin2m", "tmax2m",
                                      "prate", "cprat", "shum.2m",
                                      "EVImean", "EVIslope"
                                      )){
                    cShapeFile = list.files(path = file.path(cShapeFile_Path, 
                                                             cSpecies),
                                            pattern = paste0(bname, ".shp"), 
                                            full.names = T)
                  } else if (cVarType %in% "WindDirection"){
                    cPresLevel = gsub("[^0-9.]", "",  bname)
                    OpTime     = RemSignals[cSignalNr, "OperationOverTime"]
                    if (grepl("WindComing", OpTime)){
                      cPattern = "ComingFrom"
                    } else if (grepl("WindGoing", OpTime)){
                      cPattern = "GoingTo"
                    }
                    cShapeFile = list.files(path = file.path(cShapeFile_Path, 
                                                             cSpecies),
                                            pattern = paste0(bname, ".shp"),
                                            full.names = T)
                    # cShapeFile = grep(bname, cShapeFile, value = TRUE)
                    # cShapeFile = grep(cPresLevel, cShapeFile, value = TRUE)
                  } else if (bname %in% c("WA1000p", "WA850p", "WA925p")){
                    cPresLevel = gsub("[^0-9.]", "",  cVarType)
                    cShapeFile = list.files(path = file.path(cShapeFile_Path, 
                                                             cSpecies),
                                            pattern = paste0(".shp"), 
                                            full.names = T)
                    cShapeFile = grep(bname, cShapeFile, value = TRUE)
                    cShapeFile = grep(cPresLevel, cShapeFile, value = TRUE)
                  } 
                  
                # Read the shapefile
                # =========================================================
                  cShapeData = readOGR(dsn   = cShapeFile, 
                                       layer = gsub("\\.shp", "", basename(cShapeFile)))
                  
                # Subset to the current signal, and add to shapefile 
                # with all final signals
                # =========================================================
                  if (nrow(RemSignals) > 1){
                    ColsToAdd = colnames(cVarImp)#[1:ncol(cVarImp)]
                  }
                  if (cSignalNr == 1){
                    # ContColumn = ifelse(grepl("Europe", RemSignals[cSignalNr, "VarName"]), "PxNr_WE", "PxNr_WA")
                    cIndex = cShapeData$PixelNr == RemSignals[cSignalNr, "PixelNr"] &
                             cShapeData$WOpen   == RemSignals[cSignalNr, "WOpen"] &
                             cShapeData$WClose  == RemSignals[cSignalNr, "WClose"]
                            # !(is.na(cShapeData@data[, ContColumn]))
                    cOutputShape = subset(cShapeData, cIndex)
                    # Also add varimp columns if there are more than 1 rem signal
                    if (nrow(RemSignals) > 1){  
                      cOutputShape[1, ColsToAdd] = RemSignals[cSignalNr, ColsToAdd]
                    }
                  } else {
                    # ContColumn = ifelse(grepl("Europe", RemSignals[cSignalNr, "VarName"]), "PxNr_WE", "PxNr_WA")
                    cIndex = cShapeData$PixelNr == RemSignals[cSignalNr, "PixelNr"] &
                             cShapeData$WOpen   == RemSignals[cSignalNr, "WOpen"] &
                             cShapeData$WClose  == RemSignals[cSignalNr, "WClose"]
                             # !(is.na(cShapeData@data[, ContColumn]))
                    tmpShape = subset(cShapeData, cIndex)
                    tmpShape[1, ColsToAdd] = RemSignals[cSignalNr, ColsToAdd]
                    # Make sure binding is done appropriately by deleting rownames of shape data
                      row.names(tmpShape@data) <- NULL
                      row.names(cOutputShape@data) <- NULL
                    cOutputShape = rbind.SpatialPointsDataFrame(cOutputShape, tmpShape)
                  }
              }
              
            # Save the shapefile with all of the remaining signals
            # =============================================================
              cOutputFile = file.path(OutputDir_FinalSig, "2_Shapes", "1_All",
                                      paste0(gsub(" ", "", cSpecies), "_",
                                               cSeason, "_FinalOverview.shp"))
              dir.create(dirname(cOutputFile), showWarnings = F, recursive = T)
              writeOGR(obj = cOutputShape, 
                       dsn = cOutputFile, 
                       layer = gsub("\\.shp", "", basename(cOutputFile)), 
                       driver = "ESRI Shapefile")
              
            # CASE: ONLY 1 Remaining signals
            # =============================================================
              if (nrow(RemSignals) == 1){
                cOutputShape_NrMostImp = cOutputShape
            # CASE: More than 1 Remaining signals
            # =============================================================    
              } else {
                # Read the file with all of the model performances
                # =============================================================
                  cModelOverviewFile = file.path(VarImp_Path, cSpecies, cSeason,
                                                 "VarImp", "# =============================================================================
", "Overview_Models.csv")
                  cModelOverview = read.csv(cModelOverviewFile)
                  
                # Create subset of the remaining signals to the max. # most 
                # important signals (as requested in the input)
                # =============================================================
                  cOutputShape_NrMostImp = cOutputShape[order(cOutputShape$Importance_Mean_All, 
                                                          decreasing = T), ]
                  if (nrow(cOutputShape_NrMostImp) > NrMostImp){
                    cOutputShape_NrMostImp = cOutputShape_NrMostImp[1:NrMostImp, ]
                  }
              }
              
            # Save the shapefile with the max. # most important signals (# is given in input)
            # =============================================================
              cOutputFile = file.path(OutputDir_FinalSig, "2_Shapes", 
                                      paste0("2_Max", NrMostImp, "MostImp"), 
                                      paste0(gsub(" ", "", cSpecies), "_",
                                               cSeason, "_FinalOverview.shp"))
              dir.create(dirname(cOutputFile), showWarnings = F, recursive = T)
              writeOGR(obj = cOutputShape_NrMostImp, 
                       dsn = cOutputFile, 
                       layer = gsub("\\.shp", "", basename(cOutputFile)), 
                       driver = "ESRI Shapefile")
              
            # Get a summary with the R2, adj. R2, and pred. R2s for the 
            # model with the # most important signals
            # =============================================================
              # CASE: ONLY 1 Remaining signals
              # ===========================================================
                if (nrow(RemSignals) == 1){
                  CModelRes.4 = data.frame(X                = NA_integer_,
                                           X.Intercept.     = NA_real_,
                                           ToReplace        = NA_real_,
                                           R.2	             = NA_real_, 
                                           adjR.2           = NA_real_, 
                                           df               = NA_integer_,
                                           logLik           = NA_real_, 
                                           AICc             = NA_real_, 
                                           AllPress         = NA_real_, 
                                           AllPredRsquared  = NA_real_, 
                                           AllPressRsquared = NA_real_,
                                           model	           = NA_integer_,
                                           ModelType        = c("Best4", "TopModel"))
                  cVarName = cOverview[is.na(cOverview$RemovedPriorToVarImp), "VarName"]
                  colnames(CModelRes.4)[3] <- cVarName 
                  # Get the values
                    cFormula = formula(paste0(RespVar, " ~ ", cVarName))
                    CModelRes.4[, "X.Intercept."]    = summary(lm(cFormula, data = cData))$coefficients["(Intercept)", "Estimate"]
                    CModelRes.4[, cVarName]          = summary(lm(cFormula, data = cData))$coefficients[cVarName, "Estimate"]
                    CModelRes.4[, "R.2"]             = summary(lm(cFormula, data = cData))$r.squared
                    CModelRes.4[, "adjR.2"]          = summary(lm(cFormula, data = cData))$adj.r.squared
                    CModelRes.4[, "df"]              = summary(lm(cFormula, data = cData))$df[2]
                    CModelRes.4[, "AICc"]            = MuMIn::AICc(lm(cFormula, data = cData))
                    CModelRes.4[, "AllPredRsquared"] = model_fit_stats(lm(cFormula, data = cData))$pred.r.squared
                     
                  
              # CASE: More than 1 Remaining signals
              # =============================================================   
                } else {
                  cModelInd   = complete.cases(cModelOverview[, cOutputShape_NrMostImp$VarName[!is.na(cOutputShape_NrMostImp$VarName)]])
                  CModelRes.4 = cModelOverview[cModelInd,]
                  # change mod. R2 to adj. R2
                  cFormula = formula(paste0(RespVar, " ~ ", 
                                    paste(cOutputShape_NrMostImp$VarName, collapse = "+")))
                  CModelRes.4[, "adjR.2"] = summary(lm(cFormula, data = cData))$adj.r.squared
                }
              
            # Create subset of the remaining signals consisting of the signals 
            # from the overall top performing model
            # =============================================================
              # CASE: ONLY 1 Remaining signals
              # ===========================================================
                if (nrow(RemSignals) == 1){
                  cOutputShape_Top = cOutputShape
                  
              # CASE: More than 1 Remaining signals
              # =============================================================   
                } else {
                  cTop           = cModelOverview[1, ]
                  cTopVars       = cTop[, grep(cSeason, colnames(cTop))]
                  cTopVars       = colnames(cTopVars)[!(is.na(cTopVars))]
                  cOutputShape_Top = cOutputShape[cOutputShape$VarName %in% cTopVars, ]
                }
              
            # Save the shapefile of top model signals
            # =============================================================
              cOutputFile = file.path(OutputDir_FinalSig, "2_Shapes", "3_TopModel",
                                      paste0(gsub(" ", "", cSpecies), "_",
                                               cSeason, "_FinalOverview.shp"))
              dir.create(dirname(cOutputFile), showWarnings = F, recursive = T)
              writeOGR(obj = cOutputShape_Top, 
                       dsn = cOutputFile, 
                       layer = gsub("\\.shp", "", basename(cOutputFile)), 
                       driver = "ESRI Shapefile")
              
            # Update the summary with the adj. R2 instead of mod. R2 for the 
            # top model
            # =============================================================
              if (nrow(RemSignals) == 1){
                cModelPerf = CModelRes.4
              } else {
                CModelRes.top = cTop
                # change mod. R2 to adj. R2
                cFormula = formula(paste0(RespVar, " ~ ", 
                                  paste(cTopVars, collapse = "+")))
                CModelRes.top[, "adjR.2"] = summary(lm(cFormula, data = cData))$adj.r.squared
                cModelPerf = rbind(CModelRes.4, CModelRes.top)
                cModelPerf$ModelType = c(paste0("Best", NrMostImp), 
                                         "TopModel")
              }
              
            # Save the joined summaries to a csv file
            # =============================================================
              cOutputFile = file.path(OutputDir_FinalSig, "3_PerformanceOverview", 
                                      paste0(gsub(" ", "", cSpecies), "_",
                                               cSeason, "_PerfOverview.csv"))
              dir.create(dirname(cOutputFile), showWarnings = F, recursive = T)
              write.csv(cModelPerf, 
                        file = cOutputFile,
                        row.names = F)
          }
       }
  }
    
# Unmap the X drive used to deal with long filenames
# =============================================================================
  if (exists("currWD")){
    setwd(currWD)
    System$unmapDriveOnWindows("X:")
  }
  
# =============================================================================
# =============================================================================
# THE END
# =============================================================================
# =============================================================================
