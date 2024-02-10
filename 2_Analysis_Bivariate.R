# ============================================================================= 
# Author:           Birgen Haest, birgen.haest@protonmail.com 
# Last changed:     February 25, 2020
# =============================================================================
# Functionality:
#   Script to analyse univariate relationships to phenological dates
# =============================================================================
# =============================================================================
# Package requirements, and files to source
# =============================================================================
  ReqPackages = c("ggplot2", "ggpmisc", "ggrepel",
                  "climwin", "lubridate", 
                  "MASS", "Kendall", "zyp", "ecodist", 
                  "grid", "gridExtra", "tools", 
                  "R.utils", "binhf", 
                  "WaveletComp")
  #install.packages(ReqPackages)
  invisible(lapply(ReqPackages, require, character.only = T, quietly = T))
  source("./functions/BH_OtherFunctions.R")
  source("./functions/BH_ProgressReport_DoubleForLoop.R")

# =============================================================================
# GLOBAL INPUT SETTING 
# =============================================================================
  # Input Path of phenology dataset  
  # ===========================================================================
    InputDatabaseFolderPath = file.path("./data")
    InputDatabasePath       = file.path(InputDatabaseFolderPath,
                                        "PhenologyDataset.rds")
  
  # Output Path
  # ===========================================================================
    OutputFolderPath        = file.path("./results")
  
  # Set Variables
  # ===========================================================================
    ResponseVariables = c("MADJDAYSWS")        
    
    PredictorVariables = c("Year")

# =============================================================================
# OPTION SETTINGS - Set what you want to do 1 (or TRUE) 
# Only one at a time.
# =============================================================================
    Do_AnalyseTrendLinearity       = 0  
         
    Do_Climwin_Spatial             = 0
    Do_GatherResults_Spatial       = 1  
    
    Do_AddCandidatesToPhenDatabase = 0  
    
    Do_GetChangeOverTime           = 0  
    
# =============================================================================
# OPTION SPECIFICATIONS
# =============================================================================
  # Do_AnalyseTrendLinearity
  # ===========================================================================
    if (Do_AnalyseTrendLinearity){
      InputDatabasePath_LT  = InputDatabasePath
      OutputDir_LT          = file.path(OutputFolderPath, "1_TrendAnalysis")
      TrendSettings = list(VarToAnalyse          = "MADJDAYSWS", 
                           VarNameForPlot        = "MADJDAYSWS",
                           TimeVar               = "Year",
                           OverWritePlotFileName = list(FALSE, "MADJDAYSWS"),
                           Do_AICMethod          = TRUE,
                           AIC_Loess_Span        = 1, 
                           Data_Subset_1         = list(FALSE, "Species", c("European Pied Flycatcher")),
                           Data_Subset_2         = list(FALSE, "Year", seq(1968, 2014, 1)),
                           Do_AugDickeyFuller    = TRUE, 
                           ADF_Detrend           = list(TRUE, 1)) # Second element: 1 = linear; 2 = quadratic; 3 = cubic
    }
 
  # Do_Climwin_Spatial
  # Note: Depending on the size of the spatiotemporal weather grid you are analysing, 
  #       this step can require a lot of RAM 
  # ===========================================================================
    if (Do_Climwin_Spatial){
      PerSpeciesOrID         = 1
      SpeciesOrID_VarName    = "Species"
      DateVariable           = "RespVar"   # Variable determing the date in the bird
                                           # database. Set to "RespVar" if it is the same
                                           # as the response variable. This variable should be 
                                           # a winter-solstice based date number.
      
      DataSubset_Species     = list(1, "Species", c("yday_50"))
      DataSubset_Year        = list(1, "Year", seq(2000, 2020))
      
      Database_ClimSpat_Path = file.path("./data/6_MOD13A2_EVI_interpolated_rescaled.rds")
      
      ClimPredictorVar   = c("EVI")

      CW_cinterval       = "day"
      CW_range           = c(258, 0) # 258 go back to Jan 1st 2000 while reference day is Sep 15th; 107 go back to June
      CW_type            = "absolute"
      
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
        
      Ref_LatLon         = c(44.9072, -84.7197)  # Needs to be set when using the Wind function,
                                                 # Latlon of the location where you measured the
                                                 # biological event (e.g. phenology)
                                                 # KAPX station location (mid of GreatLakes)
      
      CW_stat            = c("slope") #mean for temp and precip within the time window; slope for EVI
                                                            #here it calculate the number of days within time window
      CW_func            = c("lin")
      CW_refday          = c(15, 9) # day, month
      
      CW_exclude         = c(14, 0) # the minimum window length
      ModelName          = "EVI_GreatLakes"          # Model name, free to choose to your personal preference 
      
      MaskPath          = file.path("./data/Masks/6_MOD13A2_EVI_interpolated_rescaled_FullMask.tif")
      
      Use_CW_BaseModel   = FALSE      # Can be used to overwrite the 
      CW_BaseModel       = "Year"     # base function used in the climwin analysis
    }
     
  # Do_GatherResults_Spatial
  # ===========================================================================
    if (Do_GatherResults_Spatial){
      # Set GatherType - only set one at a time!
      # For the climwin gathering, keep an eye on your RAM
      # =======================================================================
        GatherType = "ClimWin"
        # GatherType = "JoinModelRasters"
        # GatherType = "AICOverview"
        # GatherType = "SelectCandidatePixels"
        
      # Settings for GatherType: "ClimWin"
      # ======================================================================= 
        if (identical(GatherType, "ClimWin")){
          #InputPathMain = file.path("./results/Climwin")
          InputPathMain = file.path("./results/Climwin/MADJDAYSWS/EVI_GreatLakes")
          GatherClimWin_Settings = list(Do_MakePlots      = 1,
                                        Do_PrintBestModel = 1,
                                        Do_GetpValues     = 1,
                                        Do_PrintCombos    = 1,
                                        ClimWinInputPath  = c(file.path(InputPathMain, 
                                                                      "FullModels"),
                                                              file.path(InputPathMain, 
                                                                      "RandomModels")),
                                        ClimwinOutputPath = file.path("./results/Climwin"),
                                        pValueThreshold   = 0.3,
                                        SpeciesSubset     = list(1, "yday_50"),
                                        ReferenceType     = "Climate",
                                      
                                      # Please uncomment the relevant location-variable combo
                                        
                                      # temp
                                        # ReferenceClimRowCols = c(13, 51),
                                        # ReferenceClimExtent  = c(215.625-1.875/2, 309.375+1.875/2, 39.047-1.9047/2,  61.9033+1.9047/2))  # MinLong, MaxLong, MinLat, MaxLat ### bbox
                                        
                                      # Wind
                                        # ReferenceClimRowCols = c(11, 39),
                                        # ReferenceClimExtent  = c(215 - 2.5/2, 310+2.5/2, 37.5-2.5/2, 62.5+2.5/2))
                              
                                      # precipitation
                                        # ReferenceClimRowCols = c(13, 51),
                                        # ReferenceClimExtent  = c(215.625-1.875/2, 309.375+1.875/2, 39.047-1.9047/2,  61.9033+1.9047/2))
                                      
                                      # EVI
                                        ReferenceClimRowCols = c(38, 115),
                                        ReferenceClimExtent  = c(180, 308.8, 39.83307, 60.09163))
                                      
                                      # Aqua temp
                                      # ReferenceClimRowCols = c(115, 184),
                                      # ReferenceClimExtent  = c(267.804, 289.5896, 41.34913, 49.01273))
          # library(rasf)
          # to360(-70.4104)

                              
                                      }
      
      # Settings for GatherType: "JoinModelRasters"
      # =======================================================================
        if (identical(GatherType, "JoinModelRasters")){
          GatherModelRasters_Settings = 
            list(InputFolder = "./results/Climwin_reordered",
                 TypesToJoin = c("AICc", 
                                 "R-squares", 
                                 "Slopes",
                                 "_Open.tif",
                                 "_Close.tif",
                                 "tiffs"),
                 OutputFolder = "./results/Climwin_reordered/1_JoinResults", 
                 ClimVarList = c("air.2m",
                                 "prate.sfc",
                                 "WindModels_925p_Night_ComingFromHelgo",
                                 "WindModels_925p_Night_GoingToHelgo"), 
                 AllModels = list(
                     air2m.Models                   = c("Temp_Europe", "Temp_Africa"),         # You need to adjust this according to the names  
                     pratesfc.Models                = c("Prec_Europe", "Prec_Africa"),         # you gave to the models in the spatclimwin runs
                     Wind925pNight.Models_CF        = c("WindCF_Europe", "WindCF_Africa"),
                     Wind925pNight.Models_GT        = c("WindGT_Europe", "WindGT_Africa"))
                 )
        }    
        
      # Settings for GatherType: "AICOverview"
      # put all the variables in one table
      # =======================================================================
        if (identical(GatherType, "AICOverview")){
          GatherAIC_Settings = list(Path_Main        = "./results/Climwin_reordered/1_JoinResults",
                                    Path_Output      = "./results/Climwin_reordered/2_AICOverview", 
                                    DirsToInclude    = c("AICc", 
                                                         "R-squares", 
                                                         "Slopes", 
                                                         "WOpen",
                                                         "WClose"),
                                    RefDate           = "2017-09-15",
                                    ShapeWithLocation = "./data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp", 
                                    ShapeLayerColumn  = c("ne_10m_admin_0_countries", 
                                                          "SOVEREIGNT"), 
                                    Do_CI             = FALSE,
                                    CI_SubsetToIgnore = c(""))
        }
        
      # Settings for GatherType: "SelectCandidatePixels"
      # select the "valley" compared to the nieghbors (local minimum)
      # =======================================================================
        if (identical(GatherType, "SelectCandidatePixels")){
          GatherSCP_Settings = list(Path_Mosaics      = "./results/Climwin_reordered/1_JoinResults",
                                    Path_Output       = "./results/3_CandidatePixels", 
                                    DirsToInclude     = c("AICc"), 
                                    Path_AICcOverview = "./results/Climwin_reordered/2_AICOverview"
                                    )
        }
    }

  # Do_AddCandidatesToPhenDatabase
  # ===========================================================================
    if (Do_AddCandidatesToPhenDatabase){
      InputPath_AddCand      = "./data"
      InputFile_AddCand      = file.path(InputPath_AddCand, "PhenologyDataset.rds")
      
      OutputFile_AddCand     = file.path(InputPath_AddCand, 
                                         "PhenologyDataset_WithCandidateInfluences.rds")
      
      DatabaseClim_Path      = "./data/2_ClimateDatabase_2000to2020_Temp_AvDay.rds"
      
      SelectionFile          = file.path("./results/5_ShortlistWeatherSignals",
                                         "AllSpecies_air.2m_WestEurope.csv")
      
      RefDate                = c(1,7)      # c(day, month), e.g. June 1 -> c(1,6)
      NewVarNamePrefix       = "Eur_Temp"    
      
    }

  # Do_GetChangeOverTime
  # ===========================================================================
    if (Do_GetChangeOverTime){
      VarsToGetChangeOverTimeFor = c("WestEur_Spat37_air.2m_pixel306_Open100_Close48_mean", "WestAfr_Spat43_WindDirection_pixel107_Open142_Close22_NrDaysWithWindComingFromTarget", "WestEur_Spat37_air.2m_pixel333_Open216_Close201_mean", "WestAfr_Spat38_air.2m_pixel204_Open141_Close23_mean", "WestAfr_Spat44_WindDirection_pixel205_Open76_Close57_NrDaysWithWindGoingToTarget", "WestEur_Spat42_WindDirection_pixel128_Open146_Close16_NrDaysWithWindGoingToTarget", "WestEur_Spat39_prate.sfc_pixel59_Open109_Close81_sum", "WestEur_Spat42_WindDirection_pixel155_Open125_Close100_NrDaysWithWindGoingToTarget", "WestAfr_Spat38_air.2m_pixel284_Open192_Close178_mean", "WestEur_Spat41_WindDirection_pixel109_Open36_Close12_NrDaysWithWindComingFromTarget", "WestEur_Spat42_WindDirection_pixel128_Open95_Close48_NrDaysWithWindGoingToTarget", "WestAfr_Spat38_air.2m_pixel154_Open100_Close48_mean", "WestEur_Spat37_air.2m_pixel327_Open101_Close56_mean", "WestAfr_Spat40_prate.sfc_pixel38_Open147_Close109_sum", "WestAfr_Spat43_WindDirection_pixel148_Open131_Close83_NrDaysWithWindComingFromTarget", "WestEur_Spat37_air.2m_pixel307_Open100_Close45_mean", "WestAfr_Spat38_air.2m_pixel393_Open60_Close31_mean", "WestAfr_Spat44_WindDirection_pixel221_Open34_Close15_NrDaysWithWindGoingToTarget", "WestEur_Spat42_WindDirection_pixel128_Open121_Close39_NrDaysWithWindGoingToTarget", "WestAfr_Spat40_prate.sfc_pixel38_Open113_Close81_sum", "WestAfr_Spat38_air.2m_pixel374_Open82_Close25_mean", "WestAfr_Spat43_WindDirection_pixel181_Open140_Close10_NrDaysWithWindComingFromTarget", "WestEur_Spat42_WindDirection_pixel127_Open97_Close25_NrDaysWithWindGoingToTarget", "WestEur_Spat37_air.2m_pixel222_Open65_Close40_mean", "WestEur_Spat39_prate.sfc_pixel77_Open247_Close117_sum", "WestAfr_Spat38_air.2m_pixel374_Open135_Close28_mean")
      AddFormula  = FALSE
      PlotWidth   = 8
      PlotHeight  = 6
      PlotUnits   = "cm"
      PlotDPI     = 300
      PlotYTitle  = element_blank()   # Set to NULL to use the name of the variable
      PlotXTitle  = element_blank()   # Set to NULL to use "Year"
      PlotYBreaks = list(FALSE,      # Set first element to TRUE if you wish to
                                    # replace the breakpoints in the plot
                         seq(70, 120, 10),
                         c("Mar 1", "Mar 11", "Mar 21", "Mar 31", "Apr 10", "Apr 20"))
      # PlotYBreaks = list(TRUE,      # Set first element to TRUE if you wish to
      #                               # replace the breakpoints in the plot
      #                    seq(280, 340, 10),
      #                    c("Sep 27", "Oct 7", "Oct 17", "Oct 27", "Nov 6", "Nov 16", "Nov 26"))
      OutputPlotName = NULL     # Set to NULL to use the name of the variable
      # OutputPlotName = "filename.tiff"     # Set to NULL to use the name of the variable
    }

    
# =============================================================================
# ============================================================================= 
# ============================= END OF INPUT ==================================
# =============================================================================
# =============================================================================
  
# =============================================================================
# =============================================================================
# START PROCESSING - DO NOT CHANGE AS OF HERE
# ============================================================================= 

# =============================================================================
# Do_AnalyseTrendLinearity
# =============================================================================
  if (Do_AnalyseTrendLinearity){
    source("./functions/BH_AD_Do_TimeTrendLinearity.R")
    BH_AD_Do_TimeTrendLinearity(InputDatabasePath_LT, OutputDir_LT, TrendSettings)
  }

# =============================================================================   
# Do_Climwin_Spatial
# =============================================================================
  if (Do_Climwin_Spatial){
    source("./functions/BH_AD_Do_Climwin_Spatial.R")
    
    BH_AD_Do_Climwin_Spatial(InputDatabasePath      = InputDatabasePath, 
                             OutputFolderPath       = OutputFolderPath, 
                             ResponseVariables      = ResponseVariables, 
                             PredictorVariables     = PredictorVariables,
                             PerSpeciesOrID         = PerSpeciesOrID, 
                             SpeciesOrID_VarName    = SpeciesOrID_VarName,
                             DateVariable           = DateVariable, 
                             DataSubset_Species     = DataSubset_Species, 
                             DataSubset_Year        = DataSubset_Year, 
                             Database_ClimSpat_Path = Database_ClimSpat_Path, 
                             ClimPredictorVar       = ClimPredictorVar, 
                             CW_cinterval           = CW_cinterval, 
                             CW_range               = CW_range, 
                             CW_type                = CW_type, 
                             CW_stat                = CW_stat, 
                             CW_func                = CW_func, 
                             CW_refday              = CW_refday, 
                             CW_exclude             = CW_exclude, 
                             ModelName              = ModelName, 
                             MaskPath               = MaskPath, 
                             Ref_LatLon             = Ref_LatLon,
                             Use_CW_BaseModel       = Use_CW_BaseModel, 
                             CW_BaseModel           = CW_BaseModel, 
                             NrRandRepeats          = 5)
  }
          
   # =============================================================================   
# Do_GatherResults_Spatial
# =============================================================================
  if (Do_GatherResults_Spatial){
    source("./functions/BH_AD_Do_GatherResults_Spatial.R")
    source("./functions/BH_ProgressReport_DoubleForLoop.R")
    
    if (identical(GatherType, "ClimWin")){
      BH_AD_Do_GatherResults_Spatial(GatherType, 
                                     ResponseVariables, PredictorVariables, 
                                     GatherClimWin_Settings)
    } else if (identical(GatherType, "JoinModelRasters")){
      BH_AD_Do_GatherResults_Spatial(GatherType, 
                                     ResponseVariables, PredictorVariables, 
                                     GatherModelRasters_Settings)
    } else if (identical(GatherType, "AICOverview")){
      BH_AD_Do_GatherResults_Spatial(GatherType, 
                                     ResponseVariables, PredictorVariables, 
                                     GatherAIC_Settings)
    } else if (identical(GatherType, "SelectCandidatePixels")){
      BH_AD_Do_GatherResults_Spatial(GatherType, 
                                     ResponseVariables, PredictorVariables, 
                                     GatherSCP_Settings)
    }
  }    
  
# Do_AddCandidatesToPhenDatabase
# ===========================================================================
  if (Do_AddCandidatesToPhenDatabase){
    # Read Database_AD
    # =======================================================================
      Database_AD = readRDS(InputFile_AddCand)
    
    # Read File with selected pixels
    # =======================================================================
      SelectedPixels = read.csv(SelectionFile)
      
    # Do for each of the pixels in the selected pixels file
    # =======================================================================
      for (i in 1:nrow(SelectedPixels)){
        # Set Variables
        # ===================================================================
          ClimVar           = as.character(SelectedPixels[i, "ClimateVar"]) 
          PixelNr           = SelectedPixels[i, "PixelNr"]
          TimeOpenClose     = c(SelectedPixels[i, "WOpen"], SelectedPixels[i, "WClose"])
          species           = as.character(SelectedPixels[i, "Species"])
          OperationOverTime = as.character(SelectedPixels[i, "OperationOverTime"])
          
        # Read Spatial Climate Database
        # ===================================================================
          Database_Clim = readRDS(DatabaseClim_Path)
           
        # Subset Spatial Climate Database to desired climate variable
        # ===================================================================
          Database_Clim = Database_Clim[,,ClimVar]
          
        # Subset Spatial Climate Database to the desired pixel
        # ===================================================================
          Database_Clim = Database_Clim[PixelNr,]
          
        # Get the dates of the climate database
        # ===================================================================
          AllDates = names(Database_Clim)
          AllDates = as.Date(AllDates, format = "%Y_%m_%d_XX")
          
        # Declare variable in Database_AD to hold the new data 
        # (if it does not exist yet)
        # ===================================================================
          NewColumnName = paste0(NewVarNamePrefix, "_", 
                                 ClimVar, 
                                 "_pixel", PixelNr, 
                                 "_Open", TimeOpenClose[1], 
                                 "_Close", TimeOpenClose[2], 
                                 "_", OperationOverTime) 
          if (!(NewColumnName %in% names(Database_AD))){
            Database_AD[NewColumnName] = NA_real_
          }
        
        # Get the years for the current species
        # ===================================================================
          CurrentYears = Database_AD$Year[Database_AD$Species == species]
        
        # Get the indices of the RefDate for all years of this species
        # ===================================================================
          RefDateIndices = match(as.Date(paste0(CurrentYears, 
                                                "-", RefDate[2], 
                                                "-", RefDate[1])), AllDates) 
          
        # Get the indices for the whole time period
        # ===================================================================
          IndicesToGet = matrix(data = NA_real_, 
                                nrow = length(RefDateIndices), 
                                ncol = length(seq(RefDateIndices[1] - TimeOpenClose[1], 
                                                  RefDateIndices[1] - TimeOpenClose[2])))
          for (j in 1:length(RefDateIndices)){
            IndicesToGet[j, ] = seq(RefDateIndices[j] - TimeOpenClose[1], 
                                    RefDateIndices[j] - TimeOpenClose[2])
          }
          IndicesToGet = as.vector(t(IndicesToGet))
          
        # Get the data for the given TimePeriod 
        # ===================================================================
          ClimDataPeriod = matrix(data = Database_Clim[IndicesToGet], 
                                  nrow = length(RefDateIndices), 
                                  ncol = length(seq(RefDateIndices[1] - TimeOpenClose[1], 
                                                    RefDateIndices[1] - TimeOpenClose[2])), 
                                  byrow = T)
          
        # Adjust OperationOverTime, in case of NrDaysWithWindComingFromTarget
        # ===================================================================
          if (identical(OperationOverTime, "NrDaysWithWindComingFromTarget")){
            TargetCoords    = as.numeric(unlist(strsplit(as.character(SelectedPixels[i, "TargetCoords"]), ", *")))
            SourceCoords    = as.numeric(unlist(strsplit(as.character(SelectedPixels[i, "SourceCoords"]), ", *")))
            TargetDir       = BH_CalculateHeading(SourceCoords, TargetCoords)
            OperationOverTime = function(x, env = .GlobalEnv){
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
          }
         
        # Adjust OperationOverTime, in case of NrDaysWithWindGoingToTarget
        # ===================================================================
          if (identical(OperationOverTime, "NrDaysWithWindGoingToTarget")){
            TargetCoords    = as.numeric(unlist(strsplit(as.character(SelectedPixels[i, "TargetCoords"]), ", *")))
            SourceCoords    = as.numeric(unlist(strsplit(as.character(SelectedPixels[i, "SourceCoords"]), ", *")))
            TargetDir       = BH_CalculateHeading(TargetCoords, SourceCoords)
            OperationOverTime = function(x, env = .GlobalEnv){
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
          }
          
        # Perform the desired operation over the TimePeriod for each year
        # ===================================================================
          ValuesToAdd = apply(ClimDataPeriod, 1, OperationOverTime)
          
        # Add the values to the AD Database
        # ===================================================================
          Database_AD[Database_AD$Species == species, NewColumnName] = ValuesToAdd    
      }
      
    # Save Database_AD
    # =======================================================================
      saveRDS(Database_AD, file = OutputFile_AddCand)
      write.csv(Database_AD, 
                file = gsub("\\.rds", "\\.csv", OutputFile_AddCand))
  }

# Do_GetChangeOverTime
# ===========================================================================
  if (Do_GetChangeOverTime){
    # Read input data
    # =======================================================================
      InputData = readRDS(InputDatabasePath)
      
    # Create output subfolder
    # =======================================================================
      cOutputFolder = file.path(OutputFolderPath, "ChangeOverTime")
      dir.create(cOutputFolder, showWarnings = F, recursive = T)
      
    # Do for each of the variables
    # =======================================================================
      for (cVar in VarsToGetChangeOverTimeFor){
        # Make a plot 
        # ===================================================================
          cPlot = ggplot(data = InputData, 
                         mapping = aes(x = Year, y = InputData[, cVar])) + 
                  geom_point() + 
                  stat_smooth(method = "lm", colour = "#333333") +
                  theme_bw()
          
        # Set the Y title
        # ===================================================================
          if (is.null(PlotYTitle)){
            cPlot = cPlot + ylab(cVar)
          } else {
            cPlot = cPlot + ylab(PlotYTitle)
          }
          
        # Set the X title, if requested
        # ===================================================================
          if (!is.null(PlotXTitle)){
            cPlot = cPlot + xlab(PlotXTitle)
          }
          
        # Add the formula, if requested
        # ===================================================================
          if (AddFormula){
            cPlot = cPlot + 
                    stat_poly_eq(formula = y ~ x, 
                                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                                 parse = TRUE)
          }
          
        # Rescale the Y axis, if requested
        # ===================================================================
          if (PlotYBreaks[[1]]){
            cPlot = cPlot + 
                    scale_y_continuous(breaks = PlotYBreaks[[2]],
                                       labels = PlotYBreaks[[3]])
          }
          
        # Save the plot
        # ===================================================================
          if (is.null(OutputPlotName)){
            cOutputFile = file.path(cOutputFolder, paste0(cVar, ".tiff"))
          } else {
            cOutputFile = file.path(cOutputFolder, OutputPlotName)
          }
          ggsave(filename = cOutputFile, 
                 plot     = cPlot,
                 width    = PlotWidth,
                 height   = PlotHeight,
                 units    = PlotUnits,
                 dpi      = PlotDPI)
          
        # Print summary of lm to file
        # ===================================================================
          cSumm = summary(lm(InputData[, cVar] ~ InputData[, "Year"]))
          cOutputFile = gsub("\\.tiff", "\\.txt", cOutputFile)
          sink(cOutputFile)
          print(cSumm)
          sink()
      }
  }
        
# =============================================================================
# =============================================================================
# THE END
# =============================================================================
# =============================================================================
