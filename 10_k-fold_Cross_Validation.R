# we want to do a 21-fold cross-validation with slidingwin function to recalculate the R2 value
ReqPackages = c("lme4", "arm", "meifly", "mgcv", "MuMIn", "binhf",
               "ggplot2", "ggrepel", "ggpmisc", "R.utils", "rgdal", 
               "tidyr", "plyr", "dplyr", "BHGeometryTools", "relaimpo",
               "climwin", "lubridate", "raster")
invisible(lapply(ReqPackages, require, character.only = T, quietly = T))

RemSignalsPath     = "./results/5_FinalSignals/top5.csv"
MainOut             = "./results/8_Cross_Validation"
PhenDataPath       = file.path("./data/PhenologyDataset.rds")
climDatababasePath = "./data/ClimateDatabase"
dir.create(MainOut, showWarnings = F, recursive = T)

# Set Output Paths
# ============================================================================= 
#OutputFolderPath_Climwin = file.path(MainOut, "Climwin")
OutputFolderPath_Climwin = file.path(MainOut, "FinalPixels_21fold")

cData = readRDS(PhenDataPath) %>% 
  filter(Species=="yday_50") %>% 
  mutate(Date = as.Date(MADJDAYSWS, origin = paste0(Year,"-01-01")))

randwinSettings = list(exclude         = c(14, 0),
                       refday          = c(15, 9),
                       Ref_LatLon         = c(44.9072, -84.7197),
                       cinterval       = "day",
                       range           = c(258, 0),
                       type            = "absolute",
                       func            = "lin",
                       baselineFormula = paste0("MADJDAYSWS ~ Year"),
                       #nrRepeats       = 5,
                       nrRepeats       = 1)

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
assign("NrDaysWithWindComingFromTarget", NrDaysWithWindComingFromTarget, .GlobalEnv)

# Read the remaining pixels file
# =====================================================================
RemSignals   = read.csv(RemSignalsPath)
RemSignals[5, "ClimateVar"] = "WindDirection"
allRemModels = unique(RemSignals$Model)

for (cModel in allRemModels){
  # Determine the climate database, and adjust some randwinSettings 
  # according to the climate
  # =================================================================
  if (cModel %in% c("cprat_yday_50", "prate_yday_50", "shum.2m_yday_50")){
    cClimDatabase = file.path(climDatababasePath, 
                              "3_ClimateDatabase_2000to2020_Prec_DaySum.rds")
    randwinSettings$stat = "sum"
    randwinSettings$range = c(258,0)
    MaskPath = "./data/Masks/3_ClimateDatabase_2000to2020_Prec_DaySum_shum.2m_OceanMask.tif"
  } else if (cModel %in% c("tmin2m_yday_50", "air2m_yday_50", "tmax2m_yday_50")){
    cClimDatabase = file.path(climDatababasePath, 
                              "2_ClimateDatabase_2000to2020_Temp_AvDay.rds")
    randwinSettings$stat = "mean"
    randwinSettings$range = c(258,0)
    MaskPath = "./data/Masks/2_ClimateDatabase_2000to2020_Temp_AvDay_tmin.2m_OceanMask.tif"
  } else if (cModel %in% c("WindCF850_yday_50")){
    cClimDatabase = file.path(climDatababasePath, 
                              "4_ClimateDatabase_2000to2020_Wind_850p_Daytime.rds")
    randwinSettings$stat = "NrDaysWithWindComingFromTarget"
    randwinSettings$range = c(107,0)
    MaskPath = "./data/Masks/4_ClimateDatabase_2000to2020_Wind_850p_Daytime_WindDirection_OceanMask.tif"
  }
  
  
  # Read the climate database
  # =================================================================
  cClimateData = readRDS(cClimDatabase)
  
  # Subset to the desired variable in the climate database
  # =================================================================
  cClimVar = RemSignals$ClimateVar[RemSignals$Model %in% cModel][1]
  Database_Clim_Mat = cClimateData[,,which(dimnames(cClimateData)[[3]] == cClimVar)]
  cClimateData = cClimateData[, , cClimVar]
  
  # Read the mask raster
  # =========================================================================
  Mask_Raster = raster(MaskPath)
  Mask_Matrix = matrix(Mask_Raster, nrow = nrow(Mask_Raster), byrow = T)
  
  # Create matrices that hold the coordinates for each pixel
  # =========================================================================
  LatCoords = coordinates(Mask_Raster)[,2]
  LonCoords = coordinates(Mask_Raster)[,1]
  LatCoords_Mat = matrix(LatCoords, nrow = nrow(Mask_Raster), byrow = T)
  LonCoords_Mat = matrix(LonCoords, nrow = nrow(Mask_Raster), byrow = T)
  
  
  # Get the dates for the climate data
  # =================================================================
  climateDates    = colnames(cClimateData)
  climateDates    = gsub("_XX", "", climateDates)
  climateDates    = as.Date(climateDates, format = "%Y_%m_%d")
  climateDates    = as.character(format(climateDates, "%d/%m/%Y"))

  
  # Add randwinSettings to global Environment
  # =================================================================
  assign("randwinSettings", randwinSettings, envir = .GlobalEnv)
  
  # Create a variable to hold all models
  # =============================================================
  ClimWinResults_All = vector(mode = "list", length = length(Database_Clim_Mat))
  ClimWinRand_All    = vector(mode = "list", length = length(Database_Clim_Mat))
  
  ProgressOutputPath          = file.path(OutputFolderPath_Climwin, 
                                          cClimVar, cModel)
  OutputFolderPath_Climwin_FM = file.path(ProgressOutputPath, "FullModels") 
  OutputFolderPath_Climwin_RM = file.path(ProgressOutputPath, "RandomModels")
  dir.create(OutputFolderPath_Climwin_FM, showWarnings = F, recursive = T)
  dir.create(OutputFolderPath_Climwin_RM, showWarnings = F, recursive = T)
  
  # Create temporary file for writing progress
  # =============================================================
  tmpFile = file(paste0(ProgressOutputPath, 
                        "Progress_", cClimVar, ".txt"))
  
  # Do for each pixel in the climate image files
  # =============================================================
  pixels = unique((RemSignals %>% filter(Model == cModel))$PixelNr)
  
  #for (pixelNr in 1:nrow(Database_Clim_Mat)){    
  for (pixelNr in pixels){
    
    print(paste0("Pixel #", pixelNr, " out of ", nrow(Database_Clim_Mat)))
    
    # Check if the current pixel falls in the ocean, and if so: skip
    # =========================================================
    MaskValue = Mask_Matrix[pixelNr]
    if (is.na(MaskValue)){
      ClimWinResults_All[[pixelNr]] = NA
      ClimWinRand_All[[pixelNr]]    = NA
      next
    }
    
    # Get the values of the current pixel
    # =========================================================
    CurrentPixel = Database_Clim_Mat[pixelNr, ]
    
    # Get the location of the current pixel
    # =========================================================
    CP_LatLon = c(LatCoords_Mat[pixelNr], LonCoords_Mat[pixelNr])
    Ref_LatLon = c(44.9072, -84.7197)
    
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
    
    # Check that there are no NAs
    # =========================================================
    if (!any(is.na(CurrentPixel))){
      # Do the climwin analysis for this pixel
      # =====================================================
      ClimWinResults = climwin::slidingwin(
        xvar  = list(Climate = CurrentPixel),
        cdate = climateDates,
        bdate = cData$Date,
        baseline  = lm(as.formula(randwinSettings$baselineFormula),
                       data = cData),
        cinterval = randwinSettings$cinterval,
        range     = randwinSettings$range,
        type      = randwinSettings$type,
        refday    = randwinSettings$refday,
        stat      = randwinSettings$stat,
        func      = randwinSettings$func,
        cmissing  = FALSE,
        exclude   = randwinSettings$exclude,
        k=21)
      
      ClimwinRands = climwin::randwin(repeats   = 1,
                                      xvar      = list(Climate = CurrentPixel),
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
                                      exclude   = randwinSettings$exclude,
                                      k=21)
      
      
      # Strip non-essential data from the model to reduce output filesize
      # =====================================================
      attr(ClimWinResults[[1]]$BestModel$terms, ".Environment") <- c()
      
      # Store the results for this pixel
      # =====================================================
      ClimWinResults_All[[pixelNr]] = ClimWinResults
      ClimWinRand_All[[pixelNr]]    = ClimwinRands
      
    }
    
    # Print Progress to temp file
    # =========================================================
    Progress = pixelNr/nrow(Database_Clim_Mat)*100
    cat(Progress, file = tmpFile)
  }
  
  # Save the results for this PerVariable to a file
  # =============================================================
  saveRDS(ClimWinResults_All,
          file = file.path(OutputFolderPath_Climwin_FM, 
                           paste0(cClimVar, ".rds")))
  saveRDS(ClimWinRand_All,
          file = file.path(OutputFolderPath_Climwin_RM, 
                           paste0(cClimVar, ".rds")))
}

i = 294
climwin::plotall(dataset = ClimWinResults_All[[i]][[1]]$Dataset,
                 datasetrand = ClimWinRand_All[[i]][[1]],
                 bestmodel = ClimWinResults_All[[i]][[1]]$BestModel, 
                 bestmodeldata = ClimWinResults_All[[i]][[1]]$BestModelData)

climwin::pvalue(dataset = ClimWinResults_All[[i]][[1]]$Dataset,
       datasetrand = ClimWinRand_All[[i]][[1]], metric = "AIC", sample.size= 21)

CurrentBestSumm = summary(ClimWinResults_All[[i]][[1]]$BestModel)
CurrentBestSumm$adj.r.squared


# =============================================================================
# ============================================================================= 

Do_GatherResults_Spatial       = 1  

if (Do_GatherResults_Spatial){
  
  library(BHiO)
  
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
    InputPathMain = file.path("./results/8_Cross_Validation/Climwin/WindDirection/WindCF850_yday_50")
    GatherClimWin_Settings = list(Do_MakePlots      = 1,
                                  Do_PrintBestModel = 1,
                                  Do_GetpValues     = 1,
                                  Do_PrintCombos    = 1,
                                  ClimWinInputPath  = c(file.path(InputPathMain, 
                                                                  "FullModels"),
                                                        file.path(InputPathMain, 
                                                                  "RandomModels")),
                                  ClimwinOutputPath = file.path("./results/8_Cross_Validation/Climwin/WindDirection/"),
                                  pValueThreshold   = 0.3,
                                  #SpeciesSubset     = list(1, "yday_50"),
                                  SpeciesSubset     = 0,
                                  ReferenceType     = "Climate",
                                  
                                  # Please uncomment the relevant location-variable combo
                                  
                                  # temp
                                  # ReferenceClimRowCols = c(13, 51),
                                  # ReferenceClimExtent  = c(215.625-1.875/2, 309.375+1.875/2, 39.047-1.9047/2,  61.9033+1.9047/2))  # MinLong, MaxLong, MinLat, MaxLat ### bbox
                                  
                                  # Wind
                                  ReferenceClimRowCols = c(11, 39),
                                  ReferenceClimExtent  = c(215 - 2.5/2, 310+2.5/2, 37.5-2.5/2, 62.5+2.5/2))
                                  
                                  # precipitation
                                  # ReferenceClimRowCols = c(13, 51),
                                  # ReferenceClimExtent  = c(215.625-1.875/2, 309.375+1.875/2, 39.047-1.9047/2,  61.9033+1.9047/2))
                                  
                                  # # EVI
                                  # ReferenceClimRowCols = c(38, 115),
                                  # ReferenceClimExtent  = c(180, 308.8, 39.83307, 60.09163))
    
                                  # Aqua temp
                                  # ReferenceClimRowCols = c(115, 184),
                                  # ReferenceClimExtent  = c(267.804, 289.5896, 41.34913, 49.01273))
    # library(rasf)
    # to360(-70.4104)
    
    
  }}

ResponseVariables = c("MADJDAYSWS")        
PredictorVariables = c("Year")

# Do_GatherResults_Spatial
# =============================================================================
if (Do_GatherResults_Spatial){
  source("./functions/BH_AD_Do_GatherResults_Spatial.R")
  source("./functions/BH_ProgressReport_DoubleForLoop.R")
  
  if (identical(GatherType, "ClimWin")){
    BH_AD_Do_GatherResults_Spatial(GatherType, 
                                   ResponseVariables, PredictorVariables, 
                                   GatherClimWin_Settings)
  } 
}    

AllCombos = read.csv(file.path(ClimwinOutputPath, "Combos.csv")) %>% 
  dplyr::select(-X)

# Gather R2 results from climwin
combo = data.frame(model=NA, adj.r.squared = NA, PcVal=NA, variable=NA, PixelNr=NA)

for (i in 1:length(ClimWinResults_All)){
  if (!(is.null(ClimWinResults_All[[i]]) || is.na(ClimWinResults_All[[i]]))){
    
    CurrentBestSumm = summary(ClimWinResults_All[[i]][[1]]$BestModel)
    combo[i,2] = CurrentBestSumm$adj.r.squared
    combo[i,1] = paste0(cModel)
    combo[i,3] = climwin::pvalue(datasetrand = ClimWinRand_All[[i]][[1]], 
                                   dataset = ClimWinResults_All[[i]][[1]]$Dataset, 
                                   metric = "C", 
                                   sample.size = length(ClimWinResults_All[[i]][[1]]$BestModelData$yvar))
    combo[i,4] = strsplit(cModel, "_")[[1]][1]
    combo[i,5] = i
  }
}
combo = na.omit(combo)

combo_final = left_join(combo, AllCombos, by="PixelNr")
write.csv(combo_final, file.path(ClimwinOutputPath, "Combos_final.csv"))


