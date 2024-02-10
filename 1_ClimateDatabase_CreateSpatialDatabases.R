# ============================================================================= 
# Author:           Birgen Haest, birgen.haest@protonmail.com
# Last changed:     January  20, 2020
# ============================================================================= 
# Functionality:    
# Script to create spatial weather databases from NCEP data (using RNCEP) 
# ============================================================================= 
# ============================================================================= 
# Package requirements and clear memory
  ReqPackages = c("lubridate", "abind", "raster", "gdalUtils", "NISTunits")
  invisible(lapply(ReqPackages, library, character.only = T, quietly = T))
  
# Source necessary function to get the climate data
  source("./functions/BH_GetClimateVariable_Spatial.R")

# ============================================================================= 
# ============================================================================= 
# Choose what you want to do
# ============================================================================= 
  Do_CreateBaseDatabase   = TRUE
  Do_GetTemperatureData   = TRUE
  Do_GetPrecipitationData = TRUE
  Do_GetWindData          = TRUE
  
  Do_ConvertToRaster      = TRUE
  Do_Create_Masks         = TRUE
 
# ============================================================================= 
# Input Settings
# ============================================================================= 
  RegionToGet                 = "Europe" # Set to "Europe or "WestAfrica
  InputPath_Main              = "./data"
  InputPath_WSRefs            = "./data/WinterSolsticeTimes.csv"
    
  OutputPath_Main             = "./data"
  YearList                    = seq(2000,2020)
  OutputFile_Base             = file.path(OutputPath_Main, 
                                          "1_ClimateDatabase_Base_Spatial_2000to2020.rds")
  OutputFile_Temp_AvDay       = file.path(OutputPath_Main,
                                          "2_ClimateDatabase_2000to2020_Temp_AvDay.rds")
  OutputFile_Prec_DaySum      = file.path(OutputPath_Main,
                                          "3_ClimateDatabase_2000to2020_Prec_DaySum.rds")
  OutputFile_Wind_Daytime    = file.path(OutputPath_Main,
                                          "4_ClimateDatabase_2000to2020_Wind_850p_Daytime.rds")
 
  OutputFolder_Rasters     = file.path(OutputPath_Main, "Rasters")
  OutputFolder_Masks       = file.path(OutputPath_Main, "Masks")
  
  # General spatial settings - Uncomment Europe or West Africa
  # ===========================================================================
    StatusBar  = TRUE        # Whether or not to show progress of data download  
    #Great Lakes North
    LatMinMax  = c(39.7728386,  60.0926778)  # Set coordinates of spatial area
    LonMinMax  = c(-143.1316366, -52.0093988)  #  of interest.
    
  # Specific Settings: Do_ConvertToRaster - Uncomment the one you wish to convert
  # to a raster. Only 1 by 1.
  # ===========================================================================
    # Temperature:
      DatabaseToConvert2Ras = OutputFile_Temp_AvDay
      RasterRowCol = c(13, 51)
      RasterExtent = c(215.625-1.875/2, 309.375+1.875/2, 39.047-1.9047/2,  61.9033+1.9047/2) # MinLong, MaxLong, MinLat, MaxLat
       
    # Precipitation:
      DatabaseToConvert2Ras = OutputFile_Prec_DaySum
      RasterRowCol = c(13, 51)
      RasterExtent = c(215.625-1.875/2, 309.375+1.875/2, 39.047-1.9047/2,  61.9033+1.9047/2)

    # # Wind:
      DatabaseToConvert2Ras = OutputFile_Wind_Daytime
      RasterRowCol = c(11, 39)
      RasterExtent = c(215 - 2.5/2, 310+2.5/2, 37.5-2.5/2, 62.5+2.5/2)
    
    
  # Specific Settings: Do_Create_Masks 
  # ===========================================================================
    InputFolder_Masks = OutputFolder_Rasters
    SubsetFiles = list(TRUE, "6_MOD13A2_EVI_interpolated_rescaled.tif")
    NewMasks      = 1            # Set to 1 if this is not an update of a mask
    MaskType      = "FullMask"  # Set to "OceanMask" for masking Oceans (out - not visible)
                                 # Set to "LandMask" for masking land
                                 # Set to "FullMask" to masking nothing
    
# =============================================================================    
# =============================================================================
# Do_CreateBaseDatabase, if requested
# =============================================================================
  if (Do_CreateBaseDatabase){
    # Read Reference Database
    # ========================================================================= 
      Database_WS = read.csv(file = InputPath_WSRefs)
    
    # ========================================================================= 
    # Create the output dataframe, consisting of the variables "Date" and "Year"
    # ========================================================================= 
      Climdata_Base = data.frame(Date = seq(as.Date(paste0(YearList[1], "-01-01")), 
                                            as.Date(paste0(YearList[length(YearList)], "-12-31")), 
                                            by = 1))
      Climdata_Base$Year = year(Climdata_Base$Date)
        
    # ========================================================================= 
    # Add JDAYSWS (i.e. Julian date counting from WS previous year)
    # ========================================================================= 
      # Delete unneeded years from Database_WS
      # =======================================================================
        Database_WS = Database_WS[which((Database_WS$YYYY >= min(YearList)-1 ) &
                                        (Database_WS$YYYY < max(YearList))), ]
    
      # Add Variable to winter solstice table to hold the full date
      # =======================================================================
        Database_WS$FullDate = as.Date(paste0(as.character(Database_WS$YYYY),
                                              as.character(Database_WS$MMM),
                                              as.character(Database_WS$DD)),
                                              format = "%Y%b%d")
    
      # Add JDAYSWS to Climdata_Base
      # =======================================================================
        Climdata_Base$JDAYSWS = NA_integer_
        for (year in YearList){
          CurrentWS  = Database_WS$FullDate[Database_WS$YYYY == (year-1)]
          RowsToFill = Climdata_Base$Year == year
          Climdata_Base$JDAYSWS[RowsToFill] = Climdata_Base$Date[RowsToFill] - CurrentWS + 1  
        }
    
    # ========================================================================= 
    # Save base database to output
    # =========================================================================
      dir.create(dirname(OutputFile_Base), showWarnings = F, recursive = T)
      saveRDS(Climdata_Base, file = OutputFile_Base)
      write.csv(Climdata_Base, 
                file = gsub("\\.rds", "\\.csv", OutputFile_Base), 
                row.names = F)
  }

# =============================================================================
# Do_GetTemperatureData, if requested
# =============================================================================
  if (Do_GetTemperatureData){
    # Read spatial base database
    # ========================================================================= 
      Climdata_Base = readRDS(OutputFile_Base)
    
    # Create Mean, Max, Min Daily Temperature database 
    # (Mean, max, or min of the 4 measurements made over the day)
    # ========================================================================= 
      Climdata_Temp = BH_GetClimateVariable_Spatial(ClimateDatabase   = Climdata_Base,
                                                    DesiredVariable   = "air.2m",
                                                    LatMinMax         = LatMinMax, 
                                                    LonMinMax         = LonMinMax,
                                                    CreateNewDatabase = 1, 
                                                    status.bar        = StatusBar)
      Climdata_Temp = BH_GetClimateVariable_Spatial(ClimateDatabase   = Climdata_Temp,
                                                    DesiredVariable   = "tmax.2m",
                                                    LatMinMax         = LatMinMax, 
                                                    LonMinMax         = LonMinMax,
                                                    CreateNewDatabase = 0, 
                                                    status.bar        = StatusBar)
      Climdata_Temp = BH_GetClimateVariable_Spatial(ClimateDatabase   = Climdata_Temp,
                                                    DesiredVariable   = "tmin.2m",
                                                    LatMinMax         = LatMinMax, 
                                                    LonMinMax         = LonMinMax,
                                                    CreateNewDatabase = 0, 
                                                    status.bar        = StatusBar)
      if (RegionToGet %in% "Europe"){
        saveRDS(Climdata_Temp, file = OutputFile_Temp_AvDay)
      } else if (RegionToGet %in% "WestAfrica"){
        saveRDS(Climdata_Temp, file = OutputFile_Temp_WA_AvDay)
      }
  }    

# =============================================================================
# Do_GetPrecipitationData, if requested
# =============================================================================
  if (Do_GetPrecipitationData){
    # Read spatial base database
    # ========================================================================= 
      Climdata_Base = readRDS(OutputFile_Base)
    
    # Create Precipitation database 
    # (sum of the 4 measurements made over the day)
    # ========================================================================= 
      Climdata_Prec = BH_GetClimateVariable_Spatial(ClimateDatabase   = Climdata_Base,
                                                    DesiredVariable   = "shum.2m",
                                                    LatMinMax         = LatMinMax, 
                                                    LonMinMax         = LonMinMax,
                                                    CreateNewDatabase = 1,
                                                    DownScale         = NA,
                                                    status.bar        = StatusBar)
      Climdata_Prec = BH_GetClimateVariable_Spatial(ClimateDatabase   = Climdata_Prec,
                                                    DesiredVariable   = "cprat.sfc",
                                                    LatMinMax         = LatMinMax, 
                                                    LonMinMax         = LonMinMax,
                                                    CreateNewDatabase = 0, 
                                                    DownScale         = NA,
                                                    status.bar        = StatusBar)
      Climdata_Prec = BH_GetClimateVariable_Spatial(ClimateDatabase   = Climdata_Prec,
                                                    DesiredVariable   = "prate.sfc",
                                                    LatMinMax         = LatMinMax, 
                                                    LonMinMax         = LonMinMax,
                                                    CreateNewDatabase = 0, 
                                                    DownScale         = NA,
                                                    status.bar        = StatusBar)
      if (RegionToGet %in% "Europe"){
        saveRDS(Climdata_Prec, file = OutputFile_Prec_DaySum)
      } else if (RegionToGet %in% "WestAfrica"){
        saveRDS(Climdata_Prec, file = OutputFile_Prec_WA_DaySum)
      }
  }

# =============================================================================    
# Do_GetWindData, if requested
# =============================================================================
  if (Do_GetWindData){
    # Set variables 
    # =========================================================================  
      WindVariables = c("uwnd.1000p", "vwnd.1000p","uwnd.925p", "vwnd.925p","uwnd.850p", "vwnd.850p")
    
    # Read spatial base database
    # ========================================================================= 
      Climdata_Base = readRDS(OutputFile_Base)
      
    # Get uwind component (daytime value: 12, 18 UTC)
    # =======================================================================
      Climdata_Wind = BH_GetClimateVariable_Spatial(ClimateDatabase    = Climdata_Base, 
                                                    DesiredVariable    = WindVariables[1],
                                                    LatMinMax          = LatMinMax, 
                                                    LonMinMax          = LonMinMax,
                                                    CreateNewDatabase  = 1,
                                                    status.bar         = StatusBar,
                                                    AdditionalSettings = list(SpecificTiming = TRUE, 
                                                                              Timing         = "Daytime"))
      Climdata_Wind = BH_GetClimateVariable_Spatial(ClimateDatabase    = Climdata_Base, 
                                                    DesiredVariable    = WindVariables[3],
                                                    LatMinMax          = LatMinMax, 
                                                    LonMinMax          = LonMinMax,
                                                    CreateNewDatabase  = 1,
                                                    status.bar         = StatusBar,
                                                    AdditionalSettings = list(SpecificTiming = TRUE, 
                                                                              Timing         = "Daytime"))
      Climdata_Wind = BH_GetClimateVariable_Spatial(ClimateDatabase    = Climdata_Base, 
                                                    DesiredVariable    = WindVariables[5],
                                                    LatMinMax          = LatMinMax, 
                                                    LonMinMax          = LonMinMax,
                                                    CreateNewDatabase  = 1,
                                                    status.bar         = StatusBar,
                                                    AdditionalSettings = list(SpecificTiming = TRUE, 
                                                                              Timing         = "Daytime"))
      
    # Get vwind component (daytime value: 12, 18 UTC)
    # =======================================================================
      Climdata_Wind = BH_GetClimateVariable_Spatial(ClimateDatabase    = Climdata_Wind, 
                                                    DesiredVariable    = WindVariables[2],
                                                    LatMinMax          = LatMinMax, 
                                                    LonMinMax          = LonMinMax,
                                                    CreateNewDatabase  = 0,
                                                    status.bar         = StatusBar,
                                                    AdditionalSettings = list(SpecificTiming = TRUE, 
                                                                              Timing         = "Daytime"))
      Climdata_Wind = BH_GetClimateVariable_Spatial(ClimateDatabase    = Climdata_Wind, 
                                                    DesiredVariable    = WindVariables[4],
                                                    LatMinMax          = LatMinMax, 
                                                    LonMinMax          = LonMinMax,
                                                    CreateNewDatabase  = 0,
                                                    status.bar         = StatusBar,
                                                    AdditionalSettings = list(SpecificTiming = TRUE, 
                                                                              Timing         = "Daytime"))
      Climdata_Wind = BH_GetClimateVariable_Spatial(ClimateDatabase    = Climdata_Wind, 
                                                    DesiredVariable    = WindVariables[6],
                                                    LatMinMax          = LatMinMax, 
                                                    LonMinMax          = LonMinMax,
                                                    CreateNewDatabase  = 0,
                                                    status.bar         = StatusBar,
                                                    AdditionalSettings = list(SpecificTiming = TRUE, 
                                                                              Timing         = "Daytime"))
    }
    
    # Calculate Wind direction
    # =======================================================================
    Database_WD = array(NA_real_, dim = c(dim(Climdata_Wind)[1:2], 1))
    for (i in 1:nrow(Database_WD)){
      for (j in 1:ncol(Database_WD)){
        Database_WD[i,j, 1] = (atan2((-1*Climdata_Wind[i,j,1]), (-1*Climdata_Wind[i,j,2]))*180/pi)
      }
    }
    NegAngles = which(Database_WD < 0)
    Database_WD[NegAngles] =  Database_WD[NegAngles] + 360
      
    # Add the wind direction to the climate database 
    # =======================================================================
      Climdata_Wind = abind(Climdata_Wind, Database_WD, along = 3)
      dimnames(Climdata_Wind)[[3]][3] <- "WindDirection"
      
    # Calculate Wind assistance
    # =======================================================================
      Database_WA = array(NA_real_, dim = c(dim(Climdata_Wind)[1:2], 1))
      for (i in 1:nrow(Database_WA)){
        for (j in 1:ncol(Database_WA)){
          Database_WA[i,j,1] = NCEP.M.Groundspeed(u=Climdata_Wind[i,j,1], v=Climdata_Wind[i,j,2], direction=, p.airspeed=12)
        }
      }
      
    # Save the wind database to a file
    # =======================================================================
      if (RegionToGet %in% "Europe"){
        saveRDS(Climdata_Wind, file = OutputFile_Wind_Daytime)
      } else if (RegionToGet %in% "WestAfrica"){
        saveRDS(Climdata_Wind, file = OutputFile_Wind_Daytime)
      }
  
    
# =============================================================================
# Do_ConvertToRaster, if requested
# =============================================================================
  if (Do_ConvertToRaster){
    # Set Variables
    # =========================================================================
      RasterCrs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
      
    # Read the input climate database
    # =======================================================================
      InputDatabase = readRDS(DatabaseToConvert2Ras)

    # Get the variables in this database
    # ===================================================================
      if (length(dim(InputDatabase)) == 3){
        AllVars = dimnames(InputDatabase)[[3]]
      } else if (length(dim(InputDatabase)) == 2){
        if (grepl("WindSpeed", InputClimDatabaseFile_CTR)){
          AllVars = "WindSpeed"
        }
      }
      
      for (CurrVar in AllVars){
        # Subset the database to the current variable
        # ===============================================================
          if (CurrVar %in% "WindSpeed"){
            Database_CurrVar = InputDatabase
          } else {
            Database_CurrVar = InputDatabase[,,CurrVar]
          }

        # Convert the current database to an array of the desired size
        # ===============================================================
          Database_Arr = array(NA_real_,
                               dim = c(RasterRowCol[1],
                                       RasterRowCol[2],
                                       dim(Database_CurrVar)[2]))
          
          # j represents each day of the data (7671 days in total)
          for (j in 1:(dim(Database_CurrVar)[2])){
            Database_Arr[,,j] = Database_CurrVar[,j]
          }

        # Convert the array to a rasterbrick object
        # ===============================================================
          OutputBrick = brick(Database_Arr,
                              xmn = RasterExtent[1],
                              xmx = RasterExtent[2],
                              ymn = RasterExtent[3],
                              ymx = RasterExtent[4],
                              crs = RasterCrs)

        # Write the raster to a file
        # ===============================================================
          CurrOutputfile = file.path(OutputFolder_Rasters,
                                     paste0(gsub("\\.rds", "", basename(DatabaseToConvert2Ras)),
                                            "_", CurrVar, ".tif"))
          dir.create(dirname(CurrOutputfile), showWarnings = F, recursive = T)
          writeRaster(OutputBrick, CurrOutputfile, overwrite=TRUE)
      }
  }

# =============================================================================
# Do_Create_Masks
# =============================================================================
  if (Do_Create_Masks){
    # Set Variables
    # =========================================================================
      OutputName_Suffix = MaskType
      ShapeMaskPath     = file.path(getwd(), "data", 
                                    "ne_10m_admin_0_countries", 
                                    "ne_10m_admin_0_countries.shp")
      ShapeLayer        = "ne_10m_admin_0_countries"
      
    # Create outputfolder if it doesn`t exist
    # =======================================================================
      dir.create(OutputFolder_Masks, showWarnings = F, recursive = T)
      
    # Get list of files for which to create a mask
    # =======================================================================
      FileList = list.files(InputFolder_Masks, pattern = "*.tif$")
    
    # Subset files, if requested
    # =======================================================================
      if (SubsetFiles[[1]]){
        FileList = FileList[grepl(SubsetFiles[[2]], FileList)]
      }
      
    # Do for each of the raster files in the input folder
    # =======================================================================
      for (CurrentFile in FileList){
        # Set the input climate raster name
        # ===================================================================
          InputRaster_Path = file.path(InputFolder_Masks, CurrentFile)
         
        # Read the input raster
        # ===================================================================
          InputRaster = raster(InputRaster_Path)
          InputRaster = raster::rotate(InputRaster)
        
        # Set the output raster climate name
        # ===================================================================
          OutputRaster_Path = file.path(OutputFolder_Masks, 
                                        gsub("\\.tif", 
                                             paste0("_", OutputName_Suffix, "\\.tif"),  
                                             CurrentFile))
           
        # Case: NEW Masks
        # ===================================================================
          if (NewMasks){
            # Create a raster mask, with all values = NA
            # ===============================================================
              RasterMask = raster(InputRaster)
              RasterMask = setValues(RasterMask, rep(NA, ncell(RasterMask)))
              writeRaster(RasterMask, 
                          filename = file.path(OutputRaster_Path), overwrite=TRUE)
              
            # Set the pixel values in raster mask to 1 for the non-masked pixels
            # ===============================================================
              # OceanMask
                if (identical(MaskType, "OceanMask")){
                  tmp = gdal_rasterize(src_datasource = ShapeMaskPath, 
                                       dst_filename   = OutputRaster_Path,
                                       b              = 1,
                                       at             = TRUE,
                                       burn           = 1,
                                       l              = ShapeLayer,
                                       output_Raster = TRUE)
                }
              # LandMask
                if (identical(MaskType, "LandMask")){
                  tmp = gdal_rasterize(src_datasource = ShapeMaskPath, 
                                       dst_filename   = OutputRaster_Path,
                                       b              = 1,
                                       at             = TRUE,
                                       i              = TRUE,
                                       burn           = 1,
                                       l              = ShapeLayer,
                                       output_Raster = TRUE)
                }
              # FullMask
                if (identical(MaskType, "FullMask")){
                  tmp = raster(OutputRaster_Path)
                  tmp = setValues(tmp, rep(1, ncell(tmp)))
                  writeRaster(tmp, OutputRaster_Path, overwrite = TRUE)
                }
              
        # Case: UPDATE Masks
          } else {
            # Set the pixel values in raster mask to NA for extra masked pixels
            # ===============================================================
              # OceanMask
                if (identical(MaskType, "OceanMask")){
                  tmp = gdal_rasterize(src_datasource = ShapeMaskPath, 
                                       dst_filename   = OutputRaster_Path,
                                       b              = 1,
                                       at             = TRUE,
                                       burn           = NA,
                                       i              = TRUE,
                                       l              = ShapeLayer,
                                       output_Raster = TRUE)
                }
              # LandMask
                if (identical(MaskType, "LandMask")){
                  tmp = gdal_rasterize(src_datasource = ShapeMaskPath, 
                                       dst_filename   = OutputRaster_Path,
                                       b              = 1,
                                       at             = TRUE,
                                       burn           = NA,
                                       l              = ShapeLayer,
                                       output_Raster = TRUE)
                }
              
          }
      }
  }
         
# ============================================================================= 
# ============================================================================= 
# END OF SCRIPT
# ============================================================================= 
# ============================================================================= 
