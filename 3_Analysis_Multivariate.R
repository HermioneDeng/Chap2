# =============================================================================
# Author:           Birgen Haest, birgen.haest@protonmail.com
# Last changed:     November 26, 2019
# =============================================================================
# Functionality:
#   Script to perform multivariate analysis on phenology - climate database
# =============================================================================
# =============================================================================
# Package requirements and files to source
  ReqPackages = c("lme4", "arm", "meifly", "mgcv", "MuMIn", "binhf",
                  "ggplot2", "ggrepel", "ggpmisc")
  invisible(lapply(ReqPackages, library, character.only = T, quietly = T))
  source("/functions/BH_OtherFunctions.R")
  
# =============================================================================
# Global Input Settings
# =============================================================================
  InputDatabasePath_AD   = "./data/PhenologyDataset_WithCandidateInfluences.rds"

  OutputFolder            = "./results"
  
  RespVar                 = c("MADJDAYSWS")    

  PredictorVariables = c("WestEur_Spat37_air.2m_pixel306_Open100_Close48_mean", 
                         "WestAfr_Spat43_WindDirection_pixel107_Open142_Close22_NrDaysWithWindComingFromTarget",
                         "WestEur_Spat37_air.2m_pixel333_Open216_Close201_mean", 
                         "WestAfr_Spat38_air.2m_pixel204_Open141_Close23_mean",
                         "WestAfr_Spat44_WindDirection_pixel205_Open76_Close57_NrDaysWithWindGoingToTarget")
  
# =============================================================================
# Option Settings - Choose what you want to do
# =============================================================================
  Do_DeltaAICcIntOnly     = 1  
  Do_AnalyseCollinearity  = 0  
  
  Do_Boruta               = 0  
  Do_BestModel            = 0  
  Do_VarImpAll            = 0  
 
# =============================================================================
# Option-dependent Input Settings
# =============================================================================
  # Do_DeltaAICcIntOnly
  # ===========================================================================
    if (Do_DeltaAICcIntOnly){
      OutputFolderPath_DIO = file.path(OutputFolder, "DeltaAICcIntOnly")
    }
  
  # Do_AnalyseCollinearity
  # =========================================================================== 
    if (Do_AnalyseCollinearity){
      OutputFolderPath_CSV = file.path(OutputFolder, "CorrCSVTables")
    }

  # Do_Boruta
  # ===========================================================================
    if(Do_Boruta){
      Boruta_Settings = list(SubsetVarCond_One       = list(TRUE, "Species", "European Pied Flycatcher"),
                             SubsetVarCond_Two       = list(FALSE, "Year", seq(1980, 2014)),
                             PerSpecies              = TRUE,
                             Do_Standardise_PredVars = TRUE, 
                             MaxRuns                 = 1000,
                             doTraceSetting          = 1)
    }
  
  # Do_BestModel
  # ===========================================================================
    if (Do_BestModel){
      SubSetSpecies = list("Species", "European Pied Flycatcher")
      
      Do_CheckForRemainingTrend = 1
    }
  
  # Do_VarImpAll
  # ===========================================================================
    if (Do_VarImpAll){
      VarImpAll_Settings  = list(InputDatabasePath_AD     = InputDatabasePath_AD,
                                 OutputFolder             = OutputFolder,
                                 ResponseVariables        = RespVar,
                                 PredictorVariables       = PredictorVariables,
                                 SubsetVarCond_One        = list(TRUE, "Species", "European Pied Flycatcher"),
                                 SubsetVarCond_Two        = list(FALSE, "Year", tmp),
                                 PerSpecies               = FALSE, 
                                 Do_Standardise_PredVars  = TRUE,    # Do not do this when analysing multiple species and having PerSpecies = FALSE
                                 Do_AIC_VarImp            = TRUE,
                                 AIC_VarImp_Settings      = list(Do_CompareAllPerm       = TRUE,
                                                                 VarsToInvert            = "",
                                                                 VarsToLog10             = "",
                                                                 Interactions            = "", 
                                                                 Do_Standardise_RespVar  = FALSE,
                                                                 Do_MakePositive         = FALSE, 
                                                                 LimitModels             = list(1, "NrVariables", 4),
                                                                 MakeAddedVariablePlots  = FALSE,
                                                                 MakePredictionPlot      = TRUE,
                                                                 MakeDiagnosticsPlot     = FALSE,
                                                                 LeaveOneOutVarImp       = FALSE),  
                                 Do_LinRelaImpo           = TRUE,
                                 RelaImpo_VarImp_Settings = list(RI_type     = c("lmg"),
                                                                 BootstrapCI = FALSE), 
                                 Do_Boruta                = TRUE,
                                 Boruta_Settings          = list(MaxRuns           = 1000,
                                                                 doTraceSetting    = 1),
                                 Gather_VarImp            = TRUE)
    }
  
# =============================================================================
# =============================================================================
# START PROCESSING - DO NOT CHANGE AS OF HERE
# =============================================================================
  # Do_DeltaAICcIntOnly
  # ===========================================================================
    if (Do_DeltaAICcIntOnly){
      # Get the phenology data, and order by year
      # =======================================================================
        cData = readRDS(InputDatabasePath_AD)
        cData = cData[with(cData, order(Year)), ]
        
      # Rename response variable in data to "RespVar"
      # =======================================================================
        colnames(cData)[colnames(cData) %in% RespVar] <- "RespVar"
        
      # Get AICc for intercept only model
      # =======================================================================
        BaseModel = lm("RespVar ~ 1", data = cData)
        BaseAICc  = AICc(BaseModel)
        
      # Create dataframe to hold all aiccs
      # =======================================================================
        AllAICc = data.frame(Var             = c("Intercept Only", PredictorVariables), 
                             AICc            = NA_real_, 
                             DiffWithIntOnly = NA_real_)
        AllAICc$AICc[1] = BaseAICc
        
      # Get all AICc 
      # =======================================================================
        for (cVar in PredictorVariables){
          tmpData = cData
          colnames(tmpData)[colnames(tmpData) %in% cVar] <- "cVar"
          cModel = lm("RespVar ~ cVar", data = tmpData)
          AllAICc$AICc[AllAICc$Var %in% cVar] = AICc(cModel)
        }
        
      # Calculate Delta AICcs
      # =======================================================================
        AllAICc$DiffWithIntOnly = AllAICc$AICc - AllAICc$AICc[1]
        
      # Save to file
      # =======================================================================
        dir.create(OutputFolderPath_DIO, showWarnings = F, recursive = T)
        OutputFile = file.path(OutputFolderPath_DIO, "DeltaAICcWithIntOnly.rds")
        saveRDS(AllAICc, file = OutputFile)
        write.csv(AllAICc, file = gsub("\\.rds", "\\.csv", OutputFile))
    }
  
  # Do_AnalyseCollinearity
  # =========================================================================== 
    if (Do_AnalyseCollinearity){
      # Get the phenology data, and order by year
      # =======================================================================
        cData = readRDS(InputDatabasePath_AD)
        cData = cData[with(cData, order(Year)), ]
  
      # Set Output file name, and create outputfolder if necessary
      # =======================================================================
        OutputFileName = file.path(OutputFolderPath_CSV, 
                                   "CollinearityAnalysis.csv")
        dir.create(OutputFolderPath_CSV, showWarnings = F, recursive = T)

      # Subset Selected variables, and convert all variables to numeric, 
      # and sort columns by order  
      # =======================================================================
        VariablesToPlot = c(PredictorVariables, RespVar)
        cDataNumeric = cData[, VariablesToPlot]
        cDataNumeric = as.data.frame(apply(cDataNumeric, 2, "as.numeric"))

      # Declare variables to hold all correlation coefficients and 
      # p-values
      # =======================================================================
        CurrentCorrCoef = array(data = NA_real_, 
                                dim = c(ncol(cDataNumeric), 
                                        ncol(cDataNumeric), 
                                        3))
         dimnames(CurrentCorrCoef)[[3]] <- c("pearson", 
                                             "spearman", 
                                             "kendall")
        CurrentpValues  = CurrentCorrCoef
        
      # Calculate correlation coefficients
      # =======================================================================
        # Pearson
          CurrentCorrCoef[,,1] = cor(cDataNumeric, 
                                     method = "pearson", 
                                     use = "complete.obs")
        # Spearman
          CurrentCorrCoef[,,2] = cor(cDataNumeric, 
                                     method = "spearman", 
                                     use = "complete.obs")
        # Kendall
          CurrentCorrCoef[,,3] = cor(cDataNumeric, 
                                     method = "kendall", 
                                     use = "complete.obs")

      # Calculate p-values
      # =======================================================================
        # Pearson
        # =====================================================================
          tmppValues            = try(cor.mtest(cDataNumeric, 
                                                0.95, 
                                                method = "pearson"), 
                                      silent = TRUE)
          if (class(tmppValues) != "try-error"){
            CurrentpValues[,,1] = tmppValues[[1]]  
          }
          
        # Spearman
        # =====================================================================
          tmppValues            = try(cor.mtest(cDataNumeric, 
                                                0.95, 
                                                method = "spearman"), 
                                      silent = TRUE)
          if (class(tmppValues) != "try-error"){
            CurrentpValues[,,2] = tmppValues[[1]]  
          }
          
        # Kendall
        # =====================================================================
          tmppValues            = try(cor.mtest(cDataNumeric, 
                                                0.95, 
                                                method = "kendall"), 
                                      silent = TRUE)
          if (class(tmppValues) != "try-error"){
            CurrentpValues[,,3] = tmppValues[[1]]  
          }
      
      # Set column and row names of Corr coef and p-values arrays
      # =======================================================================
        colnames(CurrentCorrCoef) <- c(PredictorVariables, RespVar)
        rownames(CurrentCorrCoef) <- colnames(CurrentCorrCoef)  
        colnames(CurrentpValues)  <- colnames(CurrentCorrCoef)
        rownames(CurrentpValues)  <- rownames(CurrentCorrCoef)

      # Write all correlation coefficients and p-values to csv and RData file
      # =======================================================================
        BH_ADHelpers_WriteCorrelations(OutputFileName,
                                       CurrentCorrCoef,
                                       CurrentpValues)
        save(CurrentCorrCoef, CurrentpValues, 
             file = gsub("\\.csv", "\\.RData", OutputFileName))
    }
  
  # ===========================================================================
  # Do_Boruta
  # ===========================================================================
    if (Do_Boruta){
      # Update Boruta_Settings
      # =======================================================================
        Boruta_Settings$InputDatabasePath  = InputDatabasePath_AD
        Boruta_Settings$OutputFolder       = OutputFolder
        Boruta_Settings$ResponseVariables  = RespVar
        Boruta_Settings$PredictorVariables = PredictorVariables
      
      source("./functions/BH_AD_Do_Boruta.R")
      BH_AD_Do_Boruta(Boruta_Settings)
    }
  
  # ===========================================================================
  # Do_BestModel
  # ===========================================================================
    if (Do_BestModel){
      # Read database_AD
      # =======================================================================
        Database_AD = readRDS(InputDatabasePath_AD)
      
      # Set Species
      # =======================================================================
        species = SubSetSpecies[[2]]
        
      # Subset to species
      # =======================================================================
        DataToAnalyse = Database_AD[Database_AD[[SubSetSpecies[[1]]]] %in% species, ]
        
      # Check if there is only 1 response variable
      # =======================================================================
        if (length(RespVar) != 1){
          stop("Your variable ResponseVariables should contain only 1 variable 
               for the BestModel Analysis..")
        } 
        
      # Set Outputfolder
      # =======================================================================
        OutputFolder_BM = file.path(OutputFolder, "BestModelAnalysis")
        dir.create(OutputFolder_BM, showWarnings = F, recursive = T)
        
      # Case: Check for remaining trend
      # =======================================================================
        if (Do_CheckForRemainingTrend){
          # Subset the data to the data to analyse
          # ===================================================================  
            CurrentData = DataToAnalyse[c("Species", 
                                          RespVar, PredictorVariables, 
                                          "Year")]
            
          # Remove lines (i.e. years) with NA values
          # ===================================================================
            CurrentData = CurrentData[complete.cases(CurrentData), ]
          
          # Get MAD_LY
          # ===================================================================
            RespVar_LY = paste0(RespVar, "_LY")
            CurrentData_MADLY = CurrentData
            CurrentData_MADLY[RespVar_LY] = shift(CurrentData[, RespVar], 1)
            CurrentData_MADLY[1, RespVar_LY] = NA
            CurrentData_MADLY = CurrentData_MADLY[complete.cases(CurrentData_MADLY), ]
            
          # Create the full formulas for all models
          # ===================================================================
            Formula_Without        = paste(PredictorVariables, collapse = "+")
            Formula_Without        = paste0(RespVar, "~", Formula_Without)
            Formula_WithYear       = paste(c(PredictorVariables, "Year"), collapse = "+")
            Formula_WithYear       = paste0(RespVar, "~", Formula_WithYear)
            Formula_Without_MADLY  = paste(c(PredictorVariables, RespVar_LY), collapse = "+")
            Formula_Without_MADLY  = paste0(RespVar, "~", Formula_Without_MADLY)
            Formula_With_MADLY     = paste(c(PredictorVariables, "Year", RespVar_LY), collapse = "+")
            Formula_With_MADLY     = paste0(RespVar, "~", Formula_With_MADLY)
             
          # Get AIC of model without and with Year, and with MAD_LastYear
          # ===================================================================
            AICValues           = rep(NA_real_, 2)
            AICValues_withMADLY = rep(NA_real_, 4)
            AICValues[1] = AICc(lm(Formula_Without, data = CurrentData))
            AICValues[2] = AICc(lm(Formula_WithYear, data = CurrentData))
            AICValues_withMADLY[1] = AICc(lm(Formula_Without, data = CurrentData_MADLY))
            AICValues_withMADLY[2] = AICc(lm(Formula_WithYear, data = CurrentData_MADLY))
            AICValues_withMADLY[3] = AICc(lm(Formula_Without_MADLY, data = CurrentData_MADLY))
            AICValues_withMADLY[4] = AICc(lm(Formula_With_MADLY, data = CurrentData_MADLY))
            AICValues = data.frame(AICc_WithoutYear = AICValues[1], 
                                   AICc_WithYear = AICValues[2])
            AICValues_withMADLY = data.frame(AICc_WithoutYear_withoutMADLY = AICValues_withMADLY[1], 
                                             AICc_WithYear_withoutMADLY = AICValues_withMADLY[2],
                                             AICc_WithoutYear_withMADLY = AICValues_withMADLY[3], 
                                             AICc_WithYear_withMADLY = AICValues_withMADLY[4])
            AICValues_withMADLY[2,] = rank(AICValues_withMADLY)
            write.csv(AICValues, 
                      file.path(OutputFolder_BM, paste0(species, "_AICcValues.csv")))
            write.csv(AICValues_withMADLY, 
                      file.path(OutputFolder_BM, paste0(species, "_AICcValues_withMADLY.csv")))
            
          # Write model summaries to file
          # ===================================================================
            CurrentData_Scaled = CurrentData
            for (CurrVar in PredictorVariables){
              CurrentData_Scaled[CurrVar] = scale(CurrentData_Scaled[CurrVar])
            }
            capture.output(summary(lm(data = CurrentData, 
                                      formula = Formula_Without)), 
                           file = file.path(OutputFolder_BM, paste0(species, "_ModelSummary_WithoutYear_NotStand.txt")))
            capture.output(summary(lm(data = CurrentData, 
                                      formula = Formula_WithYear)), 
                           file = file.path(OutputFolder_BM, paste0(species, "_ModelSummary_WithYear_NotStand.txt")))
            capture.output(summary(lm(data = CurrentData_Scaled, 
                                      formula = Formula_Without)), 
                           file = file.path(OutputFolder_BM, paste0(species, "_ModelSummary_WithoutYear_Stand.txt")))
            capture.output(summary(lm(data = CurrentData_Scaled, 
                                      formula = Formula_WithYear)), 
                           file = file.path(OutputFolder_BM, paste0(species, "_ModelSummary_WithYear_Stand.txt")))
            capture.output(summary(lm(data = CurrentData_MADLY, 
                                      formula = Formula_Without_MADLY)), 
                           file = file.path(OutputFolder_BM, paste0(species, "_ModelSummary_WithoutYear_With", 
                                                                    RespVar_LY, "_NotStand.txt")))
            capture.output(summary(lm(data = CurrentData_MADLY, 
                                      formula = Formula_With_MADLY)), 
                           file = file.path(OutputFolder_BM, paste0(species, "_ModelSummary_WithYear_With", 
                                                                    RespVar_LY, "_NotStand.txt")))
            
          # Make plots of measured and predicted values, over time
          # ===================================================================
            DataToPlot = CurrentData
            names(DataToPlot)[which(names(DataToPlot) == RespVar)] <- "RespVar"
            DataToPlot$PredWithout = lm(Formula_Without, data = CurrentData)$fitted.values
            DataToPlot$PredWithYear = lm(Formula_WithYear, data = CurrentData)$fitted.values
            CurrentPlot = ggplot(data = DataToPlot, 
                                 mapping = aes(y = RespVar, x = Year)) + 
                          geom_point() + 
                          ylab(RespVar) + 
                          ggtitle(paste0(species, ": Measured and predicted values")) + 
                          geom_point(aes(y = PredWithout, colour = "Without Year")) + 
                          geom_line(aes(y = PredWithout, colour = "Without Year")) +
                          geom_point(aes(y = PredWithYear, colour = "With Year")) + 
                          geom_line(aes(y = PredWithYear, colour = "With Year")) + 
                          geom_text_repel(mapping = aes(label = Year))
            BH_saveGGPlotPNG(CurrentPlot, 
                             file.path(OutputFolder_BM, paste0(species, "_MeasuredVersusPredicted_OverTime.png")))
            
          # Make plots of measured versus predicted values
          # ===================================================================
            CurrentPlot = ggplot(data = DataToPlot, 
                                 mapping = aes(y = RespVar)) + 
                          geom_point(aes(x = PredWithout, colour = "Without Year")) + 
                          geom_point(aes(x = PredWithYear, colour = "With Year")) + 
                          stat_smooth(method = "lm", aes(x = PredWithout, colour = "Without Year")) +
                          stat_smooth(method = "lm", aes(x = PredWithYear, colour = "With Year")) +
                          geom_abline(mapping = aes(colour = "45 degrees line"), 
                                      intercept = 0, slope = 1) +
                          ylab(paste0("Measured", RespVar)) + 
                          ggtitle(paste0(species, ": Measured versus predicted values")) 
            BH_saveGGPlotPNG(CurrentPlot, 
                             file.path(OutputFolder_BM, paste0(species, "_MeasuredVersusPredicted.png")))
            
          # Make model evaluation plots
          # ===================================================================
            CurrentPlot = diagPlots_OLSModel(lm(data = CurrentData, formula = Formula_Without))
            BH_saveGGPlotPNG(CurrentPlot,
                             file.path(OutputFolder_BM, paste0(species, "_ModelDiag_Without.png")))
            CurrentPlot = diagPlots_OLSModel(lm(data = CurrentData, formula = Formula_WithYear))
            BH_saveGGPlotPNG(CurrentPlot,
                             file.path(OutputFolder_BM, paste0(species, "_ModelDiag_WithYear.png")))
            
          # Make plot of residuals over time
          # ===================================================================
            Resid_Without = lm(data = CurrentData, formula = Formula_Without)$residuals
            Resid_With    = lm(data = CurrentData, formula = Formula_WithYear)$residuals
            CurrentPlot = ggplot(data = DataToPlot, 
                                 mapping = aes(x = Year)) + 
                          geom_point(aes(y = Resid_Without, 
                                         colour = "Without Year")) + 
                          geom_point(aes(y = Resid_With,
                                         colour = "With Year")) + 
                          stat_smooth(method = "lm", 
                                      aes(y = Resid_Without, colour = "Without Year")) +
                          stat_smooth(method = "lm", 
                                      aes(y = Resid_With, colour = "With Year")) +
                          stat_poly_eq(formula = y ~ x, 
                                       aes(y = Resid_Without, 
                                           label = paste(..eq.label.., ..rr.label.., sep = "~~~"), 
                                           colour = "Without Year"), 
                                       parse = TRUE) + 
                          stat_poly_eq(formula = y ~ x, 
                                       aes(y = Resid_With, 
                                           label = paste(..eq.label.., ..rr.label.., sep = "~~~"), 
                                           colour = "With Year"), 
                                       vjust = 3,
                                       parse = TRUE) + 
                          ylab("Residuals") +
                          ggtitle(paste0(species, ": Residuals over time"))
            BH_saveGGPlotPNG(CurrentPlot, 
                             file.path(OutputFolder_BM, paste0(species, "_ResidualsOverTime.png")))
        }
    }
  
  # ===========================================================================
  # Do_VarImpAll
  # ===========================================================================
    if (Do_VarImpAll){
      source("./functions/BH_AD_Do_VarImpAll.R")
      BH_AD_Do_VarImpAll(VarImpAll_Settings)
    }
  
# =============================================================================
# =============================================================================
# THE END
# =============================================================================
# =============================================================================
