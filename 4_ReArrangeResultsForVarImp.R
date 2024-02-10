# =============================================================================
# Author:           Birgen Haest, birgen.haest@protonmail.com
# Last changed:     February 05, 2020
# =============================================================================
# Functionality:
#   Script to reorder the results form 2_DoSpatClimwin to be used in the 
#   multivariate step 
# =============================================================================
# ============================================================================= 
# Input Settings
# ============================================================================= 
  inputDir  = "./results/2_SpaWin"
  outputDir = "./results/3_SpaWin_Adj"

# Declare variables
# =============================================================================  
  newOrderDirs = c("BestModelSumm", "Combos", "Plots", "pValues")
  # newOrderDirs = c("BestModelSumm", "Plots", "pValues")
  
# Get list of species (directories)
# ============================================================================= 
  speciesDirsFull  = list.dirs(inputDir, full.names = T, recursive = F)
  speciesDirsShort = list.dirs(inputDir, full.names = F, recursive = F)
  
# Create outputfolder 
# =============================================================================
  dir.create(outputDir, showWarnings = F)
  
# Do for each species
# =============================================================================
  for (cSpecies in speciesDirsShort){
    # Get list of spatial grids 
    # ========================================================================= 
      spatGridDirsFull  = list.dirs(speciesDirsFull[speciesDirsShort %in% cSpecies], 
                                    full.names = T, recursive = F)
      spatGridDirsShort = list.dirs(speciesDirsFull[speciesDirsShort %in% cSpecies], 
                                    full.names = F, recursive = F)
  
    # Do for each spatial grid that was analysed
    # =========================================================================
      for (cLoc in spatGridDirsShort){
        # Set list of randomisations
        # =====================================================================
          randDirsFull  = list.dirs(spatGridDirsFull[spatGridDirsShort %in% cLoc], 
                                    full.names = T, recursive = F)
          randDirsShort = list.dirs(spatGridDirsFull[spatGridDirsShort %in% cLoc], 
                                    full.names = F, recursive = F)

        # Do for each of the seasons
        # =====================================================================    
          for (cSeason in randDirsShort){
            # Create subfolder in output
            # =================================================================
              cSeasonFull = randDirsFull[randDirsShort %in% cSeason]
              cSeasonDir = file.path(outputDir, cSeason)
              dir.create(cSeasonDir, showWarnings = F, recursive = T)
    
            # Do for each of the new order directories
            # =================================================================
              for (cNewOrder in newOrderDirs){
                # Create the output subfolder
                # =============================================================
                  cNewOrderDir = file.path(cSeasonDir, cNewOrder)
                  dir.create(cNewOrderDir, showWarnings = F)
          
                # Get list of weather variable (models)
                # =============================================================
                  wVarDirsShort = list.dirs(cSeasonFull, 
                                            full.names = F, recursive = F)
                
                # Do for each weather variable
                # =============================================================
                  for (cVar in wVarDirsShort){
                    # Create weather variable subdirectory
                    # =========================================================
                      varOut = cVar
                      cOutDir = file.path(cNewOrderDir, varOut)
                      dir.create(cOutDir, showWarnings = F)
                    
                    # Copy relevant files
                    # =========================================================
                      baseFromDir = file.path(inputDir, cSpecies, cLoc, cSeason, cVar)
                      if (cNewOrder %in% "BestModelSumm"){
                        file.copy(from      = file.path(baseFromDir, "AICc"), 
                                  to        = cOutDir,
                                  recursive = T)
                        file.copy(from      = file.path(baseFromDir, "R-squares"), 
                                  to        = cOutDir,
                                  recursive = T)
                        file.copy(from      = file.path(baseFromDir, "Slopes"), 
                                  to        = cOutDir,
                                  recursive = T)
                      } else if (cNewOrder %in% "Combos"){
                        if (file.exists(file.path(baseFromDir, "OpenCloseHistograms"))){
                          file.copy(from      = file.path(baseFromDir, "OpenCloseHistograms"),
                                    to        = cOutDir,
                                    recursive = T)
                        }
                        if (file.exists(file.path(baseFromDir, "OpenCloseRasters"))){
                          file.copy(from      = file.path(baseFromDir, "OpenCloseRasters"),
                                    to        = cOutDir,
                                    recursive = T)
                        }
                        file.copy(from      = file.path(baseFromDir, "Combos.csv"),
                                  to        = file.path(cOutDir, paste0(cSpecies, "_Combos.csv")))
                      } else if (cNewOrder %in% "Plots"){
                        if (file.exists(file.path(baseFromDir, "Climate"))){
                          tmpSpecies = basename(list.dirs(file.path(baseFromDir, "Climate"))[2])
                          file.copy(from      = file.path(baseFromDir, "Climate", tmpSpecies), 
                                    to        = cOutDir,
                                    recursive = T)
                        }
                      } else if (cNewOrder %in% "pValues"){
                        pOutDir = file.path(cOutDir, "tables")
                        dir.create(pOutDir, showWarnings = F)
                        cFile = list.files(baseFromDir, "*_pValues.csv", full.names = F)
                        file.copy(from      = file.path(baseFromDir, cFile), 
                                  to        = file.path(pOutDir, cFile))
                        pOutDir = file.path(cOutDir, "tiffs")
                        dir.create(pOutDir, showWarnings = F)
                        cFile = list.files(baseFromDir, "*_pValues.tif", full.names = F)
                        file.copy(from      = file.path(baseFromDir, cFile), 
                                  to        = file.path(pOutDir, cFile))
                      }
                  }
              }
          }
          
      }
  }
