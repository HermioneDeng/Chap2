# =============================================================================
# Author:           Birgen Haest, birgen.haest@protonmail.com
# Created:          June 02, 2020
# Last changed:     June 02, 2020
# =============================================================================
# Functionality:
#   Script to create figures for the marginal effects on population sizes 
#   of each of the selected weather signals
# =============================================================================
# =============================================================================
# Package requirements ----
# =============================================================================
  library(sjPlot)
  library(ggplot2)

# =============================================================================
# INPUT SETTING ----
# =============================================================================
  InputData_SelPixels  = "./data/SelectedWeatherSignals.csv"
  cSeason = "MADJDAYSWS"
  OutputFolderPath     = "./results/5_FinalSignals/5_MarginalEffectsPlots"
  
  if (cSeason %in% "MADJDAYSWS"){
    InputDatabasePath_AD = "./results/4_MultiVar/4_PhenDataWithCandidates/3_PhenologyDataset_MADJDAYSWS_Updated.rds"
  } else if (cSeason %in% "am_ws"){
    InputDatabasePath_AD = "./results/4_MultiVar/4_PhenDataWithCandidates/5_BrackenCavePhenology_am_ws_Updated.rds"
  }
  
  PlotWidth   = 8
  PlotHeight  = 6
  PlotUnits   = "cm"
  PlotDPI     = 300
  PlotYTitle  = element_blank()   # Set to NULL to use the name of the variable
  PlotXTitle  = element_blank()   # Set to NULL to use "Year"
  

# =============================================================================
# =============================================================================  
# START PROCESSING ----
# =============================================================================
# Read input, and create outputfolder
# =============================================================================
  Database_AD = readRDS(InputDatabasePath_AD)
  Database_SP = read.csv(InputData_SelPixels)
  dir.create(OutputFolderPath, showWarnings = F, recursive = T)

# Subset selected pixels to best only
# =============================================================================
  #Database_SP = Database_SP[Database_SP$MultiScenario %in% "StrongestOnly", ]

# Adjust settings
# =============================================================================
  # if (cSeason %in% "sm_ws"){NrYears = 23} else if (cSeason %in% "am_ws"){NrYears = 22}
  NrYears = 21
  if (cSeason %in% "MADJDAYSWS"){
    PlotYBreaks = list(TRUE,      # Set first element to TRUE if you wish to
                                # replace the breakpoints in the plot
                     seq(230, 250, 4),
                     c("Aug 9", "Aug 13", "Aug 17", "Aug 21", "Aug 25", "Aug 29"))
  } else if (cSeason %in% "am_ws"){
    # PlotYBreaks = list(TRUE,      # Set first element to TRUE if you wish to
    #                             # replace the breakpoints in the plot
    #                  seq(0, 50000, 10000),
    #                  seq(0, 50000, 10000))
    PlotYBreaks = list(TRUE,      # Set first element to TRUE if you wish to
                                  # replace the breakpoints in the plot
                       seq(280, 340, 10),
                       c("Sep 27", "Oct 7", "Oct 17", "Oct 27", "Nov 6", "Nov 16", "Nov 26"))
  }
  
# Subset Database_Sp to current season
# =========================================================================
  D#atabase_SP = Database_SP[Database_SP$Season %in% cSeason, ]
  
# Remove first year with NA value for am_ws from Database_AD
# ===========================================================================
  # if (cSeason %in% "am_ws"){Database_AD = Database_AD[-1, ]}

# Get data to plot 
# ===========================================================================
  cVarsToGet  = as.character(Database_SP$VarName[3:5])
  cDataToPlot = Database_AD[, c(cSeason, cVarsToGet)]
  
# Create full model
# ===========================================================================
  cFormula = paste0(cSeason, " ~ ", paste(cVarsToGet, collapse = " + "))
  cModel = lm(cFormula, data = cDataToPlot)

# Create effect plot for each of the signals
# ===========================================================================
  for (cVar in cVarsToGet){
    # Make base plot 
    # =======================================================================
      cPlot = plot_model(cModel, type = "pred", terms = cVar, colors = "black")
    
    # Add points
    # =======================================================================
      cPlot = cPlot + 
              geom_point(mapping = aes(x = cDataToPlot[,cVar], y = cDataToPlot[, cSeason]), 
                         data = cDataToPlot)
      
    # Set them to bw
    # =======================================================================
      cPlot = cPlot + theme_classic()
      
    # Rescale the Y axis, if requested
    # =======================================================================
      if (PlotYBreaks[[1]]){
        cPlot = cPlot + 
                scale_y_continuous(breaks = PlotYBreaks[[2]],
                                   labels = PlotYBreaks[[3]])
      }
      
    # Add or remove axis titles
    # =======================================================================
      cPlot = cPlot + xlab(PlotXTitle) + 
                      ylab(PlotYTitle) 
      
    # Remove title 
    # =======================================================================
      cPlot = cPlot + ggtitle(label = element_blank())
      
    # Save the plot
    # =======================================================================
      cOutputFile = file.path(OutputFolderPath, paste0("dPheno.dClimate_", cVar, ".pdf"))
      ggsave(filename = cOutputFile, 
             plot     = cPlot,
             width    = PlotWidth,
             height   = PlotHeight,
             units    = PlotUnits,
             dpi      = PlotDPI)
  }
  
# Print summary of spring lm to file
# ===========================================================================
  cSumm = summary(cModel)
  cOutputFile = file.path(OutputFolderPath, paste0(cSeason, "_Model_Summary.txt"))
  sink(cOutputFile)
  print(cSumm)
  sink()
  
# Get data to plot for trend contribution plots
# ===========================================================================
  cDataToPlot = Database_AD[, c(cSeason, cVarsToGet, "Year")]
  
# Create trend contribution plots for each of the signals
# ===========================================================================
  for (cVar in cVarsToGet){
    # Make base plot 
    # =======================================================================
      cPlot = ggplot(data    = cDataToPlot,
                     mapping = aes(x = Year, 
                                   y = cDataToPlot[, cSeason])) + geom_point()
    
    # Add trend contribution line and confidence interval
    # =======================================================================
      cSlope     = Database_SP$dPheno.dClimate[Database_SP$VarName %in% cVar]
      cIntercept = mean(cDataToPlot[, cSeason] - (cSlope*cDataToPlot$Year))
      y          = cDataToPlot[, cSeason]
      y.fit      = cIntercept+cSlope*cDataToPlot$Year
      n          = NrYears
      x          = cDataToPlot$Year
      se <- sqrt(sum((y - y.fit)^2) / (n - 2)) * sqrt(1 / n + (x - mean(x))^2 / sum((x - mean(x))^2))
      
      cPlot = cPlot + 
              geom_line(mapping   = aes(x = x, 
                                        y = y.fit),
                        linetype  = 2, 
                        size      = rel(1.2)) + 
              geom_ribbon(aes(ymin  = y.fit-se,
                              ymax  = y.fit+se), 
                          alpha    = 0.5)
    
    # Add overall trend
    # =======================================================================
      cPlot = cPlot + 
              stat_smooth(method = "lm", 
                          colour = "#333333")
        
    # Set them to bw
    # =======================================================================
      cPlot = cPlot + theme_classic()
      
    # Rescale the Y axis, if requested
    # =======================================================================
      if (PlotYBreaks[[1]]){
        cPlot = cPlot + 
                scale_y_continuous(breaks = PlotYBreaks[[2]],
                                   labels = PlotYBreaks[[3]],
                                   limits = c(230, 250))
      }
      
    # Add or remove axis titles
    # =======================================================================
      cPlot = cPlot + xlab(PlotXTitle) + 
                      ylab(PlotYTitle) 
      
    # Remove title 
    # =======================================================================
      cPlot = cPlot + ggtitle(label = element_blank())
      
    # Save the plot
    # =======================================================================
      cOutputFile = file.path(OutputFolderPath, paste0("TrendContr_", cVar, ".pdf"))
      ggsave(filename = cOutputFile, 
             plot     = cPlot,
             width    = PlotWidth,
             height   = PlotHeight,
             units    = PlotUnits,
             dpi      = PlotDPI)
  }

  #dClimate/dTime
  # ===========================================================================
  
  df = read.csv("./results/4_MultiVar/4_PhenDataWithCandidates/3_PhenologyDataset_MADJDAYSWS_Updated.csv")
  df2 = read.csv("./data/SelectedWeatherSignals.csv")
  
  # Find matching columns
  matching_columns <- df %>% 
    dplyr::select(Year, MADJDAYSWS, all_of(df2$VarName))
  write.csv(matching_columns, "./results/5_FinalSignals/pheno_top5_df.csv")
  
  pheno_top5_df = read.csv("./results/5_FinalSignals/pheno_top5_df.csv") %>% 
    gather(key="variable", value="value", 4:8) 
  
  for (cVar in unique(pheno_top5_df$variable)){
    print(cVar)
    print(summary(lm(value ~ Year, data = pheno_top5_df %>% filter(variable==cVar))))
    
    # Make base plot 
    # =======================================================================
    cPlot = ggplot(pheno_top5_df %>% filter(variable==cVar), aes(x=Year, y=value))+
      geom_point()+
      geom_smooth(method="lm", formula = y~x, color="black")
    
    # Set them to bw
    # =======================================================================
    cPlot = cPlot + theme_classic()
    
    # Add or remove axis titles
    # =======================================================================
    cPlot = cPlot + xlab(PlotXTitle) + 
      ylab(PlotYTitle) 
    
    # Remove title 
    # =======================================================================
    cPlot = cPlot + ggtitle(label = element_blank())
    
    # Save the plot
    # =======================================================================
    cOutputFile = file.path(OutputFolderPath, paste0("dClimate.dTime_", cVar, ".pdf"))
    ggsave(filename = cOutputFile, 
           plot     = cPlot,
           width    = PlotWidth,
           height   = PlotHeight,
           units    = PlotUnits,
           dpi      = PlotDPI)
  }
  