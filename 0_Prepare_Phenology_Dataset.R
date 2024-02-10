library(lubridate)
setwd("~/Aeiral Insectivore Research/aerial-insectivore-nsf/Projects/Chp2_WhatWhereWhen/Code_changed_2020_HaestB")

region.df = read.csv("~/Aeiral Insectivore Research/aerial-insectivore-nsf/Projects/Chp1_Bird Phenology Calculation/Data_Output/phenology_regional_weighted.csv")
region.df$Date = as.Date(region.df$weighted_yday, origin = paste0(region.df$year, "-01-01"))
region.df = subset(region.df, select=-X)

station.df = read.csv("./data/Chap1_pheno_data/station_level_pheno.csv")
station.df$Date = as.Date(station.df$value, origin = paste0(station.df$year, "-01-01"))
station.df = subset(station.df, select=-X)

Database_WS = read.csv(file = "./data/WinterSolsticeTimes.csv")
YearList = seq(2000, 2020, 1)

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
region.df$MADJDAYSWS = NA_integer_
for (year in YearList){
  CurrentWS  = Database_WS$FullDate[Database_WS$YYYY == (year-1)]
  RowsToFill = region.df$year == year
  region.df$MADJDAYSWS[RowsToFill] = region.df$Date[RowsToFill] - CurrentWS + 1  
}

names(region.df)[2] = "Species"
names(region.df)[1] = "Year"
region.df = subset(region.df, select=c(Year, Species, Station, MADJDAYSWS, total_migrants))
region.df <- region.df[, c(2,1,3)]
# ========================================================================= 
# Save base database to output
# =========================================================================
OutputPath_Main             = "./data"
OutputFile_Base             = file.path(OutputPath_Main, 
                                        "PhenologyDataset.rds")
saveRDS(region.df, file = OutputFile_Base)
write.csv(region.df, 
          file = gsub("\\.rds", "\\.csv", OutputFile_Base), 
          row.names = F)

saveRDS(region.df, "./data/PhenologyDataset_StationLevel.rds")
write.csv(region.df, "./data/PhenologyDataset_StationLevel.csv")
