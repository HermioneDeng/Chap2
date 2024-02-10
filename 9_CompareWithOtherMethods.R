#######################################################################
################## Compared to Priori assumptions #####################
#######################################################################
## 1) local variable + dynamic time  
## 2) local variable + peak roosting timing, all compared with remote variables + dynamic timing 

RemSignalsPath     = "./results/5_FinalSignals/top5.csv"
MainOut             = "./results/7_Compare_with_Assumptions"
PhenDataPath       = "./results/4_MultiVar/4_PhenDataWithCandidates/3_PhenologyDataset_MADJDAYSWS_Updated.rds"
climDatababasePath = "./data/ClimateDatabase"
climwinResultsPath = "./results/2_SpaWin/yday_50/FullArea/MADJDAYSWS"

randwinSettings = list(exclude         = c(14, 0),
                       refday          = c(15, 9),
                       nrRepeats       = 5,
                       cinterval       = "day",
                       range           = c(258, 0),
                       type            = "absolute",
                       func            = "lin",
                       baselineFormula = paste0("MADJDAYSWS ~ Year"))

# Extract env data (all days) from a few pixels around the Great Lakes

dir.create(MainOut, showWarnings = F, recursive = T)

cData = readRDS(PhenDataPath) %>% 
  mutate(Date = as.Date(MADJDAYSWS, origin = paste0(Year,"-01-01")))

# Read the remaining pixels file
# =====================================================================
RemSignals   = read.csv(RemSignalsPath)
allRemModels = unique(RemSignals$Model)
allRemModels = allRemModels[allRemModels != "WindCF850_yday_50"]

# Create a dataframe to store all measurements (static time window July-Sep)
# =============================================================
static_df = data.frame(year = seq(2000,2020))

# Read the data raster
# =========================================================================
MaskPath    = file.path("./data/Rasters/2_ClimateDatabase_2000to2020_Temp_AvDay_tmin.2m.tif") #tmin and precip have the same pixel numbers and extent, so it's okay to use the same
Mask_Raster = rast(MaskPath)
df_Mask_Raster = cbind(as.data.frame(Mask_Raster[[1]]), as.data.frame(crds(Mask_Raster[[1]])))
df_Mask_Raster$x <- df_Mask_Raster$x - 360
df_Mask_Raster <- df_Mask_Raster[, c(2,3,1)]
Mask_Raster = rast(df_Mask_Raster, type="xyz", crs="+proj=longlat +datum=WGS84", digits=6, extent=NULL)
Mask_Matrix = matrix(Mask_Raster, nrow = nrow(Mask_Raster), byrow = T)

# Get the all the local pixels numbers at the Great Lakes
# =========================================================
lake = terra::vect("./data/maps/layers/greatlakes.shp")
plot(Mask_Raster)
plot(lake, add=T)
plot(crop(terra::mask(Mask_Raster, lake), lake))
plot(lake, add=T)
allPixelNrsThisVar = terra::extract(Mask_Raster, lake, cell=TRUE)$cell # get the raster cell numbers

for (cModel in allRemModels){
  # Determine the climate database, and adjust some randwinSettings 
  # according to the climate
  # =================================================================
  if (cModel %in% c("cprat_yday_50", "prate_yday_50", "shum.2m_yday_50")){
    cClimDatabase = file.path(climDatababasePath, 
                              "3_ClimateDatabase_2000to2020_Prec_DaySum.rds")
    randwinSettings$stat = "sum"
  } else if (cModel %in% c("tmin2m_yday_50", "air2m_yday_50", "tmax2m_yday_50")){
    cClimDatabase = file.path(climDatababasePath, 
                              "2_ClimateDatabase_2000to2020_Temp_AvDay.rds")
    randwinSettings$stat = "mean"
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
    climateDates    = colnames(cClimateData)
    climateDates    = gsub("_XX", "", climateDates)
    climateDates    = as.Date(climateDates, format = "%Y_%m_%d")
    climateDates    = as.character(format(climateDates, "%d/%m/%Y"))
  
    # Read the climwin results for this variable
    # =================================================================
    cClimwinResults = readRDS(file.path(climwinResultsPath,
                                        cModel,
                                        "MADJDAYSWS/GreatLakes/FullModels",
                                        paste0("yday_50", ".rds")))
    
    # Add randwinSettings to global Environment
    # =================================================================
    assign("randwinSettings", randwinSettings, envir = .GlobalEnv)
    
    # Create a variable to hold all models
    # =============================================================
    ClimWinResults_All = vector(mode = "list", length = length(allPixelNrsThisVar))
    ClimWinRand_All    = vector(mode = "list", length = length(allPixelNrsThisVar))
    
    for (cPixel in allPixelNrsThisVar){
      
      ClimWinResults_All[[cPixel]] = NA
      ClimWinRand_All[[cPixel]]    = NA
      
      # Subset climate data to current pixel 
      # =============================================================
      climData_cPixel = cClimateData[cPixel, ]
      
      # Extract mean value July 1st - Sep 31st for each year 
      # =============================================================
      climData_df = data.frame(value = climData_cPixel, date = as.Date(climateDates, "%d/%m/%Y")) %>% 
        filter(month(date) >= 7 & month(date) <= 9) %>% 
        mutate(year = year(date)) %>% 
        group_by(year) %>% 
        summarise_at(vars(value), list(value = mean))
      names(climData_df)[2] = paste0(cModel, "_", cPixel)
      
      static_df = left_join(static_df, climData_df)
      
      # Do climwin and 5 random climwin repeats
      # =============================================================
      ClimWinResults = climwin::slidingwin(
                                  xvar  = list(Climate = climData_cPixel),
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
                                  exclude   = randwinSettings$exclude)
      
      ClimwinRands = climwin::randwin(repeats   = 5,
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

      # Strip non-essential data from the model to reduce output filesize
      # =====================================================
      attr(ClimWinResults[[1]]$BestModel$terms, ".Environment") <- c()

      # Store the results for this pixel
      # =====================================================
      ClimWinResults_All[[cPixel]] = ClimWinResults
      ClimWinRand_All[[cPixel]]    = ClimwinRands
    }
    
    # Save the results for this PerVariable to a file
    # =============================================================
    saveRDS(ClimWinResults_All,
            file = file.path(MainOut, 
                             paste0("dynamic_window_ClimWinResults_", cModel, ".rds")))
    saveRDS(ClimWinRand_All,
            file = file.path(MainOut, 
                             paste0("dynamic_window_ClimWinRand_", cModel, ".rds")))
}


write.csv(static_df,
          file.path(MainOut,
                    paste0("static_window.csv")))


# Get the results from climwin (local variable + dydnamic time-window)
# =============================================================
# ClimWinResults_All = ClimWinResults_All[!unlist(lapply(ClimWinResults_All,is.null))]

summ_df_dynamic = data.frame()

for (cModel in allRemModels){

# read in the RDS files
ClimWinResults_All = readRDS(file.path(MainOut, 
                                              paste0("dynamic_window_ClimWinResults_", cModel, ".rds")))
ClimWinRand_All = readRDS(file.path(MainOut, 
                         paste0("dynamic_window_ClimWinRand_", cModel, ".rds")))

# Gather results from climwin
dynamic = data.frame(term=NA, adj.r.squared = NA, PcVal=NA, variable=NA)

for (i in 1:length(ClimWinResults_All)){
  if (!(is.null(ClimWinResults_All[[i]]) || is.na(ClimWinResults_All[[i]]))){
    
  CurrentBestSumm = summary(ClimWinResults_All[[i]][[1]]$BestModel)
  dynamic[i,2] = CurrentBestSumm$adj.r.squared
  dynamic[i,1] = paste0(cModel, "_", i)
  dynamic[i,3] = climwin::pvalue(datasetrand = ClimWinRand_All[[i]][[1]], 
                                    dataset = ClimWinResults_All[[i]][[1]]$Dataset, 
                                    metric = "C", 
                                    sample.size = length(ClimWinResults_All[[i]][[1]]$BestModelData$yvar))
  dynamic[i,4] = strsplit(cModel, "_")[[1]][1]
        }
      }
  dynamic = na.omit(dynamic)
  summ_df_dynamic = rbind(summ_df_dynamic, dynamic)
}

hist(summ_df_dynamic$adj.r.squared)
write.csv(summ_df_dynamic,
          file.path(MainOut,
                    paste0("dynamic_window.csv")))

summ_df_dynamic = read.csv(file.path(MainOut, "dynamic_window.csv"))
  
library(hrbrthemes)
library(ggdark)
plot1 = summ_df_dynamic %>%
  ggplot(aes(x=adj.r.squared)) +
  geom_histogram(color="#e9ecef", alpha=1, position = 'identity', breaks = seq(0, 0.8, 0.05)) +
  scale_y_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30))+
  #scale_fill_viridis(discrete=TRUE) +
  theme_classic()+
  labs(fill="", x="Adjusted R2", y="Count")+
  #dark_theme_classic()+
  theme(text = element_text(size=18))


# If we want to look at the climwin plots
climwin::plotbetas(dataset = ClimWinResults_All[[i]][[1]]$Dataset)

climwin::plotbest(dataset = ClimWinResults_All[[i]][[1]]$Dataset,
         bestmodel = ClimWinResults_All[[i]][[1]]$BestModel, 
         bestmodeldata = ClimWinResults_All[[i]][[1]]$BestModelData)

climwin::plotall(dataset = ClimWinResults_All[[i]][[1]]$Dataset,
                 datasetrand = ClimWinRand_All[[i]][[1]],
                 bestmodel = ClimWinResults_All[[i]][[1]]$BestModel, 
                 bestmodeldata = ClimWinResults_All[[i]][[1]]$BestModelData)


# Linear regressions for local variable + peak roosting timing (July-Sep)
# =============================================================
static_window = read.csv(file.path(MainOut, "static_window.csv"))[2:62]

library(broom)
variable = colnames(static_window)[3:62]

pheno = read.csv("./results/5_FinalSignals/pheno_top5_df.csv") %>% 
  dplyr::select(c('Year', 'MADJDAYSWS'))

static_pheno = merge(static_window, pheno, by.x="year", by.y="Year")

summ_df_static = data.frame()
for (cVar in variable){
  print(cVar)
  regression = paste0("MADJDAYSWS ~ ", cVar, " + year")
  model = lm(as.formula(regression), data = static_pheno)
  summ_df_static = rbind(summ_df_static, 
                         cbind(tidy(model)[2,], 
                               glance(model), 
                               variable=strsplit(cVar, "_")[[1]][1]))
}

names(summ_df_static)[10] = "model.p.value"
names(summ_df_static)[9] = "model.statictic"
hist(summ_df_static$adj.r.squared)

# average across space for the static time + space approach 
static_df = data.frame()
for (i in c("tmin2m", "prate", "shum.2m")){
  df = data.frame(static_window[1], select(static_window, matches(i)) %>% rowMeans())
  colnames(df) = c("Year", i)
  df_pheno = merge(df, pheno, by="Year")
  regression = paste0("MADJDAYSWS ~ ", i, " + Year")
  model = lm(as.formula(regression), data = df_pheno)
  static_df = rbind(static_df, 
                         cbind(tidy(model)[2,], 
                               glance(model), 
                               variable=strsplit(i, "_")[[1]][1]))
}

library(ggdark)
plot2 = summ_df_static %>%
  ggplot(aes(x=adj.r.squared)) +
  geom_histogram(color="#e9ecef", alpha=1, position = 'identity', breaks = seq(0, 0.8, 0.05))+
  scale_y_continuous(breaks = seq(0, 30, by = 5), limits = c(0,30))+
  #scale_fill_viridis(discrete=TRUE) +
  theme_classic()+
  labs(fill="", x="Adjusted R2", y="Count")+
  #dark_theme_classic()+
  theme(text = element_text(size=18))

library(ggpubr)
ggarrange(plot2, plot1, ncol = 2, labels = c("A", "B"), align = "h")


RemainingSignals %>%
  ggplot(aes(x=Rsquares, fill=ClimateVar), color=ClimateVar) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', binwidth=0.025) +
  geom_vline(aes(xintercept=Rsquares, color = ClimateVar), data = RemSignals, alpha=0.7)+
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()+
  labs(fill="")

ggplot(RemainingSignals, aes(Rsquares, colour = ClimateVar)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', binwidth=0.025) +
  geom_vline(aes(xintercept=Rsquares, color = ClimateVar), data = RemSignals, alpha=0.7)+
  scale_fill_viridis(discrete=TRUE) +
  theme_classic()+
  labs(fill="")



#######################################################################
##################### Get linear model statistics #####################
#######################################################################

# Proof of concept

df = read.csv("./results/4_MultiVar/4_PhenDataWithCandidates/3_PhenologyDataset_MADJDAYSWS_Updated.csv")

lastFive = read.csv("./results/5_FinalSignals/top5.csv")[1:5,]

# Find matching columns
matching_columns <- df %>% 
  dplyr::select(Year, MADJDAYSWS, all_of(lastFive$VarName))
write.csv(matching_columns, "./results/5_FinalSignals/pheno_top5_df.csv")
pheno_top5_df = read.csv("./results/5_FinalSignals/pheno_top5_df.csv") %>% 
  gather(key="variable", value="value", 4:8) 

#PRESS - predicted residual sums of squares

PRESS <- function(linear.model) {
  #' calculate the predictive residuals
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  #' calculate the PRESS
  PRESS <- sum(pr^2)
  
  return(PRESS)
}

pred_r_squared <- function(linear.model) {
  #' Use anova() to get the sum of squares for the linear model
  lm.anova <- anova(linear.model)
  #' Calculate the total sum of squares
  tss <- sum(lm.anova$'Sum Sq')
  # Calculate the predictive R^2
  pred.r.squared <- 1-PRESS(linear.model)/(tss)
  
  return(pred.r.squared)
}

model = lm(paste0("MADJDAYSWS", " ~ ", paste(lastFive$VarName, collapse = " + ")),
           data = matching_columns)
summary(model)
pred_r_squared(model)
rockchalk::getDeltaRsquare(model)
rsq::rsq.partial(model, adj=F) #partial R^2

summary(lm(MADJDAYSWS ~ MADJDAYSWS_prate_yday_50_pixel233_Open258_Close230_sum, data = matching_columns))
pred_r_squared(lm(MADJDAYSWS ~ MADJDAYSWS_prate_yday_50_pixel233_Open258_Close230_sum, data = matching_columns))




# Use the last 3 predictors and predict onto station-level data
station.df = read.csv("./data/PhenologyDataset_StationLevel.csv")
env_station = matching_columns %>% 
  dplyr::select(-MADJDAYSWS) %>% 
  left_join(station.df %>% filter(Species=="yday_50"), by="Year")

sum.df = data.frame(Station=c(), YearLength=c(), rsquared=c(), AveMigrants=c())
for (s in unique(env_station$Station)){
  data = env_station %>% filter(Station==s)
  model.s = lm(MADJDAYSWS ~ MADJDAYSWS_EVIslope_yday_50_pixel2659_Open74_Close23_slope
               + MADJDAYSWS_tmin2m_yday_50_pixel200_Open67_Close29_mean
               + MADJDAYSWS_WindCF850_yday_50_pixel294_Open103_Close79_NrDaysWithWindComingFromTarget,
               data)
  summary(model.s)
  
  sum = list(Station = s, YearLength = nrow(data),
             rsquared = summary(model.s)$r.squared,
             AveMigrants = mean(data$total_migrants))
  sum.df = rbind(sum.df, sum)
}

# Are the distant precidctors better than the greatlake-region predictors (with the same time window and variable)?

tmin2m.pvalue = rast("./results/3_SpaWin_Adj/MADJDAYSWS/pValues/air.min_yday50_Climwin/tiffs/yday_50_pValues.tif")
plot(tmin2m.pvalue)



#######################################################################
################ Genrate files for ArcGIS Online #####################
#######################################################################

inputDir  = "./results/3_SpaWin_Adj/MADJDAYSWS/BestModelSumm"
DirsShort  = list.dirs(inputDir, full.names = F, recursive = F)

for (i in DirsShort){

rast = raster(paste0(inputDir, "/", i, "/AICc/yday_50.tif"))
## S4 method for signature 'RasterLayer'
KML(rast, paste0("./results/6_PlotKML/", i, "_AICc.kml"), col=rev(terrain.colors(255)), 
    colNA=NA, maxpixels=100000, blur=1, zip='', overwrite=FALSE)
}

#### ebird abundance map ####

inputDir  = "./data/maps/puma/purmar_abundance_seasonal_postbreeding-migration_mean_2021.tif"
DirsShort  = list.dirs(inputDir, full.names = T, recursive = F)

for (i in DirsShort){
  
  rast = raster(inputDir)
  rast_lonlat <- projectRaster(rast, crs="+proj=longlat +datum=WGS84")
  ## S4 method for signature 'RasterLayer'
  KML(rast_lonlat, paste0("./data/maps/puma/purmar_abundance_seasonal_postbreeding-migration_mean_2021.kml"), col=rev(terrain.colors(255)), 
      colNA=NA, maxpixels=80962434, blur=1, zip='', overwrite=T)
}


#### Input the tif files from k-fold results ####

inputDir  = "./results/8_Cross_Validation/Climwin/WindDirection/"
rast = raster(paste0(inputDir,"/AICc/WindDirection.tif"))
## S4 method for signature 'RasterLayer'
KML(rast, paste0(inputDir,"/AICc/WindDirection.kml"), col=rev(terrain.colors(255)), 
    colNA=NA, maxpixels=100000, blur=1, zip='', overwrite=FALSE)
