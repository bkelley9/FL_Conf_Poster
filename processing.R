library(tidyverse)
library(reshape2)
library(zoo)
library(mgcv)
library(rLakeAnalyzer)

#Load model for 2023 Chlorophyll
source("chla_conversion.R")

#List raw buoy files for WQ data, wind data, and remaining Met data
wq_files <- list.files("data/WQ")
met_nowind_files <- list.files("data/Met", pattern = "TRB")
wind_files <- list.files("data/Met", pattern = "Wind")

#Define functions to read in raw files
met_nowind_import <- function(filename){
  data <- read.table(paste("data/Met", filename, sep = "/"), sep = ",", skip = 4)
  
  year_extract <- str_extract(filename, "\\b\\d{4}\\b")
  
  names(data) <- c("Date.Time", "Record", "AirTemp", "Rhumid", "Pressure", "SolarRad", "Precip1", "Precip2")
  
  data <- data %>%
    filter(Date.Time > year_extract) %>%
    select(-contains("Precip")) %>%
    mutate(across(AirTemp:SolarRad, ~ as.numeric(.x))) %>%
    mutate(Date.Time = as.POSIXct(Date.Time),
           Date = as.Date(Date.Time))
}
wind_import <- function(filename){
  data <- read.table(paste("data/Met", filename, sep = "/"), sep = ",", skip = 4)
  
  year_extract <- str_extract(filename, "\\b\\d{4}\\b")
  
  names(data) <- c("Date.Time", "Record", "WindDir", "WindSpeed")
  
  data <- data %>%
    filter(Date.Time > year_extract) %>%
    mutate(Date.Time = as.POSIXct(Date.Time),
           Date = as.Date(Date.Time)) %>%
    select(-Date.Time)
}
wq_import <- function(filename){
  data <- read.table(paste("data/WQ", filename, sep = "/"), sep = ",", skip = 4)
  
  year_extract <- str_extract(filename, "\\b\\d{4}\\b")
  
  data <- data %>%
    select(V1, V6, V7, V8, V9, V10, V11, V12, V13)
  
  names(data) <- c("Date.Time", "Temp", "Depth", "SpCond", "PC", "Chla", "DOsat", "DO", "Turbidity")
  
  data <- data %>%
    filter(Date.Time > year_extract)
}


#Read in wind data
wind_data <- map_df(wind_files, wind_import)


#Read in remaining Met data and merge with wind data
met_data <- map_df(met_nowind_files, met_nowind_import) %>%
  merge(wind_data, by.x = c("Date", "Record"), by.y = c("Date", "Record")) %>%
  relocate(Date.Time, .before = Date) %>%
  arrange(Date.Time) %>%
  mutate(numeric_datetime = as.numeric(Date.Time),
         Time = format(Date.Time, format = "%H:%M:%S"), .after = Date.Time)

#Read in WQ data
wq_data <- map_df(wq_files, wq_import) %>%
  arrange(Date.Time) %>%
  mutate(Date.Time = as.POSIXct(Date.Time, tz = "EST"),
         numeric_datetime = as.numeric(Date.Time)) %>%
  relocate(numeric_datetime, Depth, .after = Date.Time) %>%
  mutate(PC = PC / (1 - (Temp - 25))) %>%
  pivot_longer(cols = Temp:Turbidity, names_to = "variable") %>%
  filter(!is.na(value), !is.infinite(value))

#List slope-intercept coefficients for previous year calibrations done by John
Chla_coeff <- list(
  `2014` = c(3.3554, 1.7857),
  `2015` = c(-0.5545, 1.5292),
  `2016` = c(0.0596, 0.2874),
  `2017` = c(-0.738, 1.9633),
  `2018` = c(0.8599, 1.3519),
  `2019` = c(-0.0039, 0.6065),
  `2020` = c(0.5352, 0.7428),
  `2021` = c(0.3071, 0.4272),
  `2022` = c(0.2872, 0.4713),
  `2023` = c(0, 1)
)

#List slope-intercept coefficients for previous year calibrations done by John.
#Years with (0, 1) were not calibrated
Turbidity_coeff <- list(
  `2014` = c(0, 1),
  `2015` = c(0, 1),
  `2016` = c(0, 1),
  `2017` = c(0.465, 1.0704),
  `2018` = c(0.01443, 0.6939),
  `2019` = c(0.5352, 0.7428),
  `2020` = c(4.2679, 1.3672),
  `2021` = c(3.3343, 2.6139),
  `2022` = c(2.6444, 2.56),
  `2023` = c(0, 1)
)

#Define function to apply coefficients to raw values
conv_func <- function(data, coeff){
  data %>%
    mutate(value = (value + coeff[1])/coeff[2])
}

#Convert 2014-2022 Chlorophyll data. 2023 converted using chla_model.
Chla_year_corr <- wq_data %>%
  filter(variable == "Chla") %>%
  mutate(Year = year(Date.Time)) %>%
  group_split(Year) %>%
  map2_df(Chla_coeff, conv_func) %>%
  filter(Year != "2023") %>%
  select(-Year)

#Convert Turbidity data
Turbidity_year_corr <- wq_data %>%
  filter(variable == "Turbidity") %>%
  mutate(Year = year(Date.Time)) %>%
  group_split(Year) %>%
  map2_df(Turbidity_coeff, conv_func) %>%
  select(-Year)

#Convert 2023 Chlorophyll data
Chla_2023 <- wq_data %>%
  filter(year(Date.Time) == "2023") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate(logChla = log10(Chla)) %>%
  mutate(logtotChla = predict(chla_model, newdata = .),
         Chla = 10^logtotChla) %>%
  select(-c(logChla, logtotChla)) %>%
  pivot_longer(cols = Temp:Turbidity, names_to = "variable") %>%
  filter(variable == "Chla")

pre_WQ_QAQC_data <- wq_data %>%
  filter(variable != "Turbidity", variable != "Chla") %>%
  bind_rows(Chla_year_corr) %>%
  bind_rows(Chla_2023) %>%
  bind_rows(Turbidity_year_corr)

#Define gross cutoffs for QAQC of WQ data (points below/above these thresholds are flagged)
wq_gross_cutoffs <- list(
  Chla = c(-50, 500),
  DO = c(0, 20),
  DOsat = c(0, 150),
  PC = c(-20, 20),
  SpCond = c(250, 400),
  Temp  = c(0, 30),
  Turbidity = c(-2, 10))

#Define function for rolling average/standard deviation of points.
qc_func <- function(data, params, hour_window){
  
  data %>%
    mutate(width = seq_along(numeric_datetime) - findInterval(numeric_datetime - hour_window*3600, numeric_datetime),
           rolling_mean = rollapplyr(value, width, mean),
           rolling_sd = rollapplyr(value, width, sd)) %>%
    #if value passes BOTH gross cutoff and rolling checks, Safe
    mutate(rolling_flag = ifelse(value > rolling_mean - 3*rolling_sd & value < rolling_mean + 3*rolling_sd, "Safe", "FR"),
           gross_flag = ifelse(value >= params[1] & value <= params[2], "Safe", "FG")) %>%
    unite(Flag, rolling_flag, gross_flag)
}

#Split WQ data into a list of dataframes by variable and apply QAQC functions.
#Rolling window set at 24 hours
WQ_QAQC_data <- pre_WQ_QAQC_data %>% 
  group_split(variable) %>%
  map2_df(wq_gross_cutoffs, qc_func, 24, .progress = T)

#Remove points that don't pass flagging criteria, pivot data back to wide form
#and tidy data
wq2 <- WQ_QAQC_data %>%
  filter(Flag == "Safe_Safe") %>%
  select(-c(Flag, width,  rolling_mean, rolling_sd)) %>%
  unique() %>%
  pivot_wider(names_from = variable, values_from = value)

#Tidy data for lake physics calculations
wq3 <- wq2 %>%
  mutate(Date = as.Date(Date.Time),
         Time = format(Date.Time, format = "%H:%M:%S"),
         Cast = ifelse(Time < "11:00:00", "Midnight", "Noon"), .after = Date.Time)

#Define function that will take a cast and calculate all physics variables.
#This is done by taking a cast, interpolating temperature and chla depth casts
#to 500 point curve and calculating variables.
physics_fun <- function(data){
  
  origin_data <- data %>%
    select(Depth, Temp, Chla) 
  
  thermo_data <- origin_data %>%
    select(Depth, Temp) %>%
    group_by(Depth) %>%
    mutate(Temp = mean(Temp, na.rm = T)) %>%
    unique()

  
  temp_profile <- approx(origin_data$Temp, origin_data$Depth, method = "linear", n = 500)
  temp_profile <- data.frame(Temp = temp_profile[[1]], Depth = temp_profile[[2]])
  
  chla_profile <- approx(origin_data$Chla, origin_data$Depth, method = "constant", n = 500)
  chla_profile <- data.frame(Chla = chla_profile[[1]], Depth = chla_profile[[2]])
  
  mixed_layer3 <- temp_profile %>%
    arrange(Depth) %>%
    filter((Temp - first(Temp)) < - 0.3) %>%
    slice(1)
  
  mixed_layer5 <- temp_profile %>%
    arrange(Depth) %>%
    filter((Temp - first(Temp)) < - 0.5) %>%
    slice(1)
  
  if (nrow(mixed_layer3) < 1) {
    mixed_layer3 <- tibble(Temp = NA, Depth = NA)
  } else {
    mixed_layer3
  }
  
  if (nrow(mixed_layer5) < 1) {
    mixed_layer5 <- tibble(Temp = NA, Depth = NA)
  } else {
    mixed_layer5
  }
  
  chla_layer <-  chla_profile %>%
    arrange(desc(Chla)) %>%
    slice(1)
  
  thermo_layer <- thermo.depth(wtr = thermo_data$Temp, depths = thermo_data$Depth, seasonal = T, mixed.cutoff = 8)
  
  if (is.nan(thermo_layer)== T) {
    thermo_layer <- NA
    thermo_match <- tibble(Temp = NA, Depth = NA)
  } else {
    thermo_match <- temp_profile %>%
      mutate(Tt_extract = abs(Depth - thermo_layer)) %>%
      arrange(Tt_extract) %>%
      slice(1)
  }
  
  meta_vals = meta.depths(wtr = thermo_data$Temp, depths = thermo_data$Depth, mixed.cutoff = 8)
  
  if (is.nan(meta_vals[1]) == T | sum(meta_vals > 40) > 0){
    meta_vals <- c(NA, NA)
    metaTopMatch <- tibble(Temp = NA, Depth = NA)
    metaBotMatch <- tibble(Temp = NA, Depth = NA)
  } else {
    metaTopMatch <- temp_profile %>%
      mutate(T_metaTop_extract = abs(Depth - meta_vals[1])) %>%
      arrange(T_metaTop_extract) %>%
      slice(1)
    
    metaBotMatch<- temp_profile %>%
      mutate(T_metaBot_extract = abs(Depth - meta_vals[2])) %>%
      arrange(T_metaBot_extract) %>%
      slice(1)
  }
  
  
  data <- data %>%
    mutate(h_ml3 = mixed_layer3$Depth, 
           T_ml3 = mixed_layer3$Temp,
           h_ml5 = mixed_layer5$Depth, 
           T_ml5 = mixed_layer5$Temp,
           chla_max = chla_layer$Chla,
           chla_max_depth = chla_layer$Depth,
           h_t = thermo_layer, 
           T_t = thermo_match$Temp,
           h_metaT = meta_vals[1],
           T_metaT = metaTopMatch$Temp,
           h_metaB = meta_vals[2],
           T_metaB = metaBotMatch$Temp)
}

#Apply physics function to each individual cast by splitting dataframe by Date and Cast
wq4 <- wq3 %>%
  group_split(Date, Cast) %>%
  map_df(physics_fun, .progress = T)

#met_noon grabs the met data closest to the end of a depth profile performed by the buoy
met_noon <- met_data %>%
  filter(Time > "12:30:00", Time < "12:40:00")
#Merge met_noon with met_midnight, now we have the meteorological data to pair
#with noon and midnight buoy casts
met_midnight_noon <- met_data %>%
  filter(Time > "00:30:00", Time < "00:40:00") %>%
  bind_rows(met_noon) %>%
  arrange(Date.Time) %>%
  mutate(Cast = ifelse(Time < "11:00:00", "Midnight", "Noon"), .after = Time)

#Merge met data with WQ data to calculate wind shearing, natural convection, etc.
met_wq_merge <- wq4 %>%
  mutate(Epi = ifelse(Depth < h_ml3, "Epi", "Hypo")) %>%
  group_by(Date, Cast) %>%
  mutate(Tsurf = Temp[which.min(Depth)]) %>%
  ungroup() %>%
  mutate(DCL_present = ifelse(chla_max_depth > h_ml3, "Y", "N"),
         Surf_density = (1000*(1-(Tsurf+288.9414)*(Tsurf-3.9863)^2/(508929.2*(Tsurf+68.12963))))) %>%
  group_by(Date, Cast, Epi) %>%
  mutate(Epi_temp = mean(Temp, na.rm = T),
         Epi_chla = ifelse(Epi == "Epi", mean(Chla, na.rm = T), NA)) %>%
  ungroup() %>%
  mutate(Epi_density = (1000*(1-(Epi_temp+288.9414)*(Epi_temp-3.9863)^2/(508929.2*(Epi_temp+68.12963))))) %>%
  filter(Epi == "Epi") %>%
  select(-c(Date.Time, Time, numeric_datetime, Depth, Chla, DO, DOsat, PC, SpCond, Temp, Turbidity, Epi)) %>%
  merge(met_midnight_noon, by = c("Date", "Cast")) %>%
  mutate(delta_T = Tsurf - AirTemp) %>%
  filter(month(Date) >= 6, month(Date) <= 10,
         year(Date) != "2022") %>%
  unique()

#define constants for equations from Taylor DCL paper
wndHeight <- 3.5
vonK <- 0.4
rhoAir <- 1.2
kvisc <- 1.0e-6 
Cp_h20 <- 4186 #J/kg/C
g <- 9.81 #m/s2
alpha <- 2.07e-4 #1/K
k <- 5.98e-1 #W/m/K

#Calculate reynolds number and reynolds number from natural convection
reynolds_calc <- met_wq_merge %>%
  mutate(Cd = ifelse(WindSpeed <= 5, 0.0015, 0.001),
         WindSpeed_corr = WindSpeed/(1-sqrt(Cd)/vonK*log(10/wndHeight)),
         tau = Cd*rhoAir*WindSpeed_corr^2,
         uStar = sqrt(tau/Epi_density),
         Re = uStar/kvisc) %>%
  mutate(lamda_t = 5.5*uStar / kvisc,
         Dt = (k*h_ml3) / (Surf_density*Cp_h20*lamda_t),
         Re_nc = ifelse(delta_T >= 0, (g*alpha*Dt*(delta_T))^(1/3) / kvisc, NA)) %>%
  mutate(convection = ifelse(delta_T >= 0, "On", "Off"),
         Year = year(Date) %>% as.numeric(),
         Month = month(Date) %>% as.factor(),
         DOY = yday(Date))

#Standardized dataframe for WQ variables. 2014 and 2022 have missing data at 
#beginning and end of year respectively
may15_to_oct31 <- wq4 %>%
  mutate(Epi = ifelse(Depth < h_ml3, "Epi", "Hypo"),
         Year = year(Date) %>% factor()) %>%
  filter(yday(Date) > 134, yday(Date) < 305) %>%
  group_by(Date, Cast, Epi) %>%
  mutate(Epi_temp = ifelse(Epi == "Epi", mean(Temp, na.rm = T), NA),
         Epi_chla = ifelse(Epi == "Epi", mean(Chla, na.rm = T), NA),
         Epi_DO = ifelse(Epi == "Epi", mean(DOsat, na.rm = T), NA))
