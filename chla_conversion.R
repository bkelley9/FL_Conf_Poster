library(tidyverse)
library(Metrics)
library(zoo)
library(ggpubr)

#Read 2023 Buoy WQ data
buoy_data <- read.table("data/WQ/2023 OwascoWinchInternet_PFL_Step.dat", sep = ",", skip = 4) %>%
  select(V1, V6, V7, V8, V9, V10, V11, V12, V13)

names(buoy_data) <- c("Date.Time", "Temp", "Depth", "SpCond", "PC", "Chla", "DOsat", "DO", "Turbidity")

#Define field days
field_days <- c("2023-05-16", "2023-06-01", "2023-06-13", "2023-06-27", "2023-07-11",
                "2023-07-25", "2023-08-08", "2023-08-23", "2023-09-05", "2023-09-19",
                "2023-10-03", "2023-10-17")

#Bin buoy data into discrete depths for comparison to Fluoroprobe data
#Plus some tidying of variables
buoy_data <- buoy_data %>%
  filter(grepl(paste(field_days, collapse = "|"), Date.Time)) %>%
  mutate(Date.Time = as.POSIXct(Date.Time)) %>%
  filter(format(Date.Time, "%H:%M:%S") > "02:00:00") %>%
  mutate(Date = as.Date(Date.Time),
         depth_bin = cut(Depth, breaks = c(0, 2.25, seq(3.75, 47.25, by = 1.5)))) %>%
  relocate(Date, depth_bin, Depth, .after = Date.Time) %>%
  select(-Date.Time) %>%
  group_by(Date, depth_bin) %>%
  mutate(across(Depth:Turbidity, ~ mean(.x) %>% round(4))) %>%
  ungroup() %>%
  unique()
  
#List raw FP and EXO2 field cast files
fp_files <- list.files("2023 Field Data/FPprofiles")
exo_files <- list.files("2023 Field Data/EXO2")

#Define function to read in a raw FP file, bin points, and average within bins
fp_func <- function(file){
  
  
  lower_depth_limit <- 0.5
  upper_depth_limit <- 47
  resolution <- 1.5
  
  FP_data  <- read.table(paste0("2023 Field Data/FPprofiles/", file), sep = "\t", skip = 2, 
                         col.names = c("DateTime", "Greens", "Cyano", "Diatoms", "Crypto", 
                                       "#5", "#6", "#7", "Yellow", "totChla", "Transmission",
                                       "Depth_FP", "Temperature", "GreenCells", "CyanoCells", "DiatomCells", 
                                       "CryptoCells", "#5cells", "#6cells", "#7cells", "Yellow2", 
                                       "totCellCt", "T700", "LED3", "LED4", "LED5", "LED6", "LED7", 
                                       "LED8", "Pressure", "TLED", "TSensor")) %>%
    mutate(velo = c(NA, diff(Pressure))) %>%
    mutate(accel = rollmean(velo, 3, fill = NA, align = "left")) %>%
    filter(accel > 0.01) %>%
    filter(Depth_FP > lower_depth_limit, Depth_FP < upper_depth_limit) %>%
    mutate(DateTime = as.POSIXct(DateTime, format = "%m/%d/%Y %H:%M:%S"),
           Date = as.Date(DateTime)) %>%
    mutate(depth_bin = cut(Depth_FP, breaks = c(0, 2.25, seq(3.75, 47.25, by = resolution)))) %>%
    select(Date, depth_bin, Depth_FP, Temperature, Greens, Cyano, Diatoms, Crypto, Yellow, totChla) %>%
    group_by(Date, depth_bin) %>%
    mutate(across(Depth_FP:totChla, ~ mean(.x) %>% round(4))) %>%
    ungroup() %>%
    unique()
}

#Define function to read in raw EXO2 files, bin data, and average within bins
exo_func <- function(file){
  EXO2_G <- read.csv(paste0("2023 Field Data/EXO2/", file), skip = 9, header = F) %>%
    select(V1, V2, V5, V7, V9, V11, V12, V14, V15, V17, V20)
  names(EXO2_G) <- c("Date", "Time", "Chla_lab", "Depth_lab", "ODO_sat_lab", "ODO_lab", "PSI", "SpCond_lab", "PC_lab", "Turbidity_lab", "Temp_lab")

  lower_depth_limit <- 0.5
  upper_depth_limit <- 47
  resolution <- 1.5
  
  EXO2_prof <- EXO2_G %>%
    mutate(velo = c(NA, diff(PSI))) %>%
    mutate(accel = rollmean(velo, 3, fill = NA, align = "left")) %>%
    filter(accel > 0.5) %>%
    filter(Depth_lab > lower_depth_limit, Depth_lab < upper_depth_limit) %>%
    mutate(depth_bin = cut(Depth_lab, breaks = c(0, 2.25, seq(3.75, 47.25, by = resolution))),
           Date = as.Date(Date, format = "%m/%d/%Y")) %>%
    select(Date, depth_bin, Depth_lab, Temp_lab, Chla_lab, PC_lab, ODO_sat_lab, ODO_lab, SpCond_lab, Turbidity_lab) %>%
    group_by(Date, depth_bin) %>%
    mutate(across(Depth_lab:Turbidity_lab, ~ mean(.x) %>% round(4))) %>%
    ungroup() %>%
    unique()
}

#Apply reading functions to all FP and EXO2 files
exo_data <- map_df(exo_files, exo_func)
fp_data <- map_df(fp_files, fp_func)

#Merge FP, Field EXO2, and Buoy EXO2 dataframes and create variables to be used in model
buoy_fp_merge <- buoy_data %>%
  merge(fp_data, by.x = c("Date", "depth_bin"), by.y = c("Date", "depth_bin")) %>%
  merge(exo_data, by.x = c("Date", "depth_bin"), by.y = c("Date", "depth_bin")) %>%
  filter(Depth < 30) %>%
  mutate(month = month(Date) %>% as.factor(),
         Chla = Chla * 0.6232,
         PC = PC * 0.30784 + 0.41656,
         logChla = log10(Chla),
         logtotChla = log10(totChla),
         logLABChla = log10(Chla_lab))

buoy_fp_merge %>%
  ggplot(aes(Turbidity_lab, Depth)) +
  geom_point() +
  scale_y_reverse() +
  facet_wrap(~Date, scales = "free")

#Set random seed for reproducibility and divide data into training and test sets
set.seed(123)
index <- sample(1:nrow(buoy_fp_merge), nrow(buoy_fp_merge)/2)
training_data <- buoy_fp_merge[index, ]
testing_data <- buoy_fp_merge[-index, ]

#train model
chla_model <- lm(logtotChla ~ logChla, data = training_data)