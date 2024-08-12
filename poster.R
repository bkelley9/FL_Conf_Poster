source("processing.R")


met_noon <- met_data %>%
  filter(Time > "12:30:00", Time < "12:40:00")
met_midnight_noon <- met_data %>%
  filter(Time > "00:30:00", Time < "00:40:00") %>%
  bind_rows(met_noon) %>%
  arrange(Date.Time) %>%
  mutate(Cast = ifelse(Time < "11:00:00", "Midnight", "Noon"), .after = Time)

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

#define constants
wndHeight <- 3.5
vonK <- 0.4
rhoAir <- 1.2
kvisc <- 1.0e-6 
Cp_h20 <- 4186 #J/kg/C
g <- 9.81 #m/s2
alpha <- 2.07e-4 #1/K
k <- 5.98e-1 #W/m/K


reynolds_calc <- met_wq_merge %>%
  mutate(Cd = ifelse(WindSpeed <= 5, 0.0015, 0.001),
         WindSpeed_corr = WindSpeed/(1-sqrt(Cd)/vonK*log(10/wndHeight)),
         tau = Cd*rhoAir*WindSpeed_corr^2,
         uStar = sqrt(tau/Epi_density),
         Re = uStar/kvisc) %>%
  mutate(lamda_t = 5.5*uStar / kvisc,
         Dt = (k*h_ml) / (Surf_density*Cp_h20*lamda_t),
         Re_nc = ifelse(delta_T >= 0, (g*alpha*Dt*(delta_T))^(1/3) / kvisc, NA)) %>%
  mutate(convection = ifelse(delta_T >= 0, "On", "Off"),
         Year = year(Date) %>% as.numeric(),
         Month = month(Date) %>% as.factor(),
         DOY = yday(Date))


wq5 <- wq4 %>%
  mutate(Epi = ifelse(Depth < h_ml, "Epi", "Hypo"),
         Year = year(Date) %>% factor()) %>%
  filter(yday(Date) > 134, yday(Date) < 305) %>%
  group_by(Date, Cast, Epi) %>%
  mutate(Epi_temp = ifelse(Epi == "Epi", mean(Temp, na.rm = T), NA),
         Epi_chla = ifelse(Epi == "Epi", mean(Chla, na.rm = T), NA),
         Epi_DO = ifelse(Epi == "Epi", mean(DOsat, na.rm = T), NA))


Chla_confdata <- wq_data %>%
  filter(year(Date.Time) == "2023") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate(logChla = log10(Chla))

predictions = predict(chla_model, newdata = Chla_confdata, interval = "confidence")

Chla_confdata <- Chla_confdata %>%
  mutate(Fit = 10^predictions[,1],
         Lower = 10^predictions[,2],
         Upper = 10^predictions[,3]) %>%
  select(Date.Time, Depth, Temp, Chla, Fit, Lower, Upper)

Chla_confdata %>%
  ggplot(aes(Chla, Fit)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.1) +
  labs(x = "RFUs",
       y = "\u00b5g/L") +
  theme_pubclean()

Chla_confdata %>%
  mutate(Date = as.Date(Date.Time),
         Time = format(Date.Time, format = "%H:%M:%S"),
         Cast = ifelse(Time < "11:00:00", "Midnight", "Noon"), .after = Date.Time) %>%
  filter(Date > "2023-05-24", Date < "2023-06-14", Cast == "Noon") %>%
  arrange(Date, Depth) %>%
  pivot_longer(cols = Fit:Upper) %>%
  ggplot(aes(value, Depth, color = name)) +
  scale_color_manual(values = c("red", "gray", "gray")) +
  geom_path() +
  scale_y_reverse() +
  facet_wrap(~Date, scales = "free") +
  labs(x = "Chl-a (\u00b5g/L)",
       y = "Depth",
       color = "Prediction") +
  theme_pubclean()

preRFU_WQ_QAQC_data <- wq_data %>%
  filter(variable != "Turbidity") %>%
  bind_rows(Turbidity_year_corr)

wq_gross_cutoffs <- list(
  Chla = c(-1, 500),
  DO = c(0, 20),
  DOsat = c(0, 150),
  PC = c(-20, 20),
  SpCond = c(250, 400),
  Temp  = c(0, 30),
  Turbidity = c(-2, 10))

met_gross_cutoffs <- list(
  AirTemp = c(-10, 35),
  Pressure = c(960, 1015),
  Rhumid = c(0, 100),
  SolarRad = c(0, 1200),
  WindDir = c(0, 360),
  WindSpeed = c(0, 18)
)

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

RFU_QAQC_data <- map2_df(preRFU_WQ_QAQC_data %>% group_split(variable), wq_gross_cutoffs, qc_func, 24, .progress = T)

RFU <- RFU_QAQC_data %>%
  filter(Flag == "Safe_Safe") %>%
  select(-c(Flag, width,  rolling_mean, rolling_sd)) %>%
  unique() %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  filter(year(Date.Time) == "2023", Chla < 10.061) 


met_data_Lisa <- met_data %>%
  filter(yday(Date) > 182, yday(Date) < 273) %>%
  filter(!year(Date) %in% c("2022", "2014")) %>%
  filter(Time > "06:00:00" & Time < "18:00:00") %>%
  mutate(Year = year(Date)) %>%
  group_by(Year) %>%
  mutate(AirTemp = mean(AirTemp, na.rm = T),
         SolarRad = mean(SolarRad, na.rm = T)) %>%
  select(Year, AirTemp, SolarRad) %>%
  unique()

met_2022 <- read.csv("22OwascoMet.csv") %>%
  select(1, 2, 3, 9)

names(met_2022) <- c("Date", "Time", "AirTemp", "SolarRad")

met_2022 <- met_2022 %>%
  mutate(Date = as.Date(Date),
         Year = year(Date)) %>%
  filter(yday(Date) > 182, yday(Date) < 273, year(Date) != "2014") %>%
  filter(Time > "06:00:00" & Time < "18:00:00") %>%
  group_by(Year) %>%
  mutate(AirTemp = mean(AirTemp, na.rm = T),
         SolarRad = mean(SolarRad, na.rm = T)) %>%
  select(-c(Date, Time)) %>%
  unique() %>%
  bind_rows(met_data_Lisa)

surf_temp <- wq4 %>%
  filter(yday(Date) > 182, yday(Date) < 273) %>%
  mutate(Depth < 1.2) %>%
  select(Date, Temp) %>%
  mutate(Year = year(Date)) %>%
  group_by(Year) %>%
  mutate(Surf_Temp = mean(Temp, na.rm = T),
         count = n()) %>%
  select(Year, Surf_Temp) %>%
  unique()

model_data <- surf_temp %>%
  merge(met_2022, by = "Year")

solar_model <- lm(Surf_Temp ~ SolarRad, data = model_data)
airtemp_model <- lm(Surf_Temp ~ AirTemp, data = model_data)

model_data %>%
  ggplot(aes(SolarRad, AirTemp, color = factor(Year))) +
  geom_point() +
  labs(x = "Air Temperature (\u00B0C)",
       y = "Surface Temperature (\u00B0C)",
       color = "Year") +
  theme_pubclean()

buoy_fp_merge %>%
  ggplot(aes(x = Turbidity, y = Depth)) +
  geom_point() +
  scale_y_reverse() +
  facet_wrap(~Date, scales = "free")

data <- read.csv("0308Evan.csv") %>%
  filter(grepl("FLNF", Sample_ID)) %>%
  arrange(Analyte, Sample_ID)
write.csv(data, "FLNFPonds.csv", row.names = F)
