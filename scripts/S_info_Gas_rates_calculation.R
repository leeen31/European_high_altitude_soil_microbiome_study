###
### AlpineSoils23 - Gas rates calculation
###

# produce Fig5_E_AlpineSoils23_gas_rates.tsv for gas fluxplot in Fig5E


rm(list=ls(all=TRUE)) # removes everything


# load data and libraries ####
library(tidyverse)
library(data.table)
library(ggplot2)





# load gas raw data from both measurement batches
gas_conc_1 <- read.table("data/S_info_AlpineSoils23_gas_raw_data_batch1.data", sep="\t", header=T, quote=NULL, comment='')
gas_conc_2 <- read.table("data/S_info_AlpineSoils23_gas_raw_data_batch2.data", sep="\t", header=T, quote=NULL, comment='')

# concatenate both
gas_cons <- rbind(gas_conc_1, gas_conc_2)





# remove decimals
cols_to_round <- c("H2O", "CO2", "CH4")
gas_cons[cols_to_round] <- lapply(gas_cons[cols_to_round], round)

# !!! convert Ch4 to ppm !!!
gas_cons$CH4_ppm <- gas_cons$CH4 / 1000


# load times and metadata table
times <- read.table("data/S_info_AlpSoils23_incubation_gas_times_CO2_CH4.txt",
                    sep="\t", header=TRUE, quote="", comment.char="", row.names = NULL)




### 1. Plotting raw data and filtering time ---------------------

# format time
# convert time object
gas_cons$TIME <- as.ITime(gas_cons$TIME)


# Define start and end times of the treatment period
start_time <- as.ITime("06:00:00")  # Example start time
end_time <- as.ITime("18:00:00")    # Example end time

# Filter the dataframe for the relevant time period
filtered_gas_data <- gas_cons[gas_cons$TIME >= start_time & gas_cons$TIME <= end_time, ]



# Calculate time difference in seconds
filtered_gas_data$Seconds <- as.numeric(difftime(filtered_gas_data$TIME, filtered_gas_data$TIME[1], 
                                                 units = "secs"))


#label measurements
filtered_gas_data$ID <- NA



for(i in 1:length(times$Start_1)){
  rep_start <- as.ITime(times[i, ]$Start_1)
  rep_end <- as.ITime(times[i, ]$End_1)
  
  # Assign subsite to gas data
  filtered_gas_data[filtered_gas_data$TIME > rep_start & filtered_gas_data$TIME < rep_end,]$ID <- times[i, ]$ID
}



# plot slopes to check times
# CO2
ggplot(filtered_gas_data, aes(x = Seconds, y = CO2)) +
  geom_point() + 
  geom_smooth(aes(colour = ID), method = "lm", fill = NA, fullrange = FALSE) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# CH4
ggplot(filtered_gas_data, aes(x = Seconds, y = CH4_ppm)) + 
  geom_point() + 
  geom_smooth(aes(colour = ID), method = "lm", fill = NA, fullrange = FALSE) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




### 2. Slope and Flux calculation ---------------------------------------

# clean filtered_gas_data
clean_gas_data <- filtered_gas_data %>% 
  drop_na() %>%
  select(CO2, CH4_ppm, Seconds, ID)

gas_meta <- times %>%
  select(ID, temperature, replicate, soil_volume, headspace_L )

library(dplyr)
library(tidyr)

# Constants
R <- 0.0831446261815324  # L·bar/K/mol
P_bar <- 1  # Atmospheric pressure in bar
dead_volume_L <- 0.028  # 28 mL = 0.028 L
tubing_volume_L <- ((pi * (0.2^2) * 80.25) / 1000) * 2  # cm³ to mL to L
atmospheric_conc <- c(CH4 = 1.907, CO2 = 419)  # ppm (adjust as needed)

# Clean and prepare data
clean_gas_data <- clean_gas_data %>%
  drop_na() %>%  # Remove rows with any NA
  select(ID, Seconds, CO2, CH4_ppm) %>%
  rename(CH4 = CH4_ppm)

# Merge metadata
merged <- clean_gas_data %>%
  left_join(gas_meta, by = "ID")

# Compute mols of gas in headspace
merged <- merged %>%
  mutate(
    T_K = temperature + 273.15,
    adjusted_volume_L = headspace_L + dead_volume_L + tubing_volume_L,
    mol_headspace = (P_bar * adjusted_volume_L) / (R * T_K)
  )

# Convert concentrations to μmol
merged <- merged %>%
  mutate(
    CO2_umol = CO2 * mol_headspace,
    CH4_umol = CH4 * mol_headspace
  )


# Reshape to long format
gas_long <- merged %>%
  pivot_longer(cols = c(CO2_umol, CH4_umol),
               names_to = "gas",
               values_to = "umol") %>%
  mutate(
    gas_short = sub("_umol", "", gas)
  )


### slope calculation

# Define slope function
get_slope <- function(df) {
  model <- lm(umol ~ Seconds, data = df)
  coef(model)[["Seconds"]]
}



# Compute slopes in umol/min
slopes <- gas_long %>%
  group_by(ID, gas_short) %>%
  summarise(
    slope_umol_s = get_slope(cur_data()),
    slope_umol_min = slope_umol_s * 60,
    .groups = "drop"
  )

# Add dry weight and mol_headspace per Rep
meta_summary <- merged %>%
  group_by(ID) %>%
  summarise(
    soil_fresh_weight = mean(soil_volume, na.rm = TRUE),
    mol_headspace = mean(mol_headspace, na.rm = TRUE),
    .groups = "drop"
  )


### 3. Compute final oxidation rate in μmol/g/h --------------------------

# 2. Join with metadata and calculate oxidation rate
rates <- slopes %>%
  left_join(meta_summary, by = "ID") %>%
  mutate(
    gas_rate_umol_min = slope_umol_min,
    gas_rate_umol_g_h = gas_rate_umol_min * 60 / soil_fresh_weight
  )



### 4. Export rates table --------------------------

# add site column 
rates <- rates %>%
  mutate(
    # extract the site part: 1 or 2 characters depending on prefix
    site_code = case_when(
      substr(ID, 1, 2) %in% c("C1", "C2") ~ substr(ID, 1, 2),  # handle C1/C2
      TRUE ~ substr(ID, 1, 1)                                    # everything else: first letter
    ),
    site = recode(
      site_code,
      "A"  = "SF",
      "B"  = "MF",
      "C1" = "MY1",
      "C2" = "MY2",
      "F"  = "SN", 
      "G"  = "DS", 
      "H"  = "CD",
      "J"  = "GR",
      "K"  = "BN",
      "L"  = "PM"
    )
  ) %>%
  select(-site_code)  # optional helper column

# export
#write.table(rates, "Fig5_E_AlpineSoils23_gas_rates.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



