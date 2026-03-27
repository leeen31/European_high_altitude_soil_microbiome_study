###
### AlpSoil23 - table for soil chemistry
###


rm(list=ls(all=TRUE)) # removes everything

library(formattable)
library(dplyr)

# soil metadata
soil_data <- read.csv("data/Fig1_B_AlpineSoils23_metadata_vegetation90.txt", header= TRUE, sep= "\t")



# Remove unwanted columns
soil_data_clean <- soil_data %>%
  select(c(2,6,7,9,12,13,21,23)) 

# Summarise by site: mean only (numeric columns)
soil_summary <- soil_data_clean %>%
  group_by(site) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  mutate(
    altitude = round(altitude, 0),
    TN = round(TN, 3),
    TOC = round(TOC, 2),
    C_N = round(C_N, 1),
    GWC = round(GWC, 1),
    pH = round(pH, 1), 
    vegetation = round(vegetation, 1)
  )

# Rename sites
soil_summary <- soil_summary %>%
  mutate(site = recode(site,
                       "A" = "SF",
                       "B" = "MF",
                       "C1" = "MY1",
                       "C2" = "MY2",
                       "F" = "SN",
                       "G" = "DS",
                       "H" = "CD",
                       "J" = "GR",
                       "K" = "BN",
                       "L" = "PM")) 

# add Site name and coordinates
site_info <- data.frame(
  site = c("SF", "MF", "MY1", "MY2", "SN", "DS", "CD", "GR", "BN", "PM"),
  site_long = c(
    "Saas Fee", "Mont Fort", "Moiry 1", "Moiry 2", 
    "Schilthorn", "Diablerets Scex Rouge", "Diablerets", 
    "Grächen", "Bettmerhorn", "Plaine Morte"
  ),
  coordinates_short = c(
    "46.08110, 7.91409",
    "46.08385, 7.30274",
    "46.09091, 7.59543", 
    "46.09111, 7.59537",
    "46.55827, 7.83292",
    "46.32636, 7.20837",
    "46.33825, 7.21536",
    "46.17537, 7.85005",
    "46.41149, 8.07524",
    "46.37158, 7.48559"
  ),
  Lithology = c(
    "Basic rock", "Mica shists", "Granites", "Granites", 
    "Marly shales", "Limestone", "Limestone", 
    "Gneiss", "Granites", "Limestone"
  ),
  Ecoregion = c(
    "Continental alps", "Continental alps", "Continental alps", "Continental alps", 
    "Northern pre-alps", "Northern pre-alps", "Northern pre-alps", 
    "Continental alps", "Continental alps", "Northern pre-alps"
  ),
  stringsAsFactors = FALSE
)



# join site info and reorder
soil_summary_full <- soil_summary %>%
  dplyr::left_join(site_info, by = "site") %>%
  dplyr::select(site, site_long, Lithology, Ecoregion, altitude, vegetation, pH, GWC, C_N, TOC, TN)

# Ascending order (smallest to largest)
soil_summary_full <- soil_summary_full %>%
  arrange(altitude)


# rename columns
soil_summary_clean <- soil_summary_full %>%
  rename(
    ID = site,
    Name = site_long,
    `Altitude (m.a.s.l.)` = altitude,
    `Vegetation (%)` = vegetation,
    `SWC (%)` = GWC,
    #Coordinates = coordinates_short,
    `TOC (%)` = TOC,
    `TN (%)` = TN,
    `C:N` = C_N
  )



# make nice table
formattable(
  soil_summary_clean,
  list(
    `Vegetation (%)` = color_tile("white", "forestgreen"),
    `Altitude (m.a.s.l.)` = color_tile("lightblue1", "lightblue3"),
    #Custom color rule for pH
    pH = formatter("span",
                   style = x ~ style(
                     display = "block",
                     padding = "0 4px",
                     `border-radius` = "4px",
                     `background-color` = ifelse(x < 6.5, "gold", "orange")
                   )),
    `SWC (%)` = color_tile("lightcyan", "cyan2"),
    `C:N` = color_tile("plum1", "plum3"),
    `TOC (%)` = color_tile("khaki1", "khaki3"),
    `TN (%)` = color_tile("white", "pink2")
  )
)






### calculate  means and std table

soil_summary <- soil_data_clean %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(
    across(
      where(is.numeric),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd   = ~ sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# Rename sites
soil_summary <- soil_summary %>%
  mutate(site = recode(site,
                       "A" = "SF",
                       "B" = "MF",
                       "C1" = "MY1",
                       "C2" = "MY2",
                       "F" = "SN",
                       "G" = "DS",
                       "H" = "CD",
                       "J" = "GR",
                       "K" = "BN",
                       "L" = "PM")) 


soil_summary <- soil_summary %>%
  dplyr::mutate(
    across(ends_with("_mean"), ~ round(.x, 2)),
    across(ends_with("_sd"),   ~ round(.x, 2))
  )

soil_summary_fmt <- soil_summary %>%
  dplyr::mutate(
    altitude = sprintf("%.0f ± %.0f", altitude_mean, altitude_sd),
    TN       = sprintf("%.3f ± %.3f", TN_mean, TN_sd),
    C_N       = sprintf("%.3f ± %.3f", C_N_mean, C_N_sd),
    TOC       = sprintf("%.2f ± %.2f", TOC_mean, TOC_sd),
    SWC      = sprintf("%.1f ± %.1f", GWC_mean, GWC_sd),
    pH       = sprintf("%.1f ± %.1f", pH_mean, pH_sd),
    vegetation = sprintf("%.1f ± %.1f", vegetation_mean, vegetation_sd)
  )

soil_summary_full <- soil_summary_fmt %>%
  dplyr::left_join(site_info, by = "site") %>%
  dplyr::select(
    site,
    site_long,
    Lithology, Ecoregion,
    altitude,
    vegetation,
    pH,
    SWC,
    C_N,
    TOC,
    TN
  ) %>%
  arrange(as.numeric(sub(" ±.*", "", altitude)))

soil_summary_clean <- soil_summary_full %>%
  dplyr::rename(
    ID = site,
    Name = site_long,
    `Altitude (m.a.s.l.)` = altitude,
    `Vegetation (%)` = vegetation,
    `SWC (%)` = SWC,
    `TOC (%)` = TOC,
    `TN (%)` = TN,
    `C:N` = C_N
  )

formattable(soil_summary_clean)



### calculate pH difference in MY

MY_pH <- subset(soil_data, site %in% c("C1", "C2") & !is.na(pH))

# Normality per group
by(MY_pH$pH, MY_pH$site, shapiro.test)

# Variance equality
var.test(pH ~ site, data = MY_pH)

# Welch's two-sample t-test
t.test(pH ~ site, data = MY_pH)





