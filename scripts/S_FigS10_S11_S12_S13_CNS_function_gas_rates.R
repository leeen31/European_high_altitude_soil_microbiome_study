###
### AlpineSoils - summer 23 - microbial function
### 


rm(list=ls(all=TRUE)) # removes everything


# Load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(ggordiplots)
library(ggpubr)
library(cowplot)


# read marker gene abundance table
mgene_tab_abs <- read.csv('data/Fig5_A_AlpineSoil23_contigs_functional_annotations_GreeningDB_abs.tsv', sep='\t', header = TRUE)

# read in meta data
meta_tab <- read.csv('data/Fig3_A_AlpineSoil23_contigs_metadata.tsv', sep='\t', header = TRUE)





### Identify sample columns (all ending in _d or _t)
sample_cols <- colnames(mgene_tab_abs)[grepl("_[dt]$", colnames(mgene_tab_abs))]

# pivot long
mgene_long <- mgene_tab_abs %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "Abundance"
  )

# add site info
mgene_long <- mgene_long %>%
  left_join(meta_tab %>% select(ID, soildepth, site), by = c("Sample" = "ID"))



### transforming variables
shapiro.test(meta_tab$vegetation)
meta_tab$vegetation_trans <- log1p(meta_tab$vegetation)

shapiro.test(meta_tab$GWC)
meta_tab$SWC_trans <- log1p(meta_tab$GWC)

shapiro.test(meta_tab$soil_temperature)
meta_tab$soil_temperature_scaled <- log1p(meta_tab$soil_temperature)

depth <- factor(meta_tab$soildepth)




### Supplementary figure S10 - CAZy class heatmap  ##########################################################

library(mgcv)
library(gratia)
library(ggplot2)
library(patchwork)

# read functional annotations file
ann_tab <- read.csv('data/Fig5_BCD_AlpineSoil23_contigs_functional_annotations_emapper.tsv', sep='\t', header = T)

# read coverage table
abund_tab_abs <- read.csv('data/Fig3_A_AlpineSoils23_contigs_coverage_abs.tsv', sep='\t', header = T)



# Summarize CAZy in ann
ann_clean <- ann_tab %>%
  group_by(contig_id) %>%
  filter(!(is.na(CAZy))) %>% filter(CAZy != "-", CAZy != "") %>%
  summarise(CAZy = paste(unique(CAZy), collapse = "; "))  # combine multiple annotations

# Join using the shortened ID
cov_cazy <- abund_tab_abs %>%
  left_join(ann_clean, by = c("contig_id")) %>%
  filter(!(is.na(CAZy))) %>% filter(CAZy != "-", CAZy != "") 
# keep only contigs with CAZy annotation
head(cov_cazy)


# convert to long format
cov_long <- cov_cazy %>%
  pivot_longer(
    cols = -c(contig_id, CAZy),
    names_to = "ID",
    values_to = "abundance"
  ) %>%
  left_join(meta_tab %>% select(ID, site), by = "ID")

# summarize mean abundance
cazy_expanded <- cov_long %>%
  separate_rows(CAZy, sep = ";") %>%
  separate_rows(CAZy, sep = ",") %>%
  mutate(CAZy = str_trim(CAZy)) %>%
  mutate(CAZy_class = str_extract(CAZy, "^[A-Z]+"))



# Cazy per site
cazy_site_class <- cazy_expanded %>%
  filter(!is.na(CAZy_class)) %>%
  group_by(CAZy_class, site) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE)) %>%
  mutate(log_abundance = log10(total_abundance)) %>%
  ungroup()



# rename classes
cazy_site_class <- cazy_site_class %>%
  mutate(CAZy_class_full = case_when(
    CAZy_class == "GH"  ~ "Glycoside Hydrolases (GH)",
    CAZy_class == "GT"  ~ "Glycosyltransferases (GT)",
    CAZy_class == "CE"  ~ "Carbohydrate Esterases (CE)",
    CAZy_class == "PL"  ~ "Polysaccharide Lyases (PL)",
    CAZy_class == "AA"  ~ "Auxiliary Activities (AA)",
    CAZy_class == "CBM" ~ "Carbohydrate-Binding Modules(CBM)",
    TRUE ~ CAZy_class  # fallback if something else appears
  ))


# order sites by vegetation
cazy_site_class$site <- factor(cazy_site_class$site, 
                               levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN" ,"MY2", "SF", "BN"))

# 4. Heatmap
heatmap_plot1 <- ggplot(cazy_site_class, aes(x = site, y = CAZy_class_full, fill = log_abundance)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "viridis") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12)
  ) +
  labs(x = "Site", y = "CAZy classes", fill = "Average cell number\ng⁻¹(log10)")
heatmap_plot1




### C degradation

# prep substrate data
library(readxl)

# read in reference table
substrate_cazy  <- read_excel("data/Fig5_CD_AlpSoils23_CAZy_Cdeg_substrate_reference.xlsx")


# using cazy_expand from above
cazy_long <- cazy_expanded %>%
  separate_rows(CAZy, sep = ",") %>%   # turns “CBM48,GH13” into two rows
  mutate(CAZy = trimws(CAZy))          # clean whitespace

# join with substrates
cazy_with_substrate <- cazy_long %>%
  left_join(substrate_cazy, 
            by = c("CAZy"))

cazy_plant_relevant <- cazy_with_substrate %>%
  filter(!is.na(substrate))


# summarize by site and substrate
substrate_abundance <- cazy_plant_relevant %>%
  group_by(site, CAZy, substrate) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE)) %>%
  mutate(log_abundance = log10(total_abundance)) %>%
  ungroup()




# order sites by vegetation
substrate_abundance$site <- factor(substrate_abundance$site, 
                                   levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN" ,"MY2", "SF", "BN"))

# Plot heatmap (now correctly clustered by site)
heatmap_plot2 <- ggplot(substrate_abundance, aes(x = site, y = substrate, fill = log_abundance)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "viridis") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12)
  ) +
  labs(x = "Site", y = "CAZy substrates", fill = "Average cell number\ng⁻¹(log10)")
heatmap_plot2



### exporting final Supplementary figure S10 ---------------

FigS10_final <- ggarrange(
  heatmap_plot1, heatmap_plot2,
  ncol = 1,
  labels = c("A", "B"), 
  heights = c(1, 1)
)

FigS10_final

ggsave("figures/S_FigS10_final_AlpSoils23_microbial_function.png", plot = FigS10_final, height = 5, width = 8)










### Supplemental figures S11 - Heatmap individual genes --------------------------

# use mgene_tab_abs from above


### 1. Identify sample columns (all ending in _d or _t)
sample_cols <- colnames(mgene_tab_abs)[grepl("_[dt]$", colnames(mgene_tab_abs))]

# pivot long
mgene_long <- mgene_tab_abs %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "Abundance"
  )

# add site info
mgene_long <- mgene_long %>%
  left_join(meta_tab %>% select(ID, soildepth, site), by = c("Sample" = "ID"))


# add 
fg_site <- mgene_long %>%
  filter(!is.na(site)) %>%   # remove samples without site info
  group_by(Element, Protein, site) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(log_abundance = log10(Abundance))

fg_site_clean <- fg_site %>%
  group_by(Element, Protein, site) %>%
  summarise(log_abundance = mean(log_abundance, na.rm = TRUE), .groups = "drop")


fg_matrix <- fg_site_clean %>%
  pivot_wider(
    names_from = site,
    values_from = log_abundance
  ) %>%
  distinct(Element, Protein, .keep_all = TRUE) %>%  # ensure one row per FG
  column_to_rownames("Protein") %>%
  as.matrix()

# remove too long facet names
fg_site_clean$Element_label <- fg_site_clean$Element
fg_site_clean$Element_label[fg_site_clean$Element %in% c("Metal_electron_transfer", "Other")] <- " "
fg_site_clean$Element_label <- recode(fg_site_clean$Element_label,
                                      "C cycling" = "C")

fg_site_clean$Element_label <- factor(fg_site_clean$Element_label, 
                                      levels = c("C", "CH4", "CO", "H2", "N", "S", " "))


# order sites by vegetation
fg_site_clean$site <- factor(fg_site_clean$site, 
                             levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))

# plot the heatmap
heatmap2 <- ggplot(fg_site_clean, aes(x = site, y = Protein, fill = log_abundance)) +
  geom_tile(color = "white") +
  #scale_fill_gradientn(colours = c("white", "#FDE725", "#35B779", "#31688E", "#440154"),) +
  scale_fill_viridis_c(option = "viridis") + 
  facet_grid(Element_label  ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(y) gsub("_", " ", y)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    strip.text.y = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  ) + 
  labs(x = "Site", y = "Functional marker genes", fill = "Average cell\n number g⁻¹\n (log10)")
heatmap2


### make vegetation gradient bar
n_sites <- length(levels(fg_site_clean$site))

library(scales)
n_sites <- length(levels(fg_site_clean$site))
p_triangle <- ggplot() +
  annotate("polygon",
           x = c(1, n_sites, n_sites),
           y = c(0.5, 1, 0),
           fill = scales::alpha("forestgreen", 0.5),  # transparency
           color = NA) +                      # no border (cleaner)
  
  annotate("text",
           x = mean(c(1, n_sites)),           # centered
           y = 0.5,
           label = "Vegetation extent",
           color = "black",
           size = 4,
           fontface = "bold") +
  
  scale_x_continuous(limits = c(1, n_sites), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme_void()



FigS11 <- plot_grid(
  heatmap2,
  p_triangle,
  ncol = 1,
  align = "v",
  axis = "lr",   # align left & right margins
  rel_heights = c(1, 0.04)
)
FigS11

ggsave("figures/S_FigS11_final_AlpSoils23_microbial_function.png", plot = FigS11, height = 10, width = 10)







### Supplemental figure S12 - Heatmap of marker genes - relative abundance (contig based) -----------------------------

# read marker gene abundance table
mgene_tab_rel <- read.csv('data/S_FigS11_AlpineSoil23_contigs_functional_annotations_GreeningDB_rel.tsv', sep='\t', header = TRUE)


### Heat map of Marker gene abundance

### 1. Identify sample columns (all ending in _d or _t)
sample_cols <- colnames(mgene_tab_rel)[grepl("_[dt]$", colnames(mgene_tab_rel))]

# pivot long
mgene_long_rel <- mgene_tab_rel %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "Abundance"
  )

# add site info
mgene_long_rel <- mgene_long_rel %>%
  left_join(meta_tab %>% select(ID, soildepth, site), by = c("Sample" = "ID"))

# add 
fg_site <- mgene_long_rel %>%
  filter(!is.na(site)) %>%   # remove samples without site info
  group_by(Element, Functional_Group, site) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(log_abundance = log10(Abundance))

fg_site_clean <- fg_site %>%
  group_by(Element, Functional_Group, site) %>%
  summarise(log_abundance = mean(log_abundance, na.rm = TRUE), .groups = "drop")


fg_matrix <- fg_site_clean %>%
  pivot_wider(
    names_from = site,
    values_from = log_abundance
  ) %>%
  distinct(Element, Functional_Group, .keep_all = TRUE) %>%  # ensure one row per FG
  column_to_rownames("Functional_Group") %>%
  as.matrix()

# remove too long facet names
fg_site_clean$Element_label <- fg_site_clean$Element
fg_site_clean$Element_label[fg_site_clean$Element %in% c("Metal_electron_transfer", "Other")] <- " "
fg_site_clean$Element_label <- recode(fg_site_clean$Element_label,
                                      "C cycling" = "C")

fg_site_clean$Element_label <- factor(fg_site_clean$Element_label, 
                                      levels = c("C", "CH4", "CO", "H2", "N", "S", " "))


# order sites by vegetation
fg_site_clean$site <- factor(fg_site_clean$site, 
                             levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))

# plot the heatmap
heatmap3 <- ggplot(fg_site_clean, aes(x = site, y = Functional_Group, fill = log_abundance)) +
  geom_tile(color = "white") +
  #scale_fill_gradientn(colours = c("white", "#FDE725", "#35B779", "#31688E", "#440154"),) +
  scale_fill_viridis_c(option = "viridis") + 
  facet_grid(Element_label  ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(y) gsub("_", " ", y)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    strip.text.y = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  ) + 
  labs(x = "Site", y = "Functional marker genes", fill = "Average cell\n number g⁻¹\n (log10)")

FigS12 <- plot_grid(
  heatmap3,
  p_triangle,
  ncol = 1,
  align = "v",
  axis = "lr",   # align left & right margins
  rel_heights = c(1, 0.04)
)
FigS12


ggsave("figures/S_FigS12_final_AlpSoils23_microbial_function_relabund.png", plot = FigS12, height = 6, width = 9)




### Supplemental figures S13 - GAM on N and S genes ---------------------------------------------

# marker genes based
# C, N, S, CH4, CO, H2 cyling






# GAMMs of marker gene abundance

library(mgcv)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gratia)  # for draw()
library(RColorBrewer)

## plot with stats using function

# Site levels and colors
site_levels <- c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN")
site_colors <- RColorBrewer::brewer.pal(10, "Paired")


# x-axis labels
x_labels <- c(
  vegetation_trans = "Vegetation cover",
  SWC_trans = "Soil water content",
  soil_temperature_scaled = "Soil temperature"
)


# Plotting function with stats
plot_gam_partial_effects <- function(formula_model, df,
                                     model_name,
                                     site_var = "site",
                                     site_levels = NULL,
                                     site_colors = NULL,
                                     x_labels = NULL,
                                     n_points = 200,
                                     alpha_points = 0.75) {
  
  df <- na.omit(df)
  
  gam_model <- gam(
    as.formula(formula_model),
    data = df,
    family = nb()
  )
  
  gam_sum <- summary(gam_model)
  sm <- gam_sum$s.table
  dev_exp <- round(gam_sum$dev.expl * 100, 1)
  
  print(summary(gam_model))
  
  print(appraise(gam_model))  # Check assumptions
  print(concurvity(gam_model, full = TRUE))# Check concurvity
  
  
  smooth_terms <- rownames(sm)
  smooth_terms <- smooth_terms[!grepl(site_var, smooth_terms)]
  
  vars <- gsub("^s\\(|,.*|\\)$", "", smooth_terms)
  vars <- vars[sapply(vars, function(v) is.numeric(df[[v]]))]
  
  res <- residuals(gam_model, type = "working")
  term_pred <- predict(gam_model, type = "terms", se.fit = TRUE)
  
  plots <- list()
  
  for (i in seq_along(vars)) {
    
    v <- vars[i]
    term_name <- smooth_terms[i]
    
    smooth_contrib <- term_pred$fit[, term_name]
    partial_resid <- smooth_contrib + res
    
    df_points <- data.frame(
      x = df[[v]],
      y = partial_resid,
      site = df[[site_var]]
    )
    
    grid <- seq(min(df[[v]], na.rm = TRUE),
                max(df[[v]], na.rm = TRUE),
                length.out = n_points)
    
    newdata <- df[rep(1, n_points), ]
    newdata[[v]] <- grid
    
    pred_grid <- predict(gam_model, newdata = newdata, type = "terms", se.fit = TRUE)
    
    df_smooth <- data.frame(
      x = grid,
      fit = pred_grid$fit[, term_name],
      se = pred_grid$se.fit[, term_name]
    )
    
    df_smooth$upper <- df_smooth$fit + 2 * df_smooth$se
    df_smooth$lower <- df_smooth$fit - 2 * df_smooth$se
    
    df_points$site <- factor(df_points$site, levels = site_levels)
    
    x_label <- if (!is.null(x_labels) && v %in% names(x_labels)) {
      x_labels[[v]]
    } else v
    
    edf_val <- round(sm[term_name, "edf"], 2)
    p_val <- signif(sm[term_name, "p-value"], 3)
    
    y_label <- if (i == 1) "Partial effect" else NULL
    
    p <- ggplot(df_points, aes(x = x, y = y, colour = site)) +
      geom_point(alpha = alpha_points) +
      geom_ribbon(data = df_smooth,
                  aes(x = x, ymin = lower, ymax = upper),
                  inherit.aes = FALSE,
                  fill = "grey70", alpha = 0.3) +
      geom_line(data = df_smooth,
                aes(x = x, y = fit),
                inherit.aes = FALSE,
                colour = "black", linewidth = 1.2) +
      annotate("text",
               x = -Inf, y = Inf,
               label = paste0("edf = ", edf_val, "\nP = ", p_val),
               hjust = -0.1, vjust = 1.2,
               size = 3.5) +
      scale_colour_manual(values = site_colors) +
      theme_classic() +
      labs(x = x_label, y = y_label, color = "Site") +
      theme(legend.position = "none")
    
    plots[[v]] <- p
  }
  
  combined <- ggarrange(
    plotlist = plots,
    ncol = 3,
    nrow = 1,
    legend = "none"
  )
  
  title <- paste0("Deviance explained = ", dev_exp, "%")
  
  combined <- annotate_figure(
    combined,
    top = text_grob(title, face = "bold", size = 12)
  )
  
  list(
    plot = combined,
    plots = plots,
    model = gam_model
  )
}


form = "sum_abundance ~ soildepth + s(vegetation_trans, k = 3) + s(SWC_trans, k = 3) + s(soil_temperature_scaled, k = 3) + s(site, bs = 're')"




# prep GAM data for markergenes

# sum per Element
markerg_gam <- mgene_long %>%
  group_by(Element, Sample) %>%
  summarise(sum_abundance = sum(Abundance, na.rm = TRUE)) %>%
  ungroup()


# add metadata to CAZy df
markerg_gam <- markerg_gam %>%
  left_join(meta_tab %>% select(ID, site, soildepth, vegetation_trans, SWC_trans, soil_temperature_scaled), by = c("Sample" = "ID"))





# Run GAMM for element N

N_gam_df <- markerg_gam %>% filter(Element == 'N')
N_gam_df$site <- as.factor(N_gam_df$site)

# use functions from above to plot with stats

N_gam_results <- plot_gam_partial_effects(
  form,
  df = N_gam_df,
  site_var = "site",
  site_levels = site_levels,
  site_colors = site_colors,
  x_labels = x_labels,
  alpha_points = 0.75 
)

# inspect plot
combined_plot_N <- N_gam_results$plot
combined_plot_N

# Vertical y-axis label
FigS13A <- ggdraw() +
  draw_label("N cycling", x = 0.02, y = 0.5, angle = 90,
             vjust = 0.5, fontface = "bold") +
  draw_plot(combined_plot_N, x = 0.05, y = 0, width = 0.95, height = 1)
FigS13A


# Run GAMM for element S

S_gam_df <- markerg_gam %>% filter(Element == 'S')
S_gam_df$site <- as.factor(S_gam_df$site)

# use functions from above to plot with stats

S_gam_results <- plot_gam_partial_effects(
  form,
  df = S_gam_df,
  site_var = "site",
  site_levels = site_levels,
  site_colors = site_colors,
  x_labels = x_labels,
  alpha_points = 0.75 
)

# inspect plot
combined_plot_S <- S_gam_results$plot
combined_plot_S

# Vertical y-axis label
FigS13B <- ggdraw() +
  draw_label("S cycling", x = 0.02, y = 0.5, angle = 90,
             vjust = 0.5, fontface = "bold") +
  draw_plot(combined_plot_S, x = 0.05, y = 0, width = 0.95, height = 1)
FigS13B



# Run GAMM for element C

C_gam_df <- markerg_gam %>% filter(Element == 'C cycling')
C_gam_df$site <- as.factor(C_gam_df$site)

# use functions from above to plot with stats

C_gam_results <- plot_gam_partial_effects(
  form,
  df = C_gam_df,
  site_var = "site",
  site_levels = site_levels,
  site_colors = site_colors,
  x_labels = x_labels,
  alpha_points = 0.75 
)

# inspect plot
combined_plot_C <- C_gam_results$plot
combined_plot_C

# Vertical y-axis label
FigS13C <- ggdraw() +
  draw_label("C cycling", x = 0.02, y = 0.5, angle = 90,
             vjust = 0.5, fontface = "bold") +
  draw_plot(combined_plot_C, x = 0.05, y = 0, width = 0.95, height = 1)
FigS13C



# Run GAMM for element CH4

CH4_gam_df <- markerg_gam %>% filter(Element == 'CH4')
CH4_gam_df$site <- as.factor(CH4_gam_df$site)

# use functions from above to plot with stats

CH4_gam_results <- plot_gam_partial_effects(
  form,
  df = CH4_gam_df,
  site_var = "site",
  site_levels = site_levels,
  site_colors = site_colors,
  x_labels = x_labels,
  alpha_points = 0.75 
)

# inspect plot
combined_plot_CH4 <- CH4_gam_results$plot
combined_plot_CH4

# Vertical y-axis label
FigS13D <- ggdraw() +
  draw_label("CH4 cycling", x = 0.02, y = 0.5, angle = 90,
             vjust = 0.5, fontface = "bold") +
  draw_plot(combined_plot_CH4, x = 0.05, y = 0, width = 0.95, height = 1)
FigS13D



# Run GAMM for element CO

CO_gam_df <- markerg_gam %>% filter(Element == 'CO')
CO_gam_df$site <- as.factor(CH4_gam_df$site)

# use functions from above to plot with stats

CO_gam_results <- plot_gam_partial_effects(
  form,
  df = CO_gam_df,
  site_var = "site",
  site_levels = site_levels,
  site_colors = site_colors,
  x_labels = x_labels,
  alpha_points = 0.75 
)

# inspect plot
combined_plot_CO <- CO_gam_results$plot
combined_plot_CO

# Vertical y-axis label
FigS13E <- ggdraw() +
  draw_label("CO cycling", x = 0.02, y = 0.5, angle = 90,
             vjust = 0.5, fontface = "bold") +
  draw_plot(combined_plot_CO, x = 0.05, y = 0, width = 0.95, height = 1)
FigS13E


# Run GAMM for element H2

H2_gam_df <- markerg_gam %>% filter(Element == 'H2')
H2_gam_df$site <- as.factor(H2_gam_df$site)

# use functions from above to plot with stats

H2_gam_results <- plot_gam_partial_effects(
  form,
  df = H2_gam_df,
  site_var = "site",
  site_levels = site_levels,
  site_colors = site_colors,
  x_labels = x_labels,
  alpha_points = 0.75 
)

# inspect plot
combined_plot_H2 <- H2_gam_results$plot
combined_plot_H2

# Vertical y-axis label
FigS13F <- ggdraw() +
  draw_label("H2 cycling", x = 0.02, y = 0.5, angle = 90,
             vjust = 0.5, fontface = "bold") +
  draw_plot(combined_plot_H2, x = 0.05, y = 0, width = 0.95, height = 1)
FigS13F


### exporting final Supplementary figure S13 ---------------

FigS13_final <- ggarrange(
  FigS13A, FigS13B,FigS13C, FigS13D, FigS13E, FigS13F,
  ncol = 1,
  labels = c("A", "B", "C", "D", "E", "F")
)

FigS13_final

ggsave("figures/S_FigS13_final_AlpSoils23_GAM_markergenes_NSCs.png", plot = FigS13_final, height = 17, width = 8)










