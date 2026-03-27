###
### AlpineSoils - summer 23 - Correlation of environmental factors + environmental PCA
###

rm(list=ls(all=TRUE)) # removes everything


### load data ####
data <- read.csv("data/Fig2_AB_AlpineSoil23_metadata.tsv", header= TRUE, sep= "\t")


# remove unwanted columns
data <- data[, !colnames(data) %in% c( "Shannon", "proc_euc", "OM")]






### Fig 1 A) correlation matrix --------------------------------------------------------
#https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

# load libraries
library(corrplot)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tools)
library(psych)


# Keep only numeric columns
numeric_data <- data[, sapply(data, is.numeric)]
cleaned_data <- na.omit(numeric_data)

# Remove unwanted columns
cleaned_data <- cleaned_data[, !colnames(cleaned_data) %in% c("OM", "proc_euc", "shannon")]

# Compute correlation matrix and p-values
cor_test_results <- psych::corr.test(cleaned_data, method = "spearman", adjust = "BH") # BH FDR correction
cor_matrix <- cor_test_results$r
p_matrix   <- cor_test_results$p

# Clean variable names
rownames(cor_matrix) <- toTitleCase(rownames(cor_matrix))
colnames(cor_matrix) <- toTitleCase(colnames(cor_matrix))
rownames(p_matrix) <- rownames(cor_matrix)
colnames(p_matrix) <- colnames(cor_matrix)

# Melt to long format 
cor_long <- melt(cor_matrix, varnames = c("Var1", "Var2"), value.name = "r")
p_long   <- melt(p_matrix, varnames = c("Var1", "Var2"), value.name = "p")

# Merge correlation + p-values
cor_long <- left_join(cor_long, p_long, by = c("Var1", "Var2"))

# Keep only upper triangle
cor_long <- cor_long[upper.tri(cor_matrix, diag = FALSE), ]

# Keep only significant correlations
sig_cor <- cor_long %>% filter(p < 0.05)

# Preserve order
sig_cor$Var1 <- factor(sig_cor$Var1, levels = colnames(cor_matrix))
sig_cor$Var2 <- factor(sig_cor$Var2, levels = rev(rownames(cor_matrix)))

# rename labels
sig_cor <- sig_cor %>%
  mutate(Var1 = recode(Var1,
                        "Soildepth" = "Soil depth",
                        "Vegetation" = "Vegetation cover",
                        "Soil_temperature" = "Soil temperature",
                        "GWC" = "SWC", 
                        "C_N" = "C:N ratio", 
                        "Soildepth" = "Soil depth", 
                        "mean_cells_gFW" = "Cell number",
                        "CO2_umol_g_h" = "CO₂ rate (µmol g⁻¹ h⁻¹)", 
                        "CH4_umol_g_h" = "CH₄ rate (µmol g⁻¹ h⁻¹", 
                        "NO3" = "NO₃", 
                        "NH4" = "NH₄"
                      )) %>%
  mutate(Var2 = recode(Var2,
                       "Soildepth" = "Soil depth",
                       "Vegetation" = "Vegetation cover",
                       "Soil_temperature" = "Soil temperature",
                       "GWC" = "SWC", 
                       "C_N" = "C:N ratio", 
                       "Soildepth" = "Soil depth", 
                       "mean_cells_gFW" = "Cell number",
                       "CO2_umol_g_h" = "CO₂ rate (µmol g⁻¹ h⁻¹)", 
                       "CH4_umol_g_h" = "CH₄ rate (µmol g⁻¹ h⁻¹)",
                       "NO3" = "NO₃", 
                       "NH4" = "NH₄"
  ))



# Plot only significant ones
Fig1A <- ggplot(sig_cor, aes(x = Var1, y = Var2, size = abs(r), color = r)) +
  geom_point() +
  coord_equal() +
  scale_color_distiller(palette = "PiYG", limits = c(-1, 1)) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(sig_cor$Var2))) +
  theme_minimal(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0.5,   face = "bold"),
    axis.text.x.bottom = element_blank(),
    axis.text = element_text(size = 11,  face = "bold"),
    legend.position = c(1, 0),  # adjust coordinates
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", colour = "white"),
    legend.box.background = element_blank(),
    legend.box = "horizontal",
    legend.spacing.y = unit(0.2, "cm")
  ) +
  labs(color = "Spearman ρ",
       size  = "|ρ|")
Fig1A





#### Fig 1 B) PCA -----------------------------------------------------------------


# load libraries

library(ggordiplots)
library(ggrepel)
library(vegan)
library(bestNormalize)


# load data
soil_data <- data


# remove unwanted columns
clean_data <- soil_data[, !colnames(soil_data) %in% c( "CO2_umol_g_h", "CH4_umol_g_h", "proc_euc", "OM")]


# remove all rows with NA
clean_data <- na.omit(clean_data)




### 1. transform variables

# check normal distribution
shapiro.test(clean_data$soildepth) #***
shapiro.test(clean_data$soil_temperature) #***
shapiro.test(clean_data$TN) #***
shapiro.test(clean_data$TC) #***
shapiro.test(clean_data$TOC) #***
shapiro.test(clean_data$TH) #***
shapiro.test(clean_data$altitude) #***
shapiro.test(clean_data$vegetation) #***
shapiro.test(clean_data$pH) #*
shapiro.test(clean_data$GWC) #***
shapiro.test(clean_data$C_N) #***
shapiro.test(clean_data$Ca) #***
shapiro.test(clean_data$K) #***
shapiro.test(clean_data$Mg) #***
shapiro.test(clean_data$Cl) #***
shapiro.test(clean_data$Na) #***
shapiro.test(clean_data$NH4) #***
shapiro.test(clean_data$NO3) #***
shapiro.test(clean_data$mean_cells_gFW) #***
# no normal distribution - transformation necessary 



# find transformation
numeric_cols <- clean_data %>% select(where(is.numeric))
results_list <- list()

for(col_name in names(numeric_cols)) {
  x <- numeric_cols[[col_name]]
  if(all(is.na(x))) next
  
  bn <- bestNormalize(x, allow_lambert_s = TRUE)
  
  # Extract the transformation name as a character
  best_method <- class(bn$chosen_transform)[1]
  
  results_list[[col_name]] <- data.frame(
    Variable = col_name,
    Best_Transformation = best_method,  # now a string
    Original_Min = min(x, na.rm = TRUE),
    Original_25 = quantile(x, 0.25, na.rm = TRUE),
    Original_Median = median(x, na.rm = TRUE),
    Original_75 = quantile(x, 0.75, na.rm = TRUE),
    Original_Max = max(x, na.rm = TRUE)
  )
}

normalization_results <- bind_rows(results_list)
print(normalization_results)



# implement transformations
soil_meta_trans <- soil_data

# Loop over each variable in the normalization results
for(i in seq_len(nrow(normalization_results))) {
  
  var_name <- normalization_results$Variable[i]
  method <- normalization_results$Best_Transformation[i]
  
  x <- soil_data[[var_name]]
  if(all(is.na(x))) next  # skip if variable is all NA
  
  # Fit bestNormalize allowing all methods
  bn <- bestNormalize(x, allow_lambert_s = TRUE)
  
  # Transform variable using the best method indicated in the table
  soil_meta_trans[[paste0(var_name, "_trans")]] <- predict(bn)
}


# remove untransformed and NAs
soil_meta_trans <- soil_meta_trans[,-c(5:25)]
soil_meta_trans <- na.omit(soil_meta_trans)






### 2. compute PCA


# make numeric
env_vars <- soil_meta_trans %>% select(where(is.numeric))

# Step 4: PCA
pca_result <- prcomp(env_vars, center = TRUE, scale. = TRUE)

# Step 5: Metadata for plotting
pca_df <- as.data.frame(pca_result$x)
pca_df$site <- soil_meta_trans$site
pca_df$depth <- soil_meta_trans$depth



# remove _trans from drivers
env_vars <- env_vars %>% rename_with(~ gsub("_trans$", "", .x))








### 3. Environmental driver analysis

envfit_result <- envfit(pca_result, env_vars, permutations = 999)

# Step 7: Extract significant arrows only

all_vectors <- envfit_result$vectors$arrows
all_pvals <- envfit_result$vectors$pvals

# Adjust p-values for multiple testing (FDR / Benjamini-Hochberg)
all_pvals_adj <- p.adjust(all_pvals, method = "BH")

# Filter significant variables using adjusted p-value < 0.05
sig_names <- names(all_pvals_adj)[all_pvals_adj < 0.05]

# Build dataframe with PC1, PC2 scores, adjusted p-values, and variable names
sig_env_list <- data.frame(
  label = sig_names,
  PC1 = all_vectors[sig_names, "PC1"],
  PC2 = all_vectors[sig_names, "PC2"],
  p_adj = all_pvals_adj[sig_names]
)

# Show the list
sig_env_list <- sig_env_list %>% arrange(p_adj)
print(sig_env_list)


sig_arrows <- all_vectors[sig_names, , drop=FALSE]





# Range of PCA scores (the space arrows must fit inside)
# Extract r values (correlation of each variable with the ordination)
sig_r <- envfit_result$vectors$r[sig_names]

# Scale arrows using r values
sig_arrows_scaled <- sig_arrows * sig_r

# Step 2: Scale to fit plot visually
pca_range <- apply(pca_df[, c("PC1","PC2")], 2, diff)

# further scale to fit inside the PCA plot
max_fraction <- 0.75  # 50% of the smallest PCA range
arrow_scaling <- max_fraction * min(pca_range) / max(sqrt(rowSums(sig_arrows_scaled^2)))
sig_arrows_scaled <- sig_arrows_scaled * arrow_scaling





# Step 9: Build sig_env dataframe
sig_env <- as.data.frame(sig_arrows_scaled)
colnames(sig_env) <- c("PC1","PC2")
sig_env$x <- 0
sig_env$y <- 0
sig_env$label <- rownames(sig_arrows_scaled)




# invert  axis to fit microbial NMDS
pca_df$PC1 <- -pca_df$PC1
#pca_df$PC2 <- -pca_df$PC2
#sig_env$PC1<- sig_env$PC1 * (-1)  # and drivers
sig_env$PC2<- sig_env$PC2 * (-1)  # and drivers



# rename labels
sig_env <- sig_env %>%
  mutate(label = recode(label,
                        "altitude" = "Altitude",
                        "vegetation" = "Vegetation cover",
                        "soil_temperature" = "Soil temperature",
                        "GWC" = "SWC", 
                        "C_N" = "C:N ratio", 
                        "soildepth" = "Soil depth", 
                        "mean_cells_gFW" = "Cell number", 
                        "NO3" = "NO₃", 
                        "NH4" = "NH₄"))


# Custom colors
custom_colors <- RColorBrewer::brewer.pal(10, "Paired")

# order sites by vegetation
pca_df$site = factor(pca_df$site, levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))

# move lables to look nice
sig_env$PC1[sig_env$label == 'Vegetation cover'] = sig_env$PC1[sig_env$label == 'Vegetation cover'] + 0.5
sig_env$PC1[sig_env$label == 'Soil temperature'] = sig_env$PC1[sig_env$label == 'Soil temperature'] - 1.3
sig_env$PC1[sig_env$label == 'SWC'] = sig_env$PC1[sig_env$label == 'SWC'] + 0.5
sig_env$PC1[sig_env$label == 'C:N ratio'] = sig_env$PC1[sig_env$label == 'C:N ratio'] - 1
sig_env$PC2[sig_env$label == 'C:N ratio'] = sig_env$PC2[sig_env$label == 'C:N ratio'] - 0.5
sig_env$PC1[sig_env$label == 'Cl'] = sig_env$PC1[sig_env$label == 'Cl'] + 0.2
sig_env$PC1[sig_env$label == 'Mg'] = sig_env$PC1[sig_env$label == 'Mg'] - 0.4
sig_env$PC2[sig_env$label == 'Ca'] = sig_env$PC2[sig_env$label == 'Ca'] + 0.2


# Step 10: Plot PCA with only significant arrows
Fig1B <- ggplot(pca_df, aes(x = PC1, y = PC2, color = site, shape = depth)) + # 
  geom_point(size = 4, alpha = 0.75) +
  geom_segment(data = sig_env,
               aes(x = x, y = y, xend = PC1 * 1.4, yend = PC2 * 1.4),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text(data = sig_env,aes(x = PC1 * 1.6, y = PC2 * 1.6, label = label),
                 size = 4, color = "black", inherit.aes = FALSE) +
  #geom_text_repel(data = sig_env,aes(x = PC1 * 1.1, y = PC2 * 1.1, label = label),
  #                size = 4, color = "black", inherit.aes = FALSE) +
  scale_color_manual(values = custom_colors) +
  labs(x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100,1), "%)"),
       color = "Site") +
  labs(color = "Site", shape  = "Depth") +
  scale_shape_manual(
    values = c(17, 16),  # your shapes
    labels = c("Topsoil", "Lower soil layer")
  ) +
  theme_bw()
Fig1B


### 4. PERMANOVA using site and depth


# scale exactly as PCA did
env_scaled <- scale(env_vars)

# PERMANOVA
adonis2(env_scaled ~ site, data = soil_meta_trans, method = "euclidean") # site 0.001 ***
adonis2(env_scaled ~ depth, data = soil_meta_trans, method = "euclidean") # depths 0.001 **
adonis2(env_scaled ~ site/depth, data = soil_meta_trans, method = "euclidean") # site/depths 0.001 ***


# beta dispersion
bd <- betadisper(dist(env_scaled),soil_meta_trans$site)
anova(bd) # p 0.8934
bd <- betadisper(dist(env_scaled),soil_meta_trans$depth)
anova(bd) # p 0.8605







### exporting final figure 2 ---------------
library(ggpubr)

Fig2_final <- ggarrange(Fig1A, Fig1B, ncol = 2, nrow = 1, labels = c('A', 'B'))
Fig2_final

ggsave("figures/Fig2_final_AlpSoils243_env_corr_PCA.png", plot = Fig2_final, height = 5.5, width = 13)
ggsave("figures/Fig2_final_AlpSoils243_env_corr_PCA.svg", plot = Fig2_final, height = 5.5, width = 13)













### Supplemental Figure S1A and B -----------------------------
# PCA of topsoil and lower soil layer

# use normalization from above





### 2. compute PCA for top soil 


# filter for depths
soil_meta_trans_top <- soil_meta_trans %>% filter(depth == "Top")

# remove unwanted columns
soil_meta_trans_top <- soil_meta_trans_top[, !colnames(soil_meta_trans_top) %in% c( "depth", "soildepth_trans")]


# make numeric
env_vars <- soil_meta_trans_top %>% select(where(is.numeric))

# Step 4: PCA
pca_result <- prcomp(env_vars, center = TRUE, scale. = TRUE)

# Step 5: Metadata for plotting
pca_df <- as.data.frame(pca_result$x)
pca_df$site <- soil_meta_trans_top$site



# remove _trans from drivers
env_vars <- env_vars %>% rename_with(~ gsub("_trans$", "", .x))








### 3. Environmental driver analysis

envfit_result <- envfit(pca_result, env_vars, permutations = 999)

# Step 7: Extract significant arrows only

all_vectors <- envfit_result$vectors$arrows
all_pvals <- envfit_result$vectors$pvals

# Adjust p-values for multiple testing (FDR / Benjamini-Hochberg)
all_pvals_adj <- p.adjust(all_pvals, method = "BH")

# Filter significant variables using adjusted p-value < 0.05
sig_names <- names(all_pvals_adj)[all_pvals_adj < 0.05]

# Build dataframe with PC1, PC2 scores, adjusted p-values, and variable names
sig_env_list <- data.frame(
  label = sig_names,
  PC1 = all_vectors[sig_names, "PC1"],
  PC2 = all_vectors[sig_names, "PC2"],
  p_adj = all_pvals_adj[sig_names]
)

# Show the list
sig_env_list <- sig_env_list %>% arrange(p_adj)
print(sig_env_list)


sig_arrows <- all_vectors[sig_names, , drop=FALSE]





# Range of PCA scores (the space arrows must fit inside)
# Extract r values (correlation of each variable with the ordination)
sig_r <- envfit_result$vectors$r[sig_names]

# Scale arrows using r values
sig_arrows_scaled <- sig_arrows * sig_r

# Step 2: Scale to fit plot visually
pca_range <- apply(pca_df[, c("PC1","PC2")], 2, diff)

# further scale to fit inside the PCA plot
max_fraction <- 0.75  # 50% of the smallest PCA range
arrow_scaling <- max_fraction * min(pca_range) / max(sqrt(rowSums(sig_arrows_scaled^2)))
sig_arrows_scaled <- sig_arrows_scaled * arrow_scaling





# Step 9: Build sig_env dataframe
sig_env <- as.data.frame(sig_arrows_scaled)
colnames(sig_env) <- c("PC1","PC2")
sig_env$x <- 0
sig_env$y <- 0
sig_env$label <- rownames(sig_arrows_scaled)




# invert  axis to fit microbial NMDS
pca_df$PC1 <- -pca_df$PC1
pca_df$PC2 <- pca_df$PC2
#sig_env$PC1<- sig_env$PC1 * (-1)  # and drivers
sig_env$PC2<- sig_env$PC2 * (-1)  # and drivers



# rename labels
sig_env <- sig_env %>%
  mutate(label = recode(label,
                        "altitude" = "Altitude",
                        "vegetation" = "Vegetation cover",
                        "soil_temperature" = "Soil temperature",
                        "GWC" = "SWC", 
                        "C_N" = "C:N ratio", 
                        "soildepth" = "Soil depth", 
                        "mean_cells_gFW" = "Cell number", 
                        "NO3" = "NO₃", 
                        "NH4" = "NH₄"
                        ))


# Custom colors
custom_colors <- RColorBrewer::brewer.pal(10, "Paired")

# order sites by vegetation
pca_df$site = factor(pca_df$site, levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))

# move labels to look nice


# Step 10: Plot PCA with only significant arrows
FigS1_top <- ggplot(pca_df, aes(x = PC1, y = PC2, color = site)) + # 
  geom_point(size = 4, alpha = 0.75) +
  geom_segment(data = sig_env,
               aes(x = x, y = y, xend = PC1 * 1.4, yend = PC2 * 1.4),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text(data = sig_env,aes(x = PC1 * 1.6, y = PC2 * 1.6, label = label),
            size = 4, color = "black", inherit.aes = FALSE) +
  scale_color_manual(values = custom_colors) +
  labs(x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100,1), "%)"),
       title = "Topsoil layer", 
       color = "Site") +
  theme_bw()
FigS1_top




### 4. PERMANOVA using site

# scale exactly as PCA did
env_scaled <- scale(env_vars)

# PERMANOVA
adonis2(env_scaled ~ site, data = soil_meta_trans_top, method = "euclidean") # site 0.001 ***








### 2. compute PCA for lower soil 


# filter for depths
soil_meta_trans_deep <- soil_meta_trans %>% filter(depth == "Deep")

# remove unwanted columns
soil_meta_trans_deep <- soil_meta_trans_deep[, !colnames(soil_meta_trans_deep) %in% c( "depth", "soildepth_trans")]


# make numeric
env_vars <- soil_meta_trans_deep %>% select(where(is.numeric))

# Step 4: PCA
pca_result <- prcomp(env_vars, center = TRUE, scale. = TRUE)

# Step 5: Metadata for plotting
pca_df <- as.data.frame(pca_result$x)
pca_df$site <- soil_meta_trans_deep$site



# remove _trans from drivers
env_vars <- env_vars %>% rename_with(~ gsub("_trans$", "", .x))








### 3. Environmental driver analysis

envfit_result <- envfit(pca_result, env_vars, permutations = 999)

# Step 7: Extract significant arrows only

all_vectors <- envfit_result$vectors$arrows
all_pvals <- envfit_result$vectors$pvals

# Adjust p-values for multiple testing (FDR / Benjamini-Hochberg)
all_pvals_adj <- p.adjust(all_pvals, method = "BH")

# Filter significant variables using adjusted p-value < 0.05
sig_names <- names(all_pvals_adj)[all_pvals_adj < 0.05]

# Build dataframe with PC1, PC2 scores, adjusted p-values, and variable names
sig_env_list <- data.frame(
  label = sig_names,
  PC1 = all_vectors[sig_names, "PC1"],
  PC2 = all_vectors[sig_names, "PC2"],
  p_adj = all_pvals_adj[sig_names]
)

# Show the list
sig_env_list <- sig_env_list %>% arrange(p_adj)
print(sig_env_list)


sig_arrows <- all_vectors[sig_names, , drop=FALSE]





# Range of PCA scores (the space arrows must fit inside)
# Extract r values (correlation of each variable with the ordination)
sig_r <- envfit_result$vectors$r[sig_names]

# Scale arrows using r values
sig_arrows_scaled <- sig_arrows * sig_r

# Step 2: Scale to fit plot visually
pca_range <- apply(pca_df[, c("PC1","PC2")], 2, diff)

# further scale to fit inside the PCA plot
max_fraction <- 0.75  # 50% of the smallest PCA range
arrow_scaling <- max_fraction * min(pca_range) / max(sqrt(rowSums(sig_arrows_scaled^2)))
sig_arrows_scaled <- sig_arrows_scaled * arrow_scaling





# Step 9: Build sig_env dataframe
sig_env <- as.data.frame(sig_arrows_scaled)
colnames(sig_env) <- c("PC1","PC2")
sig_env$x <- 0
sig_env$y <- 0
sig_env$label <- rownames(sig_arrows_scaled)




# invert  axis to fit microbial NMDS
pca_df$PC1 <- -pca_df$PC1
pca_df$PC2 <- pca_df$PC2
#sig_env$PC1<- sig_env$PC1 * (-1)  # and drivers
sig_env$PC2<- sig_env$PC2 * (-1)  # and drivers



# rename labels
sig_env <- sig_env %>%
  mutate(label = recode(label,
                        "altitude" = "Altitude",
                        "vegetation" = "Vegetation cover",
                        "soil_temperature" = "Soil temperature",
                        "GWC" = "SWC", 
                        "C_N" = "C:N ratio", 
                        "soildepth" = "Soil depth", 
                        "mean_cells_gFW" = "Cell number", 
                        "NO3" = "NO₃", 
                        "NH4" = "NH₄"
  ))


# Custom colors
custom_colors <- RColorBrewer::brewer.pal(10, "Paired")

# order sites by vegetation
pca_df$site = factor(pca_df$site, levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))

# move labels to look nice


# Step 10: Plot PCA with only significant arrows
FigS1_deep <- ggplot(pca_df, aes(x = PC1, y = PC2, color = site)) + # 
  geom_point(size = 4, alpha = 0.75) +
  geom_segment(data = sig_env,
               aes(x = x, y = y, xend = PC1 * 1.4, yend = PC2 * 1.4),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text(data = sig_env,aes(x = PC1 * 1.6, y = PC2 * 1.6, label = label),
            size = 4, color = "black", inherit.aes = FALSE) +
  scale_color_manual(values = custom_colors) +
  labs(x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100,1), "%)"),
       title = "Lower soil layer",
       color = "Site") +
  theme_bw()
FigS1_deep




### 4. PERMANOVA using site

# scale exactly as PCA did
env_scaled <- scale(env_vars)

# PERMANOVA
adonis2(env_scaled ~ site, data = soil_meta_trans_deep, method = "euclidean") # site 0.001 ***


### exporting final S figure S1AB ---------------
library(ggpubr)

S_FigS1_final <- ggarrange(FigS1_top, FigS1_deep, ncol = 2, nrow = 1, labels = c('A', 'B'), common.legend =  TRUE, legend = "right" )
S_FigS1_final

ggsave("figures/S_Figs1_final_AlpSoils243_PCAs.png", plot = S_FigS1_final, height = 5.5, width = 13)





