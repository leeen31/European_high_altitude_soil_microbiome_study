###
### AlpineSoils - summer 23 - microbial community, function, diversity and cell numbers
### 


rm(list=ls(all=TRUE)) # removes everything

# read in meta data
meta_tab <- read.csv('data/Fig3_A_AlpineSoil23_contigs_metadata.tsv', sep='\t', header = T)





### Fig 3 A) NMDS of contig-based community composition -------------------

# load libs
library(ggordiplots)
library(tibble)
library(reshape2)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(ggpubr)

### read in files from file prepocessing

# read coverage table
abund_tab_abs <- read.csv('data/Fig3_A_AlpineSoils23_contigs_coverage_abs.tsv', sep='\t', header = T)



### NMDS 


# Transpose: samples become rows, ASVs become columns
rownames(abund_tab_abs) = abund_tab_abs$contig_id
abund_tab_abs$contig_id = NULL


abund_tab_abs_t <- t(abund_tab_abs)
abund_mat <- as.matrix(abund_tab_abs_t)  # already numeric, samples x ASVs


# Hellinger transformation
abund_hell <- sqrt(abund_mat / rowSums(abund_mat))




# Distance matrix and NMDS
dist_mat <- vegdist(abund_hell, method = "bray")  # Bray-Curtis

nmds <- metaMDS(dist_mat)
stressplot(nmds)

# NMDS site scores
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$ID <- rownames(nmds_scores)


# Load metadata/environmental variables
env_vars <- meta_tab
env_vars$Shannon <- NULL
env_vars$proc_euc <- NULL

# Fix sample names
rownames(env_vars) <- env_vars$ID


# Keep only samples present in NMDS
env_vars <- env_vars[rownames(env_vars) %in% rownames(nmds_scores), ]

# Standardize numeric variables
env_num <- env_vars %>%
  select(where(is.numeric)) %>%
  scale() %>%
  as.data.frame()
rownames(env_num) <- rownames(env_vars)


# Environmental fitting (envfit)
envfit_res <- envfit(nmds, env_num, permutations = 999, na.rm = TRUE)



# envfit_res is your envfit result
env_vectors <- scores(envfit_res, display = "vectors")  # coordinates of arrows
env_pvals   <- envfit_res$vectors$pvals                # raw p-values
env_r2      <- envfit_res$vectors$r                    # R² values

# Combine into a data frame
env_summary <- data.frame(
  Variable = rownames(env_vectors),
  NMDS1 = env_vectors[, "NMDS1"],
  NMDS2 = env_vectors[, "NMDS2"],
  R2 = env_r2,
  P = env_pvals
)

# Add BH-adjusted p-values
env_summary$P_adj <- p.adjust(env_summary$P, method = "BH")

# Print all variables with adjusted p-values
print(env_summary[order(env_summary$P_adj), ])





#  Extract significant vectors 
sig_vectors <- envfit_res$vectors$arrows * sqrt(envfit_res$vectors$r)
sig_pvals   <- envfit_res$vectors$pvals
sig_adj     <- p.adjust(sig_pvals, method = "BH")

sig_names <- names(sig_adj)[sig_adj < 0.05]

sig_env <- as.data.frame(sig_vectors[sig_names, , drop = FALSE])
sig_env$Variable <- rownames(sig_env)



# Scale arrows to NMDS range 
arrow_scale <- 0.5 * min(apply(nmds_scores[,c("NMDS1","NMDS2")],2,diff)) /
  max(sqrt(rowSums(sig_env[,1:2]^2)))
sig_env$NMDS1 <- sig_env$NMDS1 * arrow_scale
sig_env$NMDS2 <- sig_env$NMDS2 * arrow_scale
sig_env$x <- 0
sig_env$y <- 0



# Combine NMDS and metadata for plotting
nmds_plot_df <- nmds_scores %>%
  left_join(env_vars, by = c("ID" = "ID"))

# Plot NMDS with arrows
custom_colors <- RColorBrewer::brewer.pal(10, "Paired")
stress_val <- round(nmds$stress, 2)


# order sites by vegetation
nmds_plot_df$site = factor(nmds_plot_df$site, levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))


# rename labels
sig_env <- sig_env %>%
  mutate(Variable = recode(Variable,
                        "altitude" = "Altitude",
                        "vegetation" = "Vegetation cover",
                        "C_N" = "C:N ratio", 
                        "mean_cells_gFW" = "Cell number", 
                        "CO2_umol_g_h" = "CO₂ rate (µmol g⁻¹ h⁻¹)", 
                        "CH4_umol_g_h" = "CH₄ rate (µmol g⁻¹ h⁻¹)"))


# move labels to look nice
sig_env$x[sig_env$label == 'Cell number'] = sig_env$x[sig_env$label == 'Cell number'] - 0.5
sig_env$y[sig_env$label == 'Altitude'] = sig_env$y[sig_env$label == 'Altitude'] - 0.5
sig_env$y[sig_env$label == 'CO₂ rate (µmol g⁻¹ h⁻¹)'] = sig_env$y[sig_env$label == 'CO₂ rate (µmol g⁻¹ h⁻¹)'] -2
sig_env$y[sig_env$label == 'pH'] = sig_env$y[sig_env$label == 'pH'] - 0.8
sig_env$x[sig_env$label == 'C:N ratio'] = sig_env$x[sig_env$label == 'C:N ratio'] - 0.2
sig_env$y[sig_env$label == 'Vegetation cover'] = sig_env$y[sig_env$label == 'Vegetation cover'] + 0.4



# Flip only arrows (both axes)
sig_env$NMDS1 <- -sig_env$NMDS1
sig_env$NMDS2 <- -sig_env$NMDS2

Fig3A <- ggplot(nmds_plot_df, aes(x = NMDS1, y = NMDS2, color = site, shape = depth)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_segment(data=sig_env, aes(x=x, y=y, xend=NMDS1, yend=NMDS2),
               arrow=arrow(length=unit(0.25,"cm")), color="black", inherit.aes = FALSE) +
  geom_text_repel(data=sig_env, aes(x=NMDS1*1.1, y=NMDS2*1.1, label=Variable),
                  size=4, color="black", inherit.aes = FALSE) +
  labs(x="NMDS1", y="NMDS2", color="Site", shape="Depth") +
  scale_color_manual(values = custom_colors) +
  annotate("text", x = Inf, y = -Inf, 
           label = paste("Stress =", stress_val), 
           hjust = 1.1, vjust = -0.5, size = 4) +
  labs(x = "NMDS1", y = "NMDS2", color = "Site", shape = "Depth") +
  scale_shape_manual(
    values = c(17, 16),  # your shapes
    labels = c("Topsoil", "Lower soil layer")) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14,
    ),
    plot.margin = margin(20, 10, 10, 10)
  ) +
  labs(title = "Taxonomic level") +
  theme_bw()
Fig3A





# PERMANOVA 
# Make sure samples match
meta_tab <- meta_tab %>% filter(ID %in% rownames(abund_hell))

# Reorder abund_hell to match metadata
abund_hell_sub <- abund_hell[meta_tab$ID, , drop = FALSE]

# Compute Bray-Curtis distance
dist_mat <- vegdist(abund_hell_sub, method = "bray")

# PERMANOVA: site effect
adonis2(dist_mat ~ site, data = meta_tab, permutations = 999) # site: 0.001 ***
adonis2(dist_mat ~ depth, data = meta_tab, permutations = 999) # depth: 0.729

# PERMANOVA: depth nested within site
adonis2(dist_mat ~ site / depth, data = meta_tab, permutations = 999, by = "margin") # Site/Depth  0.001 ***







### Fig 3 B) NMDS of contig-based functional composition (KO) ------------

library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tidyverse)



### read in files from file prepocessing

# use meta data from above

# read coverage table again (use newly imported table!)
abund_tab_abs <- read.csv('data/Fig3_A_AlpineSoils23_contigs_coverage_abs.tsv', sep='\t', header = T)


# read functional annotations file
ann_tab <- read.csv('data/Fig3_B_AlpineSoil23_contigs_functional_annotations_emapper.tsv', sep='\t', header = T)





### NMDS 

### Extract KEGG KOs from annotations
ann_tab <- ann_tab %>%
  mutate(KEGG_ko = str_extract(KEGG_ko, "K\\d{5}"))

# remove append from bin_id
ann_filtered <- ann_tab %>% mutate(contig_id = str_replace(gene_id, "_\\d+$", ""))


# Join abundance + metadata 
abund_long <- abund_tab_abs %>% 
  pivot_longer(-contig_id, names_to = "ID", values_to = "abundance")

# join annotations
abund_ann <- abund_long %>% 
  left_join(ann_filtered, by = c("contig_id" ))

# summarize KO abundances per sample
ko_wide <- abund_ann %>%
  filter(!is.na(KEGG_ko)) %>%              # only annotated KOs
  group_by(ID, KEGG_ko) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = KEGG_ko,
    values_from = abundance,
    values_fill = 0
  )


# Convert to matrix with sample IDs as rownames
ko_mat <- ko_wide %>%
  column_to_rownames("ID") %>%  # rownames = sample IDs
  as.matrix()                   # numeric matrix


# hellinger transform
abund_hell <- sqrt(ko_mat / rowSums(ko_mat))


# distance matrix and nmds
dist_mat <- vegdist(abund_hell, method = "bray")
nmds <- metaMDS(dist_mat)
stressplot(nmds)

nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$ID <- rownames(nmds_scores)


# metadata alignment
env_vars <- meta_tab
env_vars$Shannon <- NULL
env_vars$proc_euc <- NULL
rownames(env_vars) <- env_vars$ID

env_vars <- env_vars[rownames(env_vars) %in% rownames(nmds_scores), ]

env_num <- env_vars %>%
  select(where(is.numeric)) %>%
  scale() %>%
  as.data.frame()
rownames(env_num) <- rownames(env_vars)


# environmental variables
envfit_res <- envfit(nmds, env_num, permutations = 999, na.rm = TRUE)

env_vectors <- scores(envfit_res, display = "vectors")
env_pvals   <- envfit_res$vectors$pvals
env_r2      <- envfit_res$vectors$r

env_summary <- data.frame(
  Variable = rownames(env_vectors),
  NMDS1 = env_vectors[, "NMDS1"],
  NMDS2 = env_vectors[, "NMDS2"],
  R2 = env_r2,
  P = env_pvals
)

env_summary$P_adj <- p.adjust(env_summary$P, method = "BH")
print(env_summary[order(env_summary$P_adj), ])


# sifnificant drivers
sig_vectors <- envfit_res$vectors$arrows * sqrt(envfit_res$vectors$r)
sig_pvals   <- envfit_res$vectors$pvals
sig_adj     <- p.adjust(sig_pvals, method = "BH")

sig_names <- names(sig_adj)[sig_adj < 0.05]

sig_env <- as.data.frame(sig_vectors[sig_names, , drop = FALSE])
sig_env$Variable <- rownames(sig_env)

arrow_scale <- 0.5 * min(apply(nmds_scores[,c("NMDS1","NMDS2")],2,diff)) /
  max(sqrt(rowSums(sig_env[,1:2]^2)))

sig_env$NMDS1 <- sig_env$NMDS1 * arrow_scale
sig_env$NMDS2 <- sig_env$NMDS2 * arrow_scale
sig_env$x <- 0
sig_env$y <- 0




# plot nmds
nmds_plot_df <- nmds_scores %>% 
  left_join(env_vars, by = "ID")



# invert x axis
nmds_plot_df$NMDS1 <- -nmds_plot_df$NMDS1
nmds_plot_df$NMDS2 <- -nmds_plot_df$NMDS2

# Flip arrows
sig_env$NMDS1 <- -sig_env$NMDS1
sig_env$NMDS2 <- -sig_env$NMDS2




# order sites by vegetation
nmds_plot_df$site = factor(nmds_plot_df$site, levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))


# rename labels
sig_env <- sig_env %>%
  mutate(Variable = recode(Variable,
                           "vegetation" = "Vegetation cover",
                           "CO2_umol_g_h" = "CO₂ rate (µmol g⁻¹ h⁻¹)"))


# move labels to look nice
sig_env$x[sig_env$label == 'TOC'] = sig_env$x[sig_env$label == 'TOC'] + 1.2
sig_env$y[sig_env$label == 'pH'] = sig_env$y[sig_env$label == 'pH'] - 0.3
sig_env$x[sig_env$label == 'Vegetation cover'] = sig_env$x[sig_env$label == 'Vegetation cover'] + 1.2




custom_colors <- RColorBrewer::brewer.pal(10, "Paired")
stress_val <- round(nmds$stress, 2)

Fig3B <- ggplot(nmds_plot_df, aes(x = NMDS1, y = NMDS2, color = site, shape = depth)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_segment(data=sig_env, aes(x=x, y=y, xend=NMDS1, yend=NMDS2),
               arrow=arrow(length=unit(0.25,"cm")), color="black", inherit.aes = FALSE) +
  geom_text_repel(data=sig_env, aes(x=NMDS1*1.1, y=NMDS2*1.1, label=Variable),
                  size=4, color="black", inherit.aes = FALSE) +
  labs(x="NMDS1", y="NMDS2", color="Site", shape="Depth") +
  scale_color_manual(values = custom_colors) +
  annotate("text", x = Inf, y = -Inf, 
           label = paste("Stress =", stress_val), 
           hjust = 1.1, vjust = -0.5, size = 4) +
  scale_shape_manual(
    values = c(17, 16),  # your shapes
    labels = c( "Topsoil","Lower soil layer")) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14,
    ),
    plot.margin = margin(20, 10, 10, 10)
  ) +
  labs(title = "Functional level") +
  theme_bw()
Fig3B




### PERMANOVA
meta_tab <- meta_tab %>% filter(ID %in% rownames(abund_hell))

abund_hell_sub <- abund_hell[meta_tab$ID, , drop = FALSE]
dist_mat <- vegdist(abund_hell_sub, method = "bray")

adonis2(dist_mat ~ site, data = meta_tab, permutations = 999) # site: 0.001 ***
adonis2(dist_mat ~ depth, data = meta_tab, permutations = 999) # depth: 0.569
adonis2(dist_mat ~ site / depth, data = meta_tab, permutations = 999) # Site/Depth  0.001 ***
















### Fig 3 C) Shannon Index alpha diversity across sites (MAG based) -------------------

# use meta data from above

# read the mags abundance file
mag_abund_rel <- read.csv('data/Fig3_C_AlpineSoil23_MAGs_abund_rel.tsv', sep='\t', header = T)





### alpha diversity 

# Transpose: samples become rows, ASVs become columns
mag_abund_rel_t <- t(mag_abund_rel)

# Convert to numeric matrix correctly
mag_abund_rel_t_mat <- matrix(
  as.numeric(mag_abund_rel_t),      # flatten all values to numeric
  nrow = nrow(mag_abund_rel_t),     # keep original number of rows
  ncol = ncol(mag_abund_rel_t),     # keep original number of columns
  dimnames = dimnames(mag_abund_rel_t) # keep row and column names
)


# Compute alpha diversity
alpha_div <- data.frame(
  Sample = rownames(mag_abund_rel_t_mat),
  Shannon = vegan::diversity(mag_abund_rel_t_mat, index = "shannon"),
  Simpson = vegan::diversity(mag_abund_rel_t_mat, index = "simpson"),
  Richness = vegan::specnumber(mag_abund_rel_t_mat)
)



# add site per sample
meta_tab$Sample <- meta_tab$ID
alpha_div <- alpha_div %>% left_join(meta_tab%>% select(Sample, site, depth), by = "Sample")

alpha_div <- na.omit(alpha_div)

library(RColorBrewer)
library(ggplot2)

# order sites
alpha_div$site <- factor(alpha_div$site,
                         levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))

custom_colors <- RColorBrewer::brewer.pal(10, "Paired")
names(custom_colors) <- levels(alpha_div$site)


# plot phylum per sample
p_1 <- ggplot(alpha_div, aes(x = site, y = Shannon, color = site)) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7, color = "grey50") + 
  geom_boxplot(fill = NA) +   # outline only, no fill
  labs(x = "Site", y = "Shannon index") +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  theme(legend.position = "none")
p_1

# add vegetation gradient bar
n_sites <- length(levels(alpha_div$site))

library(scales)
n_sites <- length(levels(alpha_div$site))

p_triangle <- ggplot() +
  annotate("polygon",
           x = c(1, n_sites, n_sites),
           y = c(0.5, 1, 0),
           fill = alpha("forestgreen", 0.5),  # transparency
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

library(cowplot)
Fig3C <- plot_grid(
  p_1,
  p_triangle,
  ncol = 1,
  align = "v",
  rel_heights = c(1, 0.08)  # triangle height
)
Fig3C





### Stats:
alpha_div$site <- factor(alpha_div$site)
alpha_div$depth <- factor(alpha_div$depth)

# normality? 
shapiro.test(alpha_div$Shannon) # significant, so not normal

# Test differences across sites
kruskal.test(Shannon ~ site, data = alpha_div) # p-value = 0.0001838
# Post-hoc pairwise Wilcoxon test
pairwise.wilcox.test(alpha_div$Shannon, alpha_div$site, p.adjust.method = "BH")



# Test differences across depths
kruskal.test(Shannon ~ depth, data = alpha_div) # p-value = 0.5769










### Fig 3 D) Microbial cell numbers across sites ---------------------------

library(dplyr)
library(ggplot2)
library(RColorBrewer)


# use meta data from above

# Read the Excel file
data <- read.table('data/Fig3_D_AlpSoils23_cell_numbers_all.txt', sep='\t', header = T)


# remove all 0
data <- data %>% filter(cells.mL != 0)

data <- data %>%
  mutate(site_name = case_when(
    site == "C1" ~ "MY1",
    site == "C2" ~ "MY2",
    site == "K"  ~ "BN",
    site == "L"  ~ "PM",
    site == "A" ~ "SF",
    site == "B" ~ "MF",
    site == "F"  ~ "SN",
    site == "G"  ~ "DS",     
    site == "H" ~ "CD",
    site == "J" ~ "GR" ))









### cells/ g fresh soil weight 

data <- data %>%
  mutate(cells_gFW = cells.mL / g.fresh.soil)



# order sites
data$site_name <- factor(data$site_name,
                         levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))
custom_colors <- RColorBrewer::brewer.pal(10, "Paired")
names(custom_colors) <- levels(alpha_div$site)



# 4) plot using scale_fill_manual with the named vector
p_2 <- ggplot(data, aes(x = site_name, y = cells_gFW, colour = site_name)) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7, color = "grey50") + 
  geom_boxplot(fill = NA) +   # outline only, no fill
  labs(x = "Site", y = "log10(Cells g⁻¹)") +
  scale_colour_manual(values = custom_colors) +
  theme_bw() +
  theme(legend.position = "none") + 
  scale_y_log10()
p_2

# use vegetation gradient bar from above

Fig3D <- plot_grid(
  p_2,
  p_triangle,
  ncol = 1,
  align = "v",
  rel_heights = c(1, 0.08)  # triangle height
)
Fig3D





### Stats:
data$site <- factor(data$site)
data$depth <- factor(data$depth)

# normality? 
shapiro.test(data$cells_gFW) # significant, so not normal

# Test differences across sites
kruskal.test(cells_gFW ~ site_name, data = data) # p-value = < 2.2e-16
# Post-hoc pairwise Wilcoxon test
pairwise.wilcox.test(data$cells_gFW, data$site_name, p.adjust.method = "BH")




# Test differences across depths
kruskal.test(cells_gFW ~ depth, data = data) # p-value = 0.8853









### exporting final figure 3 ---------------
library(ggpubr)

Fig3_final <- ggarrange(
  Fig3A, Fig3B,
  Fig3C, Fig3D,
  ncol = 2, nrow = 2,
  labels = c('A', 'B', 'C', 'D'),
  common.legend = TRUE,
  legend = "right",
  heights = c(1.8, 1)
)

Fig3_final

ggsave("figures/Fig3_final_AlpSoils23_microb.png", plot = Fig3_final, height = 7, width = 10)
ggsave("figures/data/Fig3_final_AlpSoils23_microb.svg", plot = Fig3_final, height = 8, width = 10)






### Supplemental figure S3 - MAG based taxonomic and functional NMDS  ##########################################################


### S_Fig S3 A) NMDS of MAG-based community composition ----------------

# use meta data from above

# read the mags abundance file
mag_abund_abs <- read.csv('data/FigS3_A_AlpineSoil23_MAGs_abund_abs.tsv', sep='\t', header = T)

# read tax tables
mag_tax <- read.csv('data/FigS3_A_AlpineSoil23_MAGs_taxonomy_table.tsv', sep='\t', header = T)






# Transpose: samples become rows, ASVs become columns
rownames(mag_abund_abs) = mag_abund_abs$bin_id
mag_abund_abs$bin_id = NULL

mag_abund_abs_t <- t(mag_abund_abs)
mag_abund_mat <- as.matrix(mag_abund_abs_t)  # already numeric, samples x ASVs



# Hellinger transformation
abund_hell <- sqrt(mag_abund_mat / rowSums(mag_abund_mat))




# Distance matrix and NMDS
dist_mat <- vegdist(abund_hell, method = "bray")  # Bray-Curtis

nmds <- metaMDS(dist_mat)
stressplot(nmds)

# NMDS site scores
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$ID <- rownames(nmds_scores)


# Load metadata/environmental variables
env_vars <- meta_tab
env_vars$Shannon <- NULL
env_vars$proc_euc <- NULL

# Fix sample names
rownames(env_vars) <- env_vars$ID


# Keep only samples present in NMDS
env_vars <- env_vars[rownames(env_vars) %in% rownames(nmds_scores), ]

# Standardize numeric variables
env_num <- env_vars %>%
  select(where(is.numeric)) %>%
  scale() %>%
  as.data.frame()
rownames(env_num) <- rownames(env_vars)


# Environmental fitting (envfit)
envfit_res <- envfit(nmds, env_num, permutations = 999, na.rm = TRUE)



# envfit_res is your envfit result
env_vectors <- scores(envfit_res, display = "vectors")  # coordinates of arrows
env_pvals   <- envfit_res$vectors$pvals                # raw p-values
env_r2      <- envfit_res$vectors$r                    # R² values

# Combine into a data frame
env_summary <- data.frame(
  Variable = rownames(env_vectors),
  NMDS1 = env_vectors[, "NMDS1"],
  NMDS2 = env_vectors[, "NMDS2"],
  R2 = env_r2,
  P = env_pvals
)

# Add BH-adjusted p-values
env_summary$P_adj <- p.adjust(env_summary$P, method = "BH")

# Print all variables with adjusted p-values
print(env_summary[order(env_summary$P_adj), ])





#  Extract significant vectors 
sig_vectors <- envfit_res$vectors$arrows * sqrt(envfit_res$vectors$r)
sig_pvals   <- envfit_res$vectors$pvals
sig_adj     <- p.adjust(sig_pvals, method = "BH") # adjusted p value after Benjamini Hochberg

sig_names <- names(sig_adj)[sig_adj < 0.05]

sig_env <- as.data.frame(sig_vectors[sig_names, , drop = FALSE])
sig_env$Variable <- rownames(sig_env)




# Scale arrows to NMDS range 
arrow_scale <- 0.5 * min(apply(nmds_scores[,c("NMDS1","NMDS2")],2,diff)) /
  max(sqrt(rowSums(sig_env[,1:2]^2)))
sig_env$NMDS1 <- sig_env$NMDS1 * arrow_scale
sig_env$NMDS2 <- sig_env$NMDS2 * arrow_scale
sig_env$x <- 0
sig_env$y <- 0



# Combine NMDS and metadata for plotting
nmds_plot_df <- nmds_scores %>%
  left_join(env_vars, by = c("ID" = "ID"))

# Plot NMDS with arrows
custom_colors <- RColorBrewer::brewer.pal(10, "Paired")
stress_val <- round(nmds$stress, 2)

# order sites by vegetation
nmds_plot_df$site = factor(nmds_plot_df$site, levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))


# Flip only arrows (both axes)
sig_env$NMDS1 <- -sig_env$NMDS1
sig_env$NMDS2 <- -sig_env$NMDS2


FigS3A <- ggplot(nmds_plot_df, aes(x = NMDS1, y = NMDS2, color = site, shape = depth)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_segment(data=sig_env, aes(x=x, y=y, xend=NMDS1, yend=NMDS2),
               arrow=arrow(length=unit(0.25,"cm")), color="black", inherit.aes = FALSE) +
  geom_text_repel(data=sig_env, aes(x=NMDS1*1.1, y=NMDS2*1.1, label=Variable),
                  size=4, color="black", inherit.aes = FALSE) +
  labs(x="NMDS1", y="NMDS2", color="Site", shape="Depth") +
  scale_color_manual(values = custom_colors) +
  annotate("text", x = -Inf, y = -Inf, 
           label = paste("Stress =", stress_val), 
           hjust = -0.1, vjust = -0.5, size = 4) +
  labs(x = "NMDS1", y = "NMDS2", color = "Site", shape = "Depth") +
  scale_shape_manual(
    values = c(17, 16),  # your shapes
    labels = c("Topsoil","Lower soil layer")) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14,
    ),
    plot.margin = margin(20, 10, 10, 10)
  ) +
  labs(title = "Taxonomic level") +
  theme_bw()
FigS3A



# PERMANOVA 
# Make sure samples match
meta_tab <- meta_tab %>% filter(ID %in% rownames(abund_hell))

# Reorder abund_hell to match metadata
abund_hell_sub <- abund_hell[meta_tab$ID, , drop = FALSE]

# Compute Bray-Curtis distance
dist_mat <- vegdist(abund_hell_sub, method = "bray")

# PERMANOVA: site and depth effect
adonis2(dist_mat ~ site, data = meta_tab, permutations = 999) # 0.001 ***
adonis2(dist_mat ~ depth, data = meta_tab, permutations = 999) # 0.749

# PERMANOVA: depth nested within site
adonis2(dist_mat ~ site / depth, data = meta_tab, permutations = 999) # 0.001 ***




### S_Fig S3 B) NMDS of MAG-based functional composition ----------------


# use mag abundance and meta data from above 

# read the mags abundance file (use newly imported file!)
mag_abund_abs <- read.csv('data/FigS3_A_AlpineSoil23_MAGs_abund_abs.tsv', sep='\t', header = T)

# read functional annotations file
ann_tab <- read.csv('data/FigS3_B_AlpineSoil23_MAGs_functional_annotations_emapper.tsv', sep='\t', header = T)



### NMDS of KO terms 

library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tidyverse)


### Extract KEGG KOs from annotations
ann_tab <- ann_tab %>%
  mutate(KEGG_ko = str_extract(KEGG_ko, "K\\d{5}")) %>%
  filter(!is.na(KEGG_ko))


# Join abundance + metadata 
mags_long <- mag_abund_abs %>% # use absolute abundance!
  pivot_longer(-bin_id, names_to = "ID", values_to = "abundance")

# join annotations
mag_ann <- mags_long %>% 
  left_join(ann_tab, by = c("bin_id" ), relationship = "many-to-many")


# summarize KO abundances per sample
ko_wide <- mag_ann %>%
  group_by(ID, KEGG_ko) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")  %>%
  tidyr::pivot_wider(
    names_from = KEGG_ko,
    values_from = abundance,
    values_fill = 0
  )


# Convert to matrix with sample IDs as rownames
ko_mat <- ko_wide %>%
  column_to_rownames("ID") %>%  # rownames = sample IDs
  as.matrix()                   # numeric matrix


# hellinger transform
abund_hell <- sqrt(ko_mat / rowSums(ko_mat))


# distance matrix and nmds
dist_mat <- vegdist(abund_hell, method = "bray")
nmds <- metaMDS(dist_mat)
stressplot(nmds)

nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$ID <- rownames(nmds_scores)


# metadata alignment
env_vars <- meta_tab
env_vars$Shannon <- NULL
env_vars$proc_euc <- NULL
rownames(env_vars) <- env_vars$ID

env_vars <- env_vars[rownames(env_vars) %in% rownames(nmds_scores), ]

env_num <- env_vars %>%
  select(where(is.numeric)) %>%
  scale() %>%
  as.data.frame()
rownames(env_num) <- rownames(env_vars)


# environmental variables
envfit_res <- envfit(nmds, env_num, permutations = 999, na.rm = TRUE)

env_vectors <- scores(envfit_res, display = "vectors")
env_pvals   <- envfit_res$vectors$pvals
env_r2      <- envfit_res$vectors$r

env_summary <- data.frame(
  Variable = rownames(env_vectors),
  NMDS1 = env_vectors[, "NMDS1"],
  NMDS2 = env_vectors[, "NMDS2"],
  R2 = env_r2,
  P = env_pvals
)

env_summary$P_adj <- p.adjust(env_summary$P, method = "BH")
print(env_summary[order(env_summary$P_adj), ])


# sifnificant drivers
sig_vectors <- envfit_res$vectors$arrows * sqrt(envfit_res$vectors$r)
sig_pvals   <- envfit_res$vectors$pvals
sig_adj     <- p.adjust(sig_pvals, method = "BH")

sig_names <- names(sig_adj)[sig_adj < 0.05]

sig_env <- as.data.frame(sig_vectors[sig_names, , drop = FALSE])
sig_env$Variable <- rownames(sig_env)

arrow_scale <- 0.5 * min(apply(nmds_scores[,c("NMDS1","NMDS2")],2,diff)) /
  max(sqrt(rowSums(sig_env[,1:2]^2)))

sig_env$NMDS1 <- sig_env$NMDS1 * arrow_scale
sig_env$NMDS2 <- sig_env$NMDS2 * arrow_scale
sig_env$x <- 0
sig_env$y <- 0




# plot nmds
nmds_plot_df <- nmds_scores %>% 
  left_join(env_vars, by = "ID")



# invert x axis
nmds_plot_df$NMDS1 <- -nmds_plot_df$NMDS1
nmds_plot_df$NMDS2 <- -nmds_plot_df$NMDS2

# Flip arrows
sig_env$NMDS1 <- -sig_env$NMDS1
sig_env$NMDS2 <- -sig_env$NMDS2



# set colors
custom_colors <- RColorBrewer::brewer.pal(10, "Paired")
stress_val <- round(nmds$stress, 2)

# order sites by vegetation
nmds_plot_df$site = factor(nmds_plot_df$site, levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))


FigS3B <- ggplot(nmds_plot_df, aes(x = NMDS1, y = NMDS2, color = site, shape = depth)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_segment(data=sig_env, aes(x=x, y=y, xend=NMDS1, yend=NMDS2),
               arrow=arrow(length=unit(0.25,"cm")), color="black", inherit.aes = FALSE) +
  geom_text_repel(data=sig_env, aes(x=NMDS1*1.1, y=NMDS2*1.1, label=Variable),
                  size=4, color="black", inherit.aes = FALSE) +
  labs(x="NMDS1", y="NMDS2", color="Site", shape="Depth") +
  scale_color_manual(values = custom_colors) +
  annotate("text", x = -Inf, y = -Inf, 
           label = paste("Stress =", stress_val), 
           hjust = -0.1, vjust = -0.5, size = 4) +
  scale_shape_manual(
    values = c(17, 16),  # your shapes
    labels = c("Topsoil","Lower soil layer")) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14,
    ),
    plot.margin = margin(20, 10, 10, 10)
  ) +
  labs(title = "Functional level") +
  theme_bw()
FigS3B




### PERMANOVA
meta_tab <- meta_tab %>% filter(ID %in% rownames(abund_hell))

abund_hell_sub <- abund_hell[meta_tab$ID, , drop = FALSE]
dist_mat <- vegdist(abund_hell_sub, method = "bray")

adonis2(dist_mat ~ site, data = meta_tab, permutations = 999) # 0.001 ***
adonis2(dist_mat ~ depth, data = meta_tab, permutations = 999) # 0.629

adonis2(dist_mat ~ site / depth, data = meta_tab, permutations = 999) # 0.001 ***







### exporting final S figure S3 ---------------
library(ggpubr)

S_FigS3_final <- ggarrange(
  FigS3A, FigS3B,
  ncol = 2, nrow = 1,
  labels = c('A', 'B'),
  common.legend = TRUE,
  legend = "right"
)

S_FigS3_final

ggsave("figures/S_FigS3_final_AlpSoils23_NMDSs_MAG.png", plot = S_FigS3_final, height = 4, width = 9)




### Supplemental figure S5 - Functional alpha diversity (KEGG ko based) (MAG based)  ##########################################################

# use meta data from above and ko_mat


# Compute functional diversity
alpha_div_funct <- data.frame(
  Sample = rownames(ko_mat),
  Shannon = vegan::diversity(ko_mat, index = "shannon"),
  Simpson = vegan::diversity(ko_mat, index = "simpson"),
  Richness = vegan::specnumber(ko_mat)
)

# add site per sample
meta_tab$Sample <- meta_tab$ID
alpha_div_funct <- alpha_div_funct %>% left_join(meta_tab %>% select(Sample, site, depth), by = "Sample")

# Ensure consistent site order
alpha_div_funct$site <- factor(alpha_div_funct$site,
                         levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))
custom_colors <- RColorBrewer::brewer.pal(10, "Paired")


# plot
p_3 <- ggplot(alpha_div_funct, aes(x = site, y = Shannon, color = site)) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7, color = "grey50") + 
  geom_boxplot(fill = NA) +   # outline only, no fill
  labs(x = "Site", y = "Functional diverstiy\n(Shannon Index)") +
  scale_colour_manual(values = custom_colors) +
  theme_bw() +
  theme(legend.position = "none")
p_3 


# use vegetation gradient bar from above

S_FigS5 <- plot_grid(
  p_3,
  p_triangle,
  ncol = 1,
  align = "v",
  rel_heights = c(1, 0.08)  # triangle height
)
S_FigS5

ggsave("figures/S_FigS5_final_AlpSoils23_functional_diverstiy.png", plot = S_FigS6, height = 3, width = 4.5)





# stats:
alpha_div_funct$site <- factor(alpha_div_funct$site)
alpha_div_funct$depth <- factor(alpha_div_funct$depth)

# normality? 
shapiro.test(alpha_div_funct$Shannon) # significant, so not normal

# Test differences across sites
kruskal.test(Shannon ~ site, data = alpha_div_funct) # p-value = 6.146e-010

# Post-hoc pairwise Wilcoxon test
pairwise.wilcox.test(alpha_div_funct$Shannon, alpha_div_funct$site, p.adjust.method = "BH")

# Test differences across depths
kruskal.test(Shannon ~ depth, data = alpha_div_funct) # p-value = 0.3944




















### Supplemental figure S6 - abundance barplot (contig based)  ##########################################################

# use meta data from above

# read coverage table (use newly imported table!)
abund_tab_abs <- read.csv('data/Fig3_A_AlpineSoils23_contigs_coverage_abs.tsv', sep='\t', header = T)
abund_tab_rel <- read.csv('data/FigS5_B_AlpineSoils23_contigs_coverage_rel.tsv', sep='\t', header = T)


# read tax table
tax_tab <- read.csv('data/FigS5_A_AlpineSoil23_contigs_taxonomy_table.tsv', sep='\t', header = T)

# Clean taxonomy table: replace NAs with "unclassified"
tax_tab_clean <- tax_tab %>%
  select(contig_id, kingdom:species) %>%
  distinct(contig_id, .keep_all = TRUE) %>%
  mutate(across(kingdom:species, ~replace_na(.x, "Unclassified"))) %>%
  column_to_rownames("contig_id")



# A) absolute abundance 



rownames(abund_tab_abs) = abund_tab_abs$contig_id
abund_tab_abs$contig_id = NULL

abund_tab_abs_melt = reshape2::melt(as.matrix(abund_tab_abs))
colnames(abund_tab_abs_melt) = c('Contig', 'Sample', 'Abs_abundance')

#abund_tab_abs_melt$Phylum = map_chr(abund_tab_abs_melt$Contig, function(x) tax_tab_clean$phylum[rownames(tax_tab_clean) == x])
tax_tab_clean_df <- tax_tab_clean %>% tibble::rownames_to_column("Contig")

abund_tab_abs_melt <- abund_tab_abs_melt %>% left_join(tax_tab_clean_df %>% select(Contig, phylum), by = "Contig")


# sum abundance by sample, phylum
abund_tab_abs_sum = abund_tab_abs_melt %>% 
  group_by(Sample, phylum) %>% 
  summarise(Abs_abundance=sum(Abs_abundance))

# add site per sample
meta_tab$Sample <- meta_tab$ID
abund_tab_abs_sum <- abund_tab_abs_sum %>% left_join(meta_tab %>% select(Sample, site, depth), by = "Sample")


# select top 15 phyla
top_phyla_relab <- abund_tab_abs_melt %>%
  group_by(Sample, phylum) %>%
  summarise(Abs_abundance = sum(Abs_abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(
    total_abund = sum(Abs_abundance, na.rm = TRUE),
    Rel_abundance = if_else(total_abund > 0, Abs_abundance / total_abund, 0)
  ) %>%
  ungroup() %>%
  filter(!(phylum %in% c("Unclassified", "Others"))) %>%
  group_by(phylum) %>%
  summarise(mean_rel_abund = mean(Rel_abundance, na.rm = TRUE), .groups = "drop") %>%
  slice_max(mean_rel_abund, n = 9)

top_phyla_vec <- top_phyla_relab$phylum

abund_tab_abs_sum <- abund_tab_abs_sum %>%
  mutate(phylum = if_else(phylum %in% top_phyla_vec, phylum, "Others"))

# rename Phylum
abund_tab_abs_sum$Phylum <- abund_tab_abs_sum$phylum


# Ensure consistent site order
abund_tab_abs_sum$site <- factor(abund_tab_abs_sum$site,
                                 levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))

# plot phylum per sample
FigS6A <- ggplot(abund_tab_abs_sum, aes(x = Sample, y = Abs_abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(. ~ site, scales = "free_x", space = "free_x") +
  ylab("Cell number normalized abundance (cells g⁻¹)") +
  xlab("Sample") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_line()
  ) 
FigS6A






# B) relative abundance



rownames(abund_tab_rel) = abund_tab_rel$contig_id
abund_tab_rel$contig_id = NULL

abund_tab_rel_melt = reshape2::melt(as.matrix(abund_tab_rel))
colnames(abund_tab_rel_melt) = c('Contig', 'Sample', 'rel_abundance')

#abund_tab_rel_melt$Phylum = map_chr(abund_tab_rel_melt$Contig, function(x) tax_tab_clean$phylum[rownames(tax_tab_clean) == x])
tax_tab_clean_df <- tax_tab_clean %>% tibble::rownames_to_column("Contig")

abund_tab_rel_melt <- abund_tab_rel_melt %>% left_join(tax_tab_clean_df %>% select(Contig, genus), by = "Contig")


# sum abundance by sample, genus
abund_tab_rel_sum = abund_tab_rel_melt %>% 
  group_by(Sample, genus) %>% 
  summarise(rel_abundance=sum(rel_abundance))

# add site per sample
meta_tab$Sample <- meta_tab$ID
abund_tab_rel_sum <- abund_tab_rel_sum %>% left_join(meta_tab %>% select(Sample, site, depth), by = "Sample")


# select top 15 phyla
library(dplyr)
top_phyla_relab <- abund_tab_rel_melt %>%
  group_by(Sample, genus) %>%
  summarise(rel_abundance = sum(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Sample) %>%
  ungroup() %>%
  filter(!(genus %in% c("Unclassified", "Others"))) %>%
  group_by(genus) %>%
  summarise(mean_rel_abund = mean(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
  slice_max(mean_rel_abund, n = 15)

top_phyla_vec <- top_phyla_relab$genus

abund_tab_rel_sum <- abund_tab_rel_sum %>%
  mutate(genus = if_else(genus %in% top_phyla_vec, genus, "Others"))

# rename genus
abund_tab_rel_sum$Genus <- abund_tab_rel_sum$genus


# Ensure consistent site order
abund_tab_rel_sum$site <- factor(abund_tab_rel_sum$site,
                                 levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))

# plot genus per sample
FigS6B <- ggplot(abund_tab_rel_sum, aes(x = Sample, y = rel_abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(. ~ site, scales = "free_x", space = "free_x") +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_line()
  ) 
FigS6B


### exporting final S figure S6 ---------------
library(ggpubr)

S_FigS6_final <- ggarrange(
  FigS6A, FigS6B,
  ncol = 1, nrow = 2,
  labels = c('A', 'B'),
  common.legend = TRUE,
  legend = "right"
)

S_FigS6_final

ggsave("figures/S_FigS6_final_AlpSoils23_taxonomy_abs_rel.png", plot = S_FigS6_final, height = 7.5, width = 10)






