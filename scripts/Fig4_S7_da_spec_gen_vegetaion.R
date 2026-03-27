###
### AlpineSoils - summer 23 - microbial community, function, diverstiy and cell numbers
### 


rm(list=ls(all=TRUE)) # removes everything


### Fig 4 A) Differential abundance analysis -----------------------------

library(ANCOMBC)
library(dplyr)
library(ggtext)
library(ggplot2)
library(ggpubr)



# read the mags abundance file
mag_abund_raw  <- read.csv('data/Fig4_A_AlpineSoil23_MAGs_abund_raw.tsv', sep='\t', header = T)

# read tax tables
mag_tax <- read.csv('data/Fig4_A_AlpineSoil23_MAGs_taxonomy_table.tsv', sep='\t', header = T)


# read in meta data
meta_tab <- read.csv('data/Fig3_A_AlpineSoil23_contigs_metadata.tsv', sep='\t', header = T)




### Differential abundance 

# prepare raw count table 
counts <- mag_abund_raw %>% # use raw abundance 
  tibble::column_to_rownames("bin_id")


# add metadata
meta <- meta_tab %>%
  tibble::column_to_rownames("ID")

all(colnames(counts) == rownames(meta))

# transforming veg not necessary for DA




# run ANCOM-BC
out <- ancombc2(
  data = counts,        
  meta_data = meta,     
  fix_formula = "vegetation",
  p_adj_method = "BH",
  alpha = 0.05,
  global = FALSE
)



# extract results
res <- out$res

# check significant MAGs
res_flagged <- res %>%
  mutate(
    signif_vegetation = diff_robust_vegetation & passed_ss_vegetation
  )


res_flagged$bin_id <- res_flagged$taxon
res_flagged$taxon = NULL

# join tax info
res_tax <- res_flagged %>% dplyr::left_join(mag_tax, by = "bin_id")

# keep only relevant stats
res_tax2 <- res_tax %>%
  select(
    bin_id,
    lfc_vegetation, p_vegetation, q_vegetation, se_vegetation,
    diff_vegetation, diff_robust_vegetation,
    signif_vegetation,
    Phylum, Class, Order, Family, Genus
  )



# make "others"
res_tax2 <- res_tax2 %>% mutate(Family = ifelse(Family == "f__", "Others", Family))


# remove f_
res_tax2$Family <- sub("^f__", "", res_tax2$Family)
res_tax2$Phylum <- sub("^p__", "", res_tax2$Phylum)


# filter only sig flagges phyla
res_tax2 <- res_tax2 %>%
  group_by(Phylum) %>%
  filter(any(signif_vegetation)) %>%
  ungroup()


res_tax2 <- droplevels(res_tax2)

# dot plot
Fig4A <- ggplot(res_tax2, aes(x = lfc_vegetation, y = Family)) +
  ## non-significant background
  geom_errorbar(
    data = subset(res_tax2, !signif_vegetation),
    aes(
      xmin = lfc_vegetation - se_vegetation,
      xmax = lfc_vegetation + se_vegetation),
    height = 0.15, color = "grey80") +
  geom_point(
    data = subset(res_tax2, !signif_vegetation),
    color = "grey70",
    size = 2,
    position = position_jitter(width = 0, height = 0.1)) +
  ## significant foreground
  geom_errorbar(
    data = subset(res_tax2, signif_vegetation),
    aes(
      xmin = lfc_vegetation - se_vegetation,
      xmax = lfc_vegetation + se_vegetation,
      color = Phylum), height = 0.2, linewidth = 0.6) +
  geom_point(
    data = subset(res_tax2, signif_vegetation),
    aes(color = Phylum),
    size = 2.7,
    position = position_jitter(width = 0,height = 0.1)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_grid(
    Phylum ~ .,
    scales = "free_y",
    space  = "free_y",
    switch = "y") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.just = "left",
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.key.width  = unit(0.6, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.spacing.x  = unit(0.2, "lines"),
        legend.spacing.y  = unit(0.2, "lines")
  ) +
  guides(color = guide_legend(ncol = 2)) +
  labs(
    y = "Family",
    x = "Log fold change (Vegetation)",
    color = "Phylum")
Fig4A


### Fig 4 B) Shannon index vs vegetation (GAM) -------------------


# read the mags abundance file
mag_abund_rel <- read.csv('data/Fig3_C_AlpineSoil23_MAGs_abund_rel.tsv', sep='\t', header = T)

# use meta data from above



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

  
  
### plot vegetation against alpha diversity

# transform vegetation
meta_tab$vegetation_trans <- log1p(meta_tab$vegetation)
shapiro.test(meta_tab$vegetation_trans)



# add meta_data to alpha diversity df
meta_tab$Sample <- meta_tab$ID
alpha_div <- alpha_div %>% left_join(meta_tab %>% select(Sample, vegetation_trans), by = "Sample")

alpha_div <- na.omit(alpha_div)


# run GAM
library(mgcv)
gam_fit <- gam(Shannon ~ s(vegetation_trans, k = 3), data = alpha_div)
summary_gam <- summary(gam_fit)
summary_gam


# extract stats
r2_adj <- summary_gam$r.sq
dev_expl <- summary_gam$dev.expl
f_stat <- summary_gam$s.table[1, "F"]
p_val <- summary_gam$s.table[1, "p-value"]

lab <- paste0(
  #"Adj R² = ", round(r2_adj, 3), "\n",
  "Deviance explained = ", round(dev_expl * 100, 1), "%\n",
  #"F = ", round(f_stat, 2), "\n",
  "p = ", signif(p_val, 3)
)


# Ensure site order
alpha_div$site <- factor(alpha_div$site,
                         levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))

custom_colors <- RColorBrewer::brewer.pal(10, "Paired")
names(custom_colors) <- levels(alpha_div$site)


Fig4B <-ggplot(alpha_div, aes(x = vegetation_trans, y = Shannon, colour = site)) +
  geom_point() +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 3),
    color = "black",
    se = TRUE
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = lab,
    hjust = 1.1, vjust = -0.5,
    size = 3.5
  ) +
  scale_colour_manual(values = custom_colors) +
  theme_bw() +
  labs(y = "Shannon index", x = "Vegetation", color = "Site") +
  theme(legend.position = "bottom", 
        legend.box.just = "left",
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.key.width  = unit(0.6, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.spacing.x  = unit(0.2, "lines"),
        legend.spacing.y  = unit(0.2, "lines"))
Fig4B





### Fig 4 C) Specialist:Generalist MAGs vs vegetation (GAM) --------------
library(dplyr)
library(tibble)
library(tidyr)


# use meta data from above
# run alpha div section first for vegetation_transformed

# read the mags abundance file
mag_abund_abs <- read.csv('data/Fig4_C_AlpineSoil23_MAGs_abund_abs.tsv', sep='\t', header = T)


# prepare df
rownames(mag_abund_abs) <- mag_abund_abs$bin_id
mag_abund_abs$bin_id <- NULL
mag_abund_abs_t <- as.data.frame(t(mag_abund_abs))


# make presence absence table
mag_pa <- mag_abund_abs_t > 0

# aggreagte at site level
library(dplyr)
library(tidyr)
library(tibble)

mag_site <- mag_abund_abs_t %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(meta_tab, by = "ID") %>%
  group_by(site) %>%
  summarise(across(starts_with("bin_"), mean), .groups = "drop")

mag_site_mat <- mag_site %>%
  column_to_rownames("site") %>%
  as.matrix()

# normalize per MAG across site
rel_abund <- apply(mag_site_mat, 2, function(x) {
  if(sum(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  x / sum(x, na.rm = TRUE)
})

# calculate Levin's niche breadths
n_sites <- nrow(mag_site_mat)

B <- apply(rel_abund, 2, function(p) {
  if(all(p == 0)) return(NA)
  1 / (n_sites * sum(p^2, na.rm = TRUE))
})

#filter rare MAGs
mean_abund <- colMeans(mag_site_mat, na.rm = TRUE)

threshold <- quantile(mean_abund, 0.05, na.rm = TRUE)  # remove lowest 5%
keep <- mean_abund > threshold

B <- B[keep]

# define specialists and generalists
q_low  <- quantile(B, 0.1, na.rm = TRUE)
q_high <- quantile(B, 0.9, na.rm = TRUE)

sg_vec <- case_when(
  B <= q_low  ~ "specialist",
  B >= q_high ~ "generalist",
  TRUE ~ "intermediate"
)

names(sg_vec) <- names(B)


# subset presence/absence to retained MAGs
mag_pa_sub <- mag_pa[, names(sg_vec)]


# add metadata - use TRANSFORMED
mag_pa_df <- mag_pa_sub %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(meta_tab, by = "ID")


# count specialist and generalist per sample
sg_counts <- mag_pa_df %>%
  pivot_longer(cols = starts_with("bin_"),
               names_to = "MAG",
               values_to = "present") %>%
  filter(present == TRUE) %>%
  mutate(sign = sg_vec[MAG]) %>%
  group_by(ID, sign) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = sign, values_from = n, values_fill = 0) %>%
  left_join(meta_tab %>% distinct(ID, site, vegetation_trans),
            by = "ID")


# calculate ratio
sg_counts$spec_gen_ratio = sg_counts$specialist / sg_counts$generalist






# run GAM for ratio
library(mgcv)
gam_fit <- gam(spec_gen_ratio ~ s(vegetation_trans, k = 3), data = sg_counts)
summary_gam <- summary(gam_fit)
summary_gam

# extract stats
r2_adj <- summary_gam$r.sq
dev_expl <- summary_gam$dev.expl
f_stat <- summary_gam$s.table[1, "F"]
p_val <- summary_gam$s.table[1, "p-value"]

lab <- paste0(
  #"Adj R² = ", round(r2_adj, 3), "\n",
  "Deviance explained = ", round(dev_expl * 100, 1), "%\n",
  #"F = ", round(f_stat, 2), "\n",
  "p = ", signif(p_val, 3)
)


# Ensure alphabetical site order
sg_counts$site <- factor(sg_counts$site,
                         levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))
custom_colors <- RColorBrewer::brewer.pal(10, "Paired")

# plot with stats
library(ggplot2)
Fig4C <- ggplot(sg_counts, aes(x = vegetation_trans, y = spec_gen_ratio, colour = site)) +
  geom_point() +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 3),
    color = "black",
    se = TRUE
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = lab,
    hjust = 1.1, vjust = -0.5,
    size = 3.5
  ) +
  scale_color_manual(values = custom_colors) +
  theme_bw() + 
  labs(x = "Vegetation", y = "Specialist:Generalist MAG", color = "Site")  +
  theme(legend.position = "bottom", 
        legend.box.just = "left",
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.key.width  = unit(0.6, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.spacing.x  = unit(0.2, "lines"),
        legend.spacing.y  = unit(0.2, "lines"))
Fig4C




### exporting final figure 4 ---------------
library(ggpubr)

right_column <- ggarrange(
  Fig4B, Fig4C,
  ncol = 1,
  nrow = 2,
  labels = c("B", "C"),
  common.legend = TRUE,
  legend = "bottom"
)

Fig4_final <- ggarrange(
  Fig4A, right_column,
  ncol = 2,
  labels = c("A", ""),
  widths = c(1.1, 1)
)

Fig4_final

ggsave("figures/Fig4_final_AlpSoils23_abund_div_vegetation.png", plot = Fig4_final, height = 9, width = 9)
ggsave("figures/Fig4_final_AlpSoils23_abund_div_vegetation.svg", plot = Fig4_final, height = 9, width = 9)












### Supplementary Fig S7 -  Specialist:Generalist MAGs using LEVIN and OCCUPANCY --------------

library(EcolUtils)
library(spaa)
library(dplyr)
library(tibble)
library(tidyr)


# use meta data from above
# run alpha div section first for vegetation_transformed

# read the mags abundance file
mag_abund_abs <- read.csv('data/Fig4_C_AlpineSoil23_MAGs_abund_abs.tsv', sep='\t', header = T)


# prepare df
rownames(mag_abund_abs) <- mag_abund_abs$bin_id
mag_abund_abs$bin_id <- NULL
mag_abund_abs_t <- as.data.frame(t(mag_abund_abs))


# make integers
scale_factor <- 1/(min(mag_abund_abs_t[mag_abund_abs_t > 0], na.rm = TRUE))
mag_abund_int <- round(mag_abund_abs_t * scale_factor)


# make presence absence table
mag_pa <- mag_abund_abs_t > 0



## A) LEVIN INDEX
# run spaa 
sg_df <- spec.gen(mag_abund_int)
sg_df

# match with sg_df 
sg_vec <- sg_df$sign 
names(sg_vec) <- rownames(sg_df) 
mag_pa <- mag_pa[, names(sg_vec)]# Ensure column names match

names(sg_vec) <- colnames(mag_pa)


# add metadata - use TRANSFORMED
mag_pa_df <- mag_pa %>%
  as.data.frame() %>%            # convert matrix to data frame
  tibble::rownames_to_column("ID") %>%
  left_join(meta_tab, by = "ID")

# count MAGs per category
sg_counts <- mag_pa_df %>%
  pivot_longer(cols = starts_with("bin_"),
               names_to = "MAG",
               values_to = "present") %>%
  filter(present == TRUE) %>%                  # only present MAGs
  mutate(sign = sg_vec[MAG]) %>%               # attach specialist/generalist status
  group_by(ID, sign) %>%                     # count per category
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = sign, values_from = n, values_fill = 0) %>%
  left_join(meta_tab %>% distinct(ID, site, vegetation_trans),
            by = "ID")


# calculate ratio
sg_counts$spec_gen_ratio = sg_counts$SPECIALIST / sg_counts$GENERALIST



# run GAM for ratio
library(mgcv)
gam_fit <- gam(spec_gen_ratio ~ s(vegetation_trans, k = 3), data = sg_counts)
summary_gam <- summary(gam_fit)
summary_gam

# extract stats
r2_adj <- summary_gam$r.sq
dev_expl <- summary_gam$dev.expl
f_stat <- summary_gam$s.table[1, "F"]
p_val <- summary_gam$s.table[1, "p-value"]

lab <- paste0(
  #"Adj R² = ", round(r2_adj, 3), "\n",
  "Deviance explained = ", round(dev_expl * 100, 1), "%\n",
  #"F = ", round(f_stat, 2), "\n",
  "p = ", signif(p_val, 3)
)


# Ensure alphabetical site order
sg_counts$site <- factor(sg_counts$site,
                         levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))
custom_colors <- RColorBrewer::brewer.pal(10, "Paired")

# plot with stats
library(ggplot2)
FigS7A <- ggplot(sg_counts, aes(x = vegetation_trans, y = spec_gen_ratio, colour = site)) +
  geom_point() +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 3),
    color = "black",
    se = TRUE
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = lab,
    hjust = 1.1, vjust = -0.5,
    size = 3.5
  ) +
  scale_color_manual(values = custom_colors) +
  theme_bw() + 
  labs(title = "Levin's Index\n", x = "Vegetation", y = "Specialist:Generalist MAG", color = "Site")  +
  theme(legend.position = "bottom", 
        legend.box.just = "left",
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.key.width  = unit(0.6, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.spacing.x  = unit(0.2, "lines"),
        legend.spacing.y  = unit(0.2, "lines"))
FigS7A








## B) OCCUPANCY with abundance
occupancy <- colSums(mag_pa) / nrow(mag_pa)

mean_abund <- colMeans(mag_abund_abs_t)
sg_vec <- case_when(
  occupancy >= 0.9 & mean_abund >= median(mean_abund) ~ "generalist",
  occupancy <= 0.5 & mean_abund < median(mean_abund) ~ "specialist",
  TRUE ~ "intermediate"
)

names(sg_vec) <- colnames(mag_pa)


# add metadata - use TRANSFORMED
mag_pa_df <- mag_pa %>%
  as.data.frame() %>%            # convert matrix to data frame
  tibble::rownames_to_column("ID") %>%
  left_join(meta_tab, by = "ID")

# count MAGs per category
sg_counts <- mag_pa_df %>%
  pivot_longer(cols = starts_with("bin_"),
               names_to = "MAG",
               values_to = "present") %>%
  filter(present == TRUE) %>%                  # only present MAGs
  mutate(sign = sg_vec[MAG]) %>%               # attach specialist/generalist status
  group_by(ID, sign) %>%                     # count per category
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = sign, values_from = n, values_fill = 0) %>%
  left_join(meta_tab %>% distinct(ID, site, vegetation_trans),
            by = "ID")


# calculate ratio
sg_counts$spec_gen_ratio = sg_counts$specialist / sg_counts$generalist



# run GAM for ratio
library(mgcv)
gam_fit <- gam(spec_gen_ratio ~ s(vegetation_trans, k = 3), data = sg_counts)
summary_gam <- summary(gam_fit)
summary_gam

# extract stats
r2_adj <- summary_gam$r.sq
dev_expl <- summary_gam$dev.expl
f_stat <- summary_gam$s.table[1, "F"]
p_val <- summary_gam$s.table[1, "p-value"]

lab <- paste0(
  #"Adj R² = ", round(r2_adj, 3), "\n",
  "Deviance explained = ", round(dev_expl * 100, 1), "%\n",
  #"F = ", round(f_stat, 2), "\n",
  "p = ", signif(p_val, 3)
)


# Ensure alphabetical site order
sg_counts$site <- factor(sg_counts$site,
                         levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))
custom_colors <- RColorBrewer::brewer.pal(10, "Paired")

# plot with stats
library(ggplot2)
FigS7B <- ggplot(sg_counts, aes(x = vegetation_trans, y = spec_gen_ratio, colour = site)) +
  geom_point() +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 3),
    color = "black",
    se = TRUE
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = lab,
    hjust = 1.1, vjust = -0.5,
    size = 3.5
  ) +
  scale_color_manual(values = custom_colors) +
  theme_bw() + 
  labs(title= "Occupancy Thresholds\n",  x = "Vegetation", y = "Specialist:Generalist MAG", color = "Site")  +
  theme(legend.position = "bottom", 
        legend.box.just = "left",
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.key.width  = unit(0.6, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.spacing.x  = unit(0.2, "lines"),
        legend.spacing.y  = unit(0.2, "lines"))
FigS7B




## c) ABUNDANCE WEIGHTED INVERSE SIMPSON
# determine specialists and generalists
rel_abund_mag <- apply(mag_abund_abs_t, 2, function(x) {
  x / sum(x, na.rm = TRUE)
})



niche_breadth <- apply(rel_abund_mag, 2, function(p) {
  1 / sum(p^2, na.rm = TRUE)
})



n_samples <- nrow(mag_abund_abs_t)
niche_breadth_std <- (niche_breadth - 1) / (n_samples - 1)

q_low  <- quantile(niche_breadth_std, 0.25)
q_high <- quantile(niche_breadth_std, 0.75)

sg_vec <- case_when(
  niche_breadth_std <= q_low  ~ "specialist",
  niche_breadth_std >= q_high ~ "generalist",
  TRUE ~ "intermediate"
)

names(sg_vec) <- colnames(mag_pa)

# add metadata - use TRANSFORMED
mag_pa_df <- mag_pa %>%
  as.data.frame() %>%            # convert matrix to data frame
  tibble::rownames_to_column("ID") %>%
  left_join(meta_tab, by = "ID")

# count MAGs per category
sg_counts <- mag_pa_df %>%
  pivot_longer(cols = starts_with("bin_"),
               names_to = "MAG",
               values_to = "present") %>%
  filter(present == TRUE) %>%                  # only present MAGs
  mutate(sign = sg_vec[MAG]) %>%               # attach specialist/generalist status
  group_by(ID, sign) %>%                     # count per category
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = sign, values_from = n, values_fill = 0) %>%
  left_join(meta_tab %>% distinct(ID, site, vegetation_trans),
            by = "ID")


# calculate ratio
sg_counts$spec_gen_ratio = sg_counts$specialist / sg_counts$generalist



# run GAM for ratio
library(mgcv)
gam_fit <- gam(spec_gen_ratio ~ s(vegetation_trans, k = 3), data = sg_counts)
summary_gam <- summary(gam_fit)
summary_gam

# extract stats
r2_adj <- summary_gam$r.sq
dev_expl <- summary_gam$dev.expl
f_stat <- summary_gam$s.table[1, "F"]
p_val <- summary_gam$s.table[1, "p-value"]

lab <- paste0(
  #"Adj R² = ", round(r2_adj, 3), "\n",
  "Deviance explained = ", round(dev_expl * 100, 1), "%\n",
  #"F = ", round(f_stat, 2), "\n",
  "p = ", signif(p_val, 3)
)


# Ensure alphabetical site order
sg_counts$site <- factor(sg_counts$site,
                         levels = c("PM", "MF", "DS", "CD", "GR", "MY1", "SN", "MY2", "SF", "BN"))
custom_colors <- RColorBrewer::brewer.pal(10, "Paired")

# plot with stats
library(ggplot2)
FigS7C <- ggplot(sg_counts, aes(x = vegetation_trans, y = spec_gen_ratio, colour = site)) +
  geom_point() +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 3),
    color = "black",
    se = TRUE
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = lab,
    hjust = 1.1, vjust = -0.5,
    size = 3.5
  ) +
  scale_color_manual(values = custom_colors) +
  theme_bw() + 
  labs(title= "Abundance Weighted\nInverse Simpson",  x = "Vegetation", y = "Specialist:Generalist MAG", color = "Site")  +
  theme(legend.position = "bottom", 
        legend.box.just = "left",
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.key.width  = unit(0.6, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.spacing.x  = unit(0.2, "lines"),
        legend.spacing.y  = unit(0.2, "lines"))
FigS7C










### exporting final figure S7 ---------------
library(ggpubr)


FigS7_final <- ggarrange(
  FigS7A, FigS7B, FigS7C,
  ncol = 3,
  labels = c("A", "B", "C"), 
  common.legend = TRUE,
  legend = "bottom"
)

FigS7_final

ggsave("figures/FigS7_final_AlpSoils23_spec_gen.png", plot = FigS7_final, height = 4.5, width = 10)
