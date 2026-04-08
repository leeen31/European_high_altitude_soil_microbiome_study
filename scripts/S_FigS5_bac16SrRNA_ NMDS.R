###
### AlpineSoils - summer 23 - bacterial 16S rRNA gene amplicon NMDS
###

rm(list=ls(all=TRUE)) # removes everything


# load libs
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(grid)  # for arrow()



# Load data
otu_tab <- read.table('data/S_FigS5_AlpSoils23_16S_otu_table.tsv', sep='\t', header = TRUE)
rownames(otu_tab) <- otu_tab$OTU
otu_tab$OTU <- NULL

tax_tab_raw <- read.table('data/S_FigS5_AlpSoils23_16S_taxonomy.tsv', sep=',', header = TRUE)
meta_tab <- read.table('data/Fig2_AB_AlpineSoil23_metadata.tsv', sep='\t', header = TRUE)
meta_tab$ID <- gsub('-', '_', meta_tab$ID)
rownames(meta_tab) <- meta_tab$ID



# Taxonomy table processing 
tax_tab_raw <- tax_tab_raw[order(tax_tab_raw$Cluster),]
tax_tab <- tax_tab_raw %>% 
  select(Taxonomy) %>%
  separate(Taxonomy, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), ";") %>%
  as.matrix()
tax_tab <- gsub("d__|p__|c__|o__|f__|g__|s__", "", tax_tab)
rownames(tax_tab) <- rownames(otu_tab)



# Negative control filtering 
negative_mask <- !((otu_tab$NC_soil_1 > 0) | (otu_tab$NC_soil_2 > 0) |
                     (otu_tab$NC_soil_3 > 0) | (otu_tab$NC_PCR_1 > 0) |
                     (otu_tab$NC_PCR_2 > 0) | (otu_tab$NC_libprep_1 > 0))
otu_tab <- otu_tab[negative_mask, ]
tax_tab <- tax_tab[negative_mask, ]



# Remove negative control columns
nc_cols <- grep("NC_", colnames(otu_tab), value = TRUE)
otu_tab <- otu_tab[, !colnames(otu_tab) %in% nc_cols]



# If OTU IDs are still a column
if("OTU.ID" %in% colnames(otu_tab)) {
  rownames(otu_tab) <- otu_tab$OTU.ID
  otu_tab$OTU.ID <- NULL
}


# Keep only numeric columns (counts)
otu_tab <- otu_tab[, sapply(otu_tab, is.numeric)]


# Convert to numeric matrix
otu_tab <- as.matrix(otu_tab)
mode(otu_tab) <- "numeric"


# Keep only samples present in both
common_samples <- intersect(colnames(otu_tab), meta_tab$ID)
otu_tab <- otu_tab[, common_samples, drop = FALSE]
meta_tab <- meta_tab[common_samples, , drop = FALSE]

# Identify zero-sum OTUsand remove
zero_rows <- which(rowSums(otu_tab) == 0)
length(zero_rows)
otu_tab <- otu_tab[rowSums(otu_tab) > 0, ]




## NMDS ----------------------


# --- Hellinger transformation ---
abund_hell <- sqrt(otu_tab / rowSums(otu_tab))
sum(is.na(abund_hell))  # Should be 0




# aggregate OTUs by genus

# Keep only taxonomy rows that are still in otu_tab
tax_tab_filtered <- tax_tab[rownames(otu_tab), , drop = FALSE]

# Now aggregate by genus
otu_genus <- aggregate(otu_tab, by = list(tax_tab_filtered[, "Genus"]), FUN = sum)

# Set row names to genus
rownames(otu_genus) <- otu_genus$Group.1
otu_genus$Group.1 <- NULL


abund_hell <- sqrt(otu_genus / rowSums(otu_genus))


# --- Bray-Curtis distance and NMDS ---
library(parallelDist)
dist_mat <- parDist(t(abund_hell), method = "bray")




nmds <- metaMDS(dist_mat)
stressplot(nmds)

# --- NMDS site scores ---
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$ID <- rownames(nmds_scores)

# --- Prepare environmental data ---
env_vars <- meta_tab
env_vars$Shannon <- NULL
env_vars$proc_euc <- NULL
rownames(env_vars) <- env_vars$ID

# Keep only samples present in NMDS
env_vars <- env_vars[rownames(env_vars) %in% rownames(nmds_scores), ]

# Standardize numeric variables
env_num <- env_vars %>%
  select(where(is.numeric)) %>%
  scale() %>%
  as.data.frame()
rownames(env_num) <- rownames(env_vars)

# --- Environmental fitting ---
envfit_res <- envfit(nmds, env_num, permutations = 999, na.rm = TRUE)

# Extract envfit vectors and stats
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

# --- Extract significant vectors for plotting ---
sig_vectors <- envfit_res$vectors$arrows * sqrt(envfit_res$vectors$r)
sig_pvals   <- envfit_res$vectors$pvals
sig_adj     <- p.adjust(sig_pvals, method = "BH")
sig_names   <- names(sig_adj)[sig_adj < 0.05]

sig_env <- as.data.frame(sig_vectors[sig_names, , drop = FALSE])
sig_env$Variable <- rownames(sig_env)

# Scale arrows to NMDS range
arrow_scale <- 0.5 * min(apply(nmds_scores[, c("NMDS1","NMDS2")], 2, diff)) /
  max(sqrt(rowSums(sig_env[,1:2]^2)))
sig_env$NMDS1 <- sig_env$NMDS1 * arrow_scale
sig_env$NMDS2 <- sig_env$NMDS2 * arrow_scale
sig_env$x <- 0
sig_env$y <- 0

# --- Combine NMDS scores and metadata for plotting ---
nmds_plot_df <- nmds_scores %>%
  left_join(env_vars, by = "ID")

# Define colors
custom_colors <- RColorBrewer::brewer.pal(10, "Paired")
stress_val <- round(nmds$stress, 2)


# invert x axis
nmds_plot_df$NMDS1 <- -nmds_plot_df$NMDS1
nmds_plot_df$NMDS2 <- -nmds_plot_df$NMDS2


# Flip only arrows (both axes)
nmds_plot_df
sig_env$NMDS1 <- -sig_env$NMDS1
sig_env$NMDS2 <- -sig_env$NMDS2

# Order sites
nmds_plot_df$site <- factor(nmds_plot_df$site, levels = c("PM","MF","DS","CD","GR","MY1","SN","MY2","SF","BN"))

# --- Plot NMDS ---
S_FigS5 <- ggplot(nmds_plot_df, aes(x = NMDS1, y = NMDS2, color = site, shape = depth)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_segment(data = sig_env, aes(x = x, y = y, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", inherit.aes = FALSE) +
  geom_text_repel(data = sig_env, aes(x = NMDS1*1.1, y = NMDS2*1.1, label = Variable),
                  size = 4, color = "black", inherit.aes = FALSE) +
  labs(x = "NMDS1", y = "NMDS2", color = "Site", shape = "Depth") +
  scale_color_manual(values = custom_colors) +
  annotate("text", x = Inf, y = -Inf, label = paste("Stress =", stress_val),
           hjust = 1.1, vjust = -0.5, size = 4) +
  scale_shape_manual(values = c(16,17), labels = c("Lower soil layer", "Topsoil")) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14),
        plot.margin = margin(20,10,10,10)) +
  labs(title = "Bacterial taxonomic composition")
S_FigS5


ggsave("figures/S_FigS5_final_AlpSoils23_bac16srRNA_NMDS.png", plot = S_FigS5, height = 5, width = 6.5)



# --- PERMANOVA ---
# Remove taxa with zero counts (prevents NA in Hellinger)
otu_tab <- otu_tab[rowSums(otu_tab) > 0, ]

# Hellinger transformation
abund_hell <- sqrt(otu_tab / rowSums(otu_tab))

# Make sure samples match metadata
meta_tab <- meta_tab %>% filter(ID %in% colnames(abund_hell))

# Reorder abundance matrix to match metadata
abund_hell_sub <- abund_hell[, meta_tab$ID]

# Bray-Curtis distance (samples must be rows)
dist_mat <- vegdist(t(abund_hell_sub), method = "bray")

# PERMANOVA
adonis2(dist_mat ~ site, data = meta_tab, permutations = 999)

adonis2(dist_mat ~ depth, data = meta_tab, permutations = 999)

# depth nested in site
adonis2(dist_mat ~ site / depth, data = meta_tab, permutations = 999, by = "margin")


