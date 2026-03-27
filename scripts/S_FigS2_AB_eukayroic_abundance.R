###
### AlpineSoils - summer 23 - eukaryotic contigs
###


rm(list=ls(all=TRUE)) # removes everything



# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)

# coverage table
cov_abs <- read_tsv("data/Fig3_A_AlpineSoils23_contigs_coverage_abs.tsv")

# sample metadata
meta_tab <- read_tsv("data/Fig3_A_AlpineSoil23_contigs_metadata.tsv")

# tiara classification
tiara_tab <- read_tsv("data/S_FigS2_AlpineSoils23_eukaryotes_tiara.tsv")



## make kraken euk taxonomy table: 
# 1. Read the Kraken2 classification file (per contig) 
kraken_class <- read_tsv("data/S_FigS2_Alpsoil23_tax_contigs_kraken_euk.txt",
                         col_names = c("status", "contig_id", "taxid", "length", "kmer_hits"),
                         col_types = "ccicc")

# 2. Read the Kraken2 report (summary) 
kraken_report <- read_tsv("data/S_FigS2_Alpsoil23_kraken_report_euk.txt",
                          col_names = c("percent", "reads_clade", "reads_direct", "rank", "taxid", "name"),
                          col_types = cols(.default = "c")) %>%
  mutate(
    percent = as.numeric(percent),
    reads_clade = as.integer(reads_clade),
    reads_direct = as.integer(reads_direct),
    taxid = as.integer(taxid),
    name = trimws(name)
  )



# 3. Build lineage paths
# For each row, we keep track of the lineage up to that rank
lineage <- list()
current_lineage <- list()

for (i in seq_len(nrow(kraken_report))) {
  row <- kraken_report[i, ]
  rank <- row$rank
  name <- row$name
  taxid <- row$taxid
  
  # Update current lineage depending on rank
  if (rank == "R") current_lineage <- list(root = name)
  else if (rank == "K") current_lineage$kingdom <- name
  else if (rank == "P") current_lineage$phylum <- name
  else if (rank == "C") current_lineage$class <- name
  else if (rank == "O") current_lineage$order <- name
  else if (rank == "F") current_lineage$family <- name
  else if (rank == "G" | rank == "G1") current_lineage$genus <- name
  else if (rank == "S") current_lineage$species <- name
  
  # Save a snapshot of the lineage for this taxid
  lineage[[as.character(taxid)]] <- current_lineage
}

# 4. Convert lineage list to dataframe
lineage_df <- tibble(
  taxid = as.integer(names(lineage)),
  lineage = lineage
) %>%
  unnest_wider(lineage)

# 5. Merge with contigs 
tax_tab <- kraken_class %>%
  left_join(lineage_df, by = "taxid")



### clean up contigs (contamination) ------------------

tiara_euk <- tiara_tab %>%
  filter(class_fst_stage == "eukarya") %>%
  pull(sequence_id)

# Join with Kraken taxonomy
euk_tax <- tax_tab %>%
  filter(contig_id %in% tiara_euk)

# Count species to see if humans dominate
euk_tax %>%
  count(species, sort = TRUE) %>% head(20)



# Remove human contigs
contam_filtered <- tax_tab %>%
  filter(species != "Homo sapiens") %>%
  # Remove Illumina phiX sequences
  filter(!str_detect(species, "phiX")) %>%
  # Remove organelle sequences (chloroplast/mitochondria)
  filter(kingdom != "organelle")

# add domain categories
domain_map <- c(
  "Bacillati" = "Bacteria",
  "Pseudomonadati" = "Bacteria",
  "Fusobacteriati" = "Bacteria",
  "Nanobdellati" = "Bacteria",
  "Thermotogati" = "Archaea",
  "Thermoproteati" = "Archaea",
  "Methanobacteriati" = "Archaea",
  "Metazoa" = "Eukarya",
  "Fungi" = "Eukarya",
  "Viridiplantae" = "Eukarya",
  "Protista" = "Eukarya",
  "Bamfordvirae" = "Viruses",
  "Shotokuvirae" = "Viruses",
  "Heunggongvirae" = "Viruses",
  "Orthornavirae" = "Viruses"
)

contam_filtered <- contam_filtered %>%
  mutate(domain = domain_map[kingdom]) %>%
  filter(!is.na(domain)) # keep only known domains

# merge coverage and metadata
cov_long <- cov_abs %>%
  pivot_longer(-contig_id, names_to = "sample", values_to = "coverage") %>%
  left_join(contam_filtered %>% select(contig_id, kingdom, domain), by = "contig_id") %>%
  left_join(meta_tab %>% select(ID, site, depth), by = c("sample" = "ID")) %>%
  filter(!is.na(domain))

# domain level abundance per site
domain_abundance_site <- cov_long %>%
  group_by(site, domain) %>%
  summarise(total_coverage = sum(coverage, na.rm = TRUE), .groups = "drop") %>%
  arrange(site, desc(total_coverage))

# order sites
site_order <- c("PM","MF","DS","CD","GR","MY1","SN","MY2","SF","BN")
domain_abundance_site$site <- factor(domain_abundance_site$site, levels = site_order)

# order kingdoms
domain_order <- c("Bacteria","Archaea","Eukarya","Viruses")
domain_abundance_site$domain <- factor(domain_abundance_site$domain, levels = domain_order)

# plot domain composition per site
FigS2A <- ggplot(domain_abundance_site, aes(x = site, y = total_coverage, fill = domain)) +
  geom_col() +
  facet_wrap(domain ~ ., scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(x = "Site", y = "Domain abundances (Average cell number g⁻¹)", fill = "Domain")
FigS2A

# calculate ratio per site
prok_euk_ratio <- cov_long %>%
  group_by(site, domain) %>%
  summarise(total_coverage = sum(coverage), .groups = "drop") %>%
  pivot_wider(names_from = domain, values_from = total_coverage, values_fill = 0) %>%
  mutate(prok_to_euk = (Bacteria + Archaea) / Eukarya)

# order sites
site_order <- c("PM","MF","DS","CD","GR","MY1","SN","MY2","SF","BN")
prok_euk_ratio$site <- factor(prok_euk_ratio$site, levels = site_order)


FigS2B <- ggplot(prok_euk_ratio, aes(x = site, y = prok_to_euk)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(x = "Site", y = "Prokaryote:Eukaryote ratio")
FigS2B


### exporting final S figure S2AB ---------------
library(ggpubr)

S_FigS2_final <- ggarrange(FigS2A, FigS2B, ncol = 2, nrow = 1, labels = c('A', 'B'), legend = "none" )
S_FigS2_final

ggsave("figures/S_Figs2_final_AlpSoils23_eukaryotes.png", plot = S_FigS2_final, height = 5, width = 7)
