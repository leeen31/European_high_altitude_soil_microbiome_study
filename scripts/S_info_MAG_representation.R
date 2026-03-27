###
### AlpSoils23 - MAG representation
###

rm(list=ls(all=TRUE)) # removes everything


# read libraries
library(tidyverse) 
library(tidyr)
library(dplyr)



# read reads count table
counts_tab <- read.table('data/S_info_AlpineSoils23_MAGrepresentation_read_counts.tsv', sep='\t', header = TRUE)

# change samples names
# remove appendixes to sample name
counts_tab$sample <- sub("_S[0-9]+.*", "", counts_tab$sample)
counts_tab$sample <- sub(".*fa\\.", "", counts_tab$sample)


counts_tab <- counts_tab %>%
  mutate(sample = recode (sample,
                          "230713_LA_d" = "L_a_d",
                          "230713_LA_t" = "L_a_t",
                          "230713_LB_d" = "L_b_d",
                          "230713_LB_t" = "L_b_t",
                          "230713_LC_d" = "L_c_d", 
                          "230713_LC_t" = "L_c_t", 
                          "230713_LD_d" = "L_d_d",
                          "230713_LD_t" = "L_d_t",
                          "230713_LE_d" = "L_e_d",
                          "230713_LE_t" = "L_e_t"))

# Replace all '-' with '_' in column names
counts_tab$sample <- gsub("-", "_", counts_tab$sample)

# Remove negative control
counts_tab <- counts_tab %>% select(-nc1) 




# read contig-in-bins file
# Path to your concatenated MAGs FASTA
fa_file <- "S_info_AlpineSoils23_MAGrepresentation_contigs_in_bins.fa"

# Read all lines starting with ">"
contig_headers <- readLines(fa_file)
contig_headers <- contig_headers[grepl("^>", contig_headers)]

# Remove the ">" to get clean contig IDs
contig_ids_in_MAGs <- sub("^>", "", contig_headers)





# add column to counts_tab
counts_tab <- counts_tab %>%
  mutate(in_MAG = ifelse(contig %in% contig_ids_in_MAGs, "yes", "no"))






### calculate % representation

# total reads
total_reads <- sum(counts_tab$read_count)
total_reads


# reads in MAGs
reads_in_MAGs <- counts_tab %>%
  filter(in_MAG == "yes") %>%
  summarize(reads = sum(read_count)) %>%
  pull(reads)

percent_in_MAGs <- reads_in_MAGs / total_reads * 100
percent_in_MAGs



# reads in contigs
# reads mapped to any contig (exclude "unmapped")
reads_mapped <- sum(counts_tab$read_count[counts_tab$contig != "unmapped"])

# fraction of reads represented by contigs
percent_reads_covered <- reads_mapped / total_reads * 100
percent_reads_covered


