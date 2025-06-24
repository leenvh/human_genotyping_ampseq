# Load necessary libraries
library(ggplot2)
library(readr)
library(dplyr)

# Load data
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Please provide the input file path and output directory as command-line arguments.")
}
input_file <- args[1]
output_dir <- args[2]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


df <- read_tsv(input_file, col_types = cols())

# Filter for rows with valid RS IDs
df_filtered <- df %>%
  filter(RS != ".",
         SAMPLE != "unassigned",
         !grepl("^NC", SAMPLE)) %>%
  mutate(DP = as.numeric(DP))

# Add 'rs' prefix to all RS IDs
df_filtered$RS <- paste0("rs", df_filtered$RS)

# Re-factor RS to preserve order after sorting by chromosome
df_filtered$RS <- factor(df_filtered$RS, levels = unique(df_filtered$RS))

# Check if data exists
if (nrow(df_filtered) == 0) {
  warning("No SNPs with RS IDs found.")
} else {
  # Map chromosome numbers to gene names
  chrom_map <- c(
    "1" = "ACKR1",
    "4" = "Dantu",
    "11" = "HBB",
    "X" = "G6PD"
  )
  
  # Apply mapping
  df_filtered$GENE <- chrom_map[as.character(df_filtered$CHROM)]
  
  # Define chromosome order for sorting
  chrom_order <- c("1", "4", "11", "X")
  
  # Sort RS by chromosome order
  df_filtered <- df_filtered %>%
    mutate(CHROM = factor(CHROM, levels = chrom_order)) %>%
    arrange(CHROM, POS)
  
  # Re-factor RS in the order of appearance
  df_filtered$RS <- factor(df_filtered$RS, levels = unique(df_filtered$RS))
  
  # Plot
  p1<-ggplot(df_filtered, aes(x = RS, y = DP, fill = GENE)) +
    geom_boxplot() +
    labs(
      title = "Read Depth per rsID",
      x = "SNP (rsID)",
      y = "Read Depth (DP)",
      fill = "Gene"
    ) +
    theme_minimal() +
    scale_y_sqrt(breaks = c(0, 50, 100,500, 1000, 2000, 4000,8000,12000,16000)) +
    geom_hline(yintercept = 50, linetype = "dashed")+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
    # Save plot
  ggsave(filename = file.path(output_dir, "snp_read_depth_boxplot.pdf"), plot = p1, width = 10, height = 6)
}



# Graph for negative controls
NC <- df %>%
  filter(grepl("^NC", SAMPLE))


# Check if data exists
if (nrow(NC) == 0) {
  warning("Negative controls not found.")
} else {
  # Map chromosome numbers to gene names
  chrom_map <- c(
    "1" = "ACKR1",
    "4" = "Dantu",
    "11" = "HBB",
    "X" = "G6PD"
  )
  
  # Apply mapping
  NC$GENE <- chrom_map[as.character(NC$CHROM)]
  
  # Define chromosome order for sorting
  chrom_order <- c("1", "4", "11", "X")
  
  # Sort RS by chromosome order
  NC <- NC %>%
    mutate(CHROM = factor(CHROM, levels = chrom_order)) %>%
    arrange(CHROM, POS)
  
  # Plot
  p2<-ggplot(NC, aes(x = SAMPLE, y = DP, fill = GENE)) +
    geom_boxplot() +
    labs(
      title = "Read Depth per Negative control",
      y = "Read Depth (DP)",
      fill = "Gene"
    ) +
    theme_minimal() +
    scale_y_sqrt(limits= c(0,100),breaks = c(0,10,20,30,40, 50, 100)) +
    geom_hline(yintercept = 50, linetype = "dashed")+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
  ggsave(filename = file.path(output_dir, "negative_control_read_depth_boxplot.pdf"), plot = p2, width = 10, height = 6)
}

# Clean up unwanted Rplots.pdf if created
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")