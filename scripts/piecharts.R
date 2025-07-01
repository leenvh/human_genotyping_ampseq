# === Load required libraries ===
library(ggplot2)
library(ggforce)
library(dplyr)
library(tidyr)
library(readr)
library(ggh4x)

# === Parse command-line arguments ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript plot_genotype_pies.R <summary_file_ackr1_dantu_hbb> <summary_file_g6pd> <output_file_hbb_pdf> <output_file_g6pd_pdf>")
}

file_hbb <- args[1]
file_g6pd <- args[2]
output_hbb_pdf <- args[3]
output_g6pd_pdf <- args[4]

# === Plot 1: ACKR1, Dantu, HBB ===

# Load data
df <- read_tsv(file_hbb)

# Gene order
gene_order <- c("ACKR1", "Dantu", "HBB")

# Format RS labels and gene ordering
df <- df %>%
  mutate(Gene = factor(Gene, levels = gene_order)) %>%
  arrange(Gene, RS_ID) %>%
  mutate(RS_LABEL = paste0("rs", RS_ID, " (", REF, "→", ALT, ")"),
         RS_LABEL = factor(RS_LABEL, levels = unique(RS_LABEL)))

# Reshape to long format for pie chart
df_long <- df %>%
  select(RS_LABEL, Gene, District, `Ref/Ref`, `Ref/Alt`, `Alt/Alt`) %>%
  pivot_longer(cols = c(`Ref/Ref`, `Ref/Alt`, `Alt/Alt`), names_to = "Genotype", values_to = "Count") %>%
  group_by(RS_LABEL, District) %>%
  mutate(Total = sum(Count),
         prop = Count / Total,
         angle = prop * 2 * pi,
         cumulative = cumsum(angle),
         start = lag(cumulative, default = 0),
         end = cumulative) %>%
  ungroup()

# Define colors
genotype_colors_hbb <- c("Ref/Ref" = "#1b9e77", "Ref/Alt" = "#d95f02", "Alt/Alt" = "#7570b3")

# Summary labels
df_total_labels <- df %>%
  group_by(Gene, RS_LABEL) %>%
  summarise(Total = sum(Total_Samples), .groups = "drop") %>%
  mutate(District = "Total samples (n=249)")

df_long_textonly <- df_total_labels %>%
  mutate(label = paste0("n = ", Total), x = 0, y = 0) %>%
  select(Gene, RS_LABEL, District, label, x, y)

# Plot
pdf(output_hbb_pdf, width = 12, height = 25)

ggplot() +
  geom_arc_bar(data = df_long,
               aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start, end = end, fill = Genotype),
               color = NA) +
  geom_text(data = df_long_textonly,
            aes(x = x, y = y, label = label),
            size = 3.5, fontface = "bold") +
  scale_fill_manual(values = genotype_colors_hbb) +
  coord_fixed() +
  facet_grid2(Gene + RS_LABEL ~ District, switch = "y") +
  theme_void() +
  theme(strip.text.x = element_text(size = 10),
        panel.spacing = unit(1.2, "lines"),
        legend.position = "right") +
  labs(title = "Genotype Proportions per SNP and District", fill = "Genotype")

dev.off()


# === Plot 2: G6PD ===

# Load data
df <- read_tsv(file_g6pd)

# Format RS labels
df <- df %>%
  arrange(Gene, RS_ID) %>%
  mutate(RS_LABEL = paste0("rs", RS_ID, " (", REF, "→", ALT, ")"),
         RS_LABEL = factor(RS_LABEL, levels = unique(RS_LABEL)))

# Reshape to long format
df_long <- df %>%
  select(RS_LABEL, Gene, District, Male_Hemizygote, Female_Het, Female_Hom, WT) %>%
  pivot_longer(cols = c(Male_Hemizygote, Female_Het, Female_Hom, WT), names_to = "Genotype", values_to = "Count") %>%
  group_by(RS_LABEL, District) %>%
  mutate(Total = sum(Count),
         prop = Count / Total,
         angle = prop * 2 * pi,
         cumulative = cumsum(angle),
         start = lag(cumulative, default = 0),
         end = cumulative) %>%
  ungroup()

# Define colors
genotype_colors_g6pd <- c("WT" = "#1b9e77", "Male_Hemizygote" = "red", "Female_Het" = "#d95f02", "Female_Hom" = "#7570b3")

# Summary labels
df_total_labels <- df %>%
  group_by(Gene, RS_LABEL) %>%
  summarise(Total = sum(Total_Samples), .groups = "drop") %>%
  mutate(District = "Total samples (n=249)")

df_long_textonly <- df_total_labels %>%
  mutate(label = paste0("n = ", Total), x = 0, y = 0) %>%
  select(Gene, RS_LABEL, District, label, x, y)

# Plot
pdf(output_g6pd_pdf, width = 12, height = 25)

ggplot() +
  geom_arc_bar(data = df_long,
               aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start, end = end, fill = Genotype),
               color = NA) +
  geom_text(data = df_long_textonly,
            aes(x = x, y = y, label = label),
            size = 3.5, fontface = "bold") +
  scale_fill_manual(values = genotype_colors_g6pd) +
  coord_fixed() +
  facet_grid2(Gene + RS_LABEL ~ District, switch = "y") +
  theme_void() +
  theme(strip.text.x = element_text(size = 10),
        panel.spacing = unit(1.2, "lines"),
        legend.position = "right") +
  labs(title = "G6PD Genotype Proportions per SNP and District", fill = "Genotype")

dev.off()
