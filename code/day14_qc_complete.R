# ==============================================================================#
#                        FULL QC PIPELINE — DAY 14 SAMPLE                      #
#                                                                               #
# This file contains the complete, explicit QC workflow for the Day14 sample.  #
# It mirrors the Day0 QC pipeline exactly.                                     #
# Original source code only contained a comment placeholder for this section.  #
# ==============================================================================#


# ==============================================================================#
#                           READ DAY14 SAMPLE                                   #
# ==============================================================================#

# Set working directory to the Day14 sample folder
# ── UPDATE THIS PATH ──────────────────────────────────────────────────────────
setwd("/Volumes/BCP-ALL/presentations/single cell/GSE142564_RAW/Day14")
# ─────────────────────────────────────────────────────────────────────────────

# List all files in the directory to confirm data presence
list.files()

# Read 10X Genomics matrix files
# Reads: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
data_14 <- Read10X("/Volumes/BCP-ALL/presentations/single cell/GSE142564_RAW/Day14/")

# Create Seurat object from raw count matrix
Day14 <- CreateSeuratObject(counts = data_14)

# Print Seurat object summary
Day14

# Filter cells with fewer than 200 detected genes
# WHY 200: removes the most obvious empty droplets and dead cell fragments
Day14 <- subset(Day14, subset = nFeature_RNA >= 200)

# Display filtered object summary
Day14


# ==============================================================================#
#                                QUALITY CONTROL                                #
# ==============================================================================#

# Calculate percentage of mitochondrial transcripts per cell
# Mouse mitochondrial genes have the prefix ^mt- (lowercase)
# Human mitochondrial genes have the prefix ^MT- (uppercase)
Day14$mitoRatio <- PercentageFeatureSet(object = Day14, pattern = "^mt-")

# Convert mitochondrial percentage into fraction (0–1 scale)
Day14$mitoRatio <- Day14@meta.data$mitoRatio / 100

# Verify mitochondrial genes are present in this dataset
# If this returns integer(0), the dataset may use a different naming convention
grep("^mt-", rownames(Day14))

# Calculate gene complexity score:
# log10(genes) / log10(UMIs) — high = transcriptionally complex; low = debris or RBC
Day14$log10GenesPerUMI <- log10(Day14$nFeature_RNA) / log10(Day14$nCount_RNA)


# ==============================================================================#
#                           INITIAL QC VISUALIZATION                            #
# ==============================================================================#

# Visualize QC metrics BEFORE filtering
# Always look at the distribution first — thresholds should be data-driven
VlnPlot(
  Day14,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol     = 3,
  layer    = "counts"
)

# Remove low-quality cells and potential doublets
# Upper limits:
#   nFeature_RNA < 4000  — removes doublets (two cells ~ twice the genes)
#   nCount_RNA < 20000   — removes cells with anomalously high sequencing depth
#   mitoRatio < 0.05     — removes dying/damaged cells (>5% mitochondrial reads)
Day14 <- subset(
  Day14,
  nFeature_RNA < 4000 &
    nCount_RNA  < 20000 &
    mitoRatio   < 0.05
)

# Re-visualize QC metrics after filtering
VlnPlot(
  Day14,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol     = 3,
  layer    = "counts"
)


# ==============================================================================#
#                              ORGANIZE METADATA                                #
# ==============================================================================#

# Add sample identity
Day14$sample <- "Day14"

# Extract metadata dataframe
metadata_14 <- Day14@meta.data

# Store cell barcodes in metadata (they are rownames by default, not accessible after dplyr)
metadata_14$cells <- rownames(metadata_14)

# Preserve original counts and feature numbers before any renaming
metadata_14$orig_counts   <- metadata_14$nCount_RNA
metadata_14$orig_features <- metadata_14$nFeature_RNA

# Rename metadata columns for clarity
metadata_14 <- metadata_14 %>%
  dplyr::rename(
    sample = sample,
    nUMI   = nCount_RNA,
    nGene  = nFeature_RNA
  )

# Write updated metadata back into Seurat object
Day14@meta.data <- metadata_14


# ==============================================================================#
#                         NUMBER OF CELLS PER SAMPLE                            #
# ==============================================================================#

# Plot total number of cells in Day14 after filtering
metadata_14 %>%
  ggplot(aes(x = sample, fill = sample)) +
  geom_bar(width = 0.7, color = "black", show.legend = FALSE) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y     = element_text(color = "black"),
    axis.title.y    = element_text(face = "bold"),
    plot.title      = element_text(face = "bold", size = 16, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  ) +
  labs(
    title = "Number of Cells per Sample — Day 14",
    x = "",
    y = "Cell Count"
  )


# ==============================================================================#
#                           SEQUENCING DEPTH QC                                 #
# ==============================================================================#

# Plot transcript count (UMI) distribution per cell
# Log10 scale because UMI counts span several orders of magnitude
# Dashed line at nUMI = 2000 marks the suggested lower reference threshold
# The final lower threshold applied is nUMI > 3200 (below in Round 2)
metadata_14 %>%
  ggplot(aes(x = nUMI, color = sample, fill = sample)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_x_log10() +
  geom_vline(
    xintercept = 2000,
    linetype   = "dashed",
    color      = "black",
    linewidth  = 1
  ) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title    = element_blank(),
    legend.position = "top",
    axis.text       = element_text(color = "black"),
    axis.title.y    = element_text(face = "bold"),
    axis.title.x    = element_text(face = "bold"),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 16)
  ) +
  labs(
    title = "Transcript Counts per Cell — Day 14",
    x     = "nUMI (log10 scale)",
    y     = "Cell Density"
  )


# ==============================================================================#
#                         NUMBER OF GENES PER CELL                              #
# ==============================================================================#

# Plot detected gene distribution per cell
# Dashed line at nGene = 1000 marks the suggested lower reference threshold
metadata_14 %>%
  ggplot(aes(x = nGene, color = sample, fill = sample)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_x_log10() +
  geom_vline(
    xintercept = 1000,
    linetype   = "dashed",
    color      = "black",
    linewidth  = 1
  ) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title    = element_blank(),
    legend.position = "top",
    axis.text       = element_text(color = "black"),
    axis.title.y    = element_text(face = "bold"),
    axis.title.x    = element_text(face = "bold"),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 16),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  ) +
  labs(
    title = "Gene Counts per Cell — Day 14",
    x     = "nGene (log10 scale)",
    y     = "Cell Density"
  )


# ==============================================================================#
#                      MITOCHONDRIAL CONTENT DISTRIBUTION                       #
# ==============================================================================#

# Plot mitochondrial transcript ratio distribution
# Most healthy immune cells cluster below 0.05 (5%)
# A long right tail indicates poor sample quality (cell membrane damage)
metadata_14 %>%
  ggplot(aes(x = mitoRatio, color = sample, fill = sample)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_x_log10() +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text       = element_text(color = "black"),
    axis.title      = element_text(face = "bold"),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 16)
  ) +
  labs(
    title = "Distribution of Mitochondrial Ratio — Day 14",
    x     = "Mitochondrial Ratio",
    y     = "Cell Density"
  )


# ==============================================================================#
#                    RELATIONSHIP BETWEEN QC METRICS                            #
# ==============================================================================#

# Scatter plot: UMI counts vs. gene counts, coloured by mitochondrial ratio
# This is the most informative single QC figure.
#
# WHAT TO LOOK FOR:
#   - Main cloud (upper right, light grey): healthy cells — high UMI, high genes, low mito
#   - Lower-left outliers (dark colour): damaged cells — low UMI, low genes, HIGH mito
#   - Points below the trend line: low-complexity cells (RBCs, empty droplets)
#   - Points above the cloud (upper outliers): potential doublets
#
# Red lines mark the Round 2 thresholds:
#   Vertical   (nUMI = 3200): lower UMI cutoff
#   Horizontal (nGene = 1000): lower gene cutoff
metadata_14 %>%
  ggplot(aes(x = nUMI, y = nGene, color = mitoRatio)) +

  geom_point(alpha = 0.6, size = 1.5) +

  # Add regression trend line — helps identify cells deviating from expected
  stat_smooth(
    method   = "lm",
    color    = "blue",
    se       = TRUE,
    linetype = "dashed"
  ) +

  scale_color_gradient(low = "gray90", high = "black") +

  scale_x_log10() +
  scale_y_log10() +

  # Vertical line: nUMI lower threshold
  geom_vline(
    xintercept = 3200,
    linetype   = "dotted",
    color      = "red",
    linewidth  = 0.8
  ) +

  # Horizontal line: nGene lower threshold
  geom_hline(
    yintercept = 1000,
    linetype   = "dotted",
    color      = "red",
    linewidth  = 0.8
  ) +

  theme_minimal(base_size = 14) +

  theme(
    legend.position = "right",
    axis.text       = element_text(color = "black"),
    axis.title      = element_text(face = "bold"),
    strip.text      = element_text(face = "bold", size = 12)
  ) +

  labs(
    title  = "UMI vs. Gene Counts Coloured by Mitochondrial Ratio — Day 14",
    x      = "nUMI (log10 scale)",
    y      = "nGene (log10 scale)",
    color  = "Mito Ratio"
  ) +

  facet_wrap(~sample)


# ==============================================================================#
#                       GENES PER UMI COMPLEXITY SCORE                          #
# ==============================================================================#

# Plot transcriptional complexity distribution
# Threshold at 0.82 removes low-complexity cells (RBCs, empty droplets)
# Values below ~0.75 are almost always debris or transcriptionally simple cells
metadata_14 %>%
  ggplot(aes(x = log10GenesPerUMI, color = sample, fill = sample)) +

  geom_density(alpha = 0.3, linewidth = 1) +

  # Complexity threshold
  geom_vline(
    xintercept = 0.82,
    linetype   = "dashed",
    color      = "black",
    linewidth  = 1
  ) +

  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +

  theme_minimal(base_size = 14) +

  theme(
    legend.position = "top",
    axis.text       = element_text(color = "black"),
    axis.title      = element_text(face = "bold"),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 16)
  ) +

  labs(
    title = "Transcriptional Complexity Distribution — Day 14",
    x     = "log10(Genes per UMI)",
    y     = "Cell Density"
  )


# Apply final QC filtering (Round 2) based on:
#   - Complexity score: log10GenesPerUMI > 0.82
#   - Minimum gene count: nGene > 900
#   - Minimum UMI count: nUMI > 3200
#
# These lower-bound filters remove what Round 1 missed:
# low-complexity cells, cells with insufficient sequencing depth,
# and empty droplets that happen to express a few genes at high counts.
Day14 <- subset(
  x = Day14,
  subset = (
    log10GenesPerUMI > 0.82
  ) &
    (nGene > 900) &
    (nUMI > 3200)
)

# Display final dimensions after all filtering
dim(Day14)

# Final QC violin plot — distributions should now be clean and unimodal
VlnPlot(
  Day14,
  features = c("nGene", "nUMI", "mitoRatio"),
  pt.size  = 0
)


# ==============================================================================#
#                         QC COMPARISON: DAY0 vs DAY14                          #
# ==============================================================================#

# After QC-ing both samples independently, compare the key metrics.
# Large systematic differences (e.g., Day14 has 2× the median UMI of Day0)
# indicate a sequencing depth batch effect that must be corrected by integration.

cat("=== Post-QC Summary ===\n")
cat(sprintf("Day0  | Cells: %d | Median nUMI: %.0f | Median nGene: %.0f\n",
            ncol(Day0),  median(Day0$nUMI),  median(Day0$nGene)))
cat(sprintf("Day14 | Cells: %d | Median nUMI: %.0f | Median nGene: %.0f\n",
            ncol(Day14), median(Day14$nUMI), median(Day14$nGene)))

# Visual comparison of QC metrics side-by-side (requires both Day0 and Day14 to be loaded)
# Combine metadata for comparison plot
combined_meta <- bind_rows(
  Day0@meta.data,
  Day14@meta.data
)

# Side-by-side UMI distribution comparison
combined_meta %>%
  ggplot(aes(x = nUMI, color = sample, fill = sample)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  scale_x_log10() +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text       = element_text(color = "black"),
    axis.title      = element_text(face = "bold"),
    plot.title      = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "UMI Distribution Comparison: Day 0 vs. Day 14 (post-QC)",
    x     = "nUMI (log10 scale)",
    y     = "Cell Density"
  )

# Side-by-side gene count comparison
combined_meta %>%
  ggplot(aes(x = nGene, color = sample, fill = sample)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  scale_x_log10() +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text       = element_text(color = "black"),
    axis.title      = element_text(face = "bold"),
    plot.title      = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Gene Count Distribution Comparison: Day 0 vs. Day 14 (post-QC)",
    x     = "nGene (log10 scale)",
    y     = "Cell Density"
  )


# ==============================================================================#
#                           PROCEED TO SECTION: MERGE                           #
# ==============================================================================#
# After running this script, both Day0 and Day14 are QC-filtered and ready.
# Continue with code/analysis_code.txt from the "MERGE DATASETS" section.
# ==============================================================================#
