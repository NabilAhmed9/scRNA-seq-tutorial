# Code Explanation: Section-by-Section Walkthrough

> **Note:** This repository preserves the original analysis code and adds educational documentation based on the accompanying presentation.

**Source file:** `code/analysis_code.txt`  
**Language:** R  
**Analysis framework:** Seurat v5 / Monocle3 / CellChat  

---

## Section 1 — Library Loading

### Code Purpose
Loads all R packages required for the full scRNA-seq pipeline.

```r
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
```

### What the Code Does
- `SingleCellExperiment`: provides the SCE data structure used by Bioconductor tools and is a backend for Seurat object conversion
- `Seurat`: the primary analysis framework — handles data loading, QC, normalisation, integration, clustering, and annotation
- `tidyverse`: provides `dplyr`, `ggplot2`, `stringr`, and pipe operators (`%>%`) for data manipulation and visualisation
- `Matrix`: enables sparse matrix representation of count data, essential for memory efficiency given the zero-inflated nature of scRNA-seq data
- `scales`: aesthetic formatting for ggplot2 axes
- `cowplot`: publication-quality multi-panel plot layouts
- `RCurl`: URL-based data downloading (not used explicitly in the preserved code but available)

### Expected Input
None — package imports only.

### Expected Output
R session with all packages loaded and functions available.

### Interpretation
The choice to use Seurat as the central framework is deliberate — it has the most mature ecosystem, best documentation, and broadest community adoption for this type of analysis. The `Matrix` package is critical: without sparse matrix support, a 20,000 gene × 10,000 cell matrix would require ~1.6 GB in dense format vs. ~50 MB in sparse format.

---

## Section 2 — Data Loading and Initial Object Creation (Day 0)

### Code Purpose
Load 10X Genomics CellRanger output for the Day 0 sample and create a Seurat object with basic pre-filtering.

```r
setwd("/Volumes/BCP-ALL/presentations/single cell/GSE142564_RAW/Day0")
data <- Read10X("/Volumes/BCP-ALL/presentations/single cell/GSE142564_RAW/Day0/")
Day0 <- CreateSeuratObject(counts = data)
Day0 <- subset(Day0, subset = nFeature_RNA >= 200)
```

### What the Code Does
1. Sets the working directory to the Day0 data folder (note: **path is system-specific** — users must update this)
2. `Read10X()` reads three CellRanger output files simultaneously: `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`
3. `CreateSeuratObject()` wraps the count matrix into Seurat's data structure, which adds default QC metrics (`nCount_RNA` = total UMIs, `nFeature_RNA` = genes detected per cell)
4. Initial filtering removes cells with fewer than 200 detected genes — this eliminates empty droplets and debris

### Expected Input
10X CellRanger output directory containing:
- `barcodes.tsv.gz`
- `features.tsv.gz`
- `matrix.mtx.gz`

### Expected Output
A filtered `Day0` Seurat object with cells having ≥200 detected genes.

### Interpretation
The 200-gene minimum is a conservative initial filter. The downstream QC steps will refine this further. **Important:** The working directory path `/Volumes/BCP-ALL/...` is specific to the original author's machine. Users must set their own paths.

---

## Section 3 — Quality Control Metric Calculation

### Code Purpose
Calculate three core QC metrics for each cell.

```r
Day0$mitoRatio <- PercentageFeatureSet(object = Day0, pattern = "^mt-") / 100
Day0$log10GenesPerUMI <- log10(Day0$nFeature_RNA) / log10(Day0$nCount_RNA)
```

### What the Code Does
1. **Mitochondrial fraction:** `PercentageFeatureSet()` counts UMIs mapping to genes with names starting with `mt-` (mouse mitochondrial gene convention — human uses `MT-`). Dividing by 100 converts percentage to fraction.
2. **Complexity score:** `log10GenesPerUMI` captures how many unique genes are detected relative to total transcript count. A healthy cell should have a high ratio (many distinct genes per UMI).

### Expected Input
`Day0` Seurat object with raw counts

### Expected Output
Two new columns in `Day0@meta.data`: `mitoRatio` and `log10GenesPerUMI`

### Interpretation
| Metric | Low Value Means | High Value Means |
|---|---|---|
| `nCount_RNA` | Empty droplet or low-quality | Normal or potential doublet |
| `nFeature_RNA` | Debris or dying cell | Normal or potential doublet |
| `mitoRatio` | Healthy cell | Damaged/apoptotic cell |
| `log10GenesPerUMI` | Homogeneous cell (RBC, debris) | Transcriptionally complex cell |

Note: `grep("^mt-", rownames(Day0))` verifies that mitochondrial genes are present in this dataset — a critical sanity check, as some count matrices lack them.

---

## Section 4 — QC Visualisation (Pre-filtering)

### Code Purpose
Generate violin plots to visually inspect QC metric distributions before applying cutoffs.

```r
VlnPlot(Day0, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3, layer = "counts")
```

### What the Code Does
Creates a 3-panel violin plot showing the distribution of detected genes, UMI counts, and mitochondrial ratio across all cells in the Day0 sample.

### Expected Output
Three-panel violin plot saved to the active R graphics device. In RStudio, this appears in the Plots pane.

### Interpretation
Before filtering, the violin plots should show a right-skewed distribution for nCount and nFeature (few very high cells = potential doublets) and a left-skewed mitoRatio (most cells low, some high = dying cells).

---

## Section 5 — Multi-metric Cell Filtering

### Code Purpose
Apply conservative thresholds to remove low-quality cells, potential doublets, and dying cells.

```r
Day0 = subset(
  Day0,
  nFeature_RNA < 4000 &
    nCount_RNA < 20000 &
    mitoRatio < .05
)
```

### What the Code Does
Applies three simultaneous upper-bound filters:
- `nFeature_RNA < 4000`: removes potential doublets (two cells captured together have abnormally high gene counts)
- `nCount_RNA < 20000`: removes cells with aberrantly high UMI counts
- `mitoRatio < 0.05`: removes cells with >5% mitochondrial reads (damaged/dying cells)

### Expected Input
Pre-filtered `Day0` object with QC metrics

### Expected Output
`Day0` object with a reduced cell count — low-quality cells removed

### Interpretation
These thresholds should always be **dataset-specific**. The code does not apply lower bounds at this stage (those are applied later after the complexity score calculation). The 5% mitochondrial threshold is appropriate for many immune cell datasets; highly metabolic cells (e.g., cardiomyocytes) may require a higher threshold.

---

## Section 6 — Metadata Organisation

### Code Purpose
Extract, organise, and relabel the cell metadata dataframe with cleaner column names.

```r
metadata <- Day0@meta.data
metadata$cells <- rownames(metadata)
metadata$orig_counts = metadata$nCount_RNA
metadata$orig_features = metadata$nFeature_RNA
metadata <- metadata %>%
  dplyr::rename(sample = sample, nUMI = nCount_RNA, nGene = nFeature_RNA)
Day0@meta.data <- metadata
```

### What the Code Does
1. Extracts the metadata slot as a standard R dataframe
2. Preserves original column names as backup columns
3. Renames columns to more readable labels (`nUMI`, `nGene`)
4. Writes the modified dataframe back to the Seurat object

### Expected Input/Output
`Day0` Seurat object with metadata slots updated

### Interpretation
This is a housekeeping step. Preserving `orig_counts` and `orig_features` is good practice — it maintains a record of the pre-cleaning values for auditing purposes.

---

## Section 7 — Detailed QC Visualisations

### Code Purpose
Generate publication-quality diagnostic plots to validate filtering decisions.

```r
# UMI density plot with lower cutoff at 2000
metadata %>% ggplot(aes(x = nUMI, color = sample, fill = sample)) +
  geom_density(alpha = 0.3) + scale_x_log10() +
  geom_vline(xintercept = 2000, linetype = "dashed")

# Gene count density plot with lower cutoff at 1000
metadata %>% ggplot(aes(x = nGene, ...)) + geom_density() + scale_x_log10() +
  geom_vline(xintercept = 1000)

# UMI vs Gene scatter coloured by mitochondrial ratio
metadata %>% ggplot(aes(x = nUMI, y = nGene, color = mitoRatio)) +
  geom_point() + stat_smooth(method = "lm") +
  geom_vline(xintercept = 3200) + geom_hline(yintercept = 1000) +
  facet_wrap(~sample)

# Complexity score density
metadata %>% ggplot(aes(x = log10GenesPerUMI)) +
  geom_density() + geom_vline(xintercept = 0.82)
```

### What the Code Does
- Density plots show the per-sample distribution of each metric on log10 scale
- Dashed vertical lines indicate the proposed filtering thresholds
- The scatter plot coloured by `mitoRatio` reveals whether cells failing on UMI/gene counts also have high mitochondrial content — a powerful joint QC diagnostic
- `facet_wrap(~sample)` generates separate panels per sample for direct comparison

### Interpretation
In a high-quality scRNA-seq dataset, the UMI and gene count distributions should be approximately log-normal (linear on a log10 axis). The UMI vs. gene scatter should show a near-linear relationship; cells far from this line (low genes for their UMI count) are low-complexity and likely empty droplets or red blood cells.

---

## Section 8 — Final Complexity-Based Filtering

### Code Purpose
Apply the final round of QC filtering using the complexity score and refined thresholds.

```r
Day0 <- subset(
  x = Day0,
  subset = (log10GenesPerUMI > 0.82) & (nGene > 900) & (nUMI > 3200)
)
dim(Day0)
```

### What the Code Does
Applies three lower-bound filters simultaneously:
- `log10GenesPerUMI > 0.82`: removes low-complexity cells
- `nGene > 900`: minimum gene diversity threshold (stricter than initial 200-gene minimum)
- `nUMI > 3200`: minimum sequencing depth threshold

### Expected Output
Final filtered Day0 object. `dim()` reports remaining cells × genes.

### Interpretation
The complexity threshold `0.82` is a log10-scale measure. This corresponds approximately to a ratio where cells are expressing at least a baseline diversity of distinct gene transcripts. Red blood cells, empty droplets, and debris consistently fall below this threshold.

---

## Section 9 — Day 14 Sample QC (Implicit)

The code includes a section header indicating that the **identical QC workflow is applied to the Day14 sample**:

```r
# Repeat same workflow for Day14 sample
# 1. Data loading
# 2. Seurat object creation
# 3. QC metric calculation
# 4. Filtering
# 5. Metadata organization
# 6. QC visualization
```

**Note:** The explicit Day14 code is not reproduced in the source file — it is implied by this comment. Learners should duplicate all Day0 QC steps, substituting `Day14` for `Day0` and updating the input path accordingly.

---

## Section 10 — Merging Datasets

### Code Purpose
Combine Day0 and Day14 Seurat objects into a single merged object.

```r
Day0 <- RenameCells(Day0, add.cell.id = "Day0")
Day14 <- RenameCells(Day14, add.cell.id = "Day14")
merged <- merge(x = Day0, y = c(Day14))
```

### What the Code Does
1. `RenameCells()` prepends a sample prefix to each cell barcode (e.g., `ACGTACGT-1` → `Day0_ACGTACGT-1`), preventing barcode collision between samples
2. `merge()` concatenates the two Seurat objects into a single merged object, retaining all metadata and the `sample` column for downstream splitting

### Expected Output
A `merged` Seurat object containing cells from both time points, identifiable by their barcode prefix.

---

## Section 11 — BCR/TCR Gene Removal

### Code Purpose
Remove immunoglobulin (Ig) and T-cell receptor (TCR) variable genes before feature selection.

```r
merged <- merged[!grepl("^Ig[hkl]v", rownames(merged)), ]  # Ig variable
merged <- merged[!grepl("^Ig[hkl]j", rownames(merged)), ]  # Ig joining
merged <- merged[!grepl("^Ig[kl]c",  rownames(merged)), ]  # Ig constant
merged <- merged[!grepl("^Igh[adegm]", rownames(merged)),] # Heavy chain isotypes
merged <- merged[!grepl("^Tr[abdg][vjc]", rownames(merged)),] # TCR genes
```

### What the Code Does
Uses regular expression pattern matching to identify and remove rows (genes) from the count matrix that correspond to Ig/TCR receptor genes.

### Biological Rationale
In B-cell or T-cell rich datasets, Ig/TCR variable genes are clonally expressed — each B cell or T cell expresses a unique V(D)J recombination product. These genes will appear as highly variable and will drive PC1 of the PCA, causing cells to cluster by clonotype rather than transcriptional state. Removing them forces the analysis to focus on functional transcriptional programmes.

### Expected Input/Output
`merged` object with reduced gene count (Ig/TCR genes removed from rows)

---

## Section 12 — Sparse Gene Filtering

### Code Purpose
Remove genes expressed in fewer than 10 cells to reduce noise and dimensionality.

```r
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])
counts <- LayerData(merged, assay = "RNA", layer = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
dual_counts <- counts[keep_genes, ]
merged <- CreateSeuratObject(dual_counts, meta.data = merged@meta.data)
```

### What the Code Does
1. `JoinLayers()` collapses Seurat v5 layer structure into a single count matrix
2. Creates a logical (TRUE/FALSE) matrix of expressed vs. not expressed
3. Sums across cells to find genes expressed in ≥10 cells
4. Rebuilds the Seurat object with the filtered count matrix, preserving metadata

### Expected Output
`merged` object with a reduced feature (gene) count

### Interpretation
This step is critical for reducing the computational burden of downstream steps. Genes expressed in fewer than 10 cells contribute primarily noise to variable gene selection and PCA.

---

## Section 13 — Cell Cycle Scoring

### Code Purpose
Assign each cell a cell-cycle phase score (S and G2M) to enable regression of cell-cycle effects.

```r
merged <- NormalizeData(merged)
s.genes <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)
merged <- CellCycleScoring(merged, g2m.features = g2m.genes, s.features = s.genes, set.ident = TRUE)
RidgePlot(merged, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)
```

### What the Code Does
1. `NormalizeData()` applies log-normalisation to enable gene expression comparison
2. `str_to_title()` converts the built-in `cc.genes` list (stored in all-caps for human) to title case for mouse (e.g., `PCNA` → `Pcna`)
3. `CellCycleScoring()` scores each cell using Seurat's built-in S-phase and G2M-phase gene signatures
4. `RidgePlot()` visualises expression of canonical cell-cycle markers to validate the scoring

### Expected Output
Two new metadata columns: `S.Score` and `G2M.Score`; cells assigned to Phase (G1, S, or G2M)

### Interpretation
The `str_to_title()` conversion is an important mouse-vs-human data handling step. Mouse gene symbols use title case (first letter capitalised); human uses all-caps. Using the wrong case will result in the scoring failing silently.

---

## Section 14 — Variable Features, Scaling, Initial PCA

### Code Purpose
Identify highly variable genes and perform an initial PCA for cell-cycle visualisation.

```r
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged)
DimPlot(merged, reduction = "pca", group.by = "Phase", split.by = "Phase")
```

### What the Code Does
1. `FindVariableFeatures()` identifies the 2000 most variable genes using variance-stabilising transformation (VST)
2. `ScaleData()` z-score normalises each gene across cells
3. `RunPCA()` performs PCA on the variable genes
4. `DimPlot()` coloured by cell-cycle phase assesses whether phase is a dominant driver of variation

### Interpretation
This PCA is **diagnostic** — if cells cluster by phase, the S.Score and G2M.Score covariates must be regressed out in SCTransform. This step precedes the main SCTransform normalisation pipeline.

---

## Section 15 — SCTransform Normalisation (Per Sample)

### Code Purpose
Perform SCTransform normalisation independently on each sample, regressing out cell-cycle effects.

```r
split_seurat <- SplitObject(merged, split.by = "sample")
options(future.globals.maxSize = 8 * 1024^3)
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(
    split_seurat[[i]],
    vars.to.regress = c("S.Score", "G2M.Score")
  )
}
```

### What the Code Does
1. `SplitObject()` divides the merged object back into per-sample Seurat objects
2. `options(future.globals.maxSize = 8 * 1024^3)` increases the memory limit for parallel processing to 8 GB — required for large datasets
3. The `for` loop applies `SCTransform()` to each sample independently
4. `vars.to.regress = c("S.Score", "G2M.Score")` removes cell-cycle variation from the normalised matrix as a covariate in the negative binomial regression model

### Expected Output
Each `split_seurat[[i]]` object gains a new `SCT` assay slot containing pearson residuals (normalised values).

### Interpretation
SCTransform is the recommended normalisation approach for Seurat v5. It applies regularised negative binomial regression to model and remove the dependency between gene expression and sequencing depth. Running it per-sample before integration (rather than on the merged object) prevents cross-sample variation from confounding the depth correction.

---

## Section 16 — Data Integration

### Code Purpose
Integrate Day0 and Day14 datasets into a shared embedding using Seurat's anchor-based method.

```r
bm.features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)
bm.list <- PrepSCTIntegration(object.list = split_seurat, anchor.features = bm.features)
bm.anchors <- FindIntegrationAnchors(
  object.list = bm.list, normalization.method = "SCT", anchor.features = bm.features
)
remove(list = setdiff(ls(), "bm.anchors"))
bm.integrated <- IntegrateData(anchorset = bm.anchors, normalization.method = "SCT")
DefaultAssay(bm.integrated) = "integrated"
```

### What the Code Does
1. `SelectIntegrationFeatures()` identifies 3000 genes highly variable across both samples — these anchor the integration
2. `PrepSCTIntegration()` pre-processes SCT residuals for anchor finding
3. `FindIntegrationAnchors()` identifies mutual nearest-neighbour pairs across samples (cells with similar biology that should cluster together)
4. `remove(list = setdiff(ls(), "bm.anchors"))` frees memory before the computationally intensive integration step — **important for large datasets**
5. `IntegrateData()` corrects expression values using the anchor pairs, creating the `integrated` assay

### Expected Output
`bm.integrated` Seurat object with an `integrated` assay slot containing batch-corrected values

### Interpretation
After integration, cells from Day0 and Day14 with the same transcriptional identity should co-cluster. The memory clean-up step (`remove()`) is a practical necessity — integration can require >16 GB RAM on large datasets.

---

## Section 17 — Dimensionality Reduction and Clustering

### Code Purpose
Reduce dimensions (PCA → UMAP) and identify cell clusters using a graph-based approach.

```r
bm.integrated <- RunPCA(object = bm.integrated, verbose = FALSE)

for(i in c(20)) {
  ElbowPlot(bm.integrated, ndims = 50)
  bm.integrated <- RunUMAP(object = bm.integrated, dims = 1:i)
  bm.integrated <- FindNeighbors(bm.integrated, reduction = "pca", dims = 1:i)
  bm.integrated <- FindClusters(bm.integrated, resolution = c(seq(0.1, 2.0, by = 0.1)))
  print(clustree(x = bm.integrated, prefix = "integrated_snn_res.") + ggtitle(paste("clustree", i, "PCs")))
}

Idents(bm.integrated) = "integrated_snn_res.0.1"
DimPlot(bm.integrated, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(bm.integrated, reduction = "umap", split.by = "sample", label = TRUE)
DimPlot(bm.integrated, reduction = "umap", split.by = "Phase", label = TRUE)
```

### What the Code Does
1. `RunPCA()` on the integrated assay
2. `ElbowPlot()` — visualises variance explained per PC; used to choose the number of PCs
3. `RunUMAP(dims = 1:20)` — 20 PCs chosen; UMAP produces 2D embedding
4. `FindNeighbors(dims = 1:20)` — builds KNN graph in PCA space
5. `FindClusters(resolution = seq(0.1, 2.0))` — tests 20 different clustering resolutions simultaneously
6. `clustree()` visualises cluster stability across resolutions
7. `Idents() = "integrated_snn_res.0.1"` — selects resolution 0.1 as the final clustering

### Three Diagnostic DimPlots
- **Unlabelled overall UMAP** — shows cluster geography
- **Split by sample** — verifies inter-sample mixing (key integration QC)
- **Split by phase** — confirms cell-cycle effects were successfully regressed

---

## Section 18 — RNA Assay Normalisation and Marker Visualisation

### Code Purpose
Switch to RNA assay for differential expression and visualise canonical lineage marker genes.

```r
DefaultAssay(bm.integrated) = "RNA"
bm.integrated <- NormalizeData(object = bm.integrated, scale.factor = 10000, normalization.method = "LogNormalize")
bm.integrated <- ScaleData(object = bm.integrated, features = rownames(bm.integrated))

FeaturePlot(bm.integrated, features = c("Cd3e", "Cd19", "Adgre1", "Cd14"), label = TRUE, ncol = 2, order = TRUE)
DotPlot(bm.integrated, features = c("Cd3e", "Cd19", "Adgre1", "Cd14"))
bm.integrated = JoinLayers(bm.integrated)
```

### What the Code Does
1. Switches `DefaultAssay` to `RNA` — integration assay should not be used for gene expression quantification
2. Applies log-normalisation and scaling to the RNA assay
3. `FeaturePlot()` overlays individual gene expression on the UMAP — each of the four genes marks a distinct lineage
4. `DotPlot()` shows both expression level (dot size = % cells expressing) and expression magnitude (colour) per cluster
5. `JoinLayers()` collapses the RNA assay layers (required for Seurat v5 before differential expression)

### Canonical Marker Interpretation
| Gene | Cell Lineage |
|---|---|
| Cd3e | Pan-T cell |
| Cd19 | B cell |
| Adgre1 (F4/80) | Macrophage |
| Cd14 | Monocyte |

---

## Section 19 — Differential Expression (FindAllMarkers)

### Code Purpose
Identify marker genes for each cluster using a rigorous statistical test.

```r
markers <- FindAllMarkers(
  object = bm.integrated,
  logfc.threshold = 0.25,
  test.use = "negbinom",
  only.pos = TRUE,
  return.thresh = 0.05
)
markers <- markers %>% filter(p_val_adj > 0 & p_val_adj < 0.05)
```

### What the Code Does
1. `FindAllMarkers()` tests each gene for differential expression in each cluster vs. all others (one-vs-rest)
2. `logfc.threshold = 0.25` — minimum log2 fold-change; reduces computation by pre-filtering
3. `test.use = "negbinom"` — negative binomial test accounts for the overdispersed, zero-inflated nature of count data
4. `only.pos = TRUE` — returns only upregulated markers per cluster
5. `return.thresh = 0.05` — pre-filters to FDR-adjusted p-value < 0.05
6. The `filter()` step removes any markers with `p_val_adj == 0` (a Seurat artefact where identical cells produce exact p = 0)

### Heatmap Visualisation
```r
markers %>% group_by(cluster) %>% filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup() -> top10
DoHeatmap(seurat_integrated, features = top10$gene) + NoLegend()
DotPlot(seurat_integrated, features = top10$gene)
```

Selects the top 10 markers per cluster by log2FC for heatmap and DotPlot visualisation.

---

## Section 20 — T-cell Subclustering

### Code Purpose
Extract T cells from the integrated object and perform independent subclustering.

```r
T_cell = subset(seurat_integrated, subset = integrated_snn_res.0.1 == "3")
DefaultAssay(T_cell) = "integrated"
T_cell <- RunUMAP(object = T_cell, dims = 1:20)
T_cell <- FindNeighbors(T_cell, reduction = "pca", dims = 1:i)
T_cell <- FindClusters(T_cell, resolution = c(seq(0.1, 2.0, by = 0.1)))
```

### What the Code Does
1. `subset()` isolates cluster 3 (the T-cell cluster) from the integrated object
2. Runs UMAP, KNN graph construction, and clustering on this subset independently
3. `clustree()` visualises subclustering stability; resolution 0.2 is selected

### Why Subcluster?
T cells form a heterogeneous population (Th1, Th17, Treg, γδ, CD8+ CTL) that may not be separable at the whole-dataset clustering resolution. Subclustering with a finer lens reveals these subtypes.

---

## Section 21 — T-cell Marker Gene Annotation

### Code Purpose
Annotate T-cell subclusters using canonical marker genes for known T-cell subtypes.

```r
tcell_clusters <- list(
  Cytotoxic_T_cells    = c("Lef1", "Ms4a4b", "Cd8a"),
  Th17_cells           = c("Tnfsf8", "Cxcr6", "Cd4", "Il17a", "Rora"),
  Regulatory_T_cells   = c("Foxp3", "Ctla4"),
  Il17a_positive_gdT   = c("Il23r", "Trgv2"),
  Il17a_negative_gdT   = c("Birc5", "Top2a")
)
DotPlot(T_cell, features = tcell_clusters)
```

### Marker Interpretation
| Subtype | Key Markers | Function |
|---|---|---|
| Cytotoxic T cells | Cd8a, Lef1, Ms4a4b | Direct cell killing |
| Th17 cells | Cd4, Il17a, Rora, Cxcr6 | Inflammatory IL-17 secretors |
| Regulatory T cells | Foxp3, Ctla4 | Immune suppression |
| IL-17a+ γδ T cells | Il23r, Trgv2 | Innate-like T cells |
| IL-17a− γδ T cells | Birc5, Top2a | Cycling/proliferating γδ cells |

---

## Section 22 — Trajectory Inference (Monocle3)

### Code Purpose
Order T-cell subclusters along a pseudotime trajectory to reconstruct developmental relationships.

```r
library(monocle3)
library(SeuratWrappers)
DefaultAssay(T_cell) = "integrated"

bm.integrated.cds <- as.cell_data_set(T_cell, group.by = "integrated_snn_res.0.2")
bm.integrated.cds <- cluster_cells(cds = bm.integrated.cds, reduction_method = "UMAP")
bm.integrated.cds <- learn_graph(bm.integrated.cds, use_partition = TRUE)
plot_cells(cds = bm.integrated.cds, show_trajectory_graph = TRUE)
bm.integrated.cds <- order_cells(bm.integrated.cds, reduction_method = "UMAP")
plot_cells(cds = bm.integrated.cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
```

### What the Code Does
1. `SeuratWrappers::as.cell_data_set()` converts the Seurat T-cell object to Monocle3's `cell_data_set` format
2. `cluster_cells()` reclusters cells inside Monocle3 using UMAP
3. `learn_graph()` fits a principal graph through the UMAP embedding — this is the trajectory skeleton
4. `order_cells()` assigns pseudotime values; the root node (biological start point) must be chosen interactively or specified programmatically
5. The final `plot_cells()` colours cells by pseudotime — cells closer to the root have lower pseudotime values

### Interpretation
The `use_partition = TRUE` argument restricts the trajectory to each connected partition of the UMAP, preventing biologically unconnected cell states from being linked in the trajectory. Root selection is the most subjective step and should be guided by biological knowledge of the developmental starting point.

---

## Section 23 — CellChat Analysis: Day 0

### Code Purpose
Build and analyse the cell–cell communication network for the Day 0 time point.

```r
library(CellChat); library(patchwork)
options(stringsAsFactors = FALSE)

# Extract data
data.input.Day0 <- LayerData(Day0, assay = "RNA", layer = "data")
labels <- paste0("group", as.character(Idents(Day0)))
metaDay0 <- data.frame(group = as.factor(labels), row.names = colnames(data.input.Day0))

# Build CellChat object
cellchatDay0 <- createCellChat(object = data.input.Day0, meta = metaDay0, group.by = "group")

# Assign database
CellChatDB <- CellChatDB.mouse
cellchatDay0@DB <- CellChatDB

# Run pipeline
cellchatDay0 <- subsetData(cellchatDay0)
cellchatDay0 <- identifyOverExpressedGenes(cellchatDay0)
cellchatDay0 <- identifyOverExpressedInteractions(cellchatDay0)
cellchatDay0 <- computeCommunProb(cellchatDay0)
cellchatDay0 <- filterCommunication(cellchatDay0, min.cells = 10)
cellchatDay0 <- computeCommunProbPathway(cellchatDay0)
cellchatDay0 <- aggregateNet(cellchatDay0)
```

### What the Code Does
The CellChat pipeline proceeds in a defined sequence:
1. Input data is the log-normalised RNA assay (`layer = "data"`)
2. Cluster labels are prefixed with "group" to avoid issues with numeric cluster names
3. `subsetData()` keeps only genes present in both the expression matrix and CellChatDB
4. `identifyOverExpressedGenes/Interactions()` finds L-R pairs where both ligand and receptor are significantly overexpressed in at least one cell group
5. `computeCommunProb()` calculates communication probability using a mass action model
6. `filterCommunication(min.cells = 10)` removes interactions from groups with <10 cells
7. `computeCommunProbPathway()` aggregates L-R pair probabilities to the pathway level
8. `aggregateNet()` summarises the full communication network as interaction counts and weights

### Expected Output
A `cellchatDay0` object with populated `@net`, `@netP`, and `@LR` slots; circle plots showing interaction number and strength.

---

## Section 24 — CellChat Analysis: Day 14

Identical pipeline applied to the Day14 sample. Key addition:

```r
metaDay14$samples <- "Day14"
```

This adds a sample column (recommended by CellChat to suppress a warning), but is optional.

---

## Section 25 — CellChat Comparative Analysis

### Code Purpose
Merge Day0 and Day14 CellChat objects and perform comparative cell–cell communication analysis.

```r
cellchatDay0 = netAnalysis_computeCentrality(cellchatDay0)
cellchatDay14 = netAnalysis_computeCentrality(cellchatDay14)
object.list <- list(Day0 = cellchatDay0, Day14 = cellchatDay14)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
```

### Comparative Analyses Performed

| Analysis | Function | What It Shows |
|---|---|---|
| Total interactions | `compareInteractions()` | Number and strength change Day0→Day14 |
| Differential per cell type | `netVisual_diffInteraction()` | Which cell types gained/lost communication |
| Interaction heatmaps | `netVisual_heatmap()` | Sender–receiver matrices per condition |
| Signalling roles | `netAnalysis_signalingRole_scatter()` | Cells as senders vs. receivers |
| Functional pathway similarity | `computeNetSimilarityPairwise(type="functional")` | Pathways grouped by who they talk to |
| Structural pathway similarity | `computeNetSimilarityPairwise(type="structural")` | Pathways grouped by ligand-receptor structure |
| Pathway ranking | `rankNet()` | Which pathways carry the most information flow |
| Outgoing signals | `netAnalysis_signalingRole_heatmap(pattern="outgoing")` | Per-cell-type outgoing signal heatmap |
| Bubble plots | `netVisual_bubble()` | L-R pair expression per sender-receiver pair |

### Interpretation
`netAnalysis_computeCentrality()` must be run before merging — it computes network centrality metrics (betweenness, closeness, degree) that identify the most influential sender and receiver cell types in each condition. Comparing Day0 vs. Day14 communication networks reveals which signalling pathways are conserved, gained, or lost during the biological process under study.
