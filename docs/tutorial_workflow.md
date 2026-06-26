# Tutorial Workflow: Step-by-Step Guide

> This document provides a sequential walkthrough of the complete analysis, linking each code step to the conceptual knowledge from the accompanying presentation.

---

## Prerequisites

- R ≥ 4.3.0 installed
- All packages listed in `environment/requirements_or_dependencies.md` installed
- GEO data from **GSE142564** downloaded and uncompressed
- Minimum 32 GB RAM recommended (16 GB possible with memory-saving options)
- RStudio or equivalent R environment

---

## Step 1 — Environment Setup

**What happens:** Load all required R libraries.

```r
library(SingleCellExperiment); library(Seurat); library(tidyverse)
library(Matrix); library(scales); library(cowplot); library(RCurl)
```

**Presentation concept (Slide 11):** Large-scale scRNA-seq data requires specialised libraries. The `Matrix` package enables sparse matrix handling, which is the only practical way to store 20,000-gene × 10,000-cell count matrices in memory.

**Expected result:** R session loads cleanly without errors. All functions from the listed packages are available.

**Common issue:** Package version conflicts, especially between Seurat v5 and older Bioconductor packages. Use `sessionInfo()` to document your environment.

---

## Step 2 — Data Loading (Day 0)

**What happens:** Read the Day 0 10X CellRanger output into R and create a Seurat object.

```r
data <- Read10X("/path/to/GSE142564_RAW/Day0/")
Day0 <- CreateSeuratObject(counts = data)
Day0 <- subset(Day0, subset = nFeature_RNA >= 200)
```

**Presentation concept (Slides 21–23):** CellRanger produces three output files — `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` — that together define the gene × cell count matrix. The Seurat object wraps this matrix with associated metadata.

**Expected result:**
- `data` is a sparse matrix (dgCMatrix) with genes as rows and cell barcodes as columns
- `Day0` is a Seurat object; `Day0` printed to console shows dimensions (e.g., `23,000 features × 8,000 cells`)
- After `subset()`, cells with <200 genes are removed

**⚠️ User action required:** Update the file path to your local data location.

---

## Step 3 — QC Metric Calculation

**What happens:** Calculate mitochondrial ratio and complexity score for each cell.

```r
Day0$mitoRatio <- PercentageFeatureSet(object = Day0, pattern = "^mt-") / 100
Day0$log10GenesPerUMI <- log10(Day0$nFeature_RNA) / log10(Day0$nCount_RNA)
```

**Presentation concept (Slides 26–28):**
- Mitochondrial ratio >0.05 indicates membrane-compromised (dying) cells — their cytoplasmic RNA has leaked out, leaving only mitochondria-enclosed transcripts
- The log10 genes-per-UMI ratio is a proxy for transcriptional complexity; low values indicate empty droplets, red blood cells, or doublets

**Expected result:** Two new columns appear in `Day0@meta.data`: `mitoRatio` (values 0–1) and `log10GenesPerUMI` (values typically 0.7–1.0 for healthy cells).

**Sanity check:** `grep("^mt-", rownames(Day0))` should return a non-empty vector. If it returns `integer(0)`, the count matrix may not contain mitochondrial genes, or the gene naming convention differs (check if `^MT-` is needed for human data).

---

## Step 4 — Pre-filtering QC Visualisation

**What happens:** Generate violin plots of QC metrics to inspect the raw data distribution.

```r
VlnPlot(Day0, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3, layer = "counts")
```

**Presentation concept (Slides 29–32):** Visualising the distribution before filtering is essential. You cannot know the right thresholds without first seeing the data.

**Expected result:** Three violin plots showing the distribution of each metric. Look for:
- A main cell population peak in nFeature and nCount
- A long upper tail (potential doublets) in nFeature and nCount
- Most cells below 0.05 in mitoRatio, with a tail of high-mito cells

**Learner task:** Note the approximate position of the main distribution peak and the extent of the tails. These inform your filtering thresholds.

---

## Step 5 — Initial Cell Filtering

**What happens:** Remove potential doublets and dead/damaged cells.

```r
Day0 = subset(Day0, nFeature_RNA < 4000 & nCount_RNA < 20000 & mitoRatio < .05)
```

**Presentation concept (Slides 29–32):** Thresholds are dataset-specific. The rationale:
- Upper gene count limit (4000): doublets have ~2× the gene count of singlets
- Upper UMI limit (20000): captures unusually high sequencing depth (likely technical)
- Mitochondrial upper limit (0.05): removes cells with cytoplasmic RNA loss

**Expected result:** Reduced cell count. Run `dim(Day0)` before and after to quantify removed cells.

**⚠️ Note:** These thresholds were set for GSE142564. For your own data, always explore the distribution first.

---

## Step 6 — Metadata Organisation and Extended QC Visualisation

**What happens:** Clean metadata and generate detailed diagnostic plots.

```r
# Rename and organise metadata
metadata <- Day0@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>% dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA)
Day0@meta.data <- metadata
```

Then generate:
1. Cell count bar chart per sample
2. UMI density plot (log10 scale, lower cutoff at 2000)
3. Gene count density plot (log10 scale, lower cutoff at 1000)
4. Mitochondrial ratio density plot
5. UMI vs. Gene scatter coloured by mitochondrial ratio (with regression line)
6. Complexity score density plot (cutoff at 0.82)

**Presentation concept (Slides 29–32):** The UMI-vs-Gene scatter plot is the most powerful single QC figure: it simultaneously visualises three metrics and reveals the joint distribution of quality.

**Expected result:** The UMI-vs-Gene scatter should show a tight linear cloud with the majority of cells above both red threshold lines. Cells in the lower-left quadrant with dark colour (high mitoRatio) are the targets for removal.

---

## Step 7 — Final QC Filtering with Complexity Threshold

**What happens:** Apply complexity and depth thresholds to remove the final low-quality cells.

```r
Day0 <- subset(x = Day0, subset = (log10GenesPerUMI > 0.82) & (nGene > 900) & (nUMI > 3200))
dim(Day0)
VlnPlot(Day0, features = c("nGene", "nUMI", "mitoRatio"))
```

**Presentation concept (Slides 27–28):** The complexity threshold (`log10GenesPerUMI > 0.82`) is a lower-bound filter that removes transcriptionally simple cells (empty droplets, RBCs, debris) that might have passed the upper-bound filters.

**Expected result:** Further reduced cell count. The final violin plots should show clean, unimodal distributions without the artefactual tails present in the raw data.

---

## Step 8 — Repeat QC for Day 14

**What happens:** Apply the identical QC pipeline to the Day14 sample.

**Presentation concept (Slides 13–14, Batch Effects):** Consistent QC across samples is critical. Using different thresholds per sample introduces artificial differences between groups.

**Learner task:** Duplicate Steps 2–7, replacing `Day0` with `Day14` and updating the input path. Compare the final cell counts and QC metric distributions between samples.

---

## Step 9 — Merge Datasets

**What happens:** Combine Day0 and Day14 into a single merged Seurat object.

```r
Day0 <- RenameCells(Day0, add.cell.id = "Day0")
Day14 <- RenameCells(Day14, add.cell.id = "Day14")
merged <- merge(x = Day0, y = c(Day14))
```

**Presentation concept (Slides 13–14, Batch Effects):** Merging is not integrating. After `merge()`, cells will still cluster by sample identity if run through PCA directly — this is why integration (Step 15) is necessary.

**Expected result:** `merged` contains cells from both samples. Each cell barcode is prefixed (e.g., `Day0_ACGTACGT-1`). `table(merged$sample)` shows cell counts per sample.

---

## Step 10 — BCR/TCR Gene Removal and Gene Filtering

**What happens:** Remove immunoglobulin and T-cell receptor variable genes; filter sparse genes.

```r
# Remove BCR/TCR genes
merged <- merged[!grepl("^Ig[hkl]v", rownames(merged), ignore.case = FALSE), ]
# ... (5 regex patterns total)

# Remove genes expressed in <10 cells
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])
counts <- LayerData(merged, assay = "RNA", layer = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
dual_counts <- counts[keep_genes, ]
merged <- CreateSeuratObject(dual_counts, meta.data = merged@meta.data)
```

**Presentation concept (Slide 8):** In immune cell datasets, Ig/TCR variable genes are clonally expressed and will dominate PCA, causing cells to cluster by clonotype rather than functional state. Removing them is biologically justified and computationally necessary.

**Expected result:** `nrow(merged)` (gene count) decreases substantially after both filtering steps. Run `nrow(merged)` before and after to quantify.

---

## Step 11 — Cell Cycle Scoring

**What happens:** Score cells for S-phase and G2M-phase gene signatures.

```r
merged <- NormalizeData(merged)
s.genes <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)
merged <- CellCycleScoring(merged, g2m.features = g2m.genes, s.features = s.genes, set.ident = TRUE)
RidgePlot(merged, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)
```

**Presentation concept (Slide 44):** Cell cycle is a major source of transcriptional variation in proliferating cells. Without regression, cycling and non-cycling cells will cluster separately regardless of their actual cell type — masking biology.

**Expected result:** `merged@meta.data` gains `S.Score`, `G2M.Score`, and `Phase` columns. The RidgePlot shows that cycling markers (Pcna, Mki67) are specifically elevated in G2M-phase and S-phase cells.

---

## Step 12 — Variable Features, Scaling, and Diagnostic PCA

**What happens:** Identify the 2000 most variable genes, scale the data, and run PCA to check cell-cycle influence.

```r
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged)
DimPlot(merged, reduction = "pca", group.by = "Phase", split.by = "Phase")
```

**Presentation concept (Slides 45–52):** Scaling centres each gene to mean = 0 and variance = 1, making all genes contribute equally to PCA regardless of expression magnitude.

**Expected result:** If the PCA DimPlot shows cells separating by Phase (G1 cells on one side, S and G2M cells on the other), cell-cycle regression in SCTransform is necessary — which it is in this pipeline.

---

## Step 13 — SCTransform Normalisation Per Sample

**What happens:** Apply regularised negative binomial regression independently to each sample, regressing out cell-cycle effects.

```r
split_seurat <- SplitObject(merged, split.by = "sample")
options(future.globals.maxSize = 8 * 1024^3)
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("S.Score", "G2M.Score"))
}
```

**Presentation concept (Slides 54–58):** SCTransform models the relationship between gene expression and sequencing depth using a regularised negative binomial regression. This is superior to standard log-normalisation because it explicitly models count data statistics and accounts for over-dispersion.

**Expected result:** Each element of `split_seurat` gains a new `SCT` assay slot. Runtime: 5–30 minutes per sample depending on cell number and hardware.

---

## Step 14 — Anchor-based Data Integration

**What happens:** Find biologically corresponding cells across samples and integrate them into a shared embedding.

```r
bm.features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)
bm.list <- PrepSCTIntegration(object.list = split_seurat, anchor.features = bm.features)
bm.anchors <- FindIntegrationAnchors(object.list = bm.list, normalization.method = "SCT", anchor.features = bm.features)
remove(list = setdiff(ls(), "bm.anchors"))
bm.integrated <- IntegrateData(anchorset = bm.anchors, normalization.method = "SCT")
DefaultAssay(bm.integrated) = "integrated"
```

**Presentation concept (Slides 60–67):** Integration corrects for batch effects by identifying "anchor" cell pairs — cells from different samples that are biologically equivalent (mutual nearest neighbours in expression space). The before/after UMAP comparison shown on Slide 67 demonstrates the effect.

**Expected result:** `bm.integrated` Seurat object with an `integrated` assay. Runtime: 10–60 minutes depending on cell number.

**⚠️ Memory note:** The `remove()` call is essential — it frees memory before the most RAM-intensive step.

---

## Step 15 — PCA, UMAP, and Clustering

**What happens:** Reduce dimensions and identify cell communities across the integrated dataset.

```r
bm.integrated <- RunPCA(bm.integrated)
ElbowPlot(bm.integrated, ndims = 50)                        # Choose number of PCs
bm.integrated <- RunUMAP(bm.integrated, dims = 1:20)
bm.integrated <- FindNeighbors(bm.integrated, reduction = "pca", dims = 1:20)
bm.integrated <- FindClusters(bm.integrated, resolution = c(seq(0.1, 2.0, by = 0.1)))
clustree(x = bm.integrated, prefix = "integrated_snn_res.")  # Assess resolution
Idents(bm.integrated) = "integrated_snn_res.0.1"
```

**Presentation concepts (Slides 69–75):**
- ElbowPlot: the "elbow" in variance explained indicates the optimal number of PCs
- UMAP: compresses the high-dimensional PCA embedding into interpretable 2D
- clustree: visualises cluster stability across resolutions — choose the lowest resolution that yields biologically distinct clusters

**Expected result:** A UMAP plot showing discrete cell clusters. With `resolution = 0.1`, expect 5–15 clusters depending on dataset complexity.

---

## Step 16 — Integration Quality Control

**What happens:** Visually verify that integration removed batch effects.

```r
DimPlot(bm.integrated, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(bm.integrated, reduction = "umap", split.by = "sample", label = TRUE)
DimPlot(bm.integrated, reduction = "umap", split.by = "Phase", label = TRUE)
```

**Presentation concept (Slides 77–79):** Successful integration shows Day0 and Day14 cells co-occupying the same cluster positions. Clustering by sample identity indicates failed integration.

**Expected result:**
1. Labelled UMAP showing numbered clusters
2. Sample-split UMAP: Day0 and Day14 cells intermixed within each cluster (not separated)
3. Phase-split UMAP: cells distributed across the UMAP by biology, not by cell-cycle phase

**Pass/fail criterion:** If the sample-split UMAP shows Day0 and Day14 cells in completely separate UMAP regions, integration failed and must be troubleshot.

---

## Step 17 — Marker Gene Expression and RNA Normalisation

**What happens:** Switch to RNA assay, normalise for differential expression, and visualise canonical lineage markers.

```r
DefaultAssay(bm.integrated) = "RNA"
bm.integrated <- NormalizeData(bm.integrated, scale.factor = 10000)
bm.integrated <- ScaleData(bm.integrated, features = rownames(bm.integrated))
FeaturePlot(bm.integrated, features = c("Cd3e", "Cd19", "Adgre1", "Cd14"), label = TRUE, ncol = 2)
DotPlot(bm.integrated, features = c("Cd3e", "Cd19", "Adgre1", "Cd14"))
bm.integrated = JoinLayers(bm.integrated)
```

**Presentation concept (Slides 80–82):** Lineage marker genes confirm the biological identity of each cluster before running `FindAllMarkers()`. The DotPlot simultaneously shows expression level and fraction of expressing cells.

**Expected result:** FeaturePlot panels showing each marker gene expressed in a subset of UMAP clusters. Clusters with high Cd3e = T cells; Cd19 = B cells; Adgre1 = macrophages; Cd14 = monocytes.

---

## Step 18 — Differential Expression and Cluster Annotation

**What happens:** Identify top marker genes per cluster and annotate cell types.

```r
markers <- FindAllMarkers(object = bm.integrated, logfc.threshold = 0.25,
                          test.use = "negbinom", only.pos = TRUE, return.thresh = 0.05)
markers <- markers %>% filter(p_val_adj > 0 & p_val_adj < 0.05)

top10 <- markers %>% group_by(cluster) %>% filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup()
DoHeatmap(seurat_integrated, features = top10$gene) + NoLegend()
DotPlot(seurat_integrated, features = top10$gene)
```

**Presentation concept (Slides 83–84):** `FindAllMarkers()` is the key discovery step — it identifies which genes distinguish each cluster. The heatmap provides an overview; individual FeaturePlots confirm specific assignments.

**Expected result:**
- `markers` dataframe with columns: gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster
- Heatmap showing cluster-specific expression blocks
- Learner task: annotate each cluster with a cell type label based on top markers and lineage marker FeaturePlots

**Runtime:** `FindAllMarkers()` with `negbinom` is computationally intensive — expect 20–60 minutes on a standard laptop.

---

## Step 19 — T-cell Subclustering

**What happens:** Isolate the T-cell cluster and perform independent re-clustering.

```r
T_cell = subset(seurat_integrated, subset = integrated_snn_res.0.1 == "3")
DefaultAssay(T_cell) = "integrated"
T_cell <- RunUMAP(T_cell, dims = 1:20)
T_cell <- FindNeighbors(T_cell, reduction = "pca", dims = 1:20)
T_cell <- FindClusters(T_cell, resolution = c(seq(0.1, 2.0, by = 0.1)))
Idents(T_cell) = "integrated_snn_res.0.2"
```

**Presentation concept (Slide 85):** Subclustering is performed when a broad cell type (here: T cells) contains biologically distinct subtypes that require fine-grained resolution to separate.

**Expected result:** A new UMAP specific to the T-cell subset showing 3–7 subclusters at resolution 0.2.

**⚠️ Cluster number note:** The code uses `subset = integrated_snn_res.0.1 == "3"` — this assumes cluster 3 is the T-cell cluster. In your data, the T-cell cluster may have a different number. Verify using `FeaturePlot(features = "Cd3e")` before subsetting.

---

## Step 20 — T-cell Marker Visualisation and Subtype Annotation

**What happens:** Annotate T-cell subclusters using canonical subtype markers.

```r
DefaultAssay(T_cell) = "RNA"
T_cell <- NormalizeData(T_cell, scale.factor = 10000)
T_cell <- ScaleData(T_cell, features = rownames(T_cell))

tcell_clusters <- list(
  Cytotoxic_T_cells = c("Lef1", "Ms4a4b", "Cd8a"),
  Th17_cells = c("Tnfsf8", "Cxcr6", "Cd4", "Il17a", "Rora"),
  Regulatory_T_cells = c("Foxp3", "Ctla4"),
  Il17a_positive_gdT_cells = c("Il23r", "Trgv2"),
  Il17a_negative_gdT_cells = c("Birc5", "Top2a")
)
DotPlot(T_cell, features = tcell_clusters)
```

**Presentation concept (Slide 86):** T-cell subtypes are identified by canonical surface markers (Cd8a, Cd4), transcription factors (Foxp3, Rora), and cytokines (Il17a).

**Expected result:** DotPlot showing enrichment of each marker gene set in specific subclusters, enabling subtype annotation.

---

## Step 21 — Trajectory Inference with Monocle3

**What happens:** Order T-cell subclusters along a pseudotime developmental trajectory.

```r
library(monocle3); library(SeuratWrappers)
DefaultAssay(T_cell) = "integrated"
bm.integrated.cds <- as.cell_data_set(T_cell, group.by = "integrated_snn_res.0.2")
bm.integrated.cds <- cluster_cells(cds = bm.integrated.cds, reduction_method = "UMAP")
bm.integrated.cds <- learn_graph(bm.integrated.cds, use_partition = TRUE)
plot_cells(cds = bm.integrated.cds, show_trajectory_graph = TRUE)
bm.integrated.cds <- order_cells(bm.integrated.cds, reduction_method = "UMAP")
plot_cells(cds = bm.integrated.cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
```

**Presentation concept (Slides 88–91):** Trajectory inference reconstructs the continuum of transcriptional states that cells pass through during differentiation. Pseudotime is an abstract time axis — it orders cells from least to most differentiated based on their transcriptional profile.

**Expected result:**
1. UMAP with the learned principal graph overlaid (trajectory skeleton)
2. UMAP coloured by pseudotime — cells coloured from dark (early) to light (late)

**⚠️ Interactive step:** `order_cells()` opens an interactive browser window to select the root node. Choose the cell cluster representing the least differentiated state (e.g., naive T cells or precursors based on marker gene expression).

---

## Step 22 — CellChat: Day 0 Communication Network

**What happens:** Build and analyse the cell–cell communication network for Day 0 using ligand–receptor interactions.

```r
library(CellChat); library(patchwork)
# (Full pipeline: see code_explanation.md Section 23 for details)
```

**Presentation concept (Slides 92–109):** CellChat infers communication probability between cell populations based on the co-expression of known ligand–receptor pairs. The mouse CellChatDB contains hundreds of validated L-R interactions.

**Expected result:** Circle plot showing interaction number and strength between all cell groups, with edge width proportional to interaction frequency/strength.

---

## Step 23 — CellChat: Day 14 and Comparative Analysis

**What happens:** Build the Day 14 CellChat object and compare communication networks across time points.

```r
# Merge CellChat objects
cellchatDay0 = netAnalysis_computeCentrality(cellchatDay0)
cellchatDay14 = netAnalysis_computeCentrality(cellchatDay14)
cellchat <- mergeCellChat(list(Day0 = cellchatDay0, Day14 = cellchatDay14), add.names = c("Day0", "Day14"))

# Comparative visualisations
compareInteractions(cellchat, show.legend = FALSE, group = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
netVisual_heatmap(cellchat)
rankNet(cellchat, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE)
```

**Presentation concept (Slides 92–109):** The comparative analysis identifies:
- Which cell types are the dominant senders/receivers at each time point
- Which signalling pathways are gained or lost between Day0 and Day14
- Whether pathways are functionally similar (grouped by who they communicate between) or structurally similar (grouped by L-R pair family)

**Expected result:** Side-by-side comparison plots, differential interaction network, and pathway ranking bar charts showing information flow changes.

---

## Analysis Complete

At the end of this workflow, you will have:

1. ✅ Cleaned, filtered scRNA-seq data from two time points
2. ✅ Integrated dataset free of batch effects
3. ✅ Annotated cell clusters with lineage identities
4. ✅ Subclustered T-cell populations with subtype annotations
5. ✅ Pseudotime trajectory through T-cell developmental states
6. ✅ Cell–cell communication networks for each time point
7. ✅ Comparative CellChat analysis identifying rewired signalling pathways

---

## Troubleshooting Quick Reference

| Error | Likely Cause | Solution |
|---|---|---|
| `Read10X()` fails | Wrong file path or missing files | Check directory contents with `list.files()` |
| `mitoRatio` all zeros | Wrong gene prefix (human vs. mouse) | Change `^mt-` to `^MT-` for human data |
| `SCTransform()` crashes | Insufficient RAM | Reduce `nfeatures` or process fewer cells |
| Integration takes hours | Large dataset | Consider Harmony as a faster alternative |
| `as.cell_data_set()` fails | SeuratWrappers version mismatch | Reinstall SeuratWrappers from GitHub |
| CellChat empty interactions | Cluster too small | Lower `min.cells` threshold |
