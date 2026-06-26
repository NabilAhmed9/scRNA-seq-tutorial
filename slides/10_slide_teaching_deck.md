# 10-Slide Teaching Deck: Single-Cell RNA-seq Analysis

**Course:** Introduction to scRNA-seq: Concepts, Applications, and Analysis  
**Author:** Based on materials by Farag Ibrahim, University of Gothenburg, Sweden  
**Format:** Markdown teaching deck — each slide includes bullet points, code snippets, and teaching notes

---

## Slide 1 — Project Overview and Learning Goals

### Title: From Raw Counts to Cell Communication: A Complete scRNA-seq Tutorial

**Bullet points:**
- Analyse a real mouse bone marrow immune dataset (GEO: GSE142564, Day 0 vs. Day 14)
- Learn the full scRNA-seq pipeline: QC → normalisation → integration → clustering → annotation → trajectory → cell communication
- Use R with Seurat, Monocle3, and CellChat — the three dominant tools in the field
- Understand *why* each step is necessary, not just *how* to run it
- Build intuition for interpreting UMAP plots, marker heatmaps, and communication network diagrams

**Key framework:**

```
Raw 10X data → QC → SCTransform → Integration → UMAP/Clustering
     → Annotation → Subclustering → Monocle3 → CellChat
```

**Teaching note:** Begin by emphasising that this tutorial uses real data with real challenges (batch effects, cell-cycle variation, immune gene interference). Every filtering decision has a biological justification, not an arbitrary threshold. Learners should expect to spend the most time on steps 3–5 (QC, normalisation, integration) because these determine the quality of everything downstream.

---

## Slide 2 — What Single-cell RNA-seq Is About

### Title: Why Single-Cell Resolution Changes Everything

**Bullet points:**
- Bulk RNA-seq measures the **average** expression across all cells — rare populations and distinct states are invisible
- scRNA-seq profiles **each cell individually**, revealing heterogeneity within tissues
- Four transformative capabilities: (1) cell type identification, (2) rare cell detection, (3) developmental trajectory reconstruction, (4) cell-type-specific differential expression
- Modern platforms (10X Genomics Chromium) capture 10,000+ cells per sample using microfluidic GEM encapsulation
- Each cell barcode + UMI combination uniquely identifies one transcript from one cell

**Conceptual comparison:**

| Feature | Bulk RNA-seq | scRNA-seq |
|---|---|---|
| Resolution | Tissue average | Single cell |
| Heterogeneity | Hidden | Revealed |
| Rare cells | Undetectable | Detectable |
| Cost | Lower | Higher |
| Data volume | Manageable | Massive (sparse) |

**Teaching note:** Use the analogy of a smoothie vs. individual fruits. Bulk RNA-seq is like measuring the colour of a smoothie — you lose all information about individual ingredients. scRNA-seq measures each fruit separately. Emphasise that most tissues contain 5–50 distinct cell types; bulk RNA-seq masks all of this.

---

## Slide 3 — Input Data and Environment Setup

### Title: What the Analysis Starts With

**Bullet points:**
- Input: CellRanger output for each sample — three files: `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`
- These define a sparse gene × cell count matrix where each entry is a UMI count
- The `Matrix` package handles this sparse format (>90% zeros) memory-efficiently
- Seurat's `Read10X()` reads all three files simultaneously; `CreateSeuratObject()` wraps them
- Starting point: a Seurat object with `nCount_RNA` (UMIs) and `nFeature_RNA` (genes) pre-calculated

**Setup code:**

```r
library(Seurat); library(tidyverse); library(Matrix)
library(monocle3); library(CellChat)

data <- Read10X("/path/to/Day0/")
Day0 <- CreateSeuratObject(counts = data)
Day0  # Prints: 23,000 features × 8,000 cells
```

**Teaching note:** Ask learners what `nCount_RNA` and `nFeature_RNA` represent before showing the answer. These two metrics are the foundation of QC. A key insight: in 10X data, `nCount_RNA` counts **UMIs** (not reads), so PCR duplicates have already been removed by CellRanger.

---

## Slide 4 — Code Structure Overview

### Title: How the Analysis is Organised

**Bullet points:**
- The analysis has 14 distinct stages, each separated by a header block in the code
- Two main analysis modules: (A) Seurat-based single-cell analysis, (B) CellChat-based communication analysis
- Key principle: use the **integrated assay** for clustering/UMAP; use the **RNA assay** for gene expression quantification and differential expression
- Object lineage: `Day0` / `Day14` → `merged` → `split_seurat` → `bm.integrated` → `T_cell` → `bm.integrated.cds`

**Pipeline map:**

```
Day0 QC ─┐
          ├─ merge() ─── BCR/TCR filter ─── cell cycle ─── SCTransform ─── integration
Day14 QC ─┘                                                                     │
                                                                          PCA/UMAP/clustering
                                                                                 │
                                                          ┌──────────────────────┴──────────────────────┐
                                                     Annotation                                    T-cell subset
                                                          │                                              │
                                                     CellChat                                      Monocle3
```

**Teaching note:** Print the code file and annotate it together with learners before running it. Understanding the overall structure prevents the common mistake of running sections out of order (e.g., trying to use the integrated assay for `FindAllMarkers()`).

---

## Slide 5 — Data Loading and Quality Control Metrics

### Title: How to Distinguish Real Cells from Noise

**Bullet points:**
- Four QC metrics per cell: `nCount_RNA` (UMIs), `nFeature_RNA` (genes), `mitoRatio`, `log10GenesPerUMI`
- **Mitochondrial ratio:** cells with >5% mitochondrial reads have lost cytoplasmic RNA — they are dying or damaged
- **Complexity score:** `log10(genes) / log10(UMIs)` — healthy cells express many distinct genes per transcript; empty droplets do not
- Mouse mitochondrial genes use `^mt-` prefix (human: `^MT-`) — the pattern must match your organism
- Always visualise the distribution **before** choosing thresholds

**Metric calculation code:**

```r
# Mitochondrial fraction (mouse)
Day0$mitoRatio <- PercentageFeatureSet(object = Day0, pattern = "^mt-") / 100

# Complexity score
Day0$log10GenesPerUMI <- log10(Day0$nFeature_RNA) / log10(Day0$nCount_RNA)

# Visualise ALL metrics before filtering
VlnPlot(Day0, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)
```

**Teaching note:** The most common learner mistake is applying thresholds from a paper or tutorial without inspecting their own data first. Emphasise that QC thresholds are **dataset-specific**. A 5% mitochondrial cutoff is appropriate for bone marrow immune cells but would be wrong for cardiac tissue. The scatter plot of `nUMI` vs. `nGene` coloured by `mitoRatio` is the single most informative QC figure.

---

## Slide 6 — Quality Control Filtering

### Title: Removing Low-Quality Cells Systematically

**Bullet points:**
- Two filtering rounds: (1) upper bounds for doublets/debris, (2) lower bounds + complexity for empty droplets
- **Round 1 thresholds:** `nFeature_RNA < 4000`, `nCount_RNA < 20000`, `mitoRatio < 0.05`
- **Round 2 thresholds:** `log10GenesPerUMI > 0.82`, `nGene > 900`, `nUMI > 3200`
- Also remove BCR/TCR variable genes (`^Ig[hkl]v`, `^Tr[abdg][vjc]` regex patterns) to prevent clonotype-driven clustering
- Also remove genes expressed in <10 cells to reduce dimensionality and noise

**Filtering code:**

```r
# Round 1: Remove doublets and dead cells
Day0 = subset(Day0, nFeature_RNA < 4000 & nCount_RNA < 20000 & mitoRatio < .05)

# Round 2: Remove empty droplets and low-complexity cells
Day0 <- subset(x = Day0, subset = (log10GenesPerUMI > 0.82) & (nGene > 900) & (nUMI > 3200))

# Remove BCR/TCR genes (immune dataset-specific)
merged <- merged[!grepl("^Ig[hkl]v", rownames(merged), ignore.case = FALSE), ]

# Remove sparse genes (< 10 cells)
keep_genes <- Matrix::rowSums(counts > 0) >= 10
merged <- CreateSeuratObject(counts[keep_genes, ], meta.data = merged@meta.data)
```

**Teaching note:** The BCR/TCR gene removal step is frequently overlooked in published tutorials. In any B-cell or T-cell rich dataset, immunoglobulin variable genes (e.g., `Ighv1-26`, `Igkv4-55`) will be among the most variable genes in the dataset and will drive the first principal components. This causes cells to cluster by clonotype (receptor sequence) rather than functional transcriptional state. Removing them is both biologically justified and necessary for meaningful analysis.

---

## Slide 7 — Normalisation, Cell Cycle, and Integration

### Title: Making Cells from Different Samples Comparable

**Bullet points:**
- Cell cycle variation is a major confound — regress S.Score and G2M.Score covariates in SCTransform
- SCTransform (regularised negative binomial regression) is superior to log-normalisation for multi-sample data
- Run SCTransform **independently per sample** before integration — prevents cross-sample depth confounding
- Anchor-based integration identifies "biologically equivalent" cell pairs across samples as alignment anchors
- After integration: use `integrated` assay for clustering; `RNA` assay for differential expression

**SCTransform + Integration code:**

```r
# Score cell cycle (mouse gene names require str_to_title conversion)
s.genes <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)
merged <- CellCycleScoring(merged, g2m.features = g2m.genes, s.features = s.genes)

# SCTransform per sample (regress cell cycle)
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("S.Score", "G2M.Score"))
}

# Integration
bm.anchors <- FindIntegrationAnchors(object.list = bm.list,
                                     normalization.method = "SCT",
                                     anchor.features = bm.features)
bm.integrated <- IntegrateData(anchorset = bm.anchors, normalization.method = "SCT")
```

**Teaching note:** The `str_to_title()` conversion is a subtle but critical step for mouse data. Seurat's built-in `cc.genes` uses human HGNC capitalisation (all caps: `PCNA`); mouse Ensembl symbols use title case (`Pcna`). Without this conversion, cell cycle scoring will silently fail to score any genes. Integration can require 16–32+ GB RAM; the `remove()` call before `IntegrateData()` is essential for freeing memory.

---

## Slide 8 — Dimensionality Reduction, Clustering, and Annotation

### Title: Finding and Naming Cell Populations

**Bullet points:**
- PCA reduces ~3000 variable genes → 20–50 PCs; ElbowPlot guides PC number selection
- UMAP compresses PCA embeddings into 2D for visualisation; use PCA (not UMAP) for graph construction
- clustree tests 20 resolutions (0.1–2.0) simultaneously to identify stable clustering
- Resolution 0.1 selected — verified by sample-split and phase-split UMAP checks
- Annotation: canonical markers (Cd3e = T cell, Cd19 = B cell, Adgre1 = macrophage, Cd14 = monocyte)

**Clustering and annotation code:**

```r
bm.integrated <- RunPCA(bm.integrated)
bm.integrated <- RunUMAP(bm.integrated, dims = 1:20)
bm.integrated <- FindNeighbors(bm.integrated, dims = 1:20)
bm.integrated <- FindClusters(bm.integrated, resolution = seq(0.1, 2.0, 0.1))
Idents(bm.integrated) = "integrated_snn_res.0.1"

# Annotate with canonical markers
FeaturePlot(bm.integrated, features = c("Cd3e", "Cd19", "Adgre1", "Cd14"), ncol = 2)

# Find all cluster markers (negative binomial test for count data)
markers <- FindAllMarkers(bm.integrated, logfc.threshold = 0.25, test.use = "negbinom", only.pos = TRUE)
DoHeatmap(bm.integrated, features = top10$gene) + NoLegend()
```

**Teaching note:** The choice of `test.use = "negbinom"` is scientifically important. scRNA-seq counts are overdispersed relative to Poisson — they follow a negative binomial distribution. Using a Wilcoxon test (the default) is faster but ignores this structure. The negative binomial test is the biologically appropriate choice for count data. Emphasise that annotation is the most subjective step and requires domain knowledge — no algorithm fully replaces expert biological interpretation.

---

## Slide 9 — T-cell Subclustering, Trajectory, and CellChat

### Title: Advanced Downstream Analysis

**Bullet points:**
- Subclustering T cells (cluster 3) at resolution 0.2 reveals 5 distinct subtypes: CTL, Th17, Treg, IL-17a+ γδT, IL-17a− γδT
- Monocle3 trajectory orders cells in pseudotime — a proxy for developmental progression based on transcriptional distance
- Root node selection (biological start point) is an informed decision — choose the least differentiated state
- CellChat computes L-R communication probability from expression data using CellChatDB.mouse (~600 interactions)
- Day0 vs. Day14 CellChat comparison reveals gained/lost signalling pathways between time points

**T-cell subclustering + trajectory code:**

```r
# Subcluster T cells
T_cell = subset(seurat_integrated, subset = integrated_snn_res.0.1 == "3")
T_cell <- RunUMAP(T_cell, dims = 1:20); T_cell <- FindClusters(T_cell, resolution = seq(0.1, 2.0, 0.1))
Idents(T_cell) = "integrated_snn_res.0.2"

# Trajectory inference
bm.integrated.cds <- as.cell_data_set(T_cell, group.by = "integrated_snn_res.0.2")
bm.integrated.cds <- learn_graph(bm.integrated.cds, use_partition = TRUE)
bm.integrated.cds <- order_cells(bm.integrated.cds, reduction_method = "UMAP")
plot_cells(bm.integrated.cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
```

**T-cell subtype markers:**

```r
tcell_clusters <- list(
  Cytotoxic_T = c("Lef1", "Ms4a4b", "Cd8a"),
  Th17        = c("Tnfsf8", "Cxcr6", "Cd4", "Il17a", "Rora"),
  Treg        = c("Foxp3", "Ctla4"),
  gdT_IL17pos = c("Il23r", "Trgv2"),
  gdT_IL17neg = c("Birc5", "Top2a")
)
```

**Teaching note:** The inclusion of γδ T cells (identifiable by `Trgv2` and `Il23r`) is a distinctive feature of this dataset and highlights why scRNA-seq is valuable for immune profiling. These cells are rare, functionally distinct, and would be completely masked in bulk RNA-seq. The Monocle3 `order_cells()` step opens an interactive Shiny window — prepare learners for this and discuss how to biologically choose the root.

---

## Slide 10 — Summary, Limitations, and Next Steps

### Title: What We Learned and Where to Go Next

**Bullet points:**
- Complete scRNA-seq pipeline demonstrated from raw counts to cell communication networks
- Key output: annotated UMAP of immune cell types, T-cell developmental trajectory, comparative Day0 vs. Day14 communication atlas
- Limitations: (1) no automated doublet removal (DoubletFinder recommended), (2) manual root selection in trajectory is subjective, (3) CellChat communication is inferred, not experimentally validated
- Clustering is **unsupervised** — biological interpretation requires expert knowledge and literature validation
- Next steps in analysis: ligand-receptor validation (in vitro co-culture), spatial transcriptomics overlay, proteomics integration (CITE-seq)

**CellChat comparative summary code:**

```r
# Merge and compare
cellchat <- mergeCellChat(list(Day0 = cellchatDay0, Day14 = cellchatDay14), add.names = c("Day0","Day14"))

# Information flow comparison
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", stacked = FALSE, do.stat = TRUE)
gg1 + gg2

# Differential interactions per cell type
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), comparison = c(1,2), angle.x = 45)
```

**Extensions and next steps:**

| Next Analysis | Tool | What It Adds |
|---|---|---|
| Doublet removal | DoubletFinder / scDblFinder | More rigorous QC |
| Automated annotation | SingleR / Azimuth | Reference-based cell typing |
| RNA velocity | scVelo / velociraptor | Directionality in trajectories |
| Spatial context | Seurat + Visium | Where cells are in tissue |
| Multi-modal | CITE-seq / WNN | Surface protein + RNA |
| Gene regulatory networks | SCENIC | Transcription factor activity |

**Teaching note:** Close by emphasising that this pipeline represents current best practices (as of the time of the tutorial) but the field evolves rapidly. Encourage learners to check the Seurat, Monocle3, and CellChat documentation for updates. The most important skill is not memorising the code — it is understanding *why* each step exists so that learners can adapt the pipeline to new data types and biological questions.

---

## Teaching Deck Usage Notes

- **Lecture format:** Each slide is designed for approximately 8–10 minutes of teaching time (90 minutes total)
- **Code-along format:** Slides 5–9 can be run interactively, with learners executing each code block and interpreting the output
- **Assessment ideas:** Ask learners to modify QC thresholds and observe the downstream effects; ask them to identify which cluster is macrophages from the FeaturePlot
- **Prerequisites:** Basic R fluency (data frames, ggplot2), conceptual familiarity with RNA-seq

---

*End of 10-slide teaching deck*
