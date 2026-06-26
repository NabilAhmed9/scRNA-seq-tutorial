# Slide-by-Slide Knowledge Extraction

> **Source:** scRNA_seq_Presentation_pptx.pptx  
> **Author:** Farag Ibrahim, University of Gothenburg, Sweden  
> **Total slides:** 109  
> Slides with primarily visual/image content are noted where text extraction was limited.

---

## Slide 1 — Title Slide

**Title:** Introduction to Single-cell RNA-seq: Concepts, Applications, and Challenges  
**Summary:** Course title slide introducing Farag Ibrahim as Research Fellow and Bioinformatician at the University of Gothenburg.  
**Key Lesson:** This is an end-to-end tutorial integrating biology and bioinformatics for scRNA-seq.  
**Related Code:** N/A

---

## Slide 2 — Overview of Bulk RNA-Seq (Introduction)

**Title:** Overview of Bulk RNA-Seq  
**Summary:** Establishes the central dogma (DNA → RNA → Protein) as the conceptual foundation. Transcriptomics bridges genotype to phenotype.  
**Key Lesson:** RNA-seq is the primary tool for profiling the transcriptome; understanding bulk RNA-seq is a prerequisite for appreciating the need for single-cell resolution.  
**Related Code:** N/A

---

## Slide 3 — Bulk RNA-Seq: Importance

**Title:** Overview of Bulk RNA-Seq: Importance  
**Summary:** Each of the ~20,000 human genes is expressed in a tissue- and cell-type-specific manner. Tissue-specific expression is controlled by transcription factors and epigenetic mechanisms. Gene expression shifts during development, aging, stimulation, and disease.  
**Key Lesson:** Differential expression between healthy and diseased tissue reveals disease biomarkers. Epigenetic regulation (methylation, histone modification) modulates transcription.  
**Related Code:** N/A — conceptual motivation for the experiment

---

## Slide 4 — Bulk RNA-Seq: Workflow

**Title:** Overview of Bulk RNA-Seq: Workflow  
**Summary:** Visual diagram of the bulk RNA-seq wet-lab workflow (tissue dissociation → RNA extraction → library preparation → sequencing). Primarily image content.  
**Key Lesson:** The bulk workflow averages signal across the entire cell population in the input material.  
**Related Code:** N/A

---

## Slide 5 — Bulk RNA-Seq: Applications

**Title:** Applications of RNA-seq  
**Summary:** Five major application domains:  
- Differential Expression Analysis (treated vs. control)  
- Biomarker Discovery (RNA signatures for disease/prognosis)  
- Disease Research (cancer transcriptomes, oncogenes, fusion transcripts)  
- Personalised Medicine (patient-specific expression profiles)  
- Developmental Biology (temporal differentiation mapping)  
**Key Lesson:** RNA-seq is the workhorse of transcriptomics; these applications motivate upgrading to single-cell resolution.  
**Related Code:** Conceptual — motivates the entire pipeline

---

## Slide 6 — Bulk vs. scRNA-seq (Key Comparison)

**Title:** Bulk RNA-Seq vs. Single-cell RNA  
**Summary:** Side-by-side comparison:  
- **Bulk:** averages expression across millions of cells; masks heterogeneity; cannot detect rare cells  
- **scRNA-seq:** profiles each cell individually; reveals cell-type heterogeneity and rare populations; enables clustering, lineage tracing, and trajectory inference  
- **Platforms:** 10x Genomics, Drop-seq, Smart-seq3  
**Key Lesson:** Bulk RNA-seq collapses signal from thousands of distinct cell types into a single average; scRNA-seq resolves this at single-cell resolution.  
**Related Code:** Motivates entire pipeline

---

## Slide 7 — Bulk vs. scRNA-seq (Visual Comparison)

**Title:** Bulk RNA-Seq vs. Single-cell RNA (continued)  
**Summary:** Primarily a visual illustration of the resolution difference. Bulk = average; scRNA-seq = individual cell profiles.  
**Key Lesson:** Visual reinforcement of the averaging problem in bulk RNA-seq.  
**Related Code:** N/A

---

## Slide 8 — Why Single-cell RNA-seq?

**Title:** Why Single-cell RNA-seq?  
**Summary:** Four core motivations:  
1. **Cell Type Identification** — classify and annotate distinct populations using transcriptomic signatures  
2. **Rare Cell Detection** — identify low-abundance cell types (stem cells, circulating tumour cells) invisible to bulk methods  
3. **Differentiation & Dynamics** — reconstruct developmental trajectories and capture temporal gene expression changes  
4. **Differential Gene Expression** — perform cell-type-specific DEG analysis at single-cell resolution  
**Key Lesson:** Tissues are composed of heterogeneous populations; scRNA-seq is required to resolve this complexity.  
**Related Code:** Motivates the T-cell subclustering and trajectory inference sections

---

## Slide 9 — Detailed Comparison Table: Bulk vs. scRNA-seq

**Title:** Comparison with Bulk RNA-seq  
**Summary:** Structured comparison across five axes:

| Axis | Bulk RNA-seq | scRNA-seq |
|---|---|---|
| Resolution | All cells averaged | Each cell individually |
| Heterogeneity | Hidden | Revealed |
| Rare cells | Not detectable | High sensitivity |
| Cost & complexity | Lower / simpler | Higher / more complex |
| Data volume | Manageable | Massive, high-dimensional |

**Key Lesson:** scRNA-seq generates orders of magnitude more data and requires specialised computational pipelines.  
**Related Code:** Justifies the computational infrastructure choices (Seurat, SCTransform)

---

## Slide 10 — Challenges Section Header

**Title:** Challenges of scRNA-seq Analysis  
**Summary:** Section break introducing four major analytical challenges:  
1. Large-scale data volume  
2. Low sequencing depth per cell  
3. Batch effects  
4. Sequencing artefacts  
**Key Lesson:** scRNA-seq is technically demanding — understanding its failure modes is as important as running the analysis.  
**Related Code:** N/A (conceptual context for filtering and integration decisions)

---

## Slide 11 — Challenge 1: Data Volume

**Title:** Challenge 1 — Large-scale Data Volume  
**Summary:** Modern experiments profile 10,000–1,000,000+ cells; each cell has 20,000–30,000 potential gene measurements. Raw FASTQ files reach ~30 GB per sample. Computational requirements include high-RAM systems, large storage, and extended runtimes.  
**Key Lesson:** Memory management is critical. Sparse matrix representation (as used by Seurat and `Matrix::`) is essential for scalability.  
**Related Code:** `library(Matrix)` — sparse matrix handling throughout

---

## Slide 12 — Challenge 2: Low Sequencing Depth & Dropouts

**Title:** Challenge 2 — Low Sequencing Depth per Cell  
**Summary:** In droplet-based scRNA-seq, each cell is sequenced at much lower depth than bulk. This produces **dropouts** — genes with zero counts in a cell even when biologically expressed, caused by low mRNA capture efficiency.  
**Key Lesson:** Dropouts create a sparse, zero-inflated count matrix. This is a fundamental property of scRNA-seq data, not a sequencing failure. Robust statistical methods (e.g., negative binomial models) are required.  
**Related Code:** `test.use = "negbinom"` in `FindAllMarkers()` — accounts for dropout-driven zero inflation

---

## Slide 13 — Challenge 3: Batch Effects

**Title:** Challenge 3 — Batch Effects and Experimental Design  
**Summary:**  
- **Causes:** cells processed on different days/labs; reagent lot variation; different sequencing runs or instruments  
- **Detection:** PCA/UMAP where samples cluster by batch (not biology); hierarchical clustering; QC metric comparison  
- **Mitigation:** process all samples simultaneously; include biological replicates per batch; apply computational data integration  
**Key Lesson:** Batch effects are among the most common confounders in multi-sample scRNA-seq. They must be corrected before biological interpretation.  
**Related Code:** SCTransform per sample + `FindIntegrationAnchors()` + `IntegrateData()`

---

## Slide 14 — Batch Effects (Visual)

**Title:** Batch Effects and Experimental Design (continued)  
**Summary:** Primarily visual — illustrates UMAP plots before and after batch correction, showing sample-driven vs. biology-driven clustering.  
**Key Lesson:** Successful integration should yield interleaved cells from Day0 and Day14, not separate sample clusters.  
**Related Code:** `DimPlot(split.by = "sample")` — used to verify integration quality

---

## Slide 15 — Challenge 4: Sequencing Artefacts

**Title:** Challenge 4 — Sequencing Artefacts  
**Summary:** Four non-biological sources of variation:  
1. **Capture Efficiency** — not all mRNA captured → sparse counts and high dropout  
2. **Library Quality** — PCR bias, RNA degradation, low input material  
3. **Amplification Bias** — preferential amplification of certain sequences; corrected by UMIs  
4. **Batch Effects** — systematic differences from processing time, lab, or instrument  
**Key Lesson:** UMIs (Unique Molecular Identifiers) deconvolve PCR duplicates from true biological reads and are the standard unit of count in 10X data.  
**Related Code:** `nCount_RNA` (UMI counts), `nFeature_RNA` (gene counts) — primary QC metrics

---

## Slide 16 — 3′ vs Full-Length RNA Sequencing

**Title:** 3′ vs Full-Length RNA Sequencing  
**Summary:**  
| Feature | 3′ Sequencing (e.g., 10X) | Full-Length (e.g., SMART-seq) |
|---|---|---|
| Coverage | 3′ end only (~200–500 bp) | Entire transcript |
| Throughput | Very high | Lower |
| Isoform detection | No | Yes |
| Cost | Lower | Higher |
| scRNA-seq compatible | Yes (standard) | Limited |  
**Key Lesson:** 10X Genomics uses 3′ end sequencing, capturing poly-A tails. This limits isoform analysis but enables high-throughput profiling of thousands of cells.  
**Related Code:** Data loaded via `Read10X()` — 10X Genomics 3′ count matrix format

---

## Slide 17 — 10X Genomics vs SMART-seq Platforms

**Title:** 3′ vs Full-Length — Platform Comparison  
**Summary:** Platform images — 10X Genomics Chromium Controller/Chromium X vs. SMART-seq. Primarily visual.  
**Key Lesson:** The Chromium platform is the dominant droplet-based system used in high-throughput scRNA-seq.  
**Related Code:** N/A — platform context

---

## Slide 18 — 10X Genomics Microfluidics

**Title:** 10X Genomics Microfluidics  
**Summary:** Illustration of the GEM (Gel Bead in Emulsion) microfluidic system — cells and gel beads are co-encapsulated in oil droplets; each GEM contains one cell and one barcoded gel bead.  
**Key Lesson:** The barcode on each gel bead tags all mRNA from a single cell, enabling deconvolution of cell identities after pooled sequencing.  
**Related Code:** Barcodes are the column indices of the count matrix read by `Read10X()`

---

## Slide 19 — 10X Genomics Poly(dT) Capture

**Title:** 10X Genomics — Poly(dT)VN Priming  
**Summary:** The gel bead oligonucleotide contains an oligo-dT priming site that hybridises to the poly-A tail of mRNA molecules. Each transcript is tagged with a cell barcode and a UMI.  
**Key Lesson:** The structure `[16-nt Cell Barcode] + [12-nt UMI] + [Poly(dT)VN]` on the bead oligonucleotide is what generates the count matrix: each unique barcode = one cell; each unique UMI = one captured transcript.  
**Related Code:** `Read10X()` reads the barcode and feature files produced by CellRanger

---

## Slide 20 — Analysis Workflow Overview

**Title:** Analysis Workflow  
**Summary:** High-level overview of the bioinformatics pipeline from CellRanger output to biological interpretation. Primarily a workflow diagram.  
**Key Lesson:** The analysis has a canonical order: raw data → QC → normalisation → integration → dimensionality reduction → clustering → annotation → downstream.  
**Related Code:** Maps to the full code structure in `analysis_code.txt`

---

## Slide 21 — CellRanger Output: Features and Barcodes Files

**Title:** CellRanger Output Files — Features and Barcodes  
**Summary:** CellRanger produces three files per sample: `barcodes.tsv.gz` (cell barcodes), `features.tsv.gz` (gene identifiers), and `matrix.mtx.gz` (sparse count matrix in Market Exchange format).  
**Key Lesson:** These three files together define the count matrix loaded by `Read10X()`. Features = rows (genes), barcodes = columns (cells), matrix = UMI counts.  
**Related Code:** `data <- Read10X("/path/to/Day0/")`

---

## Slide 22 — CellRanger Output: Count Matrix

**Title:** CellRanger Output — Count Matrix  
**Summary:** The count matrix is sparse (most entries are zero). Market Exchange format (.mtx) stores only non-zero entries, making it memory-efficient.  
**Key Lesson:** A typical experiment produces a matrix of ~20,000 genes × thousands of cells, with >90% zeros — hence the requirement for sparse matrix libraries.  
**Related Code:** `CreateSeuratObject(counts = data)` — converts the sparse matrix into a Seurat object

---

## Slide 23 — Seurat Object Structure

**Title:** Analysis Workflow — Seurat Object  
**Summary:** The Seurat object is the central data structure. It stores: count matrices (assays), cell metadata, dimensional reductions, graphs, and cluster assignments.  
**Key Lesson:** All downstream operations (normalisation, integration, clustering) modify and extend this object. Understanding its slot structure is essential.  
**Related Code:** `Day0 <- CreateSeuratObject(counts = data)` and all subsequent `Day0@...` access

---

## Slide 24 — (Transition / Visual Slide)

**Title:** (No extractable title text)  
**Summary:** Visual transition slide — presumably showing Seurat object structure or QC overview.  
**Key Lesson:** *Unclear from source files* — likely a diagram of the Seurat object hierarchy.  
**Related Code:** N/A

---

## Slide 25 — QC: Generating the Seurat Object

**Title:** QC Step 1 — How to Generate a Seurat Object  
**Summary:** Shows the code workflow: `Read10X()` → `CreateSeuratObject()` → initial `subset()` to remove cells with fewer than 200 detected genes.  
**Key Lesson:** The minimum gene threshold (≥200 genes per cell) eliminates empty droplets and debris that were captured without a real cell.  
**Related Code:** `Day0 <- CreateSeuratObject(counts = data)` and `Day0 <- subset(Day0, subset = nFeature_RNA >= 200)`

---

## Slide 26 — QC: Mitochondrial Percentage

**Title:** QC (Cell-level) Step 2 — Mitochondrial Content  
**Summary:** `PercentageFeatureSet()` calculates the fraction of total UMIs mapping to mitochondrial genes (identified by the `^mt-` prefix for mouse). High mitochondrial content (>5%) indicates cell membrane rupture and cytoplasmic RNA loss — a marker of dying or damaged cells.  
**Key Lesson:** Mitochondrial RNA is retained in mitochondria even after cell lysis; elevated mitochondrial fraction indicates the cell's cytoplasmic RNA has leaked out, leaving only organelle-contained transcripts.  
**Related Code:** `Day0$mitoRatio <- PercentageFeatureSet(object = Day0, pattern = "^mt-") / 100`

---

## Slide 27 — QC: Genes per UMI (Complexity Score) — Part 1

**Title:** QC (Cell-level) Step 3 — Genes per UMI  
**Summary:** Introduces `log10GenesPerUMI` as a transcriptomic complexity score. Computed as `log10(nFeature_RNA) / log10(nCount_RNA)`.  
**Key Lesson:** A healthy, transcriptionally complex cell should detect many distinct genes relative to total UMI count. A low complexity score suggests a homogeneous cell (e.g., red blood cell, empty droplet, or multiplet).  
**Related Code:** `Day0$log10GenesPerUMI <- log10(Day0$nFeature_RNA) / log10(Day0$nCount_RNA)`

---

## Slide 28 — QC: Genes per UMI (Complexity Score) — Part 2

**Title:** QC (Cell-level) Step 4 — Genes per UMI (continued)  
**Summary:** Continuation of complexity score concept, likely with visual examples of low vs. high complexity cells.  
**Key Lesson:** The threshold `log10GenesPerUMI > 0.82` is the filtering cutoff applied in the code to remove low-complexity cells.  
**Related Code:** `Day0 <- subset(x = Day0, subset = (log10GenesPerUMI > 0.82))`

---

## Slide 29 — QC: Sequencing Depth (nUMI)

**Title:** QC (Cell-level) Step 5 — Sequencing Depth  
**Summary:** Density plot of UMI counts per cell (`nUMI`), plotted on a log10 scale. A vertical dashed line marks the suggested lower cutoff at nUMI = 2,000. Cells with very low UMIs are likely empty droplets.  
**Key Lesson:** The density distribution should be unimodal. A bimodal distribution or long left tail suggests a mixed population of real cells and debris. Final cutoff in code: `nUMI > 3200`.  
**Related Code:** `metadata %>% ggplot(aes(x = nUMI, ...)) + geom_density()`

---

## Slide 30 — QC: Gene Counts per Cell

**Title:** QC (Cell-level) Step 6 — Gene Counts per Cell  
**Summary:** Density plot of detected gene counts (`nGene`) per cell, log10 scaled. Suggested lower threshold shown at nGene = 1,000.  
**Key Lesson:** Cells with very few detected genes are likely empty droplets or debris. Cells with extremely high gene counts may be doublets (two cells captured in one droplet).  
**Related Code:** `metadata %>% ggplot(aes(x = nGene, ...))` and filter `nGene > 900`

---

## Slide 31 — QC: Mitochondrial Ratio

**Title:** QC (Cell-level) Step 7 — Mitochondrial Ratio  
**Summary:** Density plot of mitochondrial fraction per cell. Most healthy cells cluster at very low mitochondrial ratios (<0.05).  
**Key Lesson:** The mitochondrial ratio threshold (mitoRatio < 0.05) removes cells where cytoplasmic RNA was lost due to membrane damage — a proxy for cell death during sample preparation.  
**Related Code:** `Day0 = subset(Day0, nFeature_RNA < 4000 & nCount_RNA < 20000 & mitoRatio < .05)`

---

## Slide 32 — QC: Cells Failing Multiple Metrics

**Title:** QC (Cell-level) Step 8 — Identifying Low-Quality Cells  
**Summary:** Scatter plot of `nUMI` vs. `nGene` coloured by mitochondrial ratio. Cells failing multiple thresholds simultaneously (low UMI, low genes, high mito) are identified in the lower-left quadrant with high colour intensity.  
**Key Lesson:** QC metrics should be evaluated jointly, not in isolation. Red reference lines for nUMI = 3200 and nGene = 1000 visually delineate the filtering boundary.  
**Related Code:** `metadata %>% ggplot(aes(x = nUMI, y = nGene, color = mitoRatio)) + geom_point()`

---

## Slide 33 — (QC Visual Continuation)

**Summary:** Likely shows filtered vs. unfiltered comparison violin plots.  
**Key Lesson:** *Unclear from source files* — visual validation of filtering.  
**Related Code:** Second `VlnPlot()` call after filtering

---

## Slides 34–40 — QC Workflow Continuation, Day14, and Merging

**Summary:** These slides cover the Day14 sample QC repetition and the merging workflow. The same QC steps applied to Day0 are applied to Day14. After filtering, cell barcodes are prefixed (`RenameCells()`) and both objects are merged.  
**Key Lesson:** Consistent QC application across all samples is critical. Merging adds sample identity to each cell while preserving cell-level metadata.  
**Related Code:**  
```r
Day0 <- RenameCells(Day0, add.cell.id = "Day0")
Day14 <- RenameCells(Day14, add.cell.id = "Day14")
merged <- merge(x = Day0, y = c(Day14))
```

---

## Slides 38–40 — BCR/TCR Gene Removal and Gene Filtering

**Summary:** Before analysis, immunoglobulin (Ig) variable, joining, and constant genes, as well as T-cell receptor genes, are removed using regex patterns. Additionally, genes expressed in fewer than 10 cells are removed to reduce noise.  
**Key Lesson:** In B-cell or T-cell datasets, hypervariable Ig/TCR genes dominate the variable gene list and can distort clustering — removing them reveals biologically meaningful transcriptional programmes. Gene-level filtering reduces dimensionality and removes noise-prone sparse features.  
**Related Code:**  
```r
merged <- merged[!grepl("^Ig[hkl]v", rownames(merged), ignore.case = FALSE), ]
keep_genes <- Matrix::rowSums(nonzero) >= 10
```

---

## Slides 41–43 — Count Normalisation Concepts

**Title:** Cell Cycle: Count Normalisation  
**Summary:** Three sources of technical variation in raw UMI counts are explained:  
A. **Sequencing depth (library size)** — cells sequenced more deeply have higher raw UMI counts  
B. **Gene length** — longer genes capture more reads per transcript (less relevant for 3′ UMI data)  
C. **RNA composition** — cells expressing highly abundant genes have relatively lower counts for other genes  
**Key Lesson:** Log-normalisation to 10,000 UMIs (CPM-like) corrects for library size but not composition. SCTransform (regularised negative binomial regression) more robustly handles both depth and composition effects.  
**Related Code:** `NormalizeData(merged)` followed by `SCTransform()` per sample

---

## Slide 44 — Cell Cycle Scoring

**Title:** Cell Cycle Scoring  
**Summary:** Cell cycle phase (S, G2M) is a major source of transcriptional variation. `CellCycleScoring()` assigns each cell an S-phase and G2M-phase score based on known cell-cycle gene signatures. Setting `set.ident = TRUE` labels each cell with its predicted phase.  
**Key Lesson:** If cell cycle drives the first principal components, cells will cluster by phase rather than biology. Cell cycle scores must be regressed out during SCTransform to focus on biologically relevant variation.  
**Related Code:**  
```r
merged <- CellCycleScoring(merged, g2m.features = g2m.genes, s.features = s.genes, set.ident = TRUE)
```

---

## Slides 45–51 — Scaling Concepts

**Title:** Cell Cycle Scoring: Scaling  
**Summary:** Scaling centres each gene's expression to mean = 0 and variance = 1 across all cells, making gene contributions to PCA equal regardless of absolute expression magnitude. This is a prerequisite for PCA.  
**Key Lesson:** Unscaled data gives disproportionate weight to highly expressed genes in PCA. Scaling (z-score transformation) ensures all genes contribute equally, but note that `ScaleData()` is superseded by `SCTransform()` for the final normalisation pipeline.  
**Related Code:** `merged <- ScaleData(merged)`

---

## Slide 52 — Cell Cycle PCA Visualisation

**Title:** Cell Cycle: PCA and Visualisation  
**Summary:** After scoring and scaling, PCA is run and cells are coloured by cell-cycle phase (`group.by = "Phase"`). If the PCA separates cells by phase, regression is necessary.  
**Key Lesson:** `DimPlot(group.by = "Phase", split.by = "Phase")` visually evaluates whether cell cycle is a dominant source of variation.  
**Related Code:** `DimPlot(merged, reduction = "pca", group.by = "Phase", split.by = "Phase")`

---

## Slides 54–58 — SCTransform Normalisation

**Title:** SCTransform Normalisation  
**Summary:** SCTransform applies regularised negative binomial regression to model and remove the relationship between gene expression and sequencing depth. Cell-cycle scores (S.Score, G2M.Score) are included as covariates to be regressed out. Each sample is normalised independently before integration.  
**Key Lesson:** SCTransform is the recommended normalisation method for Seurat v5 multi-sample workflows. It replaces the `NormalizeData() → ScaleData()` pipeline for integration and produces a `SCT` assay slot in the Seurat object.  
**Related Code:**  
```r
split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("S.Score", "G2M.Score"))
```

---

## Slides 60–67 — Data Integration

**Title:** Data Integration  
**Summary:** Seurat's anchor-based integration:  
1. `SelectIntegrationFeatures()` — identifies highly variable genes common across samples  
2. `PrepSCTIntegration()` — prepares SCT-normalised objects for anchor finding  
3. `FindIntegrationAnchors()` — identifies mutual nearest-neighbour (MNN) pairs across samples as "anchors"  
4. `IntegrateData()` — uses anchors to project datasets into a shared embedding, removing batch-specific variation  
**Key Lesson:** The anchor-based approach identifies cells that are biologically similar across batches and uses them to align the datasets. The output `integrated` assay should be used for clustering; the `RNA` assay retains original counts for differential expression.  
**Related Code:**  
```r
bm.anchors <- FindIntegrationAnchors(object.list = bm.list, normalization.method = "SCT")
bm.integrated <- IntegrateData(anchorset = bm.anchors, normalization.method = "SCT")
```

---

## Slides 69–70 — PCA and Elbow Plot

**Title:** Dimensionality Reduction: PCA and Elbow Plot  
**Summary:** PCA reduces the ~2000-gene feature space to ordered principal components (PCs) capturing the largest sources of variance. The elbow plot shows variance explained per PC; the "elbow" point indicates where additional PCs contribute diminishingly.  
**Key Lesson:** The number of PCs to retain (here: 20) must be chosen carefully — too few loses biological signal; too many retains technical noise. The elbow plot guides this decision.  
**Related Code:** `ElbowPlot(bm.integrated, ndims = 50)` and `RunPCA()`

---

## Slide 71 — UMAP

**Title:** Dimensionality Reduction: UMAP  
**Summary:** UMAP (Uniform Manifold Approximation and Projection) compresses PCA embeddings (cells × 15–20 PCs) into a 2D representation while preserving local neighbourhood relationships. Similar cells end up close together in 2D space.  
**Key Lesson:** UMAP is the standard 2D visualisation for scRNA-seq. Critically, distances between UMAP clusters are not directly interpretable biologically — only local structure (within-cluster and nearest-cluster relationships) is meaningful.  
**Related Code:** `bm.integrated <- RunUMAP(object = bm.integrated, dims = 1:20)`

---

## Slide 72 — PCA vs. UMAP Comparison

**Title:** Differences Between PCA and UMAP  
**Summary:** PCA is a linear transformation preserving global variance; UMAP is a non-linear transformation preserving local neighbourhood topology. UMAP better separates discrete cell populations.  
**Key Lesson:** Always use PCA for integration and neighbourhood graph construction; use UMAP only for visualisation.  
**Related Code:** `FindNeighbors(reduction = "pca")` — uses PCA space for graph, not UMAP

---

## Slide 73 — Finding Similar Cells (KNN Graph)

**Title:** Dimensionality Reduction: Finding Similar Cells  
**Summary:** `FindNeighbors()` builds a k-nearest-neighbour (KNN) graph in PCA space. Each cell is connected to its k most similar neighbours.  
**Key Lesson:** The KNN graph is the mathematical substrate for community detection (clustering). Resolution of the graph determines the number of clusters.  
**Related Code:** `bm.integrated <- FindNeighbors(bm.integrated, reduction = "pca", dims = 1:20)`

---

## Slide 74–75 — Clustering Stability (clustree)

**Title:** Dimensionality Reduction: Clustering Stability  
**Summary:** `clustree` visualises how cluster assignments change as clustering resolution is increased from 0.1 to 2.0. A stable resolution shows clusters that are reproducible across small changes in resolution.  
**Key Lesson:** Choose the lowest resolution at which biologically distinct populations emerge as separate clusters. Overly high resolution creates artifactual splits; too low merges distinct populations.  
**Related Code:**  
```r
clustree(x = bm.integrated, prefix = "integrated_snn_res.")
Idents(bm.integrated) = "integrated_snn_res.0.1"
```

---

## Slides 77–79 — Integration QC: UMAP Visualisations

**Title:** Integration QC — Checking Clustering Results  
**Summary:** Three UMAP checks after integration:  
1. Overall cluster labels (`DimPlot(label = TRUE)`)  
2. Clusters split by sample — verifies that Day0 and Day14 cells intermix within each cluster  
3. Clusters split by cell-cycle phase — confirms that cell cycle no longer drives clustering after regression  
**Key Lesson:** A successful integration results in cells from different samples co-clustering by cell type, not sample of origin.  
**Related Code:**  
```r
DimPlot(bm.integrated, split.by = "sample", label = TRUE)
DimPlot(bm.integrated, split.by = "Phase", label = TRUE)
```

---

## Slides 80–84 — Cluster Annotation

**Title:** Cluster Annotation  
**Summary:** Cell types are identified using canonical lineage markers:  
- `Cd3e` — T cells  
- `Cd19` — B cells  
- `Adgre1` (F4/80) — macrophages  
- `Cd14` — monocytes  
`FindAllMarkers()` identifies the top differentially expressed genes per cluster using a negative binomial test. Top 10 markers per cluster are visualised as a heatmap and DotPlot.  
**Key Lesson:** Annotation is the most biologically demanding step — it requires domain knowledge to match gene signatures to cell types. Automated tools (e.g., SingleR, Azimuth) can assist but manual curation is often required.  
**Related Code:**  
```r
FeaturePlot(bm.integrated, features = c("Cd3e", "Cd19", "Adgre1", "Cd14"))
markers <- FindAllMarkers(object = bm.integrated, logfc.threshold = 0.25, test.use = "negbinom")
DoHeatmap(seurat_integrated, features = top10$gene)
```

---

## Slide 85 — T-cell Subclustering

**Title:** Subset T Cells and Recluster  
**Summary:** Cluster 3 (identified as T cells by Cd3e expression) is extracted and subjected to independent re-clustering with `RunUMAP()` + `FindNeighbors()` + `FindClusters()`. Resolution 0.2 is selected.  
**Key Lesson:** Subclustering allows fine-grained resolution of heterogeneity within a broad cell type (e.g., Th17, Tregs, γδ T cells, CTLs) that would be invisible at the whole-dataset level.  
**Related Code:**  
```r
T_cell = subset(seurat_integrated, subset = integrated_snn_res.0.1 == "3")
T_cell <- RunUMAP(T_cell, dims = 1:20)
```

---

## Slide 86 — T-cell Marker Visualisation

**Title:** Visualise Marker Genes for T-cell Clusters  
**Summary:** Defines T-cell subtype markers:  
- **Cytotoxic T cells:** Lef1, Ms4a4b, Cd8a  
- **Th17 cells:** Tnfsf8, Cxcr6, Cd4, Il17a, Rora  
- **Regulatory T cells:** Foxp3, Ctla4  
- **IL-17a+ γδ T cells:** Il23r, Trgv2  
- **IL-17a− γδ T cells:** Birc5, Top2a  
**Key Lesson:** Canonical markers remain the gold standard for T-cell subtype annotation. γδ T cells and Tregs can be distinguished from conventional αβ T cells even at single-cell resolution.  
**Related Code:** `DotPlot(T_cell, features = tcell_clusters)`

---

## Slides 88–91 — Trajectory Inference with Monocle3

**Title:** Trajectory Inference — What Is It? / Building with Monocle3  
**Summary:** Trajectory inference (also called pseudotime analysis) orders cells along a continuous developmental path based on transcriptional similarity. Monocle3 learns a principal graph through the UMAP embedding, then assigns pseudotime values — a proxy for developmental progression.  
**Key Lesson:** Pseudotime is not clock time; it represents inferred transcriptional progression. Choosing the root node (the biological starting state) is a critical, biology-informed decision.  
**Related Code:**  
```r
bm.integrated.cds <- as.cell_data_set(T_cell, group.by = "integrated_snn_res.0.2")
bm.integrated.cds <- learn_graph(bm.integrated.cds)
bm.integrated.cds <- order_cells(bm.integrated.cds)
plot_cells(cds = bm.integrated.cds, color_cells_by = "pseudotime")
```

---

## Slides 92–109 — CellChat Analysis

**Title:** Cell–Cell Communication Analysis (CellChat Package)  
**Summary:** CellChat quantifies ligand–receptor (L-R) mediated intercellular communication. The pipeline:  
1. Extract normalised expression matrix and cluster identities per sample  
2. Create CellChat objects (`createCellChat()`)  
3. Load CellChatDB.mouse ligand–receptor database  
4. Identify overexpressed L-R pairs (`identifyOverExpressedGenes/Interactions()`)  
5. Compute communication probabilities (`computeCommunProb()`)  
6. Aggregate into a network (`aggregateNet()`)  
7. Visualise with chord diagrams and bubble plots (`netVisual_circle()`, `netVisual_bubble()`)  
8. **Merge** Day0 and Day14 CellChat objects (`mergeCellChat()`)  
9. Compare: total interactions, interaction strength, differential interactions per cell type  
10. Cluster signalling pathways by functional and structural similarity  
11. Rank pathway information flow (`rankNet()`)  
12. Compare outgoing signalling patterns across conditions with heatmaps  

**Key Lessons:**  
- CellChat infers communication from expression data alone — it does not require spatial or co-culture data  
- Comparing Day0 vs. Day14 networks reveals which signalling pathways are gained, lost, or rewired during the biological process  
- `netAnalysis_computeCentrality()` identifies which cell types are dominant senders or receivers  
- Functional similarity clustering groups pathways with similar communication patterns regardless of specific L-R pairs involved  
**Related Code:** Full CellChat section (Day0 + Day14 creation, merging, comparison, and visualisation)
