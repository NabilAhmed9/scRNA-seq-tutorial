# scRNA-seq Analysis Tutorial: From Raw Counts to Cell Communication
 
**Dataset:** GSE142564 (Day 0 and Day 14 bone marrow samples, mouse)  
**Language:** R  
**Analysis Packages:** Seurat, Monocle3, CellChat  

---

## Overview

This repository provides a complete, teaching-oriented single-cell RNA sequencing (scRNA-seq) tutorial. It integrates executable R analysis code with conceptual knowledge extracted from an accompanying lecture presentation, spanning 109 slides. The tutorial covers the full analytical trajectory — from raw 10X Genomics count matrices through quality control, normalization, data integration, dimensionality reduction, cell clustering, trajectory inference, and cell–cell communication analysis.

The biological context is a **bone marrow immune profiling study** comparing two time points (Day 0 vs. Day 14), with a particular focus on T-cell heterogeneity and intercellular signalling dynamics in what appears to be a B-cell precursor acute lymphoblastic leukaemia (BCP-ALL) or related haematopoietic model system.

---

## Tutorial Goals

1. Understand the conceptual differences between bulk RNA-seq and scRNA-seq.
2. Learn how 10X Genomics droplet-based capture generates count matrices.
3. Apply rigorous, multi-metric quality control to exclude low-quality cells.
4. Perform SCTransform normalization with cell-cycle regression.
5. Integrate multi-sample data using Seurat's anchor-based integration framework.
6. Execute dimensionality reduction (PCA → UMAP) and assess clustering stability with clustree.
7. Annotate cell clusters using canonical lineage marker genes.
8. Subset and re-cluster a T-cell population; infer developmental trajectories with Monocle3.
9. Quantify and compare cell–cell communication networks across time points using CellChat.

---

## Repository Structure

```
scrnaseq-tutorial/
├── README.md                          ← This file
├── LICENSE                            ← MIT License
├── .gitignore                         ← R/RStudio ignores
├── data/
│   └── README.md                      ← Data acquisition instructions (GEO: GSE142564)
├── code/
│   └── analysis_code.txt              ← Original analysis code (preserved)
├── docs/
│   ├── slide_by_slide_knowledge.md    ← Slide-by-slide extracted knowledge
│   ├── code_explanation.md            ← Section-by-section code walkthrough
│   ├── tutorial_workflow.md           ← Full sequential tutorial workflow
│   └── github_project_summary.md     ← Project summary for sharing
├── slides/
│   └── 10_slide_teaching_deck.md      ← Condensed 10-slide teaching deck
└── environment/
    └── requirements_or_dependencies.md ← R packages and dependencies
```

---

## How to Use This Project

### For Learners

1. **Start with the concepts:** Read `docs/slide_by_slide_knowledge.md` to build foundational understanding of scRNA-seq theory.
2. **Follow the workflow:** Work through `docs/tutorial_workflow.md` step by step, running code as you go.
3. **Understand the code:** Use `docs/code_explanation.md` alongside the actual code in `code/analysis_code.txt`.
4. **Use the teaching deck:** `slides/10_slide_teaching_deck.md` is a condensed lecture guide for presenting or self-teaching.

### For Instructors

- The 10-slide deck in `slides/` can be adapted into presentation software.
- The `docs/tutorial_workflow.md` provides a structured lab exercise template.
- Each code section is independently runnable if input data are available.

---

## Expected Input Files

The analysis uses 10X Genomics CellRanger output files from GEO accession **GSE142564**:

```
GSE142564_RAW/
├── Day0/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── Day14/
    ├── barcodes.tsv.gz
    ├── features.tsv.gz
    └── matrix.mtx.gz
```

Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142564

---

## Analysis Workflow Summary

| Stage | Key Action | Primary Tool |
|---|---|---|
| 1. Data Loading | Read 10X CellRanger output | `Seurat::Read10X()` |
| 2. QC (per sample) | Filter low-quality cells | `Seurat::subset()` |
| 3. Metadata | Organise cell barcodes and metrics | `tidyverse` |
| 4. Merge | Combine Day0 + Day14 | `Seurat::merge()` |
| 5. Gene Filtering | Remove Ig/TCR genes; sparse genes | Base R / Matrix |
| 6. Cell Cycle | Score and regress S/G2M phases | `Seurat::CellCycleScoring()` |
| 7. Normalisation | SCTransform per sample | `Seurat::SCTransform()` |
| 8. Integration | Anchor-based cross-sample integration | `Seurat::IntegrateData()` |
| 9. Dim. Reduction | PCA → UMAP | `Seurat::RunPCA/RunUMAP()` |
| 10. Clustering | Leiden community detection + clustree | `Seurat::FindClusters()` + `clustree` |
| 11. Annotation | Marker genes + FindAllMarkers | `Seurat::FindAllMarkers()` |
| 12. Subclustering | T-cell subset analysis | Seurat |
| 13. Trajectory | Pseudotime ordering | `monocle3` |
| 14. Cell Comms | Ligand–receptor network comparison | `CellChat` |

---

## Who This Project Is For

- Graduate students and postdocs beginning scRNA-seq analysis
- Bioinformatics course instructors seeking a structured tutorial template
- Immunologists and cancer biologists learning computational methods
- R users transitioning from bulk RNA-seq to single-cell approaches

---

## Learning Outcomes

After completing this tutorial, learners will be able to:

- Explain why scRNA-seq reveals biology that bulk RNA-seq cannot
- Apply and interpret multi-metric QC filtering decisions
- Understand the rationale for SCTransform over log-normalisation in multi-sample settings
- Implement Seurat anchor-based integration to remove batch effects
- Interpret UMAP embeddings and clustree stability plots
- Annotate clusters using canonical marker genes and `FindAllMarkers()`
- Perform subset re-clustering and trajectory inference on a cell population of interest
- Construct and compare CellChat communication networks across experimental conditions

---



Dataset: *GEO accession GSE142564*
