# GitHub Project Summary

## Project Purpose

This repository is a **complete, teaching-oriented single-cell RNA-seq analysis tutorial** built on a real bone marrow immune profiling dataset (GEO: GSE142564). It integrates executable R analysis code with conceptual knowledge from a 109-slide lecture presentation, creating a self-contained learning resource for graduate students, postdoctoral researchers, and bioinformatics instructors.

The tutorial covers the entire scRNA-seq analytical pipeline — from raw 10X Genomics CellRanger output to advanced downstream analyses including trajectory inference (Monocle3) and cell–cell communication network comparison (CellChat).

---

## Major Components

### 1. Original Analysis Code (`code/analysis_code.txt`)
The complete R analysis script, preserved verbatim from the original author. Organised into clearly labelled sections with inline comments explaining the purpose of each block.

### 2. Slide-by-Slide Knowledge Extraction (`docs/slide_by_slide_knowledge.md`)
Structured knowledge extracted from all 109 presentation slides, covering:
- Bulk RNA-seq vs. scRNA-seq conceptual foundations
- Technical challenges (data volume, dropouts, batch effects, sequencing artefacts)
- 10X Genomics platform mechanics (GEM encapsulation, UMI barcoding)
- Full QC, normalisation, integration, and analysis rationale

### 3. Code Explanation (`docs/code_explanation.md`)
Section-by-section walkthrough of 25 code blocks, each with:
- Code purpose and what it does
- Expected inputs and outputs
- Biological or computational interpretation
- Common pitfalls and their solutions

### 4. Tutorial Workflow (`docs/tutorial_workflow.md`)
A sequential 23-step tutorial linking each code step to the corresponding presentation concept, with:
- Expected results at each step
- Pass/fail verification criteria
- Troubleshooting quick reference

### 5. 10-Slide Teaching Deck (`slides/10_slide_teaching_deck.md`)
A condensed presentation guide covering the full analysis arc — suitable for a 90-minute lecture or self-study.

### 6. Environment / Dependencies (`environment/requirements_or_dependencies.md`)
Complete R package dependency list with installation commands, version notes, and system requirements.

---

## Strengths of This Tutorial

**Scientific rigour:**
- Uses a negative binomial test for differential expression — the statistically appropriate model for count data
- BCR/TCR gene removal prevents clonotype-driven clustering artefacts in immune datasets
- Cell-cycle regression via SCTransform covariates is properly implemented
- Two rounds of QC filtering with intermediate visualisation

**Pedagogical design:**
- Every code decision is linked to a conceptual explanation from the presentation
- Diagnostic plots are embedded throughout (not just at the end) to build intuition
- Filtering thresholds are explicitly justified, not presented as magic numbers

**Practical completeness:**
- Covers data loading through advanced CellChat comparison — most tutorials stop at clustering
- Includes memory management strategies for large datasets
- Provides explicit troubleshooting guidance

**Reproducibility:**
- All thresholds and parameters are explicitly stated in the code
- Sample prefixing prevents barcode collision in multi-sample merges
- SCTransform per sample (rather than on merged object) follows current best practices

---

## What a Learner Gains

Upon completing this tutorial, a learner will have:

| Skill | Depth |
|---|---|
| scRNA-seq conceptual foundations | Thorough — from central dogma to 10X mechanics |
| Multi-metric QC decision-making | Applied — with real data and biological justification |
| SCTransform normalisation | Conceptual + practical |
| Seurat anchor-based integration | Conceptual + practical |
| UMAP interpretation | Applied + common misinterpretations addressed |
| Cluster annotation strategy | Marker-gene based, both manual and automated |
| T-cell biology | Subtype markers for 5 distinct T-cell populations |
| Trajectory inference | Pseudotime concept + Monocle3 implementation |
| CellChat L-R analysis | Full comparative pipeline between conditions |

---

## Dataset Context

**GEO accession:** GSE142564  
**Organism:** Mouse (*Mus musculus*)  
**Tissue:** Bone marrow  
**Conditions:** Day 0 and Day 14 time points  
**Disease context:** Likely BCP-ALL (B-cell precursor acute lymphoblastic leukaemia) or related haematopoietic model — inferred from the dataset path (`BCP-ALL`) and the T-cell and B-cell marker gene panel used  
**Sequencing platform:** 10X Genomics (3′ end sequencing)  
**CellRanger version:** *Not specified in source files*

---

## Suggested Citation

If you use this tutorial:

> Farag Ibrahim. *Introduction to Single-cell RNA-seq: Concepts, Applications, and Challenges*. University of Gothenburg, Sweden. GitHub tutorial repository.  
> Dataset: GEO accession GSE142564.

---

## Related Resources

| Resource | URL |
|---|---|
| Seurat documentation | https://satijalab.org/seurat/ |
| Monocle3 documentation | https://cole-trapnell-lab.github.io/monocle3/ |
| CellChat documentation | https://github.com/sqjin/CellChat |
| GEO dataset GSE142564 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142564 |
| 10X Genomics CellRanger | https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger |
| Harvard Chan Bioinformatics scRNA-seq tutorial | https://hbctraining.github.io/scRNA-seq_online/ |
