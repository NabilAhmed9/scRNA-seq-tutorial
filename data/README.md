# Data Directory

This directory is intentionally empty in the repository. Raw scRNA-seq data files are not included due to file size constraints.

## How to Obtain the Data

The analysis uses publicly available data from NCBI GEO:

**Accession:** GSE142564  
**URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142564

### Download Instructions

1. Navigate to the GEO accession page
2. Download the supplementary files: `GSE142564_RAW.tar`
3. Extract the archive:
   ```bash
   tar -xvf GSE142564_RAW.tar
   ```
4. You should find directories for each sample containing the three CellRanger output files:
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

### Update File Paths

Before running the analysis, update the file paths in `code/analysis_code.txt`:

```r
# Change this line:
setwd("/Volumes/BCP-ALL/presentations/single cell/GSE142564_RAW/Day0")
data <- Read10X("/Volumes/BCP-ALL/presentations/single cell/GSE142564_RAW/Day0/")

# To your local path, e.g.:
setwd("/your/local/path/GSE142564_RAW/Day0")
data <- Read10X("/your/local/path/GSE142564_RAW/Day0/")
```

## Expected Data Size

| File | Approximate Size |
|---|---|
| Raw FASTQ per sample | ~30 GB |
| CellRanger matrix files per sample | ~100–500 MB |
| Total for Day0 + Day14 | ~1 GB |

## Organism

*Mus musculus* (mouse) — note that mitochondrial gene prefix is `mt-` (lowercase), not `MT-` (human convention).
