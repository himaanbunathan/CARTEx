## CARTEx

**CARTEx** is an R package that calculates the CARTEx score based on precomputed gene weights from a CAR-T model of exhaustion. It provides tools for analyzing gene expression data and calculating phenotypic scores based on a PCA method — for **bulk** RNA-seq and for **single-cell** objects.

### Installation from GitHub

To install the package from GitHub, use the `devtools` package:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Single-cell HDF5 readers (.h5ad / .h5Seurat / .loom) need rhdf5:
install.packages("BiocManager")
BiocManager::install("rhdf5")

# Install CARTEx from GitHub (build_vignettes = TRUE makes vignette("cartex-single-cell") available)
devtools::install_github("himaanbunathan/CARTEx", build_vignettes = TRUE)
```

### Bulk RNA-seq

```r
library(CARTEx)
# gene_counts: normalized data.frame, genes as rows, samples as columns
scores <- cartex_score(gene_counts)
```

### Single-cell data

`cartex_sc()` scores every cell with the 200-gene CARTEx exhaustion signature.
It reads `.rds`, `.RData`, `.h5ad`, `.h5Seurat`, `.loom`, or an in-memory Seurat
object / matrix, and is built to run on a laptop — for HDF5 files it streams
only the CARTEx genes off disk instead of loading the whole matrix.

```r
library(CARTEx)

obj <- cartex_sc("pbmc.h5ad")          # or a Seurat object, .rds, .loom, ...

obj$cartex_score        # per-cell, Z-scored across cells (headline value)
obj$cartex_score_raw    # raw weighted dot-product

print(attr(obj, "cartex_report"))      # full transparency report
```

**CARTEx is only valid for CD8+ T cells**, so every cell is classified first.
By default all cells are scored and labelled (`cd8_status`); restrict to confident
CD8+ T cells with `cd8 = "restrict"`:

```r
obj <- cartex_sc("pbmc.h5ad")                    # score all, add cd8_status
cd8 <- cartex_sc("pbmc.h5ad", cd8 = "restrict")  # score only CD8+ T cells
```

Classification (`cd8_method`) uses a cell-type annotation column when available,
falling back to marker-gene gating (CD8A/CD8B+ and CD3+ and CD4−; CD3 is required
because NK cells also express CD8A). Cells it cannot resolve are labelled
`ambiguous_*` rather than guessed.

**Pseudobulk** — when cluster/cell-type info is available, score aggregated
profiles instead of single cells (sum counts → normalize → score, per the
reference CARTEx pipeline):

```r
pb <- cartex_sc("pbmc.h5ad", level = "pseudobulk", group_by = "seurat_clusters")
```

**Visualize & report** — a side-by-side UMAP (annotation vs CARTEx score) and a
self-contained HTML report:

```r
res <- cartex_sc("pbmc.h5ad")
cartex_umap(res)                                  # side-by-side UMAP
cartex_html_report(res, file = "report.html")     # HTML summary + UMAP
```

`cartex_umap()` reuses an existing UMAP (Seurat reduction, AnnData
`obsm/X_umap`, or `UMAP_1`/`UMAP_2` columns), and for a Seurat object without one
computes a UMAP with default Seurat parameters.

Every score inherits the preprocessing baked into your object, so `cartex_sc()`
prints and attaches a **transparency report** documenting:

- **Gene matching** — how many of the 200 genes matched and by which rule
  (exact symbol, cleaned symbol, or Ensembl ID), plus the unmatched list.
- **CD8 classification** — method used and the CD8+/CD8−/ambiguous composition.
- **Normalization** — value-based detection of raw counts vs log-normalized data
  (the declared label is often wrong), applying `LogNormalize` only when needed.
- **QC metadata** — detected `nFeature`/`nCount`/`percent.mt` columns and their
  ranges, with a flag for whether QC filtering looks already applied.

> The headline `cartex_score` is Z-scored across the units in one object, so it
> is comparable **within** a scored dataset, not across separately scored ones.
> Use `cartex_score_raw` for an absolute value; `cartex_score_int` is the
> rounded/capped (±4) variant from the reference pipeline.

See `vignette("cartex-single-cell")` for details.
