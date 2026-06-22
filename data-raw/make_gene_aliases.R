# data-raw/make_gene_aliases.R
#
# Builds inst/extdata/cartex_gene_map.rds: a transparent lookup that lets the
# CARTEx gene weights (named with HGNC-ish symbols, some carrying make.names
# artefacts like "IQCJ.SCHIP1") be matched against single-cell objects that may
# label genes with current symbols OR Ensembl gene IDs.
#
# Provenance of the Ensembl IDs: extracted from the var annotation of the
# CELLxGENE reference object used during development
# (Test_dataset/267e4e57-...h5ad, GENCODE annotation), supplemented with a small
# set of hand-curated symbol renames for genes whose make.names/legacy symbol no
# longer matches current references. Anyone can re-run this script to regenerate
# the table; nothing here is magic.
#
# Run from the package root:  Rscript data-raw/make_gene_aliases.R

suppressMessages({
  library(rhdf5)
})

weights_path <- "inst/extdata/weights.rds"
cartex <- readRDS(weights_path)
cartex_symbols <- names(cartex)

# --- Curated symbol renames (legacy / make.names -> current HGNC symbol) -------
# Keyed by the symbol as it appears in weights.rds.
curated_clean <- c(
  "H1F0"        = "H1-0",          # legacy symbol -> current
  "C4orf26"     = "ODAPH",         # renamed gene
  "IQCJ.SCHIP1" = "IQCJ-SCHIP1"    # make.names artefact ("." should be "-")
)

# --- Curated Ensembl IDs for the genes that do not match by symbol ------------
curated_ensembl <- c(
  "H1F0"        = "ENSG00000189060",
  "C4orf26"     = "ENSG00000174792",
  "IQCJ.SCHIP1" = "ENSG00000273611",
  "CLECL1"      = "ENSG00000184293",
  "MICALCL"     = "ENSG00000133816",
  "KIR3DX1"     = "ENSG00000275342"
)

# --- Harvest symbol -> Ensembl from the reference h5ad (if available) ----------
ref_h5ad <- Sys.getenv(
  "CARTEX_REF_H5AD",
  unset = "../Test_dataset/267e4e57-d8e0-47d7-bfed-aaee0402497d.h5ad"
)

sym2ens <- character(0)
if (file.exists(ref_h5ad)) {
  cats  <- as.character(h5read(ref_h5ad, "var/feature_name/categories"))
  codes <- as.integer(h5read(ref_h5ad, "var/feature_name/codes"))
  ens   <- as.character(h5read(ref_h5ad, "var/_index"))
  ref_syms <- cats[codes + 1L]
  keep <- !duplicated(ref_syms)
  sym2ens <- setNames(ens[keep], ref_syms[keep])
  h5closeAll()
  message("Harvested ", length(sym2ens), " symbol->Ensembl pairs from reference.")
} else {
  message("Reference h5ad not found at '", ref_h5ad,
          "'; relying on curated Ensembl IDs only.")
}

# --- Assemble the map ---------------------------------------------------------
clean_symbol <- gsub("\\.", "-", cartex_symbols)          # default: dot -> dash
ow <- intersect(names(curated_clean), cartex_symbols)
clean_symbol[match(ow, cartex_symbols)] <- curated_clean[ow]

resolve_ensembl <- function(orig, clean) {
  if (orig %in% names(curated_ensembl)) return(curated_ensembl[[orig]])
  for (key in unique(c(orig, clean))) {
    if (!is.na(sym2ens[key]) && nzchar(sym2ens[key])) return(unname(sym2ens[key]))
  }
  NA_character_
}

ensembl_id <- mapply(resolve_ensembl, cartex_symbols, clean_symbol,
                     USE.NAMES = FALSE)

gene_map <- data.frame(
  cartex_symbol = cartex_symbols,
  clean_symbol  = clean_symbol,
  ensembl_id    = ensembl_id,
  weight        = as.numeric(cartex),
  stringsAsFactors = FALSE
)

out <- "inst/extdata/cartex_gene_map.rds"
saveRDS(gene_map, out)
message("Wrote ", out, " with ", nrow(gene_map), " genes; ",
        sum(!is.na(gene_map$ensembl_id)), " have Ensembl IDs.")
print(utils::head(gene_map))
