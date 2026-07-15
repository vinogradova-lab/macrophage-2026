# Export the MSigDB C5 GO sets (the SAME data .go_msigdb() builds for go_enrich())
# to a flat CSV the Python query_gene_ontology_msigdb() reads. Guarantees the curated
# category gene lists use the exact 2026.1.Hs membership as the R Fisher enrichment.
#
# Output: reference_dbs/msigdb_go_2026.csv  (ontology, GO_id, term, description, gene)
suppressWarnings(suppressMessages(library(msigdbr)))
suppressWarnings(suppressMessages(library(dplyr)))

repo <- "/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/Vinogradova Laboratory/Henry_data processing/01_data_analysis_folders/06_macrophage/macrophage"
out_dir <- file.path(repo, "reference_dbs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(packageVersion("msigdbr") >= "26.1.0")

db_version <- NULL
one_ont <- function(ont) {
  m <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C5",
                        subcollection = paste0("GO:", ont))
  if (is.null(db_version)) db_version <<- unique(stats::na.omit(m$db_version))
  m %>%
    dplyr::transmute(
      ontology    = ont,
      # term name cleaned the same way .go_msigdb() builds t2n$name
      GO_id       = gs_exact_source,
      term        = tolower(gsub("_", " ", sub(paste0("^GO", ont, "_"), "", gs_name))),
      # gs_description = GO definition text; lets free-text queries search defs, not just names.
      description = tolower(gs_description),
      gene        = gene_symbol
    ) %>%
    dplyr::distinct()
}

all <- dplyr::bind_rows(lapply(c("BP", "CC", "MF"), one_ont))

out <- file.path(out_dir, "msigdb_go_2026.csv")
write.csv(all, out, row.names = FALSE)

cat("wrote", out, "\n")
cat("db_version:", paste(db_version, collapse = ", "), "\n")
cat("rows:", nrow(all), " GO terms:", length(unique(all$GO_id)),
    " genes:", length(unique(all$gene)), "\n")
print(table(all$ontology))
