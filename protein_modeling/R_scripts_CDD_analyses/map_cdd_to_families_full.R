#!/usr/bin/env Rscript
# Map NCBI Batch CD-Search Domain Hits to the nine reference-protein families.
# Uses base R only; preserves every original reference–homolog assignment.
#
# Usage:
# Rscript map_cdd_to_families_full.R \
#   CDD_reference_family_mapping_727.tsv \
#   NCBI_CDD_Domain_Hits_Standard.txt \
#   CDD_results
#
# Optional fifth positional argument: NCBI CDD Features text file.
# Rscript map_cdd_to_families_full.R mapping.tsv domain_hits.txt CDD_results features.txt

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 4L) {
  stop(paste(
    "Usage: Rscript map_cdd_to_families_full.R",
    "<reference_family_mapping.tsv> <domain_hits.txt> [output_prefix] [features.txt]"
  ), call. = FALSE)
}

mapping_path <- args[[1L]]
hits_path <- args[[2L]]
prefix <- if (length(args) >= 3L) args[[3L]] else "CDD_results"
features_path <- if (length(args) >= 4L) args[[4L]] else NULL

for (p in c(mapping_path, hits_path, features_path)) {
  if (!is.null(p) && !file.exists(p)) stop("File not found: ", p)
}

mapping <- read.delim(
  mapping_path, sep = "\t", header = TRUE, quote = "", comment.char = "",
  check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character"
)
required_mapping <- c("query_id", "family", "reference_accession", "role")
if (!all(required_mapping %in% names(mapping))) {
  stop("Mapping file missing columns: ",
       paste(setdiff(required_mapping, names(mapping)), collapse = ", "))
}
if (anyNA(mapping$query_id) || any(!nzchar(mapping$query_id))) {
  stop("Mapping contains blank query_id entries")
}
mapping$.source_row <- seq_len(nrow(mapping))

# NCBI Domain Hits and Features are tab-delimited, with metadata comment lines
# and a column-header line beginning with '#Query'. The latter must not be
# discarded as a comment.
read_ncbi_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header_indices <- grep("^#?Query\\t", lines)
  if (!length(header_indices)) {
    stop("Cannot find '#Query' tab-delimited header in: ", path,
         "\nCheck that this is the NCBI 'Domain Hits' or 'Features' download,",
         " not ASN/XML/JSON alignment details.")
  }
  lines <- lines[header_indices[[1L]]:length(lines)]
  lines[[1L]] <- sub("^#", "", lines[[1L]])
  lines <- lines[nzchar(trimws(lines)) & !grepl("^#", lines)]
  result <- read.delim(
    text = paste(lines, collapse = "\n"), sep = "\t", header = TRUE,
    quote = "", comment.char = "", check.names = FALSE,
    stringsAsFactors = FALSE, colClasses = "character", fill = TRUE
  )
  if (!"Query" %in% names(result)) stop("NCBI table lacks a Query column: ", path)
  result
}

# NCBI Query fields can look like
# 'Q#1 - >WP_003945816.1[WP_003945816.1]' or 'Q#1 - WP_003945816.1'.
# Extract the protein accession, omitting the optional FASTA '>' prefix
# and any bracketed defline appended by NCBI.
query_to_accession <- function(query_string) {
  query_string <- trimws(as.character(query_string))
  pat <- "^Q#[0-9]+[[:space:]]*-[[:space:]]*>?([[:alnum:]_.-]+)"
  match_list <- regmatches(query_string, regexec(pat, query_string, perl = TRUE))
  vapply(match_list, function(m) {
    if (length(m) >= 2L) m[[2L]] else NA_character_
  }, character(1L))
}

map_ncbi_data <- function(ncbi, label) {
  ncbi$query_id <- query_to_accession(ncbi$Query)
  if (anyNA(ncbi$query_id)) {
    stop(label, ": could not parse query accessions from: ",
         paste(head(unique(ncbi$Query[is.na(ncbi$query_id)]), 5L), collapse = "; "))
  }
  unexpected <- setdiff(unique(ncbi$query_id), unique(mapping$query_id))
  if (length(unexpected)) {
    stop(label, ": ", length(unexpected), " identifiers are absent from mapping: ",
         paste(head(unexpected, 10L), collapse = ", "),
         "\nCheck that the same combined FASTA was submitted and headers preserved.")
  }
  ncbi$.hit_row <- seq_len(nrow(ncbi))
  result <- merge(mapping, ncbi, by = "query_id", all.x = TRUE, sort = FALSE)
  result <- result[order(result$.source_row, result$.hit_row, na.last = TRUE), , drop = FALSE]
  rownames(result) <- NULL
  result$.source_row <- NULL
  result$.hit_row <- NULL
  result
}

write_tab <- function(x, path) {
  write.table(x, path, sep = "\t", quote = TRUE,
              row.names = FALSE, col.names = TRUE, na = "")
  message("Saved ", path, " (", nrow(x), " rows)")
}

hits <- read_ncbi_table(hits_path)
if (!"Accession" %in% names(hits)) stop("Domain Hits file lacks an Accession column")
mapped <- map_ncbi_data(hits, "Domain Hits")
mapped$CDD_status <- ifelse(is.na(mapped$Accession) | !nzchar(mapped$Accession),
                            "No reported domain hit", "Domain hit reported")
write_tab(mapped, paste0(prefix, "_mapped_domain_hits.tsv"))

annotated_ids <- unique(mapped$query_id[mapped$CDD_status == "Domain hit reported"])
missing <- mapping[!(mapping$query_id %in% annotated_ids),
                   setdiff(names(mapping), ".source_row"), drop = FALSE]
write_tab(missing, paste0(prefix, "_no_reported_hits.tsv"))

family_names <- unique(mapping$family)
family_summary <- do.call(rbind, lapply(family_names, function(fam) {
  mf <- mapping[mapping$family == fam, , drop = FALSE]
  ids <- unique(mf$query_id)
  data.frame(
    family = fam,
    reference_accession = mf$reference_accession[[1L]],
    sequence_assignments = nrow(mf),
    distinct_accessions_in_family = length(ids),
    with_reported_domain_hit = sum(ids %in% annotated_ids),
    without_reported_domain_hit = sum(!(ids %in% annotated_ids)),
    stringsAsFactors = FALSE
  )
}))
write_tab(family_summary, paste0(prefix, "_family_summary.tsv"))

if (!is.null(features_path)) {
  feats <- read_ncbi_table(features_path)
  mapped_feats <- map_ncbi_data(feats, "Features")
  write_tab(mapped_feats, paste0(prefix, "_mapped_features.tsv"))
}

message("Input mapping rows: ", nrow(mapping),
        "; unique query accessions: ", length(unique(mapping$query_id)),
        "; NCBI domain-hit rows: ", nrow(hits),
        "; unique accessions with hits: ", length(annotated_ids),
        "; unique accessions without hits: ",
        length(setdiff(unique(mapping$query_id), annotated_ids)))
message("Do not count overlapping/superfamily domain matches as distinct enzyme activities.")
