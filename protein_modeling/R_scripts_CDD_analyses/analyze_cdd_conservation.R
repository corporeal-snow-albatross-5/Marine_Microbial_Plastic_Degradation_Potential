#!/usr/bin/env Rscript
# Reproducible conserved-domain and functional-feature summary for
# NCBI Batch CD-Search dataset (nine reference-defined protein groups).
# Uses only base R. Primary annotations are user-reviewable, not inferred enzyme activity.
#
# Usage (Mac Terminal, from the folder with your CDD_results TSV files):
#   Rscript analyze_cdd_conservation_full.R \
#     CDD_results_mapped_domain_hits.tsv \
#     CDD_results_mapped_features.tsv \
#     CDD_review
#
# Outputs: marker_config.tsv; marker_by_protein.tsv; marker_summary.tsv;
# primary_domain_summary.tsv; feature_summary.tsv;
# gpx_catalytic_residue_detail.tsv; gpx_catalytic_residue_summary.tsv;
# primary_domain_conservation.pdf and .png.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript analyze_cdd_conservation_full.R <mapped_domain_hits.tsv> <mapped_features.tsv> <output_prefix>", call. = FALSE)
}

hit_file <- args[[1]]
feature_file <- args[[2]]
prefix <- args[[3]]
if (!file.exists(hit_file) || !file.exists(feature_file)) stop("Input TSV missing", call. = FALSE)

hits <- read.delim(hit_file, sep="\t", header=TRUE, quote="\"", comment.char="",
                   check.names=FALSE, stringsAsFactors=FALSE)
features <- read.delim(feature_file, sep="\t", header=TRUE, quote="\"", comment.char="",
                       check.names=FALSE, stringsAsFactors=FALSE)
required_hits <- c("family", "reference_accession", "query_id", "role", "Accession", "E-Value", "From", "To")
required_feats <- c("family", "query_id", "role", "Title", "coordinates", "complete size", "mapped size")
if (!all(required_hits %in% names(hits))) stop("Missing required columns in mapped hits: ", paste(setdiff(required_hits, names(hits)), collapse=", "))
if (!all(required_feats %in% names(features))) stop("Missing required columns in features: ", paste(setdiff(required_feats, names(features)), collapse=", "))
# Full hits should contain COG3509 in reference and all MCL homologs.
# Detect accidental use of the Standard mapped table, which suppresses these overlaps.
mcl_reference <- hits[hits$family == "MCL_PHA" & hits$role == "reference" &
                        hits$Accession == "COG3509", , drop=FALSE]
if (nrow(mcl_reference) == 0L) {
  stop("MCL_PHA reference COG3509 is absent: you likely supplied the STANDARD mapped hits. Re-run mapping with 532_unique_homologs_hitdata-2.txt (FULL).", call.=FALSE)
}

# These are EXACT CDD accessions present in the respective REFERENCE protein,
# Both COG4099 and COG3509 occur on the MCL_PHA reference and all 11
# homologs in the Full report. All those matches are NON-SPECIFIC and
# overlap, so they are alternative annotations, not independent activities.
# 'primary' is a selected reference-associated marker, NOT proof of enzymatic activity.
markers <- data.frame(
  family=c(
    "GPx_WP_003945816.1", "GPx_WP_003945816.1",
    "GPx_WP_026826575.1", "GPx_WP_026826575.1",
    "Laccase", "Laccase", "Oxidoreductase", "Oxidoreductase",
    "PETase_UUT36763.1", "PETase_UUT36763.1",
    "PETase_UUT36764.1", "PETase_UUT36764.1",
    "MCL_PHA", "MCL_PHA",
    "PHB", "PHB", "PHB", "P3HB", "P3HB", "P3HB"
  ),
  accession=c(
    "cd00340", "COG0386", "cd00340", "COG0386",
    "pfam07732", "TIGR01480", "pfam01565", "COG0277",
    "pfam12740", "COG4188", "pfam12740", "COG4188",
    "COG4099", "COG3509",
    "TIGR01840", "pfam10503", "pfam02412",
    "TIGR01840", "pfam10503", "pfam02412"
  ),
  marker=c(
    "GSH_Peroxidase", "BtuE", "GSH_Peroxidase", "BtuE",
    "Cu-oxidase_3", "copper_res_A", "FAD_binding_4", "GlcD",
    "PETase", "COG4188", "PETase", "COG4188",
    "COG4099 (predicted peptidase)", "LpqC (esterase-related)",
    "esterase_phb", "Esterase_PHB", "TSP_3 (accessory repeat)",
    "esterase_phb", "Esterase_PHB", "TSP_3 (accessory repeat)"
  ),
  marker_role=c(
    "primary", "secondary", "primary", "secondary",
    "primary", "secondary", "primary", "secondary",
    "primary", "secondary", "primary", "secondary",
    "primary", "secondary",
    "primary", "secondary", "accessory",
    "primary", "secondary", "accessory"
  ),
  stringsAsFactors=FALSE
)
if (!identical(sort(unique(hits$family)), sort(unique(markers$family)))) {
  stop("Dataset families differ from marker configuration. Review markers before analysis.", call. = FALSE)
}
write.table(markers, paste0(prefix, "_marker_config.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="")

hit_ids <- unique(hits[c("family", "reference_accession", "query_id", "role")])
if (any(duplicated(hit_ids[c("family", "query_id")]))) {
  stop("Inconsistent roles/reference metadata for the same family + query_id", call.=FALSE)
}
if (!all(hit_ids$role %in% c("reference", "marine_BLASTP_homolog"))) stop("Unexpected role values", call.=FALSE)

# One row per family x protein x annotation marker, regardless of overlapping
# CDD matches or repeated domains on one protein.
by_protein <- vector("list", nrow(markers))
by_summary <- vector("list", nrow(markers))
for (i in seq_len(nrow(markers))) {
  m <- markers[i, , drop=FALSE]
  ids <- hit_ids[hit_ids$family == m$family, , drop=FALSE]
  v <- hits[hits$family == m$family & hits$Accession == m$accession, , drop=FALSE]
  present <- ids$query_id %in% v$query_id
  best_e <- vapply(ids$query_id, function(id) {
    z <- suppressWarnings(as.numeric(v$`E-Value`[v$query_id == id]))
    z <- z[!is.na(z)]
    if (length(z) == 0L) return(NA_real_)
    min(z)
  }, numeric(1L))
  p <- data.frame(
    family=ids$family, reference_accession=ids$reference_accession,
    query_id=ids$query_id, role=ids$role,
    marker_accession=m$accession, marker=m$marker, marker_role=m$marker_role,
    has_marker=present, best_evalue=best_e,
    stringsAsFactors=FALSE
  )
  by_protein[[i]] <- p
  ref <- p$role == "reference"
  hom <- p$role == "marine_BLASTP_homolog"
  if (sum(ref) != 1L) stop("Expected exactly one reference per group: ", m$family)
  n <- sum(hom)
  k <- sum(p$has_marker[hom])
  by_summary[[i]] <- data.frame(
    family=m$family, reference_accession=ids$reference_accession[which(ref)[1]],
    marker_accession=m$accession, marker=m$marker, marker_role=m$marker_role,
    marker_found_in_reference=p$has_marker[ref],
    n_homologs=n, n_homologs_with_marker=k,
    percent_homologs_with_marker=if (n > 0) 100*k/n else NA_real_,
    stringsAsFactors=FALSE
  )
}
by_protein <- do.call(rbind, by_protein)
by_summary <- do.call(rbind, by_summary)
rownames(by_protein) <- NULL
rownames(by_summary) <- NULL
write.table(by_protein, paste0(prefix, "_marker_by_protein.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="")
write.table(by_summary, paste0(prefix, "_marker_summary.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="")

primary <- by_summary[by_summary$marker_role == "primary", , drop=FALSE]
if (!all(primary$marker_found_in_reference)) stop("A primary marker was NOT detected on its reference; revise marker config", call.=FALSE)
write.table(primary, paste0(prefix, "_primary_domain_summary.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="")


# Audit CDD-specific vs non-specific match types for every selected primary marker.
# A Full-mode overlap is not independent evidence of biochemical activity.
if ("Hit type" %in% names(hits)) {
  primary_types <- lapply(seq_len(nrow(primary)), function(i) {
    s <- primary[i,]
    h <- hits[hits$family == s$family & hits$role == "marine_BLASTP_homolog" &
              hits$Accession == s$marker_accession,,drop=FALSE]
    n_unique_specific <- length(unique(h$query_id[h$`Hit type` == "specific"]))
    n_unique_nonspecific <- length(unique(h$query_id[h$`Hit type` == "non-specific"]))
    data.frame(family=s$family, marker_accession=s$marker_accession,
      n_homologs_with_marker=s$n_homologs_with_marker,
      n_with_specific_match=n_unique_specific, n_with_non_specific_match=n_unique_nonspecific,
      n_homologs=s$n_homologs, stringsAsFactors=FALSE)
  })
  write.table(do.call(rbind, primary_types), paste0(prefix,"_match_specificity.tsv"),
              sep="\t",row.names=FALSE,quote=FALSE,na="")
}

# Catalog evidence of CDD-reported features by family. No feature annotation
# should be interpreted as absence of a motif when it was not reported.
valid_features <- features[!is.na(features$Title) & nzchar(features$Title) &
                           features$role == "marine_BLASTP_homolog", , drop=FALSE]
feature_keys <- unique(valid_features[c("family", "Title")])
feature_summary <- do.call(rbind, lapply(seq_len(nrow(feature_keys)), function(i) {
  a <- feature_keys[i, , drop=FALSE]
  x <- valid_features[valid_features$family == a$family & valid_features$Title == a$Title, , drop=FALSE]
  denom <- length(unique(hit_ids$query_id[hit_ids$family == a$family & hit_ids$role == "marine_BLASTP_homolog"]))
  data.frame(family=a$family, feature=a$Title,
             n_homologs_with_feature=length(unique(x$query_id)),
             n_homologs=denom,
             percent=100*length(unique(x$query_id))/denom,
             stringsAsFactors=FALSE)
}))
if (is.null(feature_summary)) {
  feature_summary <- data.frame(family=character(), feature=character(), n_homologs_with_feature=integer(), n_homologs=integer(), percent=numeric())
}
feature_summary <- feature_summary[order(feature_summary$family, feature_summary$feature),,drop=FALSE]
write.table(feature_summary, paste0(prefix, "_feature_summary.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="")

# For GPx only: compare ordered amino-acid identities of reference's annotated
# catalytic features with homolog's own CDD-mapped catalytic features. The
# residue NUMBERS may shift, so compare residue identities rather than numbers.
letters_of <- function(coordinates) {
  if (length(coordinates) == 0L || is.na(coordinates) || !nzchar(coordinates)) return(NA_character_)
  pieces <- strsplit(coordinates, ",", fixed=TRUE)[[1]]
  paste(gsub("[0-9[:space:]]", "", pieces), collapse="")
}
gpx_ids <- unique(hit_ids$family[grepl("^GPx_", hit_ids$family)])
residue_details <- list()
residue_summaries <- list()
for (fam in gpx_ids) {
  ids <- hit_ids[hit_ids$family == fam,,drop=FALSE]
  fc <- features[features$family == fam & !is.na(features$Title) &
                   features$Title == "catalytic residues",,drop=FALSE]
  ref <- fc[fc$role == "reference",,drop=FALSE]
  if (nrow(ref) != 1L) stop("Expected one annotated GPx reference catalytic feature: ",fam)
  ref_pattern <- letters_of(ref$coordinates[1])
  ref_n <- suppressWarnings(as.numeric(ref$`complete size`[1]))
  hom_ids <- ids$query_id[ids$role == "marine_BLASTP_homolog"]
  rows <- lapply(hom_ids, function(q) {
    rr <- fc[fc$query_id == q & fc$role == "marine_BLASTP_homolog",,drop=FALSE]
    match <- any(vapply(seq_len(nrow(rr)), function(j) {
      observed <- letters_of(rr$coordinates[j])
      valid <- !is.na(observed) && !is.na(rr$`mapped size`[j]) && !is.na(rr$`complete size`[j])
      valid && identical(observed, ref_pattern) &&
        as.numeric(rr$`mapped size`[j]) == as.numeric(rr$`complete size`[j]) &&
        as.numeric(rr$`complete size`[j]) == ref_n
    }, logical(1L)))
    data.frame(family=fam, query_id=q, reference_residue_identities=ref_pattern,
               homolog_annotated_coordinates=if (nrow(rr)) paste(rr$coordinates, collapse=";") else NA_character_,
               has_annotated_catalytic_feature=(nrow(rr)>0L),
               all_reference_residue_identities_match=match, stringsAsFactors=FALSE)
  })
  d <- do.call(rbind, rows)
  residue_details[[fam]] <- d
  residue_summaries[[fam]] <- data.frame(
    family=fam, reference_residue_identities=ref_pattern,
    n_homologs=nrow(d), n_with_cdd_catalytic_feature=sum(d$has_annotated_catalytic_feature),
    n_with_matching_residue_identities=sum(d$all_reference_residue_identities_match),
    percent_with_matching_residue_identities=100*mean(d$all_reference_residue_identities_match),
    stringsAsFactors=FALSE)
}
rd <- do.call(rbind, residue_details)
rs <- do.call(rbind, residue_summaries)
write.table(rd, paste0(prefix, "_gpx_catalytic_residue_detail.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="")
write.table(rs, paste0(prefix, "_gpx_catalytic_residue_summary.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="")

# Manuscript-ready figure showing EXACT reference-associated domain presence,
# with references excluded from denominators. Not evidence of catalytic function.
plot_primary <- function() {
  nm <- c("GPx (WP_003945816.1)", "GPx (WP_026826575.1)",
          "Laccase-associated", "Oxidoreductase", "PETase (UUT36763.1)",
          "PETase (UUT36764.1)", "MCL-PHA", "PHB", "P3HB")
  # The row order follows markers; make the order explicit for robustness.
  order_fam <- c("GPx_WP_003945816.1", "GPx_WP_026826575.1", "Laccase",
                 "Oxidoreductase", "PETase_UUT36763.1", "PETase_UUT36764.1",
                 "MCL_PHA", "PHB", "P3HB")
  s <- primary[match(order_fam, primary$family),,drop=FALSE]
  values <- rev(s$percent_homologs_with_marker)
  nlabel <- rev(sprintf("%d/%d", s$n_homologs_with_marker, s$n_homologs))
  fam <- rev(s$family)
  labs <- rev(nm)
  par(mar=c(5.4, 13, 4.2, 4), xpd=NA)
  mid <- barplot(values, names.arg=labs, horiz=TRUE, las=1, xlim=c(0,115),
                 col=ifelse(fam=="MCL_PHA", "#B35B3C", "#285474"),
                 border=NA, cex.names=0.9, cex.axis=0.95, xlab="Homologs with reference-associated CDD model (%)",
                 main="CDD domain-model retention across marine homolog groups")
  abline(v=c(25,50,75,100), col="#DFE5E9", lty="dotted")
  text(values+2, mid, labels=nlabel, adj=0, cex=0.85)
  mtext("Models selected from Full CDD hits; no biochemical activity inferred", side=1, line=3.4, cex=0.72)
}
pdf(paste0(prefix,"_primary_domain_conservation.pdf"), width=11, height=6.5, useDingbats=FALSE)
plot_primary()
dev.off()
png(paste0(prefix,"_primary_domain_conservation.png"), width=4400, height=2600, res=400)
plot_primary()
dev.off()

cat("Families analyzed:",nrow(primary),"\n")
print(primary[c("family", "marker_accession", "n_homologs_with_marker", "n_homologs", "percent_homologs_with_marker")], row.names=FALSE)
cat("\nGPx catalytic residue annotations:\n")
print(rs, row.names=FALSE)
cat("\nImportant: CDD family models and annotated residues are NOT biochemical validation.\n")
cat("Outputs written with prefix:",prefix,"\n")
