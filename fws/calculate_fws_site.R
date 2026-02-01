#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeqArray)
  library(moimix)
  library(optparse)
})

# ----------------------------
# Command-line arguments
# ----------------------------
option_list <- list(
  make_option(
    c("-d", "--workdir"),
    type = "character",
    default = ".",
    help = "Working directory"
  ),
  make_option(
    c("-g", "--gds"),
    type = "character",
    help = "Merged GDS file (e.g. merged.genotyped.gds)"
  ),
  make_option(
    c("-s", "--sites"),
    type = "character",
    help = "Sample-to-site mapping file (TSV)"
  ),
  make_option(
    c("-o", "--out"),
    type = "character",
    default = "fws_by_site.tsv",
    help = "Output file name"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ----------------------------
# Input checks
# ----------------------------
if (is.null(opt$gds) || is.null(opt$sites)) {
  stop("ERROR: Please provide both --gds and --sites")
}

gds_file   <- file.path(opt$workdir, opt$gds)
sites_file <- file.path(opt$workdir, opt$sites)
out_file   <- file.path(opt$workdir, opt$out)

if (!file.exists(gds_file)) {
  stop("ERROR: GDS file not found: ", gds_file)
}
if (!file.exists(sites_file)) {
  stop("ERROR: Site file not found: ", sites_file)
}

# ----------------------------
# Read site information
# ----------------------------
sites <- read.table(
  sites_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

required_cols <- c("sample", "site")
if (!all(required_cols %in% colnames(sites))) {
  stop("ERROR: Site file must contain columns: sample, site")
}

# ----------------------------
# Open GDS
# ----------------------------
gds <- seqOpen(gds_file)

all_samples <- seqGetData(gds, "sample.id")

# Keep only samples present in GDS
sites <- sites[sites$sample %in% all_samples, ]

if (nrow(sites) == 0) {
  stop("ERROR: No samples from site file found in GDS")
}

# ----------------------------
# Compute Fws per site
# ----------------------------
results <- list()

for (site_name in unique(sites$site)) {

  site_samples <- sites$sample[sites$site == site_name]

  if (length(site_samples) < 2) {
    message("Skipping site ", site_name, " (fewer than 2 samples)")
    next
  }

  message("Computing Fws for site: ", site_name)

  seqSetFilter(
    gds,
    sample.id = site_samples,
    verbose = FALSE
  )

  fws <- getFws(gds)

  site_df <- data.frame(
    sample = names(fws),
    site   = site_name,
    fws    = round(as.numeric(fws), 3),
    stringsAsFactors = FALSE
  )

  results[[site_name]] <- site_df

  seqResetFilter(gds)
}

seqClose(gds)

# ----------------------------
# Write output
# ----------------------------
final_df <- do.call(rbind, results)

write.table(
  final_df,
  out_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Fws calculation finished.")
message("Results written to: ", out_file)
