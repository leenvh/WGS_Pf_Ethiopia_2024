#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeqArray)
  library(optparse)
})

# ----------------------------
# Parse command-line arguments
# ----------------------------
option_list <- list(
  make_option(
    c("-d", "--workdir"),
    type = "character",
    default = ".",
    help = "Working directory containing the VCF file"
  ),
  make_option(
    c("-f", "--file"),
    type = "character",
    help = "Merged VCF file name (e.g. merged.genotyped.vcf.gz)"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ----------------------------
# Input validation
# ----------------------------
if (is.null(opt$file)) {
  stop("ERROR: Please provide a VCF file using -f")
}

vcf_file <- file.path(opt$workdir, opt$file)
gds_file <- file.path(
  opt$workdir,
  sub("\\.vcf\\.gz$", ".gds", opt$file)
)

if (!file.exists(vcf_file)) {
  stop("ERROR: VCF file does not exist: ", vcf_file)
}

# ----------------------------
# Convert VCF to GDS
# ----------------------------
if (file.exists(gds_file)) {
  message("GDS file already exists: ", gds_file)
} else {
  message("Converting VCF to GDS...")
  seqVCF2GDS(
    vcf.fn = vcf_file,
    out.fn = gds_file,
    verbose = TRUE
  )
  message("Conversion complete: ", gds_file)
}

# ----------------------------
# Quick sanity check (AD field)
# ----------------------------
message("Checking allele depth (AD) field...")
summary_ad <- seqSummary(gds_file, "annotation/format/AD")
print(summary_ad)

message("Step 1 finished successfully.")
