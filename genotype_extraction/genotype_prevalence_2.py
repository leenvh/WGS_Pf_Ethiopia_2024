import zarr
import allel
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--zarr")
parser.add_argument("--position_df")

args = parser.parse_args()

# ---------------------------
# User parameters
# ---------------------------

# Positions file (CSV or TSV)
# Must contain columns: POS
# Optional column: CHROM
#positions_file = "positions.csv"

# Default chromosome if CHROM column is missing
#target_chrom = "Pf3D7_01_v3"

# Annotation dictionary
# key = (chrom, pos, alt)
# value = (GENE, AA_change, TYPE)
ann = {}  # populate if available

# Zarr path
zarr_path = (
   args.zarr
    )

# ---------------------------
# Load positions
# ---------------------------
print("Loading target positions...")
#positions_df = pd.read_csv("/mnt/storage13/ahri/popgen_ethiopia_WGS/fst/fst_scikit_allel/hrp_status/genotype_prevalence/gametocyte_relate_high_fst.csv")
positions_df = pd.read_csv(args.position_df,sep="\t") # px1 position

if "POS" not in positions_df.columns:
    raise ValueError("Positions file must contain a POS column")

# If CHROM missing, assume all on target_chrom
if "CHROM" not in positions_df.columns:
    positions_df["CHROM"] = target_chrom

# Convert to set for fast lookup
target_sites = set(zip(positions_df["CHROM"], positions_df["POS"]))

# ---------------------------
# Load Zarr callset
# ---------------------------
print("Loading Zarr callset...")

callset = zarr.open(zarr_path, mode="r")

genotypes = allel.GenotypeArray(callset["calldata/GT"])
samples = callset["samples"][:]

chrom = callset["variants/CHROM"][:]
pos = callset["variants/POS"][:]
ref = callset["variants/REF"][:]

# Handle ALT safely (first ALT only)
alt_raw = callset["variants/ALT"][:]
alt = np.array([a[0] if len(a) > 0 else None for a in alt_raw])

# ---------------------------
# Identify matching variants
# ---------------------------
print("Matching target variants...")

mask = np.fromiter(
    ((chrom[i], pos[i]) in target_sites for i in range(len(pos))),
    dtype=bool,
    count=len(pos)
)

idxs = np.where(mask)[0]

if len(idxs) == 0:
    raise ValueError("❌ No matching variants found in the Zarr file")

print(f"Found {len(idxs)} matching variants")

# ---------------------------
# Extract genotypes
# ---------------------------
print("Extracting genotypes...")

gt = genotypes[idxs]
gt_strings = gt.to_gt()

colnames = [
    f"{chrom[i]}:{pos[i]}_{ref[i]}>{alt[i]}"
    for i in idxs
]

geno_df = pd.DataFrame(
    gt_strings.T,
    index=samples,
    columns=colnames
)

geno_df.insert(0, "sample", samples)

# ---------------------------
# Allele frequencies
# ---------------------------
print("Calculating allele frequencies...")

ac = gt.count_alleles()
af = ac.to_frequencies()

alt_af = af[:, 1] if af.shape[1] > 1 else np.zeros(len(af))

# ---------------------------
# Annotations (only for matched variants)
# ---------------------------
genes = []
aa_changes = []
types = []

for i in idxs:
    key = (chrom[i], pos[i], alt[i])
    gene, aa, vtype = ann.get(key, (None, None, None))
    genes.append(gene)
    aa_changes.append(aa)
    types.append(vtype)

# ---------------------------
# Frequency table
# ---------------------------
freq_df = pd.DataFrame({
    "CHROM": chrom[idxs],
    "POS": pos[idxs],
    "REF": ref[idxs],
    "ALT": alt[idxs],
    "REF_AF": af[:, 0],
    "ALT_AF": alt_af,
    "GENE": genes,
    "AA_change": aa_changes,
    "TYPE": types
})

# ---------------------------
# Save outputs
# ---------------------------
geno_out = "target_positions_genotypes_k13_675.csv"
freq_out = "target_positions_frequencies_k13_675.csv"
geno_df.to_csv("target_wide_format_k13_675.csv")

geno_long = geno_df.melt(
    id_vars="sample",
    var_name="variant",
    value_name="genotype"
)
# Split CHROM and rest
geno_long[["CHROM", "rest"]] = geno_long["variant"].str.split(":", expand=True)

# Split POS and alleles
geno_long[["POS", "alleles"]] = geno_long["rest"].str.split("_", expand=True)

# Split REF and ALT
geno_long[["REF", "ALT"]] = geno_long["alleles"].str.split(">", expand=True)

# Clean up
geno_long["POS"] = geno_long["POS"].astype(int)
geno_long = geno_long.drop(columns=["variant", "rest", "alleles"])


geno_long.to_csv(geno_out, index=False)
freq_df.to_csv(freq_out, index=False)

print("\n✅ Completed extracting genotypes and allele frequencies")
print(f"Genotypes saved to: {geno_out}")
print(f"Frequencies saved to: {freq_out}")
