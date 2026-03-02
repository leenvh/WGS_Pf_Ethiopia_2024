# SPDX-License-Identifier: Apache-2.0
"""
Extract genotypes for missense variants in drug-resistance genes,
filtering out sample genotypes with depth DP < min_dp (default 5).
Produces a matrix with samples as columns and variants as rows.
"""

import subprocess
import pandas as pd
import numpy as np
import os
import re
import argparse

# Define gene regions to subset (update as needed)
genes = {
    "CRT":  "Pf3D7_07_v3:403000-406000",
    "K13":  "Pf3D7_13_v3:1724817-1728000",
    "MDR1": "Pf3D7_05_v3:957500-962000",
    "DHFR": "Pf3D7_04_v3:746000-751000",
    "DHPS": "Pf3D7_08_v3:548000-552000",
    "PX1":  "Pf3D7_07_v3:889800-899213",
    "UBP1": "Pf3D7_01_v3:188400-201400",
    "AP2MU": "Pf3D7_12_v3:716700-720700",
    "AAT1": "Pf3D7_06_v3:1213100-1217350",
}

# Define gene aliases to match in BCSQ field
gene_aliases = {
    "CRT":   ["CRT"],
    "K13":   ["K13"],
    "MDR1":  ["MDR1"],
    "DHFR":  ["DHFR", "DHFR-TS"],
    "DHPS":  ["DHPS", "PPPK-DHPS"],
    "PX1":   ["PX1", "PF3D7_0720700"],
    "UBP1":  ["UBP1"],
    "AP2MU": ["AP2-MU", "AP2MU"],
    "AAT1":  ["AAT1", "PF3D7_0629500"],
}

def run_cmd(cmd):
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR running {' '.join(cmd)}:\n{result.stderr}")
        return None
    return result.stdout

def annotate_vcf(vcf, ref_fasta, gff3, out_annotated_vcf):
    """Annotate with bcftools csq and index."""
    cmd = ['bcftools', 'csq', '-p', 'a', '-f', ref_fasta, '-g', gff3,
           '-o', out_annotated_vcf, '-O', 'z', vcf]
    if run_cmd(cmd) is None:
        return False
    return run_cmd(['bcftools', 'index', out_annotated_vcf]) is not None

def query_genotypes(vcf, region, samples):
    """
    Query genotypes for all samples in a region.
    Returns a DataFrame with samples as columns and variants as rows.
    """
    # Format to get GT:DP for each sample
    fmt = '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ'
    # Add sample format for each sample
    for sample in samples:
        fmt += f'\t[%GT:%DP]'
    fmt += '\n'
    
    cmd = ['bcftools', 'query', '-f', fmt, '-r', region, '-s', ','.join(samples), vcf]
    output = run_cmd(cmd)
    
    if output is None:
        return pd.DataFrame()
    
    # Parse output
    rows = []
    for line in output.strip().split('\n'):
        if not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) < 5:
            continue
            
        chrom, pos, ref, alt, csq = fields[0:5]
        sample_data = fields[5:]
        
        # Parse BCSQ to get gene and amino acid change
        gene_found = None
        aa_change = None
        
        for eff in csq.split(','):
            parts = eff.split('|')
            if len(parts) < 6:
                continue
            consequence = parts[0].lower()
            gene_tok1 = parts[1].upper() if len(parts) > 1 else ""
            gene_tok2 = parts[2].upper() if len(parts) > 2 else ""
            
            # Check if this is a missense variant in any of our target genes
            if "missense" in consequence:
                for gene_name, aliases in gene_aliases.items():
                    alias_set = {a.upper() for a in aliases}
                    if {gene_tok1, gene_tok2} & alias_set:
                        gene_found = gene_name
                        aa_change = parts[5] if len(parts) > 5 else f"{ref}>{alt}"
                        break
                if gene_found:
                    break
        
        if not gene_found:
            continue
            
        # Create variant ID
        variant_id = f"{gene_found}_{pos}_{ref}>{alt}_{aa_change}"
        
        # Process genotypes for each sample
        genotypes = {}
        for i, sample in enumerate(samples):
            if i >= len(sample_data):
                genotypes[sample] = './.'
                continue
                
            gt_dp = sample_data[i].split(':')
            if len(gt_dp) < 2:
                genotypes[sample] = './.'
                continue
                
            gt, dp = gt_dp[0], gt_dp[1]
            
            # Filter by depth
            try:
                if dp == '.' or int(dp) < 5:
                    genotypes[sample] = './.'
                else:
                    genotypes[sample] = gt
            except ValueError:
                genotypes[sample] = './.'
        
        # Add to rows
        row = {
            'variant_id': variant_id,
            'gene': gene_found,
            'chrom': chrom,
            'pos': int(pos),
            'ref': ref,
            'alt': alt,
            'aa_change': aa_change,
            **genotypes
        }
        rows.append(row)
    
    if not rows:
        return pd.DataFrame()
    
    df = pd.DataFrame(rows)
    return df

def get_vcf_samples(vcf):
    """Return the set of sample names present in the VCF."""
    out = run_cmd(['bcftools', 'query', '-l', vcf])
    if out is None:
        return set()
    return {x.strip() for x in out.strip().splitlines() if x.strip()}

def _safe(name: str) -> str:
    """Filesystem-safe folder/file token."""
    return re.sub(r'[^A-Za-z0-9._-]+', '_', str(name)).strip('_')

def run(vcf, ref_fasta, gff3, metadata_path, outdir,
        min_dp=5, group_by=None):
    """
    Extract genotypes for missense variants.
    If group_by is provided, creates separate matrices per group.
    Otherwise creates one matrix for all samples.
    """
    os.makedirs(outdir, exist_ok=True)
    
    # Get all samples in VCF
    vcf_samples = get_vcf_samples(vcf)
    if not vcf_samples:
        raise SystemExit("Failed to read sample names from VCF")
    
    print(f"Found {len(vcf_samples)} samples in VCF")
    
    # If metadata provided, use it for grouping
    if metadata_path and group_by:
        metadata = pd.read_csv(metadata_path, sep=",", dtype=str)
        if 'sample_id' not in metadata.columns:
            raise SystemExit("Metadata must have a 'sample_id' column.")
        if group_by not in metadata.columns:
            raise SystemExit(f"Metadata does not contain '{group_by}' column.")
        
        # Filter metadata to samples in VCF
        metadata = metadata[metadata['sample_id'].isin(vcf_samples)]
        
        # Create groups
        groups = metadata.groupby(group_by)['sample_id'].apply(list).to_dict()
        
        # Process each group
        all_matrices = []
        
        for group_value, group_samples in groups.items():
            print(f"\nProcessing group: {group_by} = {group_value} (n={len(group_samples)})")
            
            # Create group directory
            group_dir = os.path.join(outdir, _safe(group_value))
            os.makedirs(group_dir, exist_ok=True)
            
            # Process all genes for this group
            group_matrix = extract_all_genes(vcf, ref_fasta, gff3, group_samples, min_dp)
            
            if not group_matrix.empty:
                # Save group matrix
                out_file = os.path.join(group_dir, "genotype_matrix.csv")
                group_matrix.to_csv(out_file, index=False)
                print(f"Saved genotype matrix for {group_value} to {out_file}")
                
                # Also save transposed version (samples as columns)
                transposed = group_matrix.set_index('variant_id').T
                transposed.to_csv(os.path.join(group_dir, "genotype_matrix_transposed.csv"))
                
                # Add group column and append to all matrices
                group_matrix['group'] = group_value
                all_matrices.append(group_matrix)
        
        # Combine all groups if we have multiple
        if len(all_matrices) > 1:
            combined = pd.concat(all_matrices, ignore_index=True)
            combined_out = os.path.join(outdir, f"genotype_matrix_all_groups.csv")
            combined.to_csv(combined_out, index=False)
            print(f"\nSaved combined matrix to {combined_out}")
    
    else:
        # No grouping - process all samples together
        print("\nProcessing all samples together...")
        all_samples = list(vcf_samples)
        matrix = extract_all_genes(vcf, ref_fasta, gff3, all_samples, min_dp)
        
        if not matrix.empty:
            out_file = os.path.join(outdir, "genotype_matrix.csv")
            matrix.to_csv(out_file, index=False)
            print(f"Saved genotype matrix to {out_file}")
            
            # Save transposed version (samples as columns, variants as rows)
            transposed = matrix.set_index('variant_id').T
            transposed_out = os.path.join(outdir, "genotype_matrix_transposed.csv")
            transposed.to_csv(transposed_out)
            print(f"Saved transposed matrix to {transposed_out}")

def extract_all_genes(vcf, ref_fasta, gff3, samples, min_dp):
    """Extract genotypes for all genes for a given set of samples."""
    all_variants = []
    
    # Create temporary subset VCF
    subset_vcf = os.path.join(os.path.dirname(vcf), f"temp_subset_{hash(tuple(samples))}.vcf.gz")
    
    # Subset VCF to these samples
    cmd = ['bcftools', 'view', '-s', ','.join(samples), '-Oz', '-o', subset_vcf, vcf]
    if run_cmd(cmd) is None:
        print("Failed to subset VCF")
        return pd.DataFrame()
    
    run_cmd(['bcftools', 'index', subset_vcf])
    
    for gene, region in genes.items():
        print(f"  Processing {gene} at {region}...")
        
        # Create gene-specific VCF
        gene_vcf = os.path.join(os.path.dirname(vcf), f"temp_{gene}_{hash(tuple(samples))}.vcf.gz")
        
        cmd = ['bcftools', 'view', '-r', region, '-Oz', '-o', gene_vcf, subset_vcf]
        if run_cmd(cmd) is None:
            continue
        
        run_cmd(['bcftools', 'index', gene_vcf])
        
        # Annotate
        annotated_vcf = gene_vcf.replace('.vcf.gz', '_csq.vcf.gz')
        if not annotate_vcf(gene_vcf, ref_fasta, gff3, annotated_vcf):
            continue
        
        # Query genotypes
        df = query_genotypes(annotated_vcf, region, samples)
        if not df.empty:
            all_variants.append(df)
        
        # Cleanup temp files
        for f in [gene_vcf, gene_vcf + '.csi', annotated_vcf, annotated_vcf + '.csi']:
            if os.path.exists(f):
                os.remove(f)
    
    # Cleanup subset VCF
    for f in [subset_vcf, subset_vcf + '.csi']:
        if os.path.exists(f):
            os.remove(f)
    
    if not all_variants:
        return pd.DataFrame()
    
    # Combine all variants
    combined = pd.concat(all_variants, ignore_index=True)
    
    # Reorder columns: metadata first, then samples
    meta_cols = ['variant_id', 'gene', 'chrom', 'pos', 'ref', 'alt', 'aa_change']
    sample_cols = [col for col in combined.columns if col not in meta_cols]
    
    combined = combined[meta_cols + sample_cols]
    
    return combined

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract genotype matrix for missense variants in drug-resistance genes")
    parser.add_argument("--vcf", required=True, help="Input VCF.gz file")
    parser.add_argument("--ref", required=True, help="Reference FASTA file")
    parser.add_argument("--gff", required=True, help="GFF3 annotation file")
    parser.add_argument("--metadata", help="Metadata TSV file with sample_id column (optional)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--min_dp", type=int, default=5, help="Minimum depth (default: 5)")
    parser.add_argument("--group_by", help="Metadata column to group by (optional)")
    
    args = parser.parse_args()
    
    run(
        vcf=args.vcf,
        ref_fasta=args.ref,
        gff3=args.gff,
        metadata_path=args.metadata,
        outdir=args.outdir,
        min_dp=args.min_dp,
        group_by=args.group_by
    )