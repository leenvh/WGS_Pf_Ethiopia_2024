
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(ggplot2)
library(ggrepel)
setwd("/mnt/storage13/ahri/popgen_ethiopia_WGS/logestic_regression/hrp_status/")
# Set working directory and number of cores
#setwd("//wsl.localhost/Ubuntu-24.04/home/metmiku/Leucodepleted/whole-genome-sequence-analysis/logestic_regression")
n_cores <- 16

# Load genotype matrix from the new generated vcf file
geno <- fread("/mnt/storage13/ahri/popgen_ethiopia_WGS/all_country_merged/WGS_Pf_Ethiopia_2024/ethiopia_vcf/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps_ethiopia.mat.bin", sep = "\t", header = TRUE, data.table = FALSE)
# geno <- fread("072025_pf_mgen_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.mat.bin", sep = "\t", header = TRUE, data.table = FALSE)
#View(geno)
metadata_hrp2_k13 <- read.csv("metadata_hrp2_k13.csv")

#first filtering for both deleted vs other and test the analysis then do for others 
hrp2_stat <- metadata_hrp2_k13%>%
  dplyr::select(Sample,hrp_deletion_call)%>%
  mutate(hrp2 = ifelse(hrp_deletion_call == "hrp2-/hrp3-" ,1,0))%>%
  dplyr::select(Sample,hrp2)%>%
  mutate(hrp2 = as.integer(hrp2))


# Mask 0.5 with NA
geno[geno == 0.5] <- NA

# Load PC covariates
pcs <- read.csv("Sep2025_pca_coordinates.csv")
colnames(pcs) <- c("sample", "PC1", "PC2", "PC3", "PC4","PC5")

# Keep only samples present in both genotype and PCs
common_samples <- intersect(colnames(geno)[4:ncol(geno)], pcs$sample)
geno <- geno[, c("chr", "pos", "ref", common_samples)]
pcs <- pcs %>% filter(sample %in% common_samples)
#write.csv(geno,"geno.csv",row.names=TRUE)

# Extract kelch13 mutation row (e.g., by annotation)
# kelch_row <- which(geno$chr == "Pf3D7_13_v3" & geno$pos == 1725133)
# kelch_geno <- as.numeric(geno[kelch_row, 4:ncol(geno)])


# Extracting hrp2 deletion samples 
common_samples <- intersect(colnames(geno)[4:ncol(geno)], hrp2_stat$Sample)

geno <- geno[, c("chr", "pos", "ref", common_samples)]
hrp2_df <- hrp2_stat %>% filter(Sample %in% common_samples)
pcs <- pcs %>% filter(sample %in% common_samples)
y <- hrp2_df$hrp2



# Filter genotype matrix to exclude kelch13 itself
geno_test <- geno 

# Parallel association test
geno_rows <- split(geno_test, seq(nrow(geno_test)))

#View(geno_rows)

# With PC components
results <- lapply(geno_rows, function(row) {
  snp2 <- as.numeric(row[1, 4:ncol(row)])
  interaction_term <- y * snp2
  if (sum(interaction_term, na.rm = TRUE) > 1) {
    print("passed interaction")
    chr <- row["chr"]
    pos <- row["pos"]
    ref <- row["ref"]
    df <- data.frame(
      y = y,
      snp2 = snp2,
      PC1 = pcs$PC1,
      PC2 = pcs$PC2,
      PC3 = pcs$PC3,
      PC4 = pcs$PC4,
      PC5 = pcs$PC5
    )
    df <- df[complete.cases(df), ]
    if (nrow(df) > 10) {
      print("exceeds sample number")
      
      # Check if snp2 has variation
      if (length(unique(df$snp2[!is.na(df$snp2)])) < 2) {
      print("snp2 has no variation, skipping")
      return(NULL)
      }
      
      # Fit model with error handling
      fit <- tryCatch({
        glm(y ~ snp2 + PC1 + PC2 + PC3 + PC4 + PC5, data = df, family = "binomial") 
      }, error = function(e) {
        print(paste("Model fitting error:", e$message))
        return(NULL)
      })
      
      if (is.null(fit)) return(NULL)
      
      coef_table <- coef(summary(fit))
      
      # Check if snp2 coefficient exists
      if (!"snp2" %in% rownames(coef_table)) {
        print("snp2 coefficient not in model, skipping")
        return(NULL)
      }
      
      return(data.frame(
        chr = chr,
        pos = pos,
        ref = ref,
        beta_snp2 = coef_table["snp2", "Estimate"],
        p_snp2 = coef_table["snp2", "Pr(>|z|)"],
        n = nrow(df)
      ))
    }
  }
  return(NULL)
})

results_df <- do.call(rbind, results)
results_df$p_adj <- p.adjust(results_df$p_snp2, method="BH")
write.csv(results_df, "/mnt/storage13/ahri/popgen_ethiopia_WGS/all_country_merged/WGS_Pf_Ethiopia_2024/logestic_regression/hrp_status/both_deleted_vs_wildinteraction_results_.csv", row.names = FALSE)

system("bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%BCSQ\\n' /mnt/storage13/ahri/popgen_ethiopia_WGS/072025_pf_mgen_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz > annotations.tsv")

make_anno <- function(annotations) {
  anno <- read.table(annotations, sep = "\t", header = FALSE)
  
  # Split the 5th column by '|' and extract gene names
  anno_split <- anno %>%
    separate(V5, into = c("field1", "gene_name1", "gene_name2"), 
             sep = "\\|", extra = "drop", fill = "right")
  
  # Create combined gene name
  anno_final <- anno_split %>%
    mutate(Gene_name = paste0(gene_name1, " (", gene_name2, ")")) %>%
    select(V1, V2, V3, V4, Gene_name)
  
  # Rename columns
  colnames(anno_final) <- c("chr", "pos", "ref", "alt", "Gene_name")
  
  # Write to file
  write.table(anno_final, "annotations_final.tsv", 
              sep = "\t", row.names = FALSE, quote = FALSE)
}
make_anno("annotations.tsv")

annotations_final <- read.csv("annotations_final.tsv", sep = "\t")
results_df_ann <- left_join(results_df, annotations_final, 
                                        by = c("chr" = "chr", "pos" = "pos"))
write.csv(results_df_ann, "/mnt/storage13/ahri/popgen_ethiopia_WGS/all_country_merged/WGS_Pf_Ethiopia_2024/logestic_regression/hrp_status/both_deleted_vs_wildinteraction_results_ann.csv", row.names = FALSE)

View(results_df_ann)


results_df_ann <- read.csv("hrp2_deletion_interaction_results_ann.csv")
# used for plotting the interaction
data<-results_df_ann

plot_data <- data %>%
  mutate(
    logp = -log10(p_adj),                # avoid log(0)
    chr = factor(chr, levels = paste0("Pf3D7_", sprintf("%02d", 1:14), "_v3")),
    pos_mb = pos / 1e6                            # convert to Mb
  )

plot_data <- plot_data %>%
  mutate(significant = ifelse(p_adj < 0.05, "Significant", "Not significant"))

plot_data <- plot_data %>%
  mutate(
    Gene_name_clean = sub("\\s*\\(.*\\)", "", Gene_name)  # keep only name before '('
  ) %>%
  filter(chr!="NA")


plot<-ggplot(plot_data, aes(x = pos_mb, y = logp, color = significant)) +
  geom_point(alpha = 0.8, size = 1) +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  geom_text_repel(
    data = subset(
      plot_data,
      p_adj < 0.05 &
        !is.na(Gene_name) &
        Gene_name_clean != "NA" &
        Gene_name_clean != "NA (NA)" &
        !grepl("^PF3D7", Gene_name)  &         # exclude PF3D7 genes
        !grepl("^\\s*$", Gene_name)            # exclude empty names
    ) %>%
      transform(Gene_name = sub("\\s*\\(.*\\)", "", Gene_name)),  # keep only name before '('
    aes(label = Gene_name),
    size = 3,
    max.overlaps = 20,
    min.segment.length = 0,
    segment.color = "grey60"
  ) +
  facet_wrap(~ chr, nrow = 2, scales = "free_x") +
  scale_color_manual(
    values = c("Significant" = "red", "Not significant" = "grey70")
  ) +
  theme_minimal() +
  labs(
    x = "Chromosome position (Mb)",
    y = expression(-log[10](p[adj])),
    color = ""
  ) +
  theme(legend.position = "none")
plot
ggsave("/mnt/storage13/ahri/popgen_ethiopia_WGS/all_country_merged/WGS_Pf_Ethiopia_2024/logestic_regression/hrp_status/both_deleted_vs_wild_interaction_manhatten.pdf", plot,
       width = 30, height = 18, units = "cm", dpi = 600)

