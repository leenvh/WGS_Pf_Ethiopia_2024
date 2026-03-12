# the command can be run 
# Rscript code_ehh.R ../selection/ethio_eritrea/hap_chr13_Ethiopia ../selection/ethio_eritrea/snp.info.inp.Ethiopia
library(rehh)
library(R.utils)
#install.packages("vcfR" , "data.table")
library(vcfR)
library(data.table)
library(qqman)
library(R.utils)
args = commandArgs(trailingOnly=TRUE)
 hap_chrhap_chrhap_chr_pop <- data2haplohh(hap_file= args[1] ,map_file = args[2],
                                recode.allele = FALSE,
                                chr.name = 13, min_perc_geno.hap = 0.8 , min_perc_geno.mrk = 0.8)
  
  res_k13 <- calc_ehh(hap_chrhap_chrhap_chr_pop, 
                mrk = "13_1725133", 
                include_nhaplo = TRUE,
                include_zero_values = TRUE,
                discard_integration_at_border = FALSE , polarized = TRUE   )

  res_675 <- calc_ehh(hap_chrhap_chrhap_chr_pop,
                mrk = "13_1724974",
                include_nhaplo = TRUE,
                include_zero_values = TRUE,
                discard_integration_at_border = FALSE , polarized = TRUE   )
jpeg("k13_ehh_plot_new_675.jpeg", width = 480, height = 480 , unit = "px")
plot(res_675, xlim = c(1.6E+6, 1.85E+6) )
dev.off()

jpeg("k13_ehh_plot_new.jpeg", width = 480, height = 480 , unit = "px")
plot(res_k13, xlim = c(1.6E+6, 1.85E+6) )
dev.off()

jpeg("k13_furcation_plot_new.jpeg", width = 480, height = 480 , unit = "px")
f <- calc_furcation(hap_chrhap_chrhap_chr_pop, mrk = "13_1725133")
plot(f,lwd = 1,
hap.names = hap.names(hap_chrhap_chrhap_chr_pop),
# no place for a legend inside the plot
legend.xy.coords = "none", xlim = c(1.7E+6, 1.8E+6) )
dev.off()