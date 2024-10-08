# Rscript /Users/xx20081/git/heRmes/scripts/gwas_qc.R -h

Rscript /Users/xx20081/git/heRmes/scripts/gwas_qc.R \
  --gwas '/Users/xx20081/Desktop/meta.all.allcause_death.autosomes.tsv' \
  --gwas_chr 'chr' \
  --gwas_bp 'bp' \
  --gwas_ea 'ea' \
  --gwas_oa 'oa' \
  --gwas_eaf 'eaf' \
  --gwas_beta 'beta' \
  --gwas_se 'se' \
  --gwas_p 'p' \
  --gwas_n 'n' \
  --gwas_info 'info' \
  --ref '/Users/xx20081/Documents/local_data/genome_reference/hrc_37/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz' \
  --ref_chr '#CHROM' \
  --ref_bp 'POS' \
  --ref_ea 'ALT' \
  --ref_oa 'REF' \
  --ref_eaf 'AF' \
  --freq_diff 0.2 \
  --out '/Users/xx20081/Desktop/qc_tests' \
  --no_indel_alleles

