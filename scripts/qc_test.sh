# Rscript /Users/xx20081/git/heRmes/scripts/gwas_qc.R -h

Rscript /Users/xx20081/git/heRmes/scripts/gwas_qc.R \
  --gwas '/Users/xx20081/Downloads/bmi_knowledge_portal_subset.tsv' \
  --gwas_chr 'chromosome' \
  --gwas_bp 'position' \
  --gwas_ea 'alt' \
  --gwas_oa 'reference' \
  --gwas_eaf 'freq' \
  --gwas_beta 'beta' \
  --gwas_se 'stdErr' \
  --gwas_p 'pValue' \
  --gwas_n 'n' \
  --ref '/Users/xx20081/Documents/local_data/genome_reference/hrc_37/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz' \
  --ref_id 'ID' \
  --ref_chr '#CHROM' \
  --ref_bp 'POS' \
  --ref_ea 'ALT' \
  --ref_oa 'REF' \
  --ref_eaf 'AF' \
  --freq_diff 0.2 \
  --out '/Users/xx20081/Desktop/qc_tests' \
  --no_indel_alleles


Rscript /Users/xx20081/git/heRmes/scripts/gwas_qc.R \
  --gwas '/Users/xx20081/Downloads/sample-gwas-expanded.csv' \
  --gwas_chr 'CHR' \
  --gwas_bp 'BP' \
  --gwas_ea 'EA' \
  --gwas_oa 'OA' \
  --gwas_eaf 'EUR_EAF' \
  --gwas_beta 'BETA' \
  --gwas_se 'SE' \
  --gwas_p 'P' \
  --gwas_n 'N' \
  --ref '/Users/xx20081/Documents/local_data/genome_reference/hrc_37/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz' \
  --ref_id 'ID' \
  --ref_chr '#CHROM' \
  --ref_bp 'POS' \
  --ref_ea 'ALT' \
  --ref_oa 'REF' \
  --ref_eaf 'AF' \
  --freq_diff 0.2 \
  --out '/Users/xx20081/Desktop/qc_tests' \
  --no_indel_alleles

