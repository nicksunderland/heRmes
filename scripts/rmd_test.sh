rmd_file="/Users/xx20081/git/heRmes/scripts/gwas_qc.Rmd"
gwas_qc_dir="/Users/xx20081/Desktop/qc_tests"
Rscript -e "rmarkdown::render('$rmd_file', params = list(fig_dir = '$gwas_qc_dir'), output_dir = '$gwas_qc_dir')"
