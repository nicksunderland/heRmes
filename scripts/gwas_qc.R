#!/usr/bin/env Rscript
#
#  Nick Sunderland <nicholas.sunderland@bristol.ac.uk>
#
#  GWAS QC script
#
#  Command line options: run this script with the -h/--help flag to see options
#


#=============================================================================
# required packages
#=============================================================================
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(cli))
suppressPackageStartupMessages(library(data.table))


#=============================================================================
# parse cli arguments
#=============================================================================
cli_h2("Parsing arguments")

parser <- ArgumentParser()
# basic controls
parser$add_argument("-v", "--verbose",    action="store_true",  default=TRUE,     help="Print extra output [default]")
parser$add_argument("-q", "--quietly",    action="store_false", dest="verbose",   help="Print little output")
# file paths and file handling
parser$add_argument("-g", "--gwas",       action="store",       type="character", help="GWAS file path")
parser$add_argument("-r", "--ref",        action="store",       type="character", help="Reference file path")
parser$add_argument("-o", "--output",     action="store",       type="character", help="Output file path")
parser$add_argument("-ow", "--overwrite", action="store_true",  default=FALSE,    help="Overwrite output file if present")
# GWAS file columns
parser$add_argument("-chr", "--gwas_chr",  action="store", type="character", help="Chromosome column name in the GWAS file", default = "chr")
parser$add_argument("-bp",  "--gwas_bp",   action="store", type="character", help="Base position column name in the GWAS file", default = "bp")
parser$add_argument("-ea",  "--gwas_ea",   action="store", type="character", help="Effect allele column name in the GWAS file", default = "ea")
parser$add_argument("-oa",  "--gwas_oa",   action="store", type="character", help="Other allele column name in the GWAS file", default = "oa")
parser$add_argument("-eaf", "--gwas_eaf",  action="store", type="character", help="Effect allele frequency column name in the GWAS file", default = "eaf")
parser$add_argument("-beta","--gwas_beta", action="store", type="character", help="Beta column name in the GWAS file", default = "beta")
parser$add_argument("-se",  "--gwas_se",   action="store", type="character", help="Standard error column name in the GWAS file", default = "se")
parser$add_argument("-p",   "--gwas_p",    action="store", type="character", help="P-value column name in the GWAS file", default = "p")
parser$add_argument("-n",   "--gwas_n",    action="store", type="character", help="Sample size column name in the GWAS file", default = "n")
parser$add_argument("-info","--gwas_info", action="store", type="character", help="Imputation information score column name in the GWAS file", default = "info")
# Reference file columns
parser$add_argument("-r_chr", "--ref_chr", action="store", type="character", help="Reference chromosome column name in the GWAS file", default = "chr")
parser$add_argument("-r_bp",  "--ref_bp",  action="store", type="character", help="Reference base position column name in the GWAS file", default = "bp")
parser$add_argument("-r_ea",  "--ref_ea",  action="store", type="character", help="Reference effect allele column name in the GWAS file", default = "ea")
parser$add_argument("-r_oa",  "--ref_oa",  action="store", type="character", help="Reference other allele column name in the GWAS file", default = "oa")
parser$add_argument("-r_eaf", "--ref_eaf", action="store", type="character", help="Reference effect allele frequency column name in the GWAS file", default = "eaf")
# QC parameters
parser$add_argument("-gc", "--genomic_control", action="store_true", default=FALSE, help="Apply genomic control for adjusting study results [default=FALSE]")
parser$add_argument("-indel", "--indel_alleles", action="store_true", default=TRUE, help="Reatin indel alleles [default=TRUE]")
# parse CLI arguments
args <- parser$parse_args()

# print the arguments requested to the console
for (i in seq_along(args)) {
  cli_text("{.strong {names(args)[i]}}:          {.val {args[[i]]}}")
}


#=============================================================================
# start QC run
#=============================================================================
cli_h1("Running GWAS quality control")


#=============================================================================
# basic check files
#=============================================================================
cli_progress_step("checking file paths and columns names")

# Check GWAS file path provided and that we can read from in
if (is.null(args$gwas) || !file.exists(args$gwas) || file.access(args$gwas, 4) != 0) {
  cli_abort("GWAS file path `{.file {args$gwas}}` does not exist or is not readable")
}

# Check GWAS file column names all present
gwas_header <- fread(args$gwas, nrows = 0)
gwas_cols <- sapply(args[c("gwas_chr","gwas_bp","gwas_ea","gwas_oa","gwas_eaf","gwas_beta","gwas_se","gwas_p","gwas_n")], function(x) x %in% names(gwas_header))
if (!all(gwas_cols)) {
  for (i in seq_along(gwas_cols)) {
    if (gwas_cols[[i]]) {
      cli_alert_success("{names(gwas_cols)[i]}: {args[names(gwas_cols)[i]]}")
    } else {
      cli_alert_danger("{names(gwas_cols)[i]}: {args[names(gwas_cols)[i]]}")
    }
  }
  cli_abort("Columns not all found in GWAS file columns [{names(gwas_header)}]")
}

# Check reference file path provided and that we can read from in
if (is.null(args$ref) || !file.exists(args$ref) || file.access(args$ref, 4) != 0) {
  cli_abort("Reference file path `{.file {args$reference}}` does not exist or is not readable")
}

# Check GWAS file column names all present
ref_header <- fread(args$ref, nrows = 0)
ref_cols <- sapply(args[c("ref_chr","ref_bp","ref_ea","ref_oa","ref_eaf")], function(x) x %in% names(ref_header))
if (!all(ref_cols)) {
  for (i in seq_along(ref_cols)) {
    if (ref_cols[[i]]) {
      cli_alert_success("{names(ref_cols)[i]}: {args[names(ref_cols)[i]]}")
    } else {
      cli_alert_danger("{names(ref_cols)[i]}: {args[names(ref_cols)[i]]}")
    }
  }
  cli_abort("Columns not all found in reference file columns [{names(ref_header)}]")
}

# Check output file path provided and that we can write to it / overwrite it
if (is.null(args$output) || (file.access(args$output, 2) != 0 && args$overwrite)) {
  cli_abort("Output file path `{.file {args$output}}` not provided or does not have write permission")
} else if (file.exists(args$output) && !args$overwrite) {
  cli_abort("Output file path `{.file {args$output}}` already exists but overwrite flag is {.val {args$overwrite}}")
}


#=============================================================================
# read data
#=============================================================================
cli_progress_step("Reading GWAS")
Sys.sleep(1)
cli_progress_step("Reading reference")
Sys.sleep(1)


#=============================================================================
# check column data
#=============================================================================
cli_progress_step("Pre-QC summary")
Sys.sleep(1)


#=============================================================================
# data fixes
#=============================================================================
cli_progress_step("Data standardisation")
Sys.sleep(1)


#=============================================================================
# apply filter
#=============================================================================
cli_progress_step("Cleaning data")
Sys.sleep(1)


#=============================================================================
# harmonise alleles
#=============================================================================
cli_progress_step("Reference harmonisation")
Sys.sleep(1)


#=============================================================================
# allele frequency analysis
#=============================================================================
cli_progress_step("Analysing allele frequency")
Sys.sleep(1)


#=============================================================================
# PZ plot
#=============================================================================
cli_progress_step("Creating PZ plot")
Sys.sleep(1)


#=============================================================================
# population stratification analysis 1
#=============================================================================
cli_progress_step("Calculating lambda-GC")
Sys.sleep(1)


#=============================================================================
# population stratification analysis 2
#=============================================================================
cli_progress_step("Running LDSC")
Sys.sleep(1)


#=============================================================================
# save
#=============================================================================
cli_progress_step("File saved to {.file {args$output}}")
Sys.sleep(1)


#=============================================================================
# finished
#=============================================================================
cli_process_done()

