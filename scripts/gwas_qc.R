#!/usr/bin/env Rscript
#
#  Nick Sunderland <nicholas.sunderland@bristol.ac.uk>
#
#  GWAS QC script
#
#  Command line options= run this script with the -h/--help flag to see options
#  Some good resources here: https://github.com/swvanderlaan/MetaGWASToolKit/tree/master/SCRIPTS
#


#=============================================================================
# required packages
#=============================================================================
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(cli))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(stats))
#   devtools::install_github("mglev1n/ldscr")
suppressPackageStartupMessages(library(ldscr))


#=============================================================================
# parse cli arguments
#=============================================================================
cli_h1("QC script - active arguments")

parser <- ArgumentParser()
# basic controls
parser$add_argument("-v", "--verbose",    action="store_true",  default=TRUE,     help="Print extra output [default]")
parser$add_argument("-q", "--quietly",    action="store_false", dest="verbose",   help="Print little output")
# file paths and file handling
parser$add_argument("-g", "--gwas",       action="store",       type="character", help="GWAS file path")
parser$add_argument("-r", "--ref",        action="store",       type="character", help="Reference file path")
parser$add_argument("-o", "--output",     action="store",       type="character", help="Output directory path")
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
parser$add_argument("-info","--gwas_info", action="store", type="character", help="Imputation information score column name in the GWAS file", default = NULL)
# Reference file columns
parser$add_argument("-r_id",  "--ref_id",  action="store", type="character", help="Reference rsid column name in the REF file", default = "chr")
parser$add_argument("-r_chr", "--ref_chr", action="store", type="character", help="Reference chromosome column name in the REF file", default = "chr")
parser$add_argument("-r_bp",  "--ref_bp",  action="store", type="character", help="Reference base position column name in the REF file", default = "bp")
parser$add_argument("-r_ea",  "--ref_ea",  action="store", type="character", help="Reference effect allele column name in the REF file", default = "ea")
parser$add_argument("-r_oa",  "--ref_oa",  action="store", type="character", help="Reference other allele column name in the REF file", default = "oa")
parser$add_argument("-r_eaf", "--ref_eaf", action="store", type="character", help="Reference effect allele frequency column name in the REF file", default = "eaf")
# QC parameters
parser$add_argument("-adj",   "--adjustment", action="store", type="character", default=NULL, choices = c("lambda", "ldsc"), help="Apply genomic control p-value & standard error adjustment using either the GC lambda ('lambda') or the LD score regression intercept ('ldsc') [default=OFF]")
parser$add_argument("-noind", "--no_indel_alleles", action="store_true", default=FALSE, help="Remove indel alleles [default=FALSE]")
parser$add_argument("-fd",    "--freq_diff", action="store", type="numeric", default=0.2, help="Retain variants with GWAS:REF allele frequency different less than ... [default=0.2]")
parser$add_argument("-it",    "--info_thresh", action="store", type="numeric", default=NULL, help="Retain variants with info score greater than ... [default=NULL (off)]")
parser$add_argument("-an",    "--ancestry", action="store", type="character", default="EUR", choices = c("AFR", "AMR", "CSA", "EAS", "EUR", "MID"), help="Ancestry code [default=EUR]")

# parse CLI arguments
args <- parser$parse_args()

# print out the arguments used to the console
for (i in seq_along(args)) {
  cli_text("{.strong {names(args)[i]}} = {.val {args[[i]]}}")
}


#=============================================================================
# start QC run
#=============================================================================
cli_h1("Running GWAS quality control")

# summary table to be iteratively filled
summary <- data.table()

#=============================================================================
# basic check files
#=============================================================================
cli_progress_step("checking file paths and columns names")

# Check GWAS file path provided and that we can read from it
if (is.null(args$gwas) || !file.exists(args$gwas) || file.access(args$gwas, 4) != 0) {
  cli_abort("GWAS file path `{.file {args$gwas}}` does not exist or is not readable")
}

# Check GWAS file column names all present
gwas_header <- fread(args$gwas, nrows = 0)
gwas_cols <- sapply(args[c("gwas_chr","gwas_bp","gwas_ea","gwas_oa","gwas_eaf","gwas_beta","gwas_se","gwas_p","gwas_n")], function(x) x %in% names(gwas_header))
if (!all(gwas_cols)) {
  for (i in seq_along(gwas_cols)) {
    if (gwas_cols[[i]]) {
      cli_alert_success("{names(gwas_cols)[i]} = {args[names(gwas_cols)[i]]}")
    } else {
      cli_alert_danger("{names(gwas_cols)[i]} = {args[names(gwas_cols)[i]]}")
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
ref_cols <- sapply(args[c("ref_id", "ref_chr","ref_bp","ref_ea","ref_oa","ref_eaf")], function(x) x %in% names(ref_header))
if (!all(ref_cols)) {
  for (i in seq_along(ref_cols)) {
    if (ref_cols[[i]]) {
      cli_alert_success("{names(ref_cols)[i]} = {args[names(ref_cols)[i]]}")
    } else {
      cli_alert_danger("{names(ref_cols)[i]} = {args[names(ref_cols)[i]]}")
    }
  }
  cli_abort("Columns not all found in reference file columns [{names(ref_header)}]")
}

# Check output directory exists
if (is.null(args$output) || !dir.exists(args$output)) {
  cli_abort("Output directory path `{.file {args$output}}` not provided or does not exist")
}


#=============================================================================
# read GWAS and reference data
#=============================================================================
cli_progress_step("reading GWAS")
gwas <- fread(args$gwas)
cli_progress_step("reading reference")
ref <- fread(args$ref)
cli_process_done()

#=============================================================================
# check column data missingness
#=============================================================================
cli_h1("Checking column data")
cli_progress_step("checking column data")
gwas_cols <- unlist(args[ names(args)[grepl("^gwas_", names(args))] ])
ref_cols  <- unlist(args[ names(args)[grepl("^ref_", names(args))] ])

# total and NA counts
num_na_summary <- gwas[, lapply(.SD, function(col) c(num = .N, num_na = sum(is.na(col)), pct_na = sprintf("%.1f%%", 100 * (sum(is.na(col)) / .N)))), .SDcols = gwas_cols]

# add to summary table
summary <- rbind(summary,
                 data.table(column   = gwas_cols,
                            std_name = names(gwas_cols),
                            type     = sapply(gwas[, mget(gwas_cols)], typeof),
                            num      = as.integer(sapply(num_na_summary, `[[`, 1)),
                            num_na   = as.integer(sapply(num_na_summary, `[[`, 2)),
                            pct_na   = sapply(num_na_summary, `[[`, 3)))

# report to console
cli_process_done()
print(summary)


#=============================================================================
# check input data validity
#=============================================================================
cli_h1("Checking data validity")
cli_progress_step("checking data validity")

# per column functions that return true if that row is valid
col_fn_list <- list(
  gwas_chr  = function(x) x %in% 1:26,
  gwas_bp   = function(x) is.integer(x) & x > 0,
  gwas_ea   = function(x) grepl("^[ACTG]+$|^[DI]$", x),
  gwas_oa   = function(x) grepl("^[ACTG]+$|^[DI]$", x),
  gwas_eaf  = function(x) is.numeric(x) & x >= 0 & x <= 1,
  gwas_beta = function(x) is.numeric(x) & is.finite(x),
  gwas_se   = function(x) is.numeric(x) & x > 0 & is.finite(x),
  gwas_p    = function(x) is.numeric(x) & x >= 0 & x <=1,
  gwas_n    = function(x) is.integer(x) & x > 0,
  gwas_info = function(x) is.numeric(x) & x >= 0 & x <=1,
  ref_chr   = function(x) x %in% 1:26,
  ref_bp    = function(x) is.integer(x) & x > 0,
  ref_ea    = function(x) grepl("^[ACTG]+$|^[DI]$", x),
  ref_oa    = function(x) grepl("^[ACTG]+$|^[DI]$", x),
  ref_eaf   = function(x) is.numeric(x) & x >= 0 & x <= 1
)

# apply the functions to the corresponding columns and create summary table
valid_summary <- gwas[, Map(function(fn, col) {
  valid <- sum(fn(col), na.rm = TRUE)
  valid_pct <- sprintf("%.1f%%", 100 * (valid / .N))
  c(valid = valid, valid_pct = valid_pct)
}, col_fn_list[names(gwas_cols)], .SD), .SDcols = gwas_cols]

# add to main summary table
summary <- cbind(summary,
                 data.table(valid     = as.integer(sapply(valid_summary, `[[`, 1)),
                            valid_pct = sapply(valid_summary, `[[`, 2)))

# print summary to console as well as the validity functions used
cli_process_done()
print(summary)
cli_div(theme = list(span.var = list(color = "gray"), span.code = list(color = "gray")))
cli_h2("Column validation functions:")
for (i in seq_along(gwas_cols)) {
  cli_text("{.field {gwas_cols[i]}:} {.code {deparse(col_fn_list[[names(gwas_cols)[i]]])[2]}}")
}


#=============================================================================
# attempt recoding / data fixes
#=============================================================================
cli_h1("Formatting data")
cli_progress_step("formatting data")

# (re)code chromosome column as integer
gwas[, (args$gwas_chr) := fcase(as.integer(get(args$gwas_chr)) %in% 1:26, as.integer(get(args$gwas_chr)),
                                grepl("(?i)^X$",   get(args$gwas_chr)),         23L,
                                grepl("(?i)^Y$",   get(args$gwas_chr)),         24L,
                                grepl("(?i)^(PAR|XY|YX)$", get(args$gwas_chr)), 25L,
                                grepl("(?i)^MT$",  get(args$gwas_chr)),         26L,
                                default = NA_integer_)]

# parse base position and n-sample columns to integer
integer_cols <- c(args$gwas_bp, args$gwas_n)
gwas[, (integer_cols) := lapply(.SD, as.integer), .SDcols = integer_cols]

# parse frequency, beta, se, p, info columns to numeric
numeric_cols <- c(args$gwas_eaf, args$gwas_beta, args$gwas_se, args$gwas_p, args$gwas_info)
gwas[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]

# remove allele frequencies <0 or >1
gwas[get(args$gwas_eaf) < 0 | get(args$gwas_eaf) > 1, (args$gwas_eaf) := NA_real_]

# remove infinite betas
gwas[is.infinite(get(args$gwas_beta)), (args$gwas_beta) := NA_real_]

# remove infinite, zero, or negative standard errors
gwas[is.infinite(get(args$gwas_se)) | get(args$gwas_se) <= 0, (args$gwas_se) := NA_real_]

# recode pvalue=0 to minimum machine precision and remove p>1 or p<0
gwas[, (args$gwas_p) := fcase(get(args$gwas_p) == 0, .Machine$double.xmin,
                              get(args$gwas_p) > 0 & get(args$gwas_p) <= 1, get(args$gwas_p),
                              default = NA_real_)]

# remove info scores >1 or <0
if (!is.null(args$info_thresh) & !is.null(args$gwas_info)) {
  gwas[get(args$gwas_info) < 0 | get(args$gwas_info) > 1 | get(args$gwas_info) < args$info_thresh, (args$gwas_info) := NA_real_]
} else if (is.null(args$info_thresh) & !is.null(args$gwas_info)) {
  cli_alert_info("info_thresh not provided, skipping info score filtering and setting info score to 1 for all variants")
  gwas[, (args$gwas_info) := 1]
}

# alleles must be characters
gwas[, c(args$gwas_ea, args$gwas_oa) := lapply(.SD, as.character), .SDcols = c(args$gwas_ea, args$gwas_oa)]

# find indels and and number to summary, remove indels if flag specified
indel_idx <- which(gwas[, .(nchar(get(args$gwas_ea)) > 1 | nchar(get(args$gwas_oa)) > 1 | grepl("(?i)^[DI]$", get(args$gwas_ea)) | grepl("(?i)^[DI]$", get(args$gwas_oa)))]$V1)
if (args$no_indel_alleles) {
  gwas[indel_idx, c(args$gwas_ea, args$gwas_oa) := NA_character_]
  cli_alert_warning("removing {length(indel_idx)} indel alleles")
}
summary[grep("^gwas_(ea|oa)$", std_name), indels := length(indel_idx)]

# ensure alleles are upper case ACTG or D/I
gwas[, c(args$gwas_ea, args$gwas_oa) := lapply(.SD, function(x) ifelse(grepl("(?i)^[ACTG]+$|^[DI]$", x), toupper(x), NA_character_)), .SDcols = c(args$gwas_ea, args$gwas_oa)]

# summarise the counts during this recoding process
postfix_na_summary <- gwas[, lapply(.SD, function(col) c(postfix_valid = sum(!is.na(col)), pct_na = sprintf("%.1f%%", 100 * (sum(!is.na(col)) / .N)))), .SDcols = gwas_cols]

# add to overall summary table
summary <- cbind(summary,
                 data.table(postfix_valid     = as.integer(sapply(postfix_na_summary, `[[`, 1)),
                            postfix_valid_pct = sapply(postfix_na_summary, `[[`, 2)))

# remove invalid rows (those containing NAs)
gwas <- gwas[ stats::complete.cases(gwas[, mget(gwas_cols)]) ]

# print updated summary to console
print(summary)

# exit if no variants remain
if (nrow(gwas) == 0) {
  cli_abort("no variants remain after data validation filters - check data and summary above")
} else {
  cli_process_done()
}

#=============================================================================
# formatting reference
#=============================================================================
col_fn_list <- list(
  ref_id    = as.character,
  ref_chr   = as.character,
  ref_bp    = as.integer,
  ref_ea    = as.character,
  ref_oa    = as.character,
  ref_eaf   = as.numeric
)
ref[ , (ref_cols)  := Map(function(fn, col) fn(col), col_fn_list[names(ref_cols) ], .SD), .SDcols = ref_cols ]


#=============================================================================
# harmonise alleles
#=============================================================================
cli_h1("Harmonising to reference")
cli_progress_step("harmonising data")

# currently the harmonising function is in my package, but probably best to supply as a stand alone script
h <- genepi.utils::harmonise_gwas(gwas, ref, join = "chr:bp", action = 2,
                                  gmap = c(chr = args$gwas_chr, bp = args$gwas_bp, ea = args$gwas_ea, oa = args$gwas_oa, eaf = args$gwas_eaf, beta = args$gwas_beta),
                                  rmap = c(chr = args$ref_chr, bp = args$ref_bp, ea = args$ref_ea, oa = args$ref_oa, eaf = args$ref_eaf))

# harmonisation summary, add to main summary
summary[, `:=`(harmonised = nrow(h), harmonised_pct = sprintf("%.1f%%", 100*(nrow(h)/num[1])))]

# print updated summary to console
print(summary)

# exit if no variants remain
if (nrow(h) == 0) {
  cli_warn("no variants remain after harmonisation - check data and summary above")
}
cli_process_done()


#=============================================================================
# allele frequency analysis
#=============================================================================
cli_h1("Anaylsing GWAS:REF allele frequency")
cli_progress_step("analysing allele frequency difference")

# absolute frequency difference cohort vs reference
h[, freq_diff := abs(eaf - eaf_ref)]

# add frequency difference counts to summary table
diff_col <- paste0("freq_diff_lt", args$freq_diff)
summary[, c(diff_col, paste0(diff_col, "_pct")) := .(sum(h$freq_diff < args$freq_diff), sprintf("%.1f%%", 100*(sum(h$freq_diff < args$freq_diff)/num[1])))]

# print summary to console
cli_process_done()
print(summary)

# plot the frequency differences
eaf_plot <- ggplot(h[freq_diff > args$freq_diff], aes(x = eaf_ref, y = eaf, color = freq_diff)) +
  geom_point(size = 3) +
  geom_point(data = h[freq_diff > args$freq_diff & strand_flip == TRUE], shape = 17, color = "red", size = 3) +
  geom_abline(slope = 1, intercept =  args$freq_diff, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = -1 * args$freq_diff, linetype = "dashed", color = "red") +
  viridis::scale_colour_viridis(option = "mako") +
  labs(x = "Reference allele frequency", y = "Cohort allele frequency",
       caption = "*red triangles strand flip") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none")

# save plot
cli_progress_step("plotting allele frequency difference")
grDevices::png(file.path(args$output, "eaf_plot.png"), width = 600, height = 600)
print(eaf_plot)
invisible(dev.off())
cli_process_done()

# filter out the frequency differences over the provided threshold
h <- h[freq_diff < args$freq_diff]

# extract the data from the harmonisation table (don't need the reference columns now, except for the rsid)
setnames(h,
         old = c(paste0(ref_cols[["ref_id"]], "_ref"), "chr","bp","ea","oa","eaf","beta", gwas_cols[["gwas_se"]], gwas_cols[["gwas_p"]], gwas_cols[["gwas_n"]]),
         new = c(ref_cols[["ref_id"]],
                 gwas_cols[["gwas_chr"]],
                 gwas_cols[["gwas_bp"]],
                 gwas_cols[["gwas_ea"]],
                 gwas_cols[["gwas_oa"]],
                 gwas_cols[["gwas_eaf"]],
                 gwas_cols[["gwas_beta"]],
                 gwas_cols[["gwas_se"]],
                 gwas_cols[["gwas_p"]],
                 gwas_cols[["gwas_n"]]))


#=============================================================================
# PZ plot
#=============================================================================
cli_h1("Assessing for analytical issues")

# generate PZ plot
pz_plot <- h[, `:=`(observed = -log10(get(gwas_cols[["gwas_p"]])),
                    expected = -log10(2*pnorm(abs(get(gwas_cols[["gwas_beta"]]) / get(gwas_cols[["gwas_se"]])), lower.tail=FALSE)))] |>
  ggplot(aes(x = expected, y = observed)) +
  geom_point(size = 0.5, color="darkblue") +
  geom_abline(slope=1, intercept=0, color="darkred", linetype = "dotted") +
  labs(x     = "P.ztest (-log10)",
       y     = "P (-log10)",
       color = NULL) +
  theme_minimal(base_size = 18) +
  theme(aspect.ratio = 1)

# save plot
cli_progress_step("plotting PZ plot")
grDevices::png(file.path(args$output, "pz_plot.png"), width = 600, height = 600)
print(pz_plot)
invisible(dev.off())
cli_process_done()


#=============================================================================
# population stratification analysis 1
#=============================================================================
cli_h1("Assessing for population stratification issues")
cli_progress_step("running LDSC")

# get Z score and run LDSC
h[, z := get(gwas_cols[["gwas_beta"]]) / get(gwas_cols[["gwas_se"]])]

# try to run LDSC
tryCatch({
  ldsc_res <- ldsc_h2(munged_sumstats = h[, list(SNP = get(ref_cols[["ref_id"]]), A1=get(gwas_cols[["gwas_ea"]]), A2=get(gwas_cols[["gwas_oa"]]), Z=z, N=get(gwas_cols[["gwas_n"]]))],
                      ancestry        = args$ancestry)
},
error = function(e) {
  cli_warn(paste0("LDSC failed: \n", e))
  ldsc_res <<- list(intercept = NA_real_)
})


#=============================================================================
# population stratification analysis 2
#=============================================================================
cli_progress_step("calculating lambda-GC")

# data for the QQ plot
h[, `:=`(chisq          = qchisq(log(get(gwas_cols[["gwas_p"]])), df = 1, lower.tail = FALSE, log.p = TRUE))] # https://stackoverflow.com/questions/40144267/calculating-miniscule-numbers-for-chi-squared-distribution-numerical-precisio
h[, `:=`(lambda         = median(chisq) / qchisq(0.5, 1))]
h[, `:=`(adj_chisq_gc   = chisq/lambda)]
h[, `:=`(adj_p_gc       = pchisq(adj_chisq_gc, 1, lower.tail=FALSE))]
h[, `:=`(adj_se_gc      = get(gwas_cols[["gwas_se"]]) * sqrt(lambda[1]))]
h[, `:=`(adj_chisq_ldsc = chisq/ldsc_res$intercept)]
h[, `:=`(adj_p_ldsc     = pchisq(adj_chisq_ldsc, 1, lower.tail=FALSE))]
h[, `:=`(adj_se_ldsc    = get(gwas_cols[["gwas_se"]]) * sqrt(ldsc_res$intercept))]

qq_data <- rbind(
  # unadjusted
  h[, .(stat     = "Unadjusted",
        value    = NA_real_,
        observed = -log10(sort(get(gwas_cols[["gwas_p"]]), decreasing=FALSE, na.last=TRUE)),
        expected = -log10(ppoints(.N)),
        clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:.N, shape2 = .N:1)),
        cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:.N, shape2 = .N:1)))],
  # GC lambda adjusted
  h[, .(stat     = "GC lambda",
        value    = lambda[1],
        observed = -log10(sort(adj_p_gc, decreasing=FALSE, na.last=TRUE)),
        expected = -log10(ppoints(.N)),
        clower   = NA_real_,
        cupper   = NA_real_)],
  # LDSC intercept adjusted
  h[, .(stat     = "LDSC intercept",
        value    = ldsc_res$intercept,
        observed = -log10(sort(adj_p_ldsc, decreasing=FALSE, na.last=TRUE)),
        expected = -log10(ppoints(.N)),
        clower   = NA_real_,
        cupper   = NA_real_)]
)
setorder(qq_data, expected)

# generate QQ-plot
# axis labels
log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

# lambda labels
labels <- qq_data[stat != "Unadjusted",
                  list(value    = value[1],
                       expected = 3.0,
                       observed = .GRP-1 + 0.25), by = "stat"]
labels[, label := sprintf("%s = %.3f", stat, value)]

# plot
# colors
point_colors <- if (is.null(args$adjustment)) {
                  c(Unadjusted = "blue", `GC lambda` = "lightgray", `LDSC intercept` = "darkgray")
                } else if (args$adjustment=="lambda") {
                  c(Unadjusted = "lightgray", `GC lambda` = "blue", `LDSC intercept` = "darkgray")
                } else if (args$adjustment=="ldsc") {
                  c(Unadjusted = "lightgray", `GC lambda` = "darkgray", `LDSC intercept` = "blue")
                }

qq_plot <- qq_data |>
  ggplot(aes(x = expected, y = observed, color = stat)) +
  geom_point(size = 0.5) +
  geom_ribbon(aes(x = expected, ymin = clower, ymax = cupper), alpha = 0.1, color="transparent") +
  geom_abline(slope=1, intercept=0, color="darkred", linetype = "dotted") +
  geom_text(data = labels, aes(x=expected, y=observed, label=label), hjust = 0, color="black", show.legend = FALSE) +
  scale_color_manual(values = point_colors) +
  labs(x = log10Pe,
       y = log10Po,
       color = "Adjustment") +
  theme_minimal() +
  theme(aspect.ratio    = 1,
        legend.position = "top")

# save plot
cli_progress_step("plotting QQ plot")
grDevices::png(file.path(args$output, "qq_plot.png"), width = 600, height = 600)
print(qq_plot)
invisible(dev.off())
cli_process_done()


#=============================================================================
# extract then save the clean GWAS and summary table
#=============================================================================

# if using adjustment with lambda GC or LDSC intercept, switch the P and SE columns
if (args$adjustment=="lambda") {
  h[, c(gwas_cols[["gwas_se"]], gwas_cols[["gwas_p"]]) := .(adj_se_gc, adj_p_gc)]
} else if (args$adjustment=="ldsc") {
  h[, c(gwas_cols[["gwas_se"]], gwas_cols[["gwas_p"]]) := .(adj_se_ldsc, adj_p_ldsc)]
}

# extract wanted columns
gwas <- h[, .SD, .SDcols = c(ref_cols[["ref_id"]], gwas_cols)]

# save gwas
sans_ext <- sub("([^.]+)\\.[[:alnum:]]+$", "\\1", sub("[.](gz|bz2|xz)$", "", basename(args$gwas)))
out_path <- file.path(args$output, paste0(sans_ext, '_clean.tsv.gz'))
cli_progress_step("saving clean GWAS file to {.file {out_path}}")
fwrite(gwas, out_path, sep = "\t")

# save summary
log_path <- file.path(args$output, paste0(sans_ext, '_log.tsv'))
cli_progress_step("saving log file to {.file {log_path}}")
fwrite(summary, log_path, sep = "\t")


#=============================================================================
# finished
#=============================================================================
cli_process_done()







#=============================================================================
# testing
#=============================================================================
if (FALSE) {
  args = list()
  args$gwas = '/Users/xx20081/Documents/local_data/hermes_incidence/raw/Pheno5-DCM_EUR/FORMAT-METAL_Pheno5-DCM_EUR.tsv.gz'
  args$ref ='/Users/xx20081/Documents/local_data/genome_reference/hrc_37/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz'
  args$gwas_beta= "A1_beta"
  args$gwas_bp= "pos_b37"
  args$gwas_chr= "chr"
  args$gwas_ea= "A1"
  args$gwas_eaf= "A1_freq"
  args$gwas_info= NULL
  args$gwas_n= "N_total"
  args$gwas_oa= "A2"
  args$gwas_p= "pval"
  args$gwas_se= "se"
  args$output= "/Users/xx20081/Desktop/qc_tests"
  args$overwrite= TRUE
  args$ref_id= "ID"
  args$ref_bp= "POS"
  args$ref_chr= "#CHROM"
  args$ref_ea= "ALT"
  args$ref_eaf= "AF"
  args$ref_oa= "REF"
  args$verbose= TRUE
  args$freq_diff= 0.2
  args$adjustment="ldsc"
  args$no_indel_alleles= FALSE
  args$info_thresh= NULL
  args$ancestry="EUR"
}
