library(data.table)

# data dir
data_dir <- file.path(Sys.getenv("DATA_DIR"), "ukbb_81499_20241114")

# read the data
withdraw  <- fread(file.path(data_dir, "withdraw81499_55_20230915.txt"))
demog     <- fread(file.path(data_dir, "data_participant.tsv"))
gp_clinical<-fread(file.path(data_dir, "dummy_gp_clinical.tsv.gz"))
gp_medication<-fread(file.path(data_dir, "dummy_gp_medication.tsv.gz"))

# read the codes
codes <- rbindlist(list(
  fread(system.file("extdata", "hermes_3_codes", "hermes_3_codes.tsv", package = "heRmes")),
  fread(system.file("extdata", "hermes_furosemide_codes", "hermes_furosemide_codes.tsv", package = "heRmes"))
))
self_reported_codes <- list(
  list(name = "Heart Failure",                      code = "1076", code_type = "ukbb_self_reported_illness"),
  list(name = "Myocardial infarction",              code = "1075", code_type = "ukbb_self_reported_illness"),
  list(name = "Hypertrophic cardiomyopathy",        code = "1588", code_type = "ukbb_self_reported_illness"),
  list(name = "Coronary artery bypass grafting",    code = "1095", code_type = "ukbb_self_reported_procedure"),
  list(name = "Percutaneous coronary intervention", code = "1070", code_type = "ukbb_self_reported_procedure")
)
codes <- rbind(codes,
               data.table(Concept     = paste0(sapply(self_reported_codes, function(x) x$name), " Self Reported"),
                          Code        = sapply(self_reported_codes, function(x) x$code),
                          Source      = sapply(self_reported_codes, function(x) x$code_type),
                          Description = sapply(self_reported_codes, function(x) x$name)))
codes[, `:=`(code      = sub("^'(.*?)'$", "\\1", Code),
             code_type = fcase(Source=="ICD10", "icd10",
                               Source=="ICD9",  "icd9",
                               Source=="OPCS4", "opcs4",
                               Source=="READ2", "read2",
                               Source=="CTV3",  "read3",
                               Source=="BNF",   "bnf",
                               Source=="DMD",   "dmd",
                               Source=="ukbb_self_reported_illness", "ukbb_self_reported_illness",
                               Source=="ukbb_self_reported_procedure", "ukbb_self_reported_procedure"))]
codes <- codes[!is.na(code_type)]

# long data of codes
cohort <- demog[, list(eid     = eid,
                       age     = as.integer(`21022-0.0`),
                       sex     = factor(`31-0.0`, levels = 0:1, labels = c("female", "male")))]

# gp clinical
gp_clinical[read_2 == "", read_2 := NA_character_]
gp_clinical[read_3 == "", read_3 := NA_character_]
gp_clinical[, date := as.Date(date)]
gp_clinical <- data.table::melt(gp_clinical,
                                id.vars = c("eid", "date"),
                                measure.vars  = c("read_2", "read_3"),
                                variable.name = "code_type",
                                value.name = "code",
                                na.rm = TRUE)
gp_clinical[, code_type := data.table::fcase(code_type == "read_2", "read2",
                                             code_type == "read_3", "read3")]

# gp medication
gp_medication[read_2 == "", read_2 := NA_character_]
gp_medication[bnf_code == "", bnf_code := NA_character_]
gp_medication[dmd_code == "", dmd_code := NA_character_]
gp_medication[, date := as.Date(date)]
code_cols <- c("read_2","bnf_code","dmd_code")
gp_medication[, (code_cols) := lapply(.SD, as.character), .SDcol = code_cols]
gp_medication <- data.table::melt(gp_medication,
                                  id.vars = c("eid", "date"),
                                  measure.vars  = c("read_2", "bnf_code", "dmd_code"),
                                  variable.name = "code_type",
                                  value.name = "code",
                                  na.rm = TRUE)
gp_medication[, code_type := data.table::fcase(code_type == "read_2",   "read2",
                                               code_type == "bnf_code", "bnf",
                                               code_type == "dmd_code", "dmd")]


