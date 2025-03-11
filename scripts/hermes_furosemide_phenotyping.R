library(data.table)
library(bit64)
library(readxl)

# data dir
data_dir <- file.path(Sys.getenv("DATA_DIR"), "ukbb_81499_20241114")

# read all of the medication data
gp_medication <- fread(file.path(data_dir, "dummy_gp_medication.tsv.gz"))

# filters
hf_regex       <- "(?i)((heart|cardiac|ventric.*?).*(failure|insuffi.*))|entresto|sacubitril"
diuretic_regex <- "(?i)furosemide|butmetanide"

# read the code mapping file, combine and filter
# all_lkps_maps_v4.xlsx << https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=592
mapping_file <- "/Users/xx20081/Downloads/primarycare_codings/all_lkps_maps_v4.xlsx"
codes <- rbindlist(list(
  # heart_failure
  icd9     = as.data.table(read_xlsx(mapping_file, sheet="icd9_lkp",          col_types = "text"))[, .(concept="Heart Failure", code=paste0("'",ICD9,"'"), desc=DESCRIPTION_ICD9)][grepl(hf_regex, desc)],
  icd10    = as.data.table(read_xlsx(mapping_file, sheet="icd10_lkp",         col_types = "text"))[, .(concept="Heart Failure", code=paste0("'",ALT_CODE,"'"), desc=DESCRIPTION)][grepl(hf_regex, desc)],
  read2    = as.data.table(read_xlsx(mapping_file, sheet="read_v2_lkp",       col_types = "text"))[, .(concept="Heart Failure", code=paste0("'",read_code,"'"), desc=term_description)][grepl(hf_regex, desc)],
  read3    = as.data.table(read_xlsx(mapping_file, sheet="read_ctv3_lkp",     col_types = "text"))[, .(concept="Heart Failure", code=paste0("'",read_code,"'"), desc=term_description)][grepl(hf_regex, desc)],
  read_med = as.data.table(read_xlsx(mapping_file, sheet="read_v2_drugs_lkp", col_types = "text"))[, .(concept="Heart Failure", code=paste0("'",read_code,"'"), desc=term_description)][grepl(hf_regex, desc)],
  bnf      = as.data.table(read_xlsx(mapping_file, sheet="bnf_lkp",           col_types = "text"))[, .(concept="Heart Failure", code=paste0("'",BNF_Presentation_Code,"'"), desc=BNF_Presentation)][grepl(hf_regex, desc)],
  dmd      = as.data.table(read_xlsx(mapping_file, sheet="dmd_lkp",           col_types = "text"))[, .(concept="Heart Failure", code=paste0("'",concept_id,"'"), desc=term)][grepl(hf_regex, desc)],
  read_med = gp_medication[read_2!="",       .(concept="Heart Failure", code=paste0("'",read_2,"'"),   desc=drug_name)][grepl(hf_regex, desc)],
  bnf      = gp_medication[bnf_code!="",     .(concept="Heart Failure", code=paste0("'",bnf_code,"'"), desc=drug_name)][grepl(hf_regex, desc)],
  dmd      = gp_medication[!is.na(dmd_code), .(concept="Heart Failure", code=paste0("'",dmd_code,"'"), desc=drug_name)][grepl(hf_regex, desc)],
  # diruetic
  read_med = as.data.table(read_xlsx(mapping_file, sheet="read_v2_drugs_lkp", col_types = "text"))[, .(concept="Loop Diuretic", code=paste0("'",read_code,"'"), desc=term_description)][grepl(diuretic_regex, desc)],
  bnf      = as.data.table(read_xlsx(mapping_file, sheet="bnf_lkp",           col_types = "text"))[, .(concept="Loop Diuretic", code=paste0("'",BNF_Presentation_Code,"'"), desc=BNF_Presentation)][grepl(diuretic_regex, desc)],
  dmd      = as.data.table(read_xlsx(mapping_file, sheet="dmd_lkp",           col_types = "text"))[, .(concept="Loop Diuretic", code=paste0("'",concept_id,"'"), desc=term)][grepl(diuretic_regex, desc)],
  read_med = gp_medication[read_2!="",       .(concept="Loop Diuretic", code=paste0("'",read_2,"'"),   desc=drug_name)][grepl(diuretic_regex, desc)],
  bnf      = gp_medication[bnf_code!="",     .(concept="Loop Diuretic", code=paste0("'",bnf_code,"'"), desc=drug_name)][grepl(diuretic_regex, desc)],
  dmd      = gp_medication[!is.na(dmd_code), .(concept="Loop Diuretic", code=paste0("'",dmd_code,"'"), desc=drug_name)][grepl(diuretic_regex, desc)]
), idcol="code_type")
codes <- unique(codes, by=c("code_type","concept","code"))

# write out for manual review
fwrite(codes[,.(concept,code,code_type,desc)], "~/Desktop/furosemide_codes_for_review.tsv", sep="\t")

# do manual review in excel and save codes to /inst/extdata/hermes_furosemide_codes/hermes_furosemide_codes.tsv

# end code curration

# now run hermes_furosemide_phenotyping_ukbb_rap.ipynb on the RAP





