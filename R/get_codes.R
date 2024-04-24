# Silence R CMD check
globalVariables(c("description", "code", "name", "phenotype_id", "V1", "V2", "V3", "ReadCode", "Description", "conceptId", "term", "pheno_codes"),
                package = "heRmes")

#' @title get phenotype codes
#' @param pheno_id a string, a valid phenotype id (see `get_phenotypes()`)
#' @return a data.table
#' @export
#'
get_codes <- function(pheno_id) {

  # get the phenotypes in the local library
  pheno_file <- system.file("extdata", "ukhdr_phenotypes", paste0(pheno_id, "_codes.tsv"), package = "heRmes")

  # read the file
  dat <- data.table::fread(pheno_file)
  dat[, code := as.character(code)]

  # return
  return(dat)

}


#' @title get icd10cm codes
#' @return a data.table
#' @export
#'
get_icd10cm <- function() {

  code_file <- system.file("extdata", "icd10cm", "icd10cm_april_2024", "icd10cm-codes-April-2024.txt", package = "heRmes")
  codes <- data.table::fread(code_file, header = FALSE, sep = NULL)
  codes[, V1 := sub("\\s+", "#,#", V1)]
  codes[, c("code", "description") := data.table::tstrsplit(V1, "#,#", fixed = TRUE)]
  codes[, V1 := NULL]
  return(codes)

}


#' @title get icd9cm codes
#' @return a data.table
#' @export
#'
get_icd9cm <- function() {

  code_file <- system.file("extdata", "icd9cm", "CMS32_DESC_LONG_DX.txt", package = "heRmes")
  codes <- data.table::fread(code_file, header = FALSE, sep = NULL)
  codes$V1 <- iconv(codes$V1, from = "ISO-8859-1", to = "UTF-8")
  codes[, V1 := sub("\\s+", "#,#", V1)]
  codes[, c("code", "description") := data.table::tstrsplit(V1, "#,#", fixed = TRUE)]
  codes[, V1 := NULL]
  return(codes)

}


#' @title get icd10 codes
#' @return a data.table
#' @export
#'
get_icd10 <- function() {

  code_file <- system.file("extdata", "icd10", "icd10_ed4_20120401.txt", package = "heRmes")
  codes <- data.table::fread(code_file, header = FALSE, skip = 1)
  codes <- codes[, list(code = V1, description = V2)]
  return(codes)

}


#' @title get SNOMED CT codes
#' @param version a string, either "int" (international) or "uk"
#' @return a data.table
#' @export
#'
get_snomed <- function(version = "int") {

  version <- match.arg(version, choices = c("uk", "int"))
  code_file <- system.file("extdata", "snomed", paste0("sct2_description_", version, ".txt"), package = "heRmes")
  codes <- data.table::fread(code_file, quote="")
  codes <- codes[, list(code = as.character(conceptId), description = term)]
  return(codes)

}

#' @title get read codes v3
#' @return a data.table
#' @export
#'
get_readcodes_v3 <- function() {

  code_file <- system.file("extdata", "read_codes", "read_codes_v3_terms.txt", package = "heRmes")
  codes <- data.table::fread(code_file, sep = "|", quote="")
  codes <- codes[, list(code = V1, description = V3)]
  return(codes)

}

#' @title get read codes v2
#' @return a data.table
#' @export
#'
get_readcodes_v2 <- function() {

  code_file <- system.file("extdata", "read_codes", "read_codes_v2_terms.txt", package = "heRmes")
  codes <- data.table::fread(code_file, sep = ",")
  codes <- codes[, list(code = ReadCode, description = Description)]
  return(codes)

}


#' @title get opcs codes
#' @return a data.table
#' @export
#'
get_opcs <- function() {

  code_file <- system.file("extdata", "opcs", "opcs_nov2022_v1.0.txt", package = "heRmes")
  codes <- data.table::fread(code_file, header = FALSE, sep = NULL)
  codes[, V1 := sub("\\s+", "#,#", V1)]
  codes[, c("code", "description") := data.table::tstrsplit(V1, "#,#", fixed = TRUE)]
  codes[, V1 := NULL]
  return(codes)

}



#' @title get read codes
#' @param regexes a list of regular expressions
#' @param terms a list of strings to search for
#' @param sources a vector of strings, one or more of c("icd9cm", "icd10cm", "icd10", "snomed", "opcs", "readcodes_v2", "readcodes_v3")
#' @return a data.table
#' @export
#'
search_codes <- function(regexes = list("heart\\s*failure",
                                  "ventric[ulare]+\\s*(dysfunction|failure)"),
                   terms = list("Transplantation of heart and lung"),
                   sources = c("icd9cm", "icd10cm", "icd10", "snomed", "opcs", "readcodes_v2", "readcodes_v3")) {

  sources <- match.arg(sources, choices = c("icd9cm", "icd10cm", "icd10", "snomed", "opcs", "readcodes_v2", "readcodes_v3"), several.ok = TRUE)
  codes <- list()
  if ("icd9cm" %in% sources)       codes <- c(codes, list(icd9cm = get_icd9cm()))
  if ("icd10cm" %in% sources)      codes <- c(codes, list(icd10cm = get_icd10cm()))
  if ("icd10" %in% sources)        codes <- c(codes, list(icd10 = get_icd10()))
  if ("snomed" %in% sources)       codes <- c(codes, list(snomed = get_snomed()))
  if ("opcs" %in% sources)         codes <- c(codes, list(opcs = get_opcs()))
  if ("readcodes_v2" %in% sources) codes <- c(codes, list(readcodes_v2 = get_readcodes_v2()))
  if ("readcodes_v3" %in% sources) codes <- c(codes, list(readcodes_v3 = get_readcodes_v3()))
  codes   <- data.table::rbindlist(codes, idcol = "source")
  codes   <- unique(codes)
  regexes <- paste0(regexes, collapse = "|")
  codes   <- codes[grepl(regexes, description) |
                     grepl(regexes, code) |
                     sapply(description, function(x) x %in% terms) |
                     sapply(code, function(x) x %in% terms), ]

  return(codes)

}

#' @title Get phenotypes
#' @return a character vector of potential phenotypes
#' @export
#'
get_phenotypes <- function() {

  # get the phenotypes in the local library
  pheno_lib <- list.files(system.file("extdata", "ukhdr_phenotypes", package = "heRmes"), pattern = ".yaml$", full.names = TRUE)

  # get the names
  phenos <- lapply(pheno_lib, function(file) {

    # get the meta data
    meta <- yaml::read_yaml(file)

    # extract id and name
    d <- list()
    d["name"] <- meta$name
    d["id"] <- meta$phenotype_id

    # return
    return(d)

  })

  n <- sapply(phenos, function(x) x[["name"]])
  phenos <- sapply(phenos, function(x) x[["id"]])
  names(phenos) <- n

  return(phenos)
}
