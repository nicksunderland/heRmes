# Silence R CMD check
globalVariables(c("description", "code", "name", "phenotype_id", "V1", "V2", "V3", "ReadCode", "Description", "conceptId", "term", "pheno_codes"),
                package = "heRmes")

#' @title get phenotype codes
#' @param pheno_id a string, a valid phenotype id (see `get_phenotypes()`)
#' @param strip a regular expression to strip from the codes, NULL for no stripping
#' @return a data.table
#' @export
#'
get_codes <- function(pheno_id, strip = "\\.") {

  # get the phenotypes in the local library
  pheno_file <- system.file("extdata", "ukhdr_phenotypes", paste0(pheno_id, "_codes.tsv"), package = "heRmes")

  # check
  stopifnot("phenotype file not found" = file.exists(pheno_file))

  # read the file
  dat <- data.table::fread(pheno_file)
  dat[, code := as.character(code)]

  if (!is.null(strip)) {
    dat[, code := gsub(strip, "", code)]
  }

  # return
  return(dat)

}

#' @title get phenotype meta-data
#' @param pheno_id a string, a valid phenotype id (see `get_phenotypes()`)
#' @return a data.table
#' @export
#'
get_metadata <- function(pheno_id) {

  # get the phenotypes in the local library
  meta_file <- system.file("extdata", "ukhdr_phenotypes", paste0(pheno_id, ".yaml"), package = "heRmes")

  # check
  if(!file.exists(meta_file)) stop("meta data file not found for ID: `", pheno_id, "`")

  # read the file
  meta <- yaml::read_yaml(meta_file)

  # return
  return(meta)

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
