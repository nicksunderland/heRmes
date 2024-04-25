#' @title Update code library
#' @description
#' Update the internal code library from the Concept Library API
#' @param search_terms a vector of strings to search the concept library phenotypes for
#' @param ids a vector of phenotypes IDs to pull
#' @param UKHDR_UN UK-HDR website username (if pulling unpublished phenotypes)
#' @param UKHDR_PW UK-HDR website password (if pulling unpublished phenotypes)
#' @return boolean, success status
#' @export
#'
#' @import ConceptLibraryClient
#' @import data.table
update_library <- function(search_terms = c("heart failure", "cardiomyopathy", "myocardial infarction"), ids = NULL, UKHDR_UN = NULL, UKHDR_PW = NULL) {

  # connect to API
  public_flag <- ifelse(is.null(UKHDR_UN) | is.null(UKHDR_PW), TRUE, FALSE)
  client <- ConceptLibraryClient::Connection$new(username = UKHDR_UN,
                                                 password = UKHDR_PW,
                                                 public   = public_flag)

  # get phenotypes
  search_results <- lapply(search_terms, function(x) client$phenotypes$get(search = x)[, c("name", "phenotype_id")]) |>
    data.table::rbindlist()
  pheno_ids <- c()
  if (length(search_results) > 0) {
    pheno_ids <- search_results[grepl(paste0(search_terms, collapse = "|"), name, ignore.case = TRUE), phenotype_id]
  }
  pheno_ids <- c(pheno_ids, ids)
  names(pheno_ids) <- pheno_ids

  for (i in seq_len(length(pheno_ids))) {

    id <- pheno_ids[[i]]

    cat("[i] reading phenotype id:", id)

    pheno_file <- file.path(system.file("extdata", "ukhdr_phenotypes", package = "heRmes"), paste0(id, ".yaml"))

    if (!file.exists(pheno_file)) {
      cat(" - downloading")

      client$phenotypes$save_definition_file(pheno_file, id)
      codes <- client$phenotypes$get_codelist(id) |> data.table::as.data.table()
      data.table::fwrite(codes, sub(".yaml$", "_codes.tsv", pheno_file), sep = "\t")

      cat(" - finished\n")

    } else {

      cat(" - skipping, already exists\n")

    }

  }

}


