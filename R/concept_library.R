#' @title Update code library
#' @description
#' Update the internal code library from the Concept Library API
#' @return boolean, success status
#' @export
#'
#' @import ConceptLibraryClient
#' @import data.table
update_library <- function() {

  # connect to API
  client <- ConceptLibraryClient::Connection$new(public = TRUE)

  # get phenotypes
  search_terms <- c("heart failure", "cardiomyopathy", "myocardial infarction")
  search_results <- lapply(search_terms, function(x) client$phenotypes$get(search = x)[, c("name", "phenotype_id")]) |>
    data.table::rbindlist()
  pheno_ids <- search_results[grepl(paste0(search_terms, collapse = "|"), name, ignore.case = TRUE), phenotype_id]
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


