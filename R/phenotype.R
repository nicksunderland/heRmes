#' @title Phenotype
#' @param ids a vector of characters IDs
#' @param codes a vector of diagnosis codes
#' @param name a string, a name for the phenotype
#' @param include a string or list of strings, a valid phenotype id (see `get_phenotypes()`)
#' @param exclude a string or list of strings, a valid phenotype id (see `get_phenotypes()`)
#' @return a data.table
#' @export
#'
phenotype <- function(ids, codes, name, include, exclude = NULL) {

  # checks
  stopifnot("ids and codes should be the same length" = length(ids)==length(codes))

  # create the data.table
  dat <- data.table::data.table(id = ids, code = codes)

  # get the pheno codes
  include_codes <- lapply(include, function(x) get_codes(x)) |> data.table::rbindlist()

  # get the exclusion codes
  if (!is.null(exclude)) {
    exclude_codes <- lapply(exclude, function(x) get_codes(x)) |> data.table::rbindlist()
  } else {
    exclude_codes <- get_codes(include)[rep(FALSE, nrow(pheno_codes)), ]
  }

  # parse the input codes
  dat[include_codes, include := TRUE, on = "code"]
  dat[is.na(include), include := FALSE]
  dat[exclude_codes, exclude := TRUE, on = "code"]
  dat[is.na(exclude), exclude := FALSE]

  # compute the result
  res <- dat[, list(include = any(include), exclude = any(exclude)), by = "id"]
  res[, (name) := include & !exclude]

  # return
  return(res)
}
