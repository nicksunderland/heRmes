# Silence R CMD check
globalVariables(c("code_type", "i.code_type", "col_name", "none"), package = "heRmes")

#' @title Phenotype
#' @param x a data.frame like object or file path readable by data.table::fread
#' @param id_col a string, the name of the id column
#' @param code_cols a list of named strings, the diagnosis and procedure code column(s). List elements must be one or more of
#' `icd9`, `icd10`, `opcs`, `readv2`, `readv3`
#' @param name a string, a name for the phenotype
#' @param gsub a list of characters of length 3, pre-processing of phenotype codes and/or codes in `x`. E.g. c(".", "", "both"). NULL to turn off.
#' Element 1: a character, regex pattern.
#' Element 2: a character, replacement string.
#' Element 3: a character vector, one or more of: 'x' (apply to codes in `x` only), 'pheno' (apply to codes in pheno files only), 'both' (apply to both `x` and all pheno files), or a phenotype_id in `include` or `exclude` (apply to that phenotype).
#' @param ... other parameters passed to gsub
#' @param include a string or list of strings, a valid phenotype id (see `get_phenotypes()`)
#' @param include_multi a string, how to deal with multiple include phenotypes, either 'all' (all required) or 'any' (all required)
#' @param exclude a string or list of strings, a valid phenotype id (see `get_phenotypes()`)
#' @param exclude_multi a string, how to deal with multiple exclude phenotypes, either 'all' (all required) or 'any' (all required)
#' @param verbose a logical, print progress
#' @return a data.table
#' @export
#'
phenotype <- function(x,
                      id_col,
                      include,
                      exclude = NULL,
                      code_cols = list("ICD9 codes"    = NULL,
                                       "ICD10 codes"   = NULL,
                                       "OPCS4 codes"   = NULL,
                                       "Read codes v2" = NULL,
                                       "Read codes v3" = NULL),
                      gsub = NULL,
                      name = "overall",
                      include_multi = "any",
                      exclude_multi = "any",
                      verbose = TRUE,
                      ...) {

  if (verbose) cat("Phenotyping...\n")

  # read if needed
  if (inherits(x, "data.frame")) {

    dat <- data.table::as.data.table(x)

  } else if (file.exists(x)) {

    if (verbose) cat("[i] reading file\n")
    dat <- data.table::fread(x, select = c(id_col, unname(unlist(code_cols))))

  } else {

    stop("`x` must be a valid data.frame like object or a valid file path")

  }

  if (verbose) cat("[i] processing", nrow(x), "records\n")

  # checks
  stopifnot("`code_cols` must contain at least 1 non-NULL value" = !all(sapply(code_cols, is.null)))
  stopifnot("Invalid `code_cols` list item detected. Options: `ICD9 codes`, `ICD10 codes`, `OPCS4 codes`, `Read codes v2`, `Read codes v3`" = all(names(code_cols) %in% c("ICD9 codes", "ICD10 codes", "OPCS4 codes", "Read codes v2", "Read codes v3")))
  include_multi <- match.arg(include_multi, choices = c("any", "all"))
  exclude_multi <- match.arg(exclude_multi, choices = c("any", "all"))
  if (!is.null(gsub)) {
    gsub[[3]] <- match.arg(gsub[[3]], choices = c("x", "pheno", "both", include, exclude))
    if (gsub[[1]] == ".") warning("Match all regex `.` used, did you want to match a literal period? If so, use `\\\\.`")
  }

  # pivot longer
  if (verbose) cat("[i] pivoting data longer\n")
  code_types <- data.table::rbindlist(
    lapply(names(code_cols), function(x) {
      data.frame(code_type = x, col_name = code_cols[[x]])
    })
  )
  dat <- data.table::melt(dat,
                          id.vars       = id_col,
                          measure.vars  = unlist(code_cols),
                          variable.name = "col_name",
                          value.name    = "code",
                          na.rm         = TRUE)
  dat <- dat[code != "", ]
  dat[code_types, code_type := i.code_type, on = "col_name"]
  dat[, col_name := NULL]

  # gsub codes if needed
  if (!is.null(gsub) && any(gsub[[3]] %in% c("x", "both"))) {
    if (verbose) cat("[i] cleaning input codes with regex [", gsub[[1]], "], replacement [", gsub[[2]], "]\n")
    dat[, code := gsub(gsub[[1]], gsub[[2]], code, ...)]
  }

  # get the inclusion codes
  if (verbose) cat("[i] getting inclusion phenotype codes from PhenoID(s)", paste0(include, collapse = ", "), "\n")
  include_codes <- lapply(include, function(x) get_codes(x))
  if (is.null(names(include))) {
    names(include_codes) <- include
  }

  # get the exclusion codes
  if (!is.null(exclude)) {
    if (verbose) cat("[i] getting exclusion phenotype codes from PhenoID(s)", paste0(exclude, collapse = ", "), "\n")
    exclude_codes <- lapply(exclude, function(x) get_codes(x))
    if (is.null(names(exclude))) {
      names(exclude_codes) <- exclude
    }
  } else {
    exclude_codes <- list()
  }

  # parse the include codes
  for (i in seq_along(include_codes)) {

    if (verbose) cat("[i] assessing phenotype", include[[i]], "\n")

    # gsub codes if needed
    if (!is.null(gsub) && any(gsub[[3]] %in% c("pheno", "both", include[[i]]))) {
      if (verbose) cat("[i] cleaning phenotype codes with regex [", gsub[[1]], "], replacement [", gsub[[2]], "]\n")
      include_codes[[i]][, code := gsub(gsub[[1]], gsub[[2]], code, ...)]
      include_codes[[i]] <- include_codes[[i]][code != "", ]
    }
    dat[include_codes[[i]], names(include_codes)[i] := TRUE, on = c("code" = "code", "code_type" = "coding_system.name")]
    dat[is.na(get(names(include_codes)[i])), names(include_codes)[i] := FALSE]

  }

  # parse the exclude codes
  for (i in seq_along(exclude_codes)) {

    if (verbose) cat("[i] assessing phenotype", exclude[[i]], "\n")

    # gsub codes if needed
    if (!is.null(exclude) && !is.null(gsub) && any(gsub[[3]] %in% c("pheno", "both", exclude[[i]]))) {
      if (verbose) cat("[i] cleaning phenotype codes with regex [", gsub[[1]], "], replacement [", gsub[[2]], "]\n")
      exclude_codes[[i]][, code := gsub(gsub[[1]], gsub[[2]], code, ...)]
      exclude_codes[[i]] <- exclude_codes[[i]][code != "", ]
    }
    dat[exclude_codes[[i]], names(exclude_codes)[i] := TRUE, on = c("code" = "code", "code_type" = "coding_system.name")]
    dat[is.na(get(names(exclude_codes)[i])), names(exclude_codes)[i] := FALSE]

  }

  # summarise / compute the result
  if (verbose) cat("[i] summarising phenotyping of participants\n")
  res <- dat[, lapply(mget(c(names(include_codes), names(exclude_codes))), any), by = id_col]

  # no phenotype match
  res[, none := rowSums(.SD)==0, .SDcols = names(res)[-1]]

  # inclusions
  if (length(include_codes) > 1) {
    count <- ifelse(include_multi == "any", 0, length(include_codes) - 1)
    res[, include := rowSums(.SD) > count, .SDcols = names(include_codes)]
  } else {
    res[, include := get(names(include_codes))]
  }

  # exclusions
  if (length(names(exclude_codes)) == 0) {
    res[, exclude := FALSE]
  } else if (length(names(exclude_codes)) == 1) {
    res[, exclude := get(names(exclude_codes))]
  } else {
    count <- ifelse(exclude_multi == "any", 0, length(exclude_codes) - 1)
    res[, exclude := rowSums(.SD) > count, .SDcols = names(exclude_codes)]
  }

  # final phenotype
  res[, (name) := include & !exclude]

  # return
  if (verbose) cat("[i] finished\n")
  return(res)
}
