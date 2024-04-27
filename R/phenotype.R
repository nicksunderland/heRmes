# Silence R CMD check
globalVariables(c("code_var"), package = "heRmes")

#' @title Phenotype
#' @param x a data.frame like object
#' @param id_col a string, the name of the id column
#' @param code_col a single string, or vector of strings, the diagnosis code column(s)
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
#' @return a data.table
#' @export
#'
phenotype <- function(x, id_col, code_col, include, exclude = NULL, gsub = NULL, name = "overall", include_multi = "any", exclude_multi = "any", ...) {

  # checks
  include_multi <- match.arg(include_multi, choices = c("any", "all"))
  exclude_multi <- match.arg(exclude_multi, choices = c("any", "all"))
  if (!is.null(gsub)) {
    gsub[[3]] <- match.arg(gsub[[3]], choices = c("x", "pheno", "both", include, exclude))
    if (gsub[[1]] == ".") warning("Match all regex `.` used, did you want to match a literal period? If so, use `\\\\.`")
  }

  # create the data.table
  dat <- data.table::as.data.table(x)

  # remove other data
  dat[, names(dat)[!names(dat) %in% c(id_col, code_col)] := NULL]

  # pivot longer if needed
  if (length(code_col) > 1) {
    dat <- data.table::melt(dat, id.vars = id_col, variable.name = "code_var", value.name = "code", na.rm = TRUE)
    dat[, code_var := NULL]
    dat <- dat[code != "", ]
  }

  # gsub codes if needed
  if (!is.null(gsub) && any(gsub[[3]] %in% c("x", "both"))) {
    dat[, code := gsub(gsub[[1]], gsub[[2]], code, ...)]
  }

  # get the inclusion codes
  include_codes <- lapply(include, function(x) get_codes(x))
  if (is.null(names(include))) {
    names(include_codes) <- include
  }

  # get the exclusion codes
  if (!is.null(exclude)) {
    exclude_codes <- lapply(exclude, function(x) get_codes(x))
    if (is.null(names(exclude))) {
      names(exclude_codes) <- exclude
    }
  } else {
    exclude_codes <- list()
  }

  # parse the include codes
  for (i in seq_along(include_codes)) {

    # gsub codes if needed
    if (!is.null(gsub) && any(gsub[[3]] %in% c("pheno", "both", include[i]))) {
      include_codes[[i]][, code := gsub(gsub[[1]], gsub[[2]], code, ...)]
      include_codes[[i]] <- include_codes[[i]][code != "", ]
    }
    dat[include_codes[[i]], names(include_codes)[i] := TRUE, on = "code"]
    dat[is.na(get(names(include_codes)[i])), names(include_codes)[i] := FALSE]

  }

  # parse the exclude codes
  for (i in seq_along(exclude_codes)) {

    # gsub codes if needed
    if (!is.null(exclude) && !is.null(gsub) && any(gsub[[3]] %in% c("pheno", "both", exclude[i]))) {
      exclude_codes[[i]][, code := gsub(gsub[[1]], gsub[[2]], code, ...)]
      exclude_codes[[i]] <- exclude_codes[[i]][code != "", ]
    }
    dat[exclude_codes[[i]], names(exclude_codes)[i] := TRUE, on = "code"]
    dat[is.na(get(names(exclude_codes)[i])), names(exclude_codes)[i] := FALSE]

  }

  # summarise / compute the result
  res <- dat[, lapply(mget(c(names(include_codes), names(exclude_codes))), any), by = id_col]

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
  return(res)
}
