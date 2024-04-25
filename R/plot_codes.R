# Silence R CMD check
globalVariables(c("group", "coding_system.name", "code"),
                package = "heRmes")


#' @title Plot phenotype code overlap
#' @param pheno_ids a list of strings, valid phenotype IDs
#' @param types a vector of strings, one or more of c("ICD10 codes", "ICD9 codes", "OPCS4 codes", "Read codes v2", "Med codes", "SNOMED  CT codes", "OXMIS codes")
#' @return a plot
#' @export
#'
plot_code_overlap <- function(pheno_ids = c("PH_HF_HERMES_3.0", "PH25", "PH129"),
                              types = c("ICD10 codes", "ICD9 codes", "OPCS4 codes", "Read codes v2",
                                        "Med codes", "SNOMED  CT codes", "OXMIS codes")) {

  stopifnot("pheno_ids must be length > 1" = length(pheno_ids) > 1)

  # gather the data
  dat <- lapply(pheno_ids, function(id) {

    # get the meta data
    meta <- get_metadata(id)

    # get the codes
    codes <- get_codes(id)[, c("code", "coding_system.name")]
    codes[, group := paste0(meta$name, " (", meta$phenotype_id, ")")]

  }) |> data.table::rbindlist()

  # list to take the plots
  code_types <- unique(dat$coding_system.name)
  code_types <- code_types[code_types %in% types]
  plot_list <- list()

  # plot each type of code separately
  for (i in seq_along(code_types)) {

    d <- dat[coding_system.name == code_types[[i]], ]
    d <- split(d, by = "group")
    d <- lapply(d, function(x) unique(x[, code]))

    # plot
    p <- plot(eulerr::euler(d, shape = "ellipse"),
              quantities = TRUE,
              labels     = FALSE,
              main       = list(label = code_types[[i]], fontsize = 8, font = 2),
              legend     = list(fontsize = 8))
    plot_list[[ code_types[[i]] ]] <- p

  }

  # combine
  if (length(plot_list) > 1) {
    p0 <- ggpubr::ggarrange(plotlist = plot_list, ncol = 2, nrow = ceiling(length(plot_list)/2))
  } else {
    p0 <- plot_list[[1]]
  }

  # return
  return(p0)
}




