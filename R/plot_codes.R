# Silence R CMD check
globalVariables(c("group", "coding_system.name", "code"),
                package = "heRmes")

plot_code_overlap <- function(pheno_ids = NULL) {

  stopifnot("pheno_ids must be length > 1" = length(pheno_ids) > 1 | is.null(pheno_ids))

  # get the phenotypes in the local library
  pheno_lib <- list.files(system.file("extdata", "ukhdr_phenotypes", package = "heRmes"), pattern = ".yaml$", full.names = TRUE)

  # filter based on the input ids
  if (!is.null(pheno_ids)) {
    pheno_lib <- pheno_lib[grepl(paste0(pheno_ids, collapse = "|"), pheno_lib)]
  }

  # gather the data
  dat <- lapply(pheno_lib, function(file) {

    # get the meta data
    meta <- yaml::read_yaml(file)
    name <- paste0(meta$phenotype_id, " - ", meta$name, " (", strsplit(meta$author, ",", fixed = TRUE)[[1]][[1]], " et al.)")

    # get the codes
    codes <- data.table::fread(sub(".yaml$", "_codes.tsv", file), select = c("code", "coding_system.name"))
    codes[, group := name]

  }) |> data.table::rbindlist()

  # list to take the plots
  code_types <- unique(dat$coding_system.name)
  plot_list <- list()

  for (i in seq_along(code_types)) {

    d <- dat[coding_system.name == code_types[[i]], ]
    d <- split(d, by = "group", drop = TRUE)
    d <- lapply(d, function(x) x[, code])

    # plot
    if (length(d) > 1) {
      p <- UpSetR::upset(UpSetR::fromList(d), order.by = "freq")
      plot_list <- c(plot_list, list(p))
    }

  }

  p0 <- ggpubr::ggarrange(plotlist = plot_list, ncol = 3)

}
#
#   size <- 10
#   mat <- matrix(sample(c(0,1), size ^ 2, replace = TRUE, prob = c(0.9, 0.1)), nrow = size, ncol = size)
#   codes <- LETTERS[1:size]
#   code_colours <- sample(c("red", "green"), size, replace = TRUE)
#   colnames(mat) = rownames(mat) = LETTERS[1:size]
#   network <- igraph::graph_from_adjacency_matrix(mat, diag = FALSE)
#   igraph::V(network)$color <- code_colours
#   plot(network)
#
# }



