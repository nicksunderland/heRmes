#' @title filter_data_dict
#'
#' @param dict_path, str, path to the dataset.data_dictionary.csv
#' @param codes_str, list, list of lists representing UKBB column name, table entity, and search strategy list(name=, entity=, search=). 
#'   name must be a valid column name in the data_dictionary, entity a valid entity in the entity dictionary, and search either "matches"
#'   for exact matches, or starts with to match cases of multiple instances (repeated measures usually)
#'
#' @returns a filtered subset of the data_dictionary 
#'
filter_data_dict <- function(dict_path, codes_struc) {
    
    data_dict <- fread(dict_path)
    
    d <- lapply(codes_struc, function(x) {
        
        d0 <- data.table()
        if (x$search=="matches") {
            d0 <- data_dict[entity==x$entity & name==x$name]
        } else if (x$search=="startswith") {
            d0 <- data_dict[entity==x$entity & grepl(paste0("^", x$name), name)]
        }
        
        if (nrow(d0)==0) {
            cat(glue("Code [{x$name}] not found in data dictionary\n"))
            stop("Code not found error")
        }
        
        d0
        
    }) |> rbindlist(idcol = "item")
    
    return(d)
}



#' @title extract_data
#'
#' @param dataset, str, a valid dataset id - format "{projectid}:{recordid}" 
#' @param fields, str, vector of UK-BB format column names e.g. p31
#' @param entity, str, string of length one - the entity to extract from e.g. participants
#' @param output, str, the base name for the output file, no extension
#'
#' @returns NULL side effect is starting a table-exporter job which outputs the file to /hermes3_data directory in the RAP
#'
#' @importFrom glue glue
extract_data <- function(dataset, fields, entity, output) {
    
    field_str <- paste0('-ifield_names="', fields, '"', collapse=" ") 
    
    cmd <- glue::glue(
      "dx run table-exporter ",
      "-idataset_or_cohort_or_dashboard={dataset} ",
      "-ioutput={output_name} ",
      "-ioutput_format=TSV ",
      "-iheader_style=FIELD-NAME ",
      "-icoding_option=RAW ",
      "{field_str} ",
      "-ientity={entity} ",
      "--destination {output_dir}/"
    )    

    o <- system(cmd, intern = TRUE)
    cat(o, sep = "\n")
}