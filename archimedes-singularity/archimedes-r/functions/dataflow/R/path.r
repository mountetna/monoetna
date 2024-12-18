input_path <- function(
  key # a string or index
) {

  if (is.character(key) && ! key %in% names(snakemake@input)) {
    stop("No input keyed ", key, " exists in the snakemake object for this job.")
  }
  
  path <- snakemake@input[[key]]

  if (!file.exists(path)) {
    stop("The input file for key ", key, " does not exist for this job.")
  }

  path
}

input_single_var <- function(
  key,
  type
) {
  scan(
    input_path(key),
    what = type,
    quiet = TRUE)
}

input_str <- function(key) {
  input_single_var(key, character())
}

input_int <- function(key) {
  input_single_var(key, integer())
}

input_num <- function(key) {
  input_single_var(key, numeric())
}

input_bool <- function(key) {
  # Checks is input, as string, is a TRUE-like string.
  # Else FALSE.
  str <- input_single_var(key, character())
  tolower(str) %in% c("yes", "y", "true", "t", "1")
}

input_json <- function(key) {
  jsonlite::fromJSON(input_path(key), simplifyMatrix = FALSE)
}

output_path <- function(
  key # a string or index
) {

  if (is.character(key) && ! key %in% names(snakemake@output)) {
    stop("No output keyed ", key, " exists in the snakemake object for this job.")
  }
  
  snakemake@output[[key]]
}

output_string <- function(
  data,
  key
) {
  if (is.character(data)) {
    data <- paste0("\"",data,"\"")
  }

  write(data, output_path(key))
}

output_var <- function(
  data,
  key
) {

  write(data, output_path(key))
}

output_json <- function(data, key) {
  json_str <- jsonlite::toJSON(
    data,
    null = "null",
    na = "string",
    factor = "string"
  )
  write(json_str, output_path(key))
}