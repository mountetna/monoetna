input_path <- function(
  name,
  inputs_env = Sys.getenv(),
  inputs_dir = NULL
) {
  
  if (is.null(inputs_dir)) {
    inputs_dir <- inputs_env['INPUTS_DIR']
  }

  path <- file.path(inputs_dir, name);

  if (!file.exists(path)) {
    stop("Input key ", name, " does not exist for this cell")
  }
  
  path
};

input_single_var <- function(
  name,
  inputs_env = Sys.getenv(),
  inputs_dir = NULL,
  type
) {
  scan(
    input_path(name, inputs_env, inputs_dir),
    what = type,
    quiet = TRUE)
};

input_str <- function(name, inputs_env = Sys.getenv(), inputs_dir = NULL) {
  input_single_var(name, inputs_env, inputs_dir, character())
};

input_int <- function(name, inputs_env = Sys.getenv(), inputs_dir = NULL) {
  input_single_var(name, inputs_env, inputs_dir, integer())
};

input_num <- function(name, inputs_env = Sys.getenv(), inputs_dir = NULL) {
  input_single_var(name, inputs_env, inputs_dir, numeric())
};

input_bool <- function(name, inputs_env = Sys.getenv(), inputs_dir = NULL) {
  # Checks is input, as string, is a TRUE-like string.
  # Else FALSE.
  str <- input_single_var(name, inputs_env, inputs_dir, character())
  tolower(str) %in% c("yes", "y", "true", "t", "1")
};

output_path <- function(
  name,
  outputs_env = Sys.getenv(),
  outputs_dir = NULL
) {
  
  if (is.null(outputs_dir)) {
    outputs_dir <- outputs_env['OUTPUTS_DIR']
  }

  path <- file.path(outputs_dir, name);

  if (!is.na(outputs_env['ENFORCE_OUTPUTS_EXIST']) && outputs_env['ENFORCE_OUTPUTS_EXIST']==TRUE) {
    if (!file.exists(path)) {
      stop("Output key ", name, " does not exist for this cell.")
    }
  }

  path
};

output_var <- function(
  data,
  name,
  outputs_env = Sys.getenv(),
  outputs_dir = NULL
) {
  if (is.character(data)) {
    data <- paste0("\"",data,"\"")
  }

  write(data, output_path(name, outputs_env, outputs_dir))
}
