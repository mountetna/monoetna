input_path <- function(
  name,
  inputsEnv = Sys.getenv(),
  inputsDir = NULL
) {
  
  if (is.null(inputsDir)) {
    inputsDir <- inputsEnv['INPUTS_DIR']
  }

  path <- file.path(inputsDir, name);

  if (!file.exists(path)) {
    stop("Input key ", name, " does not exist for this cell")
  }
  
  path
};

input_single_var <- function(
  name,
  inputsEnv = Sys.getenv(),
  inputsDir = NULL,
  type
) {
  scan(
    input_path(name, inputsEnv, inputsDir),
    what = type,
    quiet = TRUE)
};

input_str <- function(name, inputsEnv = Sys.getenv(), inputsDir = NULL) {
  input_single_var(name, inputsEnv, inputsDir, "character")
};

input_int <- function(name, inputsEnv = Sys.getenv(), inputsDir = NULL) {
  input_single_var(name, inputsEnv, inputsDir, "integer")
};

input_num <- function(name, inputsEnv = Sys.getenv(), inputsDir = NULL) {
  input_single_var(name, inputsEnv, inputsDir, "numeric")
};

input_bool <- function(name, inputsEnv = Sys.getenv(), inputsDir = NULL) {
  input_single_var(name, inputsEnv, inputsDir, "logical")
};

output_path <- function(
  name,
  outputsEnv = Sys.getenv(),
  outputsDir = NULL
) {
  
  if (is.null(outputsDir)) {
    outputsDir <- inputsEnv$OUTPUTS_DIR
  }

  path <- file.path(inputsDir, name);

  if (!is.na(outputsEnv['ENFORCE_OUTPUTS_EXIST']) && outputsEnv['ENFORCE_OUTPUTS_EXIST']) {
    stop("Output key ", name, " does not exist for this cell.")
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
