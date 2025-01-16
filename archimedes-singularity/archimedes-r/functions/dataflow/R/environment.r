token <- function() {
    Sys.getenv('TOKEN')
}
# In future, convert from the user's standard token to a longer-lived but read only, and this project only, version.

magma_host <- function() {
    Sys.getenv('MAGMA_HOST')
}

project_name <- function() {
    Sys.getenv('PROJECT_NAME')
}