#' Set up your magma environment and authentication 
#' @param token Single string. Your personal token from \url{https://janus.ucsf.edu}. 
#' 
#' When not explicitly given, you will be prompted to input it via the console.
#' @param url Single string. The url of the production, staging, or development version of magma that you you would like to target.
#' See \code{\link{authentication-and-environments}} for more information.
#' @param opts A named list of curl options and the values to give them (ex: \code{list(followlocation = FALSE, othersetting = 42)}).
#' Generally, most users can ignore this input, but it can be useful for adjusting proxy settings for particular development environment setup.
#' @return A list with three components: token, url, and opts.
#' @details This function compiles a list, from the given inputs, of the information needed by other \code{magmaR} functions
#' to properly route and authenticate a call to magma.
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     
#'     # THE DEFAULT:
#'     # When run in this way, it will ask you to give your token.
#'     # And the resulting $url will be the standard, production, magma url.
#'     prod <- magmaRset()
#'     print(prod)
#'          
#'     # TARGET = staging:
#'     # Give the proper url.
#'     # Again, because we are not providing our token to the call, it will ask.
#'     stage <- magmaRset(url = "https://magma-stage.ucsf.edu")
#'     print(stage)
#'     
#'     # We can also give additional curl options to the 'opts' input:
#'     prod_opts <- magmaRset(token = prod$token,
#'         opts = list(proxyport = 1234))
#'     print(prod_opts)
#'     
#'     # Now we can retrieve data with...
#'     retrieve(
#'         target = prod,
#'         projectName = "example",
#'         modelName = "rna_seq",
#'         recordNames = "all",
#'         attributeNames = "all",
#'         filter = "")
#' }
#' 
magmaRset <- function(
    token = NULL,
    url = "https://magma.ucsf.edu",
    opts = list(followlocation=FALSE)) {
    
    if (is.null(token)) {
        token <- readline(prompt = "Enter your Janus TOKEN (without quotes):")
    }
    
    url <- gsub("/$", "", url)
    
    if (is.null(opts$followlocation)) {
        message("Curl option 'followlocation = FALSE' added for the purpose of proper authenication error handling.")
        opts$followlocation <- FALSE
    }
    
    list(
        token = token,
        url = url,
        opts = opts
    )
}
