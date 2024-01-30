#' Starting out and working with different Etna environments
#' @name authentication-and-environments
#' 
#' @section Authorization via 'token's:
#' Access to magma via magmaR is authenticated via 'token's which users can obtain from Janus.
#' A valid token must be provided to magmaR before any calls to magma can be successfully performed.
#' To provide this token, users should obtain their token from Janus, then provide it to magmaR functions using the \code{\link{magmaRset}} function.
#' See the function's own documentation and other functions' examples for further details and usage code.
#' 
#' @section Non-production Environments:
#' The code base relies on 2 different ecosystems of all of its components for purposes of
#' "development" of new tools and features, and
#' "production", the release version which most users see.
#' We won't get into the details too much more here, but users with access to the development environment can access it magmaR.
#' 
#' To do so, users should provide the \code{url} of their alternative version of magma to magmaR functions using the \code{\link{magmaRset}} function.
#' If proxy or other curl-request settings need to be adjusted, users can provide these via the \code{opts} input of this same \code{\link{magmaRset}} function.
#' At a minimum, you will likely want to use \code{opts = list(ssl_verifyhost = FALSE, ssl_verifypeer = FALSE)}.
#' 
#' @seealso
#' \code{\link{magmaRset}}
NULL
