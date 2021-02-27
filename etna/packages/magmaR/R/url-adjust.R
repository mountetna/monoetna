#' Working with different Etna environments
#' @name magma-environments
#' @section The different Etna environments:
#' The Mount Etna code base relies on 3 different ecosystems of all of its components for purposes of
#' "development" of new tools and features,
#' "staging" of updates and data prior to release, and 
#' "production", the release version which most users see.
#' We won't get into the details too much more here, but the different Etna environments are all accessible via different URLs. 
#' 
#' @section How magmaR decides which version of magma to target:
#' magmaR uses two methods to decide which version of magma to target for a given function call: \itemize{
#' \item All functions accept a \code{url.base} input which allows a user to directly provide the url for the version of magma that they wish to target.
#' \item When no \code{url.base} input is given, magmaR will initialize and utilize a \code{.MAGMAR_URL} global variable with which it sets its default \code{url.base}.
#' Unless adjusted by a user, this will be filled in with the url of production magma.
#' }
#' 
#' @section How to use this machinery to adjust your magma environment:
#' 
#' For advanced users with access to the staging and/or development versions of magma, switching can be achieved by:
#' 
#' \strong{1. For each individual call to a magmaR function:}
#' 
#' Add the input \code{url.base="<environment-magma-url>"} (with the actual url of course).
#' 
#' OR
#' 
#' \strong{2. To adjust the default behavior & set up all future calls to target a particular environment:}
#' 
#' Run \code{.MAGMAR_URL<-"<environment-magma-url>"} (with the actual url of course) a single time.
#' Afterwards, any future magmaR call which are not given an explicit \code{url.base} input will target this \code{.MAGMAR_URL}.
NULL
