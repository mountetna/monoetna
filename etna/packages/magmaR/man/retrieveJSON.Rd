% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/retrieve.R
\name{retrieveJSON}
\alias{retrieveJSON}
\title{Download data from magma as a json, and convert to a list}
\usage{
retrieveJSON(
  projectName,
  modelName,
  recordNames = "all",
  attributeNames = "all",
  filter = "",
  page = NULL,
  pageSize = 10,
  token = .get_TOKEN(),
  ...
)
}
\arguments{
\item{projectName}{Single string. The name of the project you would like to interact with. For options, see \code{\link{retrieveProjects}}.}

\item{modelName}{Single string. The name of the subset data structure within the project, which are referred to as 'model's in magma, to interact with.
For options, see \code{\link{retrieveModels}} or https://timur.ucsf.edu/<projectName>/map.}

\item{recordNames}{Single string or string vector. Which particular sample/tube/etc. records to grab data for.
Options are "all" or any combination of individual record names. To retrieve individual options, see \code{\link{retrieveIds}}.}

\item{attributeNames}{Single string or string vector. Which features of the data to obtain.
Options are "all" or any combination of individual attribute names. To retrieve individual options, see \code{\link{retrieveAttributes}}.}

\item{filter}{String. Potential filter of the data.
Example: "<targetAttributeName>~GYN" to filter to records where <targetAttributeName> contains "GYN".

Refer to \url{https://mountetna.github.io/magma.html#retrieve} for more details.}

\item{page}{Integers. For retrieving just a portion of the data, sets slice/page size, which is equivalent to the a number of rows, and which slice to get.}

\item{pageSize}{Integers. For retrieving just a portion of the data, sets slice/page size, which is equivalent to the a number of rows, and which slice to get.}

\item{token}{Single string. Your personal token from \url{https://janus.ucsf.edu}. 

When not explicitly given in a function call: you will be prompted to input your token, one time.
This user provided token will then be stored as a hidden variable, \code{.MAGMAR_TOKEN}, in the global R environment,
and all future magmaR calls without a \code{token} explicitly provided will turn to this \code{.MAGMAR_TOKEN}.}

\item{...}{Additional parameters passed along to internal `.retrieve()`, `.query()`, or `.update()` functions,
for troubleshooting or advanced-user purposes only: \itemize{
\item \code{request.only} (Logical) & \code{json.params.only} (Logical) which stop the function before it performs any curl requests and instead outputs the values that would have been sent to magma in, either of two formats.
\item \code{verbose} (Logical) sets whether to report the status of the curl request after it is performed.
\item \code{url.base} (String) used to direct towards production versus staging versus development versions of magma. See \code{\link{magma-environments}}
}}
}
\value{
A list
}
\description{
Download data from magma as a json, and convert to a list
}
\details{
This function makes a call to magma/retrieve with \code{format = "json"}.
Then, it converts the json output into a list which is more compatible with R.
}
\examples{

if (interactive()) {
    # Running like this will ask for input of your janus token one time.
    retrieveJSON(
        projectName = "example",
        modelName = "rna_seq",
        recordNames = "all",
        attributeNames = "all",
        filter = "")
}

}
\seealso{
\code{\link{retrieve}} for similar functionality, but where the call to magma/retrieve will be made with \code{format = "tsv"} and the output is a dataframe.

\code{\link{retrieveMatrix}} for matrix data-targeted utilization of this current \code{retreiveJSON} function, followed by automated restructuring of the return into a matrix format.

\url{https://mountetna.github.io/magma.html#retrieve} for documentation of the underlying magma/retrieve function.

\code{\link{retrieveProjects}} for exploring options for the \code{projectName} input.

\code{\link{retrieveModels}}, \code{\link{retrieveIds}}, and \code{\link{retrieveAttributes}} for exploring options for the \code{modelName}, \code{recordNames}, and \code{attributeNames} inputs, respectively.
}
