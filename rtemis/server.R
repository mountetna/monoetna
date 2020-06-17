library(Rook)
library(rjson)

source('./parse_args.R')
source('./diff_exp.R')
source('./wilcox.R')
source('./matrix.R')
source('./pca.R')
source('./sd.R')
source('./beeswarm.R')
source('./normal.R')
source('./density.R')

bad_request = function(msg) {
  list(
    status = 422,
    headers = list(
      'Content-Type' = 'application/json'
    ),
    body = toJSON(
      list( error = msg )
    )
  )
}

function_name = function(func) {
  return(paste("func", func, sep="_"))
}

has_func = function(input) {
  return("func" %in% names(input)
    && length(input$func) == 1
    && grep("^\\w+$", input$func)
    && exists(function_name(input$func)))
}

has_args = function(input) {
  return("args" %in% names(input) && class(input$args) == "list")
}

rook = function(env) {
  # what DO we want to do?
  request = Request$new(env)
  body = rawToChar(request$body()$read())
  input = fromJSON(
    body,
    simplify=FALSE
  )

  if (!has_func(input)) return(bad_request("Invalid function specified"))

  if (!has_args(input)) return(bad_request("Invalid arguments specified"))

  func = get(function_name(input$func))

  output = do.call(func,parse_args(input$args))

  list(
    status = 200,
    headers = list(
      'Content-Type' = 'application/json'
    ),
    body = toJSON(
      list(
        output=transform_outputs(output)
      )
    )
  )
}
