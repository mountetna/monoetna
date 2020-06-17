# Rtemis

Rtemis is a calculation service in R - the inputs are in
JSON format and look like so:

    {
      "func": "function_name",
      "args": [
        1, 2, 3, 4
      ]
    }

"function_name" is an endpoint function for Rtemis - this is
a normal R function that takes an arbitrary number of arguments
and returns a list() as output. The list() is transformed to
JSON and returned.

The argument list consists of numbers, strings, etc. The usual
Etna JSON types of Matrix and Vector (i.e. { matrix: { rows,
row_names, col_names } } and { vector: [ { label, value }, {
label, value } ] }) are also acceptable and will be transformed
into data.frame and vector R types respectively.

# Apache config

Rtemis is built with Rook, which can be run using rApache. The crucial line in the Location handler for the service is:

    SetHandler r-handler
    RFileEval /rtemis/path/server.R:Rook::Server$call(rook)
