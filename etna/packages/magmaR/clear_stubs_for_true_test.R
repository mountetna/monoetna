### This script is for running a live check of the package based on real calls to magma

# Update the token in '.Renviron' file, *manually*.  Adjust url too wanting to test in production / dev.
# This code just opens the file.
# It should be `TOKEN="<your_token>"`.
usethis::edit_r_environ(scope = "project")

### !!! RESTART R for the above to take effect. ###

# Remove previous casettes 
casette_path <- file.path("tests","fixtures")
casettes <- file.path(casette_path, list.files(casette_path))
file.remove(casettes, path = casette_path)
