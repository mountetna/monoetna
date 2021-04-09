### This code will be run before all test-____.R files.

### Setup relating to vcr
library("vcr")

# vcr_dir <- "../fixtures"
# if (!nzchar(Sys.getenv("TOKEN"))) {
#     if (dir.exists(vcr_dir)) {
#         Sys.setenv("TOKEN" = "foobar")
#     } else {
#         stop("No API key nor cassettes, tests cannot be run.",
#             call. = FALSE)
#     }
# }

invisible(vcr::vcr_configure(
    filter_sensitive_data = list("<<<my_token>>>" = Sys.getenv('TOKEN')),
    dir = "../fixtures"
))
vcr::check_cassette_names()
