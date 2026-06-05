### This code will be run before all test-____.R files.

### Setup relating to vcr
library("vcr")

TOKEN <- .get_sysenv_or_mock("TOKEN")
URL <- magmaRset("")$url

invisible(vcr::vcr_configure(
    filter_sensitive_data = list("<<<my_token>>>" = TOKEN),
    dir = "../fixtures"
))
