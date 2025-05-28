library("testthat")
if (!requireNamespace('vcr')) {
    "Skipping tests as the vcr webmockr package is integral to testing yet unavailable."
} else {
    test_check("magmaR")
}
