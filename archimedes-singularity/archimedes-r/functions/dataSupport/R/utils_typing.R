.is_base <- function(x, FUN, na.ok=TRUE) {
    as <- suppressWarnings(FUN(x))
    if (na.ok) {
        !is.na(as) | is.na(x)
    } else {
        !is.na(as)
    }
}

.is_numeric <- function(x, na.ok=TRUE) {
    # TRUE = if x val is number, string that can be made a number, or, if na.ok, NA
    .is_base(x, as.numeric, na.ok)
}

.is_logical <- function(x, na.ok=TRUE) {
    # TRUE = if x val is logical, string that can be made a logical, or, if na.ok, NA
    .is_base(x, as.logical, na.ok)
}
