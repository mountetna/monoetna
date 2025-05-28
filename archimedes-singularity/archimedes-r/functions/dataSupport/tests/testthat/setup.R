library(jsonlite)
df <- read.csv("data_frame.csv", row.names = 1)
df2 <- as.data.frame(fromJSON("data_frame.json", simplifyDataFrame = FALSE, simplifyVector = TRUE))
cons <- fromJSON("group_selection.json", simplifyDataFrame = FALSE, simplifyVector = FALSE)
cons2 <- fromJSON("group_selection2.json", simplifyDataFrame = FALSE, simplifyVector = FALSE)

test_vector <- c(1:5, "a", "b", "c", "t", "T", "true", "false", "True", "False")
test_vector_na <- c(1:5, "a", "b", "c", "t", "T", "true", "false", "True", "False", NA)


