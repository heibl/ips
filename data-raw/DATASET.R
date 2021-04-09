## code to prepare `log_list` and `operator_list` datasets goes here

log_list <- read.csv("data-raw/log_list.csv", stringsAsFactors = FALSE)
usethis::use_data(log_list, overwrite = TRUE)

operator_list <- read.csv("data-raw/operator_list.csv", stringsAsFactors = FALSE)
usethis::use_data(operator_list, overwrite = TRUE)
